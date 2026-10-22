"""Generate source-linked findings from measured segment outputs, without phenotype inference."""
from __future__ import annotations

import copy
import itertools
from collections import Counter
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Phylo._io import read as read_tree
from Bio.Phylo.BaseTree import Tree
from Bio.Seq import Seq

from shav.common import ROOT, SEGMENTS, Row, read_fasta, read_json, read_tsv, write_json, write_tsv
from shav.phylogeny import is_strong, prune_with_support, root_tree, split_set, support, tree_plot

EUROPE = {"France", "Germany", "Switzerland", "Austria", "Belgium", "Netherlands", "Ireland", "Italy", "Spain", "Portugal", "Luxembourg", "Liechtenstein", "Jersey"}


def european_summary(tree: Tree, metadata: dict[str, Row]) -> dict[str, Any]:
    europe = {a for a, r in metadata.items() if r["group"] == "SHAV" and r["country"] in EUROPE}
    if len(europe) < 2:
        return {"tips": len(europe), "status": "insufficient_data"}
    mrca = tree.common_ancestor(sorted(europe))
    all_members = {str(t.name) for t in mrca.get_terminals()}
    others = sorted(all_members - europe)
    groups: list[dict[str, Any]] = []
    covered: set[str] = set()
    # Preorder yields maximal supported Europe-exclusive groups first.
    for node in tree.get_nonterminals(order="preorder"):
        members = {str(t.name) for t in node.get_terminals()}
        if len(members) >= 2 and members <= europe and is_strong(node) and not members <= covered:
            sh, uf = support(node)
            groups.append({"accessions": sorted(members), "sh_alrt": sh, "ufboot": uf,
                           "countries": sorted({metadata[a]["country"] for a in members})})
            covered |= members
    projects: dict[str, Any] = {}
    for name, project in [("EU1_German_project", "PRJEB123174"), ("EU2_German_project", "PRJEB125115")]:
        members = {a for a, row in metadata.items() if project in row["bioproject"].split(";")}
        if len(members) < 2:
            projects[name] = {"accessions": sorted(members), "status": "fewer_than_two_tips"}
            continue
        ancestor = tree.common_ancestor(sorted(members))
        descendants = {str(t.name) for t in ancestor.get_terminals()}
        sh, uf = support(ancestor)
        projects[name] = {"project": project, "accessions": sorted(members), "mrca_sh_alrt": sh, "mrca_ufboot": uf,
                          "other_mrca_members": sorted(descendants - members)}
    earlier = {a for a in europe if "PRJEB125115" not in metadata[a]["bioproject"].split(";")}
    if len(earlier) >= 2:
        ancestor = tree.common_ancestor(sorted(earlier))
        sh, uf = support(ancestor)
        projects["other_European_tip_set"] = {"accessions": sorted(earlier), "mrca_sh_alrt": sh, "mrca_ufboot": uf,
            "other_mrca_members": sorted({str(t.name) for t in ancestor.get_terminals()} - earlier),
            "interpretation": "Earlier European accessions in this snapshot; not a geography-only assignment of future sequences to EU1."}
    return {"tips": len(europe), "europe_exclusive_monophyly": not others,
            "other_mrca_members": others, "maximal_supported_europe_only_groups": groups,
            "europe_tips_outside_supported_europe_only_groups": sorted(europe - covered), "project_labels": projects}


def concordance() -> dict[str, Any]:
    selected: dict[str, dict[str, str]] = {}
    trees: dict[str, Tree] = {}
    strict_trees: dict[str, Tree] = {}
    for segment in SEGMENTS:
        work = ROOT / "results" / segment
        rows = read_tsv(work / "metadata.tsv")
        lengths = {a: len(s.replace("-", "")) for a, s in read_fasta(work / "backbone.fasta").items()}
        identities: dict[str, str] = {}
        for row in sorted(rows, key=lambda r: (-lengths[r["accession"]], r["accession"])):
            if not row["isolate"]:
                continue
            # Conflicting geography or collection dates remain different units.
            key = "|".join(row[k] for k in ("group", "isolate", "country", "date"))
            identities.setdefault(key, row["accession"])
        selected[segment] = identities
        trees[segment] = read_tree(work / "primary.treefile", "newick")
        strict_trees[segment] = read_tree(work / "strict.treefile", "newick")
    all_units = set.intersection(*(set(value) for value in selected.values()))
    units = {unit for unit in all_units if unit.startswith("SHAV|")}

    def unit_splits(tree: Tree, segment: str, strong: bool) -> set[frozenset[str]]:
        mapping = {selected[segment][key]: key for key in units}
        output: set[frozenset[str]] = set()
        for node in tree.get_nonterminals():
            if strong and not is_strong(node):
                continue
            side = {mapping[str(t.name)] for t in node.get_terminals() if str(t.name) in mapping}
            other = units - side
            if min(len(side), len(other)) >= 2:
                chosen = side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other
                output.add(frozenset(chosen))
        return output

    supported = {s: unit_splits(trees[s], s, True) & unit_splits(strict_trees[s], s, True) for s in SEGMENTS}
    comparisons: list[dict[str, Any]] = []
    for a, b in itertools.combinations(SEGMENTS, 2):
        conflicts = []
        for sa, sb in itertools.product(supported[a], supported[b]):
            # Four nonempty intersections define incompatible unrooted splits.
            if all((sa & sb, sa & (units - sb), (units - sa) & sb, (units - sa) & (units - sb))):
                conflicts.append({"split_a": sorted(sa), "split_b": sorted(sb)})
        comparisons.append({"segment_a": a, "segment_b": b, "stable_supported_splits_a": len(supported[a]),
                            "stable_supported_splits_b": len(supported[b]), "incompatible_splits": conflicts})
    return {"matched_units": sorted(units), "units": len(units), "all_matched_units_including_context": len(all_units),
            "accession_mapping": {s: {key: selected[s][key] for key in sorted(units)} for s in SEGMENTS},
            "comparisons": comparisons, "interpretation": "This screen is restricted to explicitly matched SHAV units. A candidate is supported discordance, not proof of reassortment, donor, direction, date, or phenotypic effect. Both primary and complete-case trees must support the conflicting splits. Unmatched or conflicting source identities are excluded."}


def variants(segment: str) -> list[Row]:
    work = ROOT / "results" / segment
    aligned = read_fasta(work / "export.aligned.fasta")
    reference = next(iter(read_fasta(work / "reference.fasta").values()))
    metadata = {r["accession"]: r for r in read_tsv(work / "metadata.tsv")}
    rows: list[Row] = []
    for accession, sequence in aligned.items():
        if metadata[accession]["group"] != "SHAV":
            continue
        for start in range(0, len(reference) - 2, 3):
            ref, query = reference[start:start + 3], sequence[start:start + 3]
            if ref == query or set(ref + query) - set("ACGT"):
                continue
            before, after = str(Seq(ref).translate()), str(Seq(query).translate())
            rows.append({"segment": segment, "accession": accession, "cds_codon_1based": str(start // 3 + 1),
                         "reference_codon": ref, "query_codon": query, "reference_aa": before, "query_aa": after,
                         "effect_category": "synonymous" if before == after else "nonsynonymous",
                         "phenotype": "not inferred"})
    return rows


def support_conclusion(summaries: list[dict[str, Any]]) -> str:
    weak: list[str] = []
    for item in summaries:
        for key, label in [("other_European_tip_set", "earlier European group"), ("EU2_German_project", "EU2")]:
            group = item["europe"]["project_labels"].get(key, {})
            if "mrca_sh_alrt" not in group:
                weak.append(f"{item['segment']} {label} (insufficient labelled tips)")
            elif group["mrca_sh_alrt"] < 80 or group["mrca_ufboot"] < 95:
                weak.append(f"{item['segment']} {label}")
    return ("Groups below the prespecified joint strong-support thresholds or without enough labelled tips: " +
            (", ".join(weak) if weak else "none") + ". Grouping is not a transmission inference.")


def human_conclusion(summary: dict[str, Any], manifest: dict[str, Any]) -> str:
    accessions = summary["explicit_human_records"]
    observed = ("Explicitly human-labelled SHAV records in this snapshot: " + ", ".join(accessions) + ".") if accessions else "No explicitly human-labelled SHAV nucleotide sequence was identified in this snapshot."
    checks = "; ".join(f"{label}: {len(check['uids'])} results (checked {check['checked_at']})"
                       for label, check in sorted(manifest["human_checks"].items()))
    return (f"{observed} Recorded NCBI keyword sensitivity checks: {checks}. "
            f"{summary['focal_missing_host']} focal records lack a host qualifier, so global absence cannot be established. "
            "Source inspection on 2026-10-22 identified no case accession in the ECDC or Inserm announcements checked that day. "
            "That dated inspection is not a claim about later announcements. IHUMI_Shav000001 (PZ860852–PZ860854) was deposited as cattle in the reviewed snapshot.")


def main() -> None:
    output = ROOT / "reports"; output.mkdir(exist_ok=True)
    metadata = read_tsv(ROOT / "data/public/metadata.tsv")
    manifest = read_json(ROOT / "data/public/manifest.json")
    focal = [r for r in metadata if r["group"] == "SHAV"]
    included = [r for r in focal if r["inclusion"] == "include"]
    summaries: list[dict[str, Any]] = []
    cluster_rows: list[Row] = []
    variation: list[Row] = []
    for segment in SEGMENTS:
        work = ROOT / "results" / segment
        measured = read_json(work / "analysis.json")
        rows = {r["accession"]: r for r in read_tsv(work / "metadata.tsv")}
        tree = root_tree(read_tree(work / "primary.treefile", "newick"), rows)
        measured["europe"] = european_summary(tree, rows)
        summaries.append(measured)
        for n, clade in enumerate(measured["europe"]["maximal_supported_europe_only_groups"], 1):
            cluster_rows.append({"segment": segment, "group": str(n), "tips": str(len(clade["accessions"])),
                                 "sh_alrt": str(clade["sh_alrt"]), "ufboot": str(clade["ufboot"]),
                                 "accessions": ";".join(clade["accessions"]), "countries": ";".join(clade["countries"])})
        # Focus figure uses the same inferred topology with context tips pruned.
        focused, provenance = prune_with_support(tree, {a for a, row in rows.items() if row["group"] == "SHAV"})
        write_tsv(work / "focused_support_provenance.tsv", provenance)
        tree_plot(focused, rows, segment, ROOT / "figures" / f"tree_{segment}_SHAV")
        variation.extend(variants(segment))
    concord = concordance()
    write_json(output / "segment_concordance.json", concord)
    write_tsv(output / "europe_clusters.tsv", cluster_rows, ["segment", "group", "tips", "sh_alrt", "ufboot", "accessions", "countries"])
    write_tsv(output / "observed_coding_variants.tsv", variation,
              ["segment", "accession", "cds_codon_1based", "reference_codon", "query_codon", "reference_aa", "query_aa", "effect_category", "phenotype"])
    summary = {"retrieved_at": manifest["retrieved_at"], "focal_records": len(focal), "canonical_focal_records": len(included),
               "explicit_human_records": [r["accession"] for r in focal if "human" in r["host"].lower() or "homo sapiens" in r["host"].lower()],
               "focal_missing_host": sum(not r["host"] for r in focal), "segments": summaries,
               "read_catalog_runs": len(read_tsv(ROOT / "data/public/read_catalog.tsv")), "country_record_counts": dict(Counter(r["country"] or "unknown" for r in included))}
    write_json(output / "summary.json", summary)
    acknowledgments: list[Row] = []
    for record in SeqIO.parse(ROOT / "data/public/records.gb", "genbank"):
        for citation in record.annotations.get("references", []):
            acknowledgments.append({"accession": str(record.id), "authors": str(citation.authors), "title": str(citation.title),
                                   "journal": str(citation.journal), "pubmed_id": str(citation.pubmed_id)})
    write_tsv(output / "sequence_acknowledgments.tsv", acknowledgments)
    lines = ["# Shamonda virus phylogenomics: baseline findings", "", f"Public snapshot retrieved {manifest['retrieved_at']}. Research analysis, not a clinical diagnostic or transmission inference.", "",
             "## Data coverage", "", f"The inventory contains {len(focal)} SHAV nucleotide records and {len(included)} canonical focal segment records, plus curated Simbu context and a Bunyamwera outgroup. These are segment/accession counts, not independent animals or genomes. The linked raw-read catalog contains {summary['read_catalog_runs']} runs; reads were not downloaded.", "",
             "Every inclusion, duplicate, original annotation and source query is auditable in [the manifest](../data/public/manifest.json), [metadata](../data/public/metadata.tsv) and [sequence acknowledgments](sequence_acknowledgments.tsv). Unknown metadata remain unknown.", "",
             "## Do European sequences cluster together?", "", "Europe is tested from deposited country metadata; EU1/EU2 project labels are kept separately from tree conclusions. A Europe-exclusive cluster requires that the MRCA contains no non-European or context tips.", "",
             "| Segment | Focal / total backbone tips | Coding sites | Europe-exclusive monophyly | Supported Europe-only groups | Model |", "|---|---:|---:|---|---:|---|"]
    for item in summaries:
        euro = item["europe"]
        lines.append(f"| {item['segment']} | {item['focal_backbone_tips']} / {item['backbone_tips']} | {item['primary_alignment_nt']} | {'Yes' if euro['europe_exclusive_monophyly'] else 'No'} | {len(euro['maximal_supported_europe_only_groups'])} | {item['best_fit_model']} |")
    lines += ["", "Groups use SH-aLRT ≥80 and UFBoot ≥95. Full-tree support is restored by unrooted split after rooting. Focused figures use the componentwise minimum across full-tree edges collapsed into each focal split, including zero for missing support; these are conservative annotations, not reestimated focal-tree supports. Each segment has focused_support_provenance.tsv. Separate Europe-only groups do not establish separate transmission chains. Sparse sampling and unresolved branches limit interpretation.", ""]
    lines += ["Project labels and the earlier European tip set are compared using measured MRCA support. Labels do not establish separate transmission chains:", "",
              "| Segment | Earlier European group, SH-aLRT / UFBoot | EU2, SH-aLRT / UFBoot |", "|---|---:|---:|"]
    for item in summaries:
        labels = item["europe"]["project_labels"]
        earlier, second = labels["other_European_tip_set"], labels["EU2_German_project"]
        lines.append(f"| {item['segment']} | {earlier['mrca_sh_alrt']:g} / {earlier['mrca_ufboot']:g} | {second['mrca_sh_alrt']:g} / {second['mrca_ufboot']:g} |")
    lines += ["", support_conclusion(summaries), ""]
    for item in summaries:
        seg, euro, sensitive = item["segment"], item["europe"], item["sensitivity"]
        lines += [f"### {seg} segment", "", f"![{seg} SHAV tree](../figures/tree_{seg}_SHAV.png)", "",
                  f"[Full context tree](../figures/tree_{seg}.pdf), [interactive dataset](../auspice/shamonda_{seg}.json), [support table](../results/{seg}/clades.tsv), [callable-site distances](../results/{seg}/distances.tsv).", "",
                  f"European focal tips: {euro['tips']}. Other members within their shared ancestor: {', '.join(euro['other_mrca_members']) or 'none'}. Primary supported unrooted splits retained in the complete-case alignment: {sensitive['primary_strong_splits_present_strict']} of {sensitive['primary_strong_splits']}. This comparison measures topology sensitivity, not an accuracy percentage.", "",
                  f"Primary sites: {item['primary_alignment_nt']}; complete-case sites: {item['strict_alignment_nt']}. Partial-CDS supplementary windows: {len(item['partial_windows'])}; their coordinates, included records and alignments are recorded in prepared.json. The backbone uses the primary CDS, so S describes N, not the full overlapping N/NSs sequence.", ""]
    conflicts = sum(len(item["incompatible_splits"]) for item in concord["comparisons"])
    lines += ["## Segment histories and evolutionary timing", "", f"The concordance analysis matches {concord['units']} explicit SHAV source-identity units across all three backbones and detects {conflicts} supported incompatible split pairs that survive complete-case masking. This conservative screen does not establish absence of historical reassortment. [The complete comparison](segment_concordance.json) retains exact accession mappings. Historical context viruses are not counted as new SHAV reassortment candidates. Discordance does not identify a donor, timing, direction or phenotypic change.", "",
              "Root-to-tip collection-date diagnostics include 100 seeded date permutations and an alternative midpoint root. These diagnostics are phylogenetically non-independent and sensitive to repeated isolates, missing dates and uneven sampling. No calibrated evolutionary rate, dated origin or migration route is reported. See each segment's analysis.json.", "",
              "## Is a human-case sequence publicly available?", "", human_conclusion(summary, manifest), "",
              "## Disease, transmission and countermeasures", "", "The [canonical evidence table](../docs/evidence.tsv) separates direct SHAV observations, related-virus evidence and hypotheses. It records the newly reported human infection, animal disease associations, gene functions and current countermeasure evidence. WOAH reports no commercially available SHAV vaccine. Evidence does not establish that an observed mutation predicts virulence, human susceptibility, transmissibility or countermeasure failure. [Observed coding differences](observed_coding_variants.tsv) are descriptive and carry no phenotype score.", "",
              "## Methods and reproducibility", "", "Independent protein-guided MAFFT alignments were backtranslated to nucleotides. The near-complete threshold is 90% of the upper-half median focal CDS length; up to 5% ambiguous nucleotides are allowed. Unannotated focal records use explicitly flagged reference-guided ORF recovery with ≥70% aligned amino-acid identity and ≥50 matched residues. Invalid translations are excluded with reasons. The primary mask retains codons with at most 50% missing sequence; the complete-case sensitivity excludes every incomplete codon. IQ-TREE ModelFinder selects the model with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates, seed 159150. Main trees use the Bunyamwera outgroup; midpoint rooting is a sensitivity comparison. Neither root establishes the geographic origin of SHAV.", "",
              "Run `make ingest`, `make build`, `make report` and `make validate`. View locally with `make view`. New cases use `make place SEQUENCES=case.fasta METADATA=case.tsv`; see [placement](../docs/placement.md). Private inputs and reports remain ignored by Git.", ""]
    (output / "findings.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()

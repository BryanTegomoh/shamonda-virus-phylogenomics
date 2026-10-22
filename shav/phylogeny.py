"""Coding-aware segment phylogenies, uncertainty checks, and local Auspice export."""
from __future__ import annotations

import argparse
import copy
import itertools
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Phylo._io import read as read_tree, write as write_tree
from Bio.Phylo._utils import draw as draw_tree
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import linregress

from shav.common import ROOT, SEGMENTS, Row, read_fasta, read_json, read_tsv, run, sha256, write_fasta, write_json, write_tsv


@dataclass(frozen=True)
class CodingSequence:
    accession: str
    nucleotide: str
    protein: str
    location: str
    product: str
    strand: int
    codon_start: int


def extract_cds(record: SeqRecord, segment: str) -> CodingSequence:
    """Extract the annotated primary CDS, respecting strand and partial codon start."""
    candidates = [f for f in record.features if f.type == "CDS" and f.location is not None]
    if segment == "S":
        named = [f for f in candidates if any(
            term in " ".join(f.qualifiers.get("product", []) + f.qualifiers.get("gene", [])).lower()
            for term in ("nucleocapsid", "nucleoprotein")
        )]
        if named:
            candidates = named
    if not candidates:
        raise ValueError("No annotated primary CDS")
    feature = max(candidates, key=lambda f: len(f))
    if feature.location is None:
        raise ValueError("Primary CDS has no location")
    offset = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
    if offset not in (0, 1, 2):
        raise ValueError("Invalid annotated codon_start")
    nucleotide = str(feature.extract(record.seq)).upper()[offset:]
    nucleotide = nucleotide[:len(nucleotide) // 3 * 3]
    if len(nucleotide) < 3:
        raise ValueError("Empty CDS")
    if set(nucleotide) - set("ACGTRYSWKMBDHVN"):
        raise ValueError("Invalid nucleotide alphabet")
    protein = str(Seq(nucleotide).translate())
    if protein.endswith("*"):
        nucleotide, protein = nucleotide[:-3], protein[:-1]
    if "*" in protein:
        raise ValueError("Internal stop codon in primary CDS")
    return CodingSequence(str(record.id), nucleotide, protein, str(feature.location),
                          " ".join(feature.qualifiers.get("product", [])),
                          int(feature.location.strand or 1), offset + 1)


def backtranslate(proteins: dict[str, str], coding: dict[str, CodingSequence]) -> dict[str, str]:
    aligned: dict[str, str] = {}
    lengths = {len(v) for v in proteins.values()}
    if len(lengths) != 1:
        raise ValueError("Protein alignment has inconsistent lengths")
    for accession, amino in proteins.items():
        cds = coding[accession]
        if amino.replace("-", "") != cds.protein:
            raise ValueError(f"Protein alignment changed sequence: {accession}")
        offset, codons = 0, []
        for residue in amino:
            if residue == "-":
                codons.append("---")
            else:
                codons.append(cds.nucleotide[offset:offset + 3])
                offset += 3
        if offset != len(cds.nucleotide):
            raise ValueError(f"Backtranslation failed: {accession}")
        aligned[accession] = "".join(codons)
    return aligned


def infer_unannotated_cds(record: SeqRecord, template: CodingSequence) -> CodingSequence:
    """Recover a homologous ORF analytically, without treating it as a deposited annotation.

    Only focal SHAV records use this fallback. Both strands and all reading frames
    compete against an annotated SHAV protein; a low-identity match is rejected.
    Coordinates are zero-based half-open on the oriented sequence, stated explicitly.
    """
    aligner = PairwiseAligner(mode="local", match_score=2, mismatch_score=-2,
                              open_gap_score=-8, extend_gap_score=-.5)
    best: tuple[float, str, int, int, int, int] | None = None
    sequence = str(record.seq).upper()
    for strand, oriented in [(1, sequence), (-1, str(Seq(sequence).reverse_complement()))]:
        for frame in range(3):
            nucleotides = oriented[frame:]
            peptide = str(Seq(nucleotides[:len(nucleotides) // 3 * 3]).translate())
            aa_offset = 0
            for orf in peptide.split("*"):
                if len(orf) >= 30:
                    if aligner.score(template.protein, orf) <= 0:
                        aa_offset += len(orf) + 1
                        continue
                    alignment = aligner.align(template.protein, orf)[0]
                    starts, ends = alignment.aligned[1][0], alignment.aligned[1][-1]
                    start, end = int(starts[0]), int(ends[1])
                    pairs = [(int(a), int(b)) for a, b in zip(alignment.indices[0], alignment.indices[1], strict=True) if a >= 0 and b >= 0]
                    identity = sum(template.protein[a] == orf[b] for a, b in pairs) / len(pairs) if pairs else 0
                    if identity >= .70 and (best is None or float(alignment.score) > best[0]):
                        nt_start, nt_end = frame + (aa_offset + start) * 3, frame + (aa_offset + end) * 3
                        best = (float(alignment.score), oriented[nt_start:nt_end], strand, nt_start, nt_end, len(pairs))
                aa_offset += len(orf) + 1
    if best is None:
        raise ValueError("No credible SHAV primary-CDS homolog in six reading frames")
    _, nucleotide, strand, start, end, matched = best
    if matched < 50:
        raise ValueError("Too little protein homology to recover unannotated CDS")
    protein = str(Seq(nucleotide).translate())
    return CodingSequence(str(record.id), nucleotide, protein,
        f"inferred_from={template.accession};oriented_0based_halfopen=[{start}:{end}];strand={strand}",
        "reference-guided primary-CDS homolog (not a submitted annotation)", strand, 1)


def mask_codons(aligned: dict[str, str], maximum_missing: float) -> tuple[dict[str, str], list[int]]:
    if not aligned or len({len(s) for s in aligned.values()}) != 1:
        raise ValueError("Invalid coding alignment")
    length = len(next(iter(aligned.values())))
    if length % 3:
        raise ValueError("Alignment is not codon aligned")
    columns: list[int] = []
    for start in range(0, length, 3):
        missing = sum(any(n not in "ACGT" for n in s[start:start + 3]) for s in aligned.values())
        if missing / len(aligned) <= maximum_missing:
            columns.extend(range(start, start + 3))
    if len(columns) < 150:
        raise ValueError("Insufficient comparable coding sites after masking")
    return {key: "".join(s[i] for i in columns) for key, s in aligned.items()}, columns


def distance(a: str, b: str) -> tuple[int, int, float | None]:
    if len(a) != len(b):
        raise ValueError("Distance requires aligned sequences of equal length")
    pairs = [(x, y) for x, y in zip(a, b, strict=True) if x in "ACGT" and y in "ACGT"]
    differences = sum(x != y for x, y in pairs)
    return differences, len(pairs), differences / len(pairs) if pairs else None


def support(clade: Clade) -> tuple[float, float]:
    """IQ-TREE labels are SH-aLRT/UFBoot when both tests are enabled."""
    label = str(clade.name or "")
    if "/" not in label:
        return 0.0, 0.0
    try:
        sh, uf = label.split("/")[:2]
        return float(sh), float(uf)
    except ValueError:
        return 0.0, 0.0


def is_strong(clade: Clade) -> bool:
    sh, uf = support(clade)
    return sh >= 80 and uf >= 95


def split_set(tree: Tree, universe: set[str] | None = None, strong_only: bool = False) -> set[frozenset[str]]:
    tips = {str(t.name) for t in tree.get_terminals()}
    universe = universe or tips
    result: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals():
        if strong_only and not is_strong(clade):
            continue
        part = {str(t.name) for t in clade.get_terminals()} & universe
        other = universe - part
        if min(len(part), len(other)) < 2:
            continue
        # Canonicalize an unrooted split, including equally sized partitions.
        side = min((part, other), key=lambda s: (len(s), sorted(s)))
        result.add(frozenset(side))
    return result


def root_tree(tree: Tree, metadata: dict[str, Row]) -> Tree:
    result = copy.deepcopy(tree)
    outgroups = [tip for tip in result.get_terminals()
                 if "bunyamwera" in metadata[str(tip.name)]["group"].lower()]
    if outgroups:
        result.root_with_outgroup(outgroups[0])
    else:
        result.root_at_midpoint()
    return result


def build_tree(alignment: Path, prefix: Path, settings: dict[str, Any]) -> Path:
    output = Path(str(prefix) + ".treefile")
    signature = {"alignment_sha256": sha256(alignment), "settings": settings}
    signature_path = Path(str(prefix) + ".signature.json")
    if output.exists() and signature_path.exists() and read_json(signature_path) == signature:
        return output
    run(["iqtree", "-s", str(alignment), "-m", "MFP", "-B", str(settings["bootstrap"]),
         "--alrt", str(settings["shalrt"]), "--seed", str(settings["seed"]),
         "-T", str(settings["threads"]), "--prefix", str(prefix), "--keep-ident", "--redo"],
        log=Path(str(prefix) + ".console.log"))
    if not output.exists():
        raise RuntimeError(f"IQ-TREE did not produce {output}")
    write_json(signature_path, signature)
    return output


def prepare(segment: str) -> None:
    settings = read_json(ROOT / "config/analysis.json")
    work = ROOT / "results" / segment
    work.mkdir(parents=True, exist_ok=True)
    all_metadata = read_tsv(ROOT / "data/public/metadata.tsv")
    metadata = {r["accession"]: r for r in all_metadata if r["segment"] == segment and r["inclusion"] == "include"}
    records = {str(r.id): r for r in SeqIO.parse(ROOT / "data/public/records.gb", "genbank")}
    template_id = {"L": "NC_018463.1", "M": "NC_018467.1", "S": "NC_018464.1"}[segment]
    template = extract_cds(records[template_id], segment)
    coding: dict[str, CodingSequence] = {}
    qc: list[Row] = []
    for accession, row in metadata.items():
        try:
            try:
                cds = extract_cds(records[accession], segment)
            except ValueError as exc:
                if str(exc) != "No annotated primary CDS" or row["group"] != "SHAV":
                    raise
                cds = infer_unannotated_cds(records[accession], template)
            ambiguity = sum(x not in "ACGT" for x in cds.nucleotide) / len(cds.nucleotide)
            if ambiguity > float(settings["maximum_ambiguity_fraction"]):
                raise ValueError("Excess ambiguous nucleotides in primary CDS")
            coding[accession] = cds
            qc.append({"accession": accession, "status": "eligible", "reason": "", "length": str(len(cds.nucleotide)),
                       "location": cds.location, "strand": str(cds.strand), "codon_start": str(cds.codon_start), "product": cds.product})
        except ValueError as exc:
            qc.append({"accession": accession, "status": "excluded", "reason": str(exc), "length": "", "location": "",
                       "strand": "", "codon_start": "", "product": ""})
    focal_lengths = [len(c.nucleotide) for a, c in coding.items() if metadata[a]["group"] == "SHAV"]
    if not focal_lengths:
        raise RuntimeError(f"No eligible focal CDS for {segment}")
    # Robust upper-half median avoids setting a complete-CDS threshold from short amplicons.
    median_length = statistics.median(sorted(focal_lengths)[len(focal_lengths) // 2:])
    minimum = median_length * float(settings["minimum_cds_fraction"])
    backbone = {a for a, cds in coding.items() if len(cds.nucleotide) >= minimum}
    for row in qc:
        if row["accession"] in coding:
            row["status"] = "backbone" if row["accession"] in backbone else "partial"
    if len(backbone) < 4:
        raise RuntimeError(f"Fewer than four backbone sequences for {segment}")
    write_tsv(work / "coding_qc.tsv", qc)
    write_fasta(work / "cds.fasta", ((a, c.nucleotide) for a, c in sorted(coding.items())))
    write_fasta(work / "proteins.fasta", ((a, c.protein) for a, c in sorted(coding.items())))
    run(["mafft", "--auto", "--thread", str(settings["threads"]), str(work / "proteins.fasta")],
        stdout=work / "proteins.aligned.fasta", log=work / "mafft.log")
    aligned = backtranslate(read_fasta(work / "proteins.aligned.fasta"), coding)
    write_fasta(work / "all.coding.aligned.fasta", sorted(aligned.items()))
    backbone_aligned = {a: aligned[a] for a in sorted(backbone)}
    primary, primary_columns = mask_codons(backbone_aligned, float(settings["primary_maximum_codon_missing_fraction"]))
    strict, strict_columns = mask_codons(backbone_aligned, float(settings["sensitivity_maximum_codon_missing_fraction"]))
    write_fasta(work / "backbone.fasta", sorted(primary.items()))
    write_fasta(work / "strict.fasta", sorted(strict.items()))
    write_tsv(work / "metadata.tsv", [metadata[a] for a in sorted(backbone)])
    # Group partial CDS by its actual aligned span. Do not fabricate a common region
    # across nonoverlapping amplicons. Each window explicitly reports its coordinates.
    windows: dict[tuple[int, int], list[str]] = {}
    for accession in sorted(set(coding) - backbone):
        covered = [i for i, base in enumerate(aligned[accession]) if base in "ACGT"]
        if len(covered) < int(settings["minimum_partial_overlap_nt"]):
            continue
        start, end = covered[0] // 3 * 3, (covered[-1] // 3 + 1) * 3
        windows.setdefault((start, end), []).append(accession)
    partial_outputs: list[dict[str, Any]] = []
    for number, ((start, end), focal_partial) in enumerate(sorted(windows.items()), 1):
        subset = {a: sequence[start:end] for a, sequence in aligned.items()
                  if sum(c in "ACGT" for c in sequence[start:end]) >= .9 * (end - start)}
        if len(subset) < 4:
            continue
        window, columns = mask_codons(subset, .5)
        name = f"partial_{number}"
        write_fasta(work / f"{name}.fasta", sorted(window.items()))
        partial_outputs.append({"name": name, "alignment_start_1based": start + 1, "alignment_end_1based": end,
                                "retained_nt": len(columns), "partial_accessions": focal_partial, "tips": len(subset)})
    write_json(work / "prepared.json", {"segment": segment, "backbone_accessions": sorted(backbone),
               "primary_columns": primary_columns, "strict_columns": strict_columns,
               "minimum_cds_nt": minimum, "partial_windows": partial_outputs,
               "input_sha256": {"genbank": sha256(ROOT / "data/public/records.gb"), "metadata": sha256(ROOT / "data/public/metadata.tsv")}})


def infer(segment: str) -> None:
    work = ROOT / "results" / segment
    settings = read_json(ROOT / "config/analysis.json")
    prepared = read_json(work / "prepared.json")
    build_tree(work / "backbone.fasta", work / "primary", settings)
    build_tree(work / "strict.fasta", work / "strict", settings)
    for window in prepared["partial_windows"]:
        build_tree(work / f"{window['name']}.fasta", work / window["name"], settings)
    write_json(work / "inferred.json", {"segment": segment, "settings": settings, "prepared_sha256": sha256(work / "prepared.json")})


def temporal_test(tree: Tree, metadata: dict[str, Row], settings: dict[str, Any]) -> dict[str, Any]:
    """Descriptive root-to-tip randomization, with repeat isolate labels collapsed."""
    rows: list[tuple[str, float, float]] = []
    seen: set[str] = set()
    for tip in tree.get_terminals():
        row = metadata[str(tip.name)]
        if row["group"] != "SHAV" or not row["date"]:
            continue
        isolate = row["isolate"]
        if not isolate or isolate in seen:
            continue
        date = decimal_date(row["date"])
        if date is None:
            continue
        seen.add(isolate)
        rows.append((str(tip.name), date, float(tree.distance(tip))))
    if len(rows) < 5 or len({row[1] for row in rows}) < 3:
        return {"status": "insufficient_metadata", "dated_unique_isolates": len(rows), "rate_estimate": None,
                "reason": "Need at least five distinct labelled isolates and three collection times for this diagnostic; no dated tree or evolutionary rate inferred."}
    x = np.array([r[1] for r in rows]); y = np.array([r[2] for r in rows])
    fit = linregress(x, y)
    generator = np.random.default_rng(int(settings["seed"]))
    observed = float(fit.rvalue) ** 2
    randomized = [float(linregress(generator.permutation(x), y).rvalue) ** 2
                  for _ in range(int(settings["temporal_permutations"]))]
    p = (1 + sum(value >= observed for value in randomized)) / (len(randomized) + 1)
    return {"status": "diagnostic_only", "dated_unique_isolates": len(rows), "r_squared": observed,
            "permutation_p": p, "positive_slope": float(fit.slope) > 0, "permutations": len(randomized),
            "rate_estimate": None, "reason": "Root-to-tip correlation is descriptive, phylogenetically non-independent, and sensitive to sampling/root. A rate or dated origin requires a separate calibrated clock analysis; neither is reported."}


def decimal_date(value: str) -> float | None:
    from datetime import date
    try:
        chunks = value.split("-")
        year = int(chunks[0])
        if len(chunks) == 1:
            return year + .5
        month = int(chunks[1])
        day = int(chunks[2]) if len(chunks) > 2 else 15
        actual = date(year, month, day)
        return year + (actual - date(year, 1, 1)).days / (date(year + 1, 1, 1) - date(year, 1, 1)).days
    except (ValueError, OverflowError):
        return None


def clade_rows(tree: Tree, metadata: dict[str, Row]) -> list[Row]:
    result: list[Row] = []
    for i, clade in enumerate(tree.get_nonterminals()):
        members = [str(t.name) for t in clade.get_terminals()]
        focal = [name for name in members if metadata[name]["group"] == "SHAV"]
        if len(focal) < 2:
            continue
        sh, uf = support(clade)
        result.append({"node": f"node_{i}", "sh_alrt": f"{sh:g}", "ufboot": f"{uf:g}",
                       "strong": str(is_strong(clade)).lower(), "focal_tips": str(len(focal)), "total_tips": str(len(members)),
                       "countries": ";".join(sorted({metadata[n]["country"] or "unknown" for n in focal})),
                       "hosts": ";".join(sorted({metadata[n]["host"] or "unknown" for n in focal})),
                       "accessions": ";".join(members)})
    return result


def tree_plot(tree: Tree, metadata: dict[str, Row], segment: str, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colours = {"France": "#a33a32", "Germany": "#bc7635", "Switzerland": "#b39043", "Japan": "#376b8a",
               "South Africa": "#497d58", "Nigeria": "#73578e"}
    fig, ax = plt.subplots(figsize=(13, max(7, len(tree.get_terminals()) * .19)))
    tree.ladderize()
    focal_only = all(metadata[str(t.name)]["group"] == "SHAV" for t in tree.get_terminals())
    def tip_label(clade: Clade) -> str | None:
        if not clade.is_terminal():
            return None
        row = metadata[str(clade.name)]
        suffix = " [EU2]" if "PRJEB125115" in row["bioproject"].split(";") else (" [EU1]" if "PRJEB123174" in row["bioproject"].split(";") else "")
        group = "" if focal_only else f" | {row['group']}"
        return f"{clade.name}{group} | {row['country'] or 'unknown'}{suffix}"
    draw_tree(tree, axes=ax, do_show=False,
               label_func=tip_label,
               branch_labels=lambda c: "•" if is_strong(c) else None,
               label_colors=lambda name: colours.get(metadata.get(str(name).split(" | ")[0], {}).get("country", ""), "#333333"))
    ax.set_title(f"Shamonda virus {segment}: primary coding-sequence ML tree\nDots: SH-aLRT ≥80 and UFBoot ≥95; branch length: substitutions/site", loc="left", fontsize=12)
    ax.set_ylabel("")
    ax.set_yticks([])
    maximum_depth = max(float(tree.distance(t)) for t in tree.get_terminals())
    ax.set_xlim(-maximum_depth * .02, maximum_depth * 1.65)
    ax.set_xlabel("Substitutions per nucleotide site")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def export(segment: str) -> None:
    work = ROOT / "results" / segment
    prepared = read_json(work / "prepared.json")
    metadata = {r["accession"]: r for r in read_tsv(work / "metadata.tsv")}
    tree = read_tree(work / "primary.treefile", "newick")
    strict = read_tree(work / "strict.treefile", "newick")
    rooted = root_tree(tree, metadata)
    midpoint = copy.deepcopy(tree); midpoint.root_at_midpoint()
    primary_splits, strict_splits = split_set(tree), split_set(strict)
    strong = split_set(tree, strong_only=True)
    sensitivity = {"primary_splits": len(primary_splits), "strict_splits": len(strict_splits),
                   "shared_splits": len(primary_splits & strict_splits),
                   "primary_strong_splits": len(strong), "primary_strong_splits_present_strict": len(strong & strict_splits),
                   "rooting_primary": "Bunyamwera outgroup" if any("bunyamwera" in r["group"].lower() for r in metadata.values()) else "midpoint",
                   "alternative_root": "midpoint", "rooting_changes_unrooted_splits": split_set(rooted) != split_set(midpoint)}
    # Stable names preserve a separate support table while allowing Augur ancestral inference.
    rows = clade_rows(rooted, metadata)
    write_tsv(work / "clades.tsv", rows, ["node", "sh_alrt", "ufboot", "strong", "focal_tips", "total_tips", "countries", "hosts", "accessions"])
    for i, node in enumerate(rooted.get_nonterminals()):
        node.name = f"NODE_{i:05d}"; node.confidence = None
    write_tree(rooted, work / "tree.nwk", "newick")
    aligned = read_fasta(work / "backbone.fasta")
    pairs: list[Row] = []
    for a, b in itertools.combinations(sorted(aligned), 2):
        mismatches, callable_sites, fraction = distance(aligned[a], aligned[b])
        pairs.append({"accession_a": a, "accession_b": b, "mismatches": str(mismatches),
                      "callable_sites": str(callable_sites), "p_distance": f"{fraction:.8f}" if fraction is not None else ""})
    write_tsv(work / "distances.tsv", pairs)
    # Export in ungapped reference-CDS coordinates, preserving the full alignment
    # separately for inference and distances. These coordinate systems are explicit.
    focal = [a for a in sorted(metadata) if metadata[a]["group"] == "SHAV"]
    preferred = {"L": "NC_018463.1", "M": "NC_018467.1", "S": "NC_018464.1"}[segment]
    ref = preferred if preferred in aligned else max(focal, key=lambda a: sum(c in "ACGT" for c in aligned[a]))
    # Use the original full CDS alignment, not the inference mask, for Nextclade.
    full = read_fasta(work / "all.coding.aligned.fasta")
    export_columns = [i for i, base in enumerate(full[ref]) if base != "-"]
    export_aligned = {a: "".join(full[a][i] for i in export_columns) for a in aligned}
    reference = export_aligned[ref]
    if set(reference) - set("ACGT"):
        raise ValueError("Nextclade reference must be unambiguous")
    write_fasta(work / "export.aligned.fasta", sorted(export_aligned.items()))
    write_fasta(work / "reference.fasta", [(ref, reference)])
    run(["augur", "ancestral", "--tree", str(work / "tree.nwk"), "--alignment", str(work / "export.aligned.fasta"),
         "--root-sequence", str(work / "reference.fasta"), "--keep-ambiguous", "--keep-overhangs", "--seed", "159150",
         "--output-node-data", str(work / "ancestral.json"), "--output-sequences", str(work / "ancestral.fasta")], log=work / "ancestral.log")
    node_data = read_json(work / "ancestral.json")
    # Derive divergence from the IQ-TREE branches. No clock or inferred dates.
    for clade in rooted.find_clades():
        name = str(clade.name)
        node_data["nodes"].setdefault(name, {})["branch_length"] = float(clade.branch_length or 0.0)
    write_json(work / "ancestral.json", node_data)
    dataset = ROOT / "nextclade" / segment
    dataset.mkdir(parents=True, exist_ok=True)
    config = {"title": f"Shamonda virus {segment} coding-sequence phylogeny", "maintainers": [{"name": "Bryan Tegomoh"}],
              "colorings": [{"key": key, "title": title, "type": "categorical"} for key, title in
                            [("group", "Virus"), ("country", "Country"), ("host", "Host"), ("date", "Collection date")]],
              "filters": ["group", "country", "host"], "panels": ["tree"],
              "display_defaults": {"distance_measure": "div", "color_by": "group"}}
    write_json(work / "auspice_config.json", config)
    output = ROOT / "auspice" / f"shamonda_{segment}.json"
    run(["augur", "export", "v2", "--tree", str(work / "tree.nwk"), "--metadata", str(work / "metadata.tsv"),
         "--metadata-id-columns", "accession", "--node-data", str(work / "ancestral.json"),
         "--auspice-config", str(work / "auspice_config.json"), "--include-root-sequence-inline", "--output", str(output)], log=work / "export.log")
    data = read_json(output)
    # Augur export preserves public metadata only. Embed annotation for the primary CDS.
    gene = {"L": "RdRp", "M": "M_polyprotein", "S": "N"}[segment]
    data["meta"]["genome_annotations"] = {"nuc": {"start": 1, "end": len(reference), "strand": "+"},
        gene: {"start": 1, "end": len(reference), "strand": "+", "type": "CDS"}}
    data["meta"]["description"] = (f"Independent {segment} segment analysis using its primary CDS. Branches show genetic divergence, not time. "
        "Public sampling is sparse and uneven. Genetic clustering does not establish transmission, virulence, or source attribution. "
        "Full methods, sequence acknowledgments and support tables are in the private repository.")
    write_json(output, data)
    write_fasta(dataset / "reference.fasta", [(ref, reference)])
    write_json(dataset / "tree.json", data)
    (dataset / "genome_annotation.gff3").write_text(f"##gff-version 3\n{ref}\tGenBank\tCDS\t1\t{len(reference)}\t.\t+\t0\tID={gene};Name={gene}\n")
    write_json(dataset / "pathogen.json", {"schemaVersion": "3.0.0", "attributes": {"name": f"SHAV {segment} primary CDS"},
        "files": {"reference": "reference.fasta", "pathogenJson": "pathogen.json", "genomeAnnotation": "genome_annotation.gff3", "treeJson": "tree.json"}})
    (dataset / "README.md").write_text(f"# SHAV {segment} local dataset\n\nReference: {ref}, primary CDS {gene}.\n\nCoordinates are primary-CDS coordinates, not whole-segment coordinates.\nUse through `make place`; lineage labels are not assigned.\n")
    temporal = temporal_test(root_tree(tree, metadata), metadata, read_json(ROOT / "config/analysis.json"))
    temporal["midpoint_sensitivity"] = temporal_test(midpoint, metadata, read_json(ROOT / "config/analysis.json"))
    model = "unknown"
    for line in (work / "primary.iqtree").read_text().splitlines():
        if "Best-fit model according to BIC:" in line:
            model = line.split(":", 1)[1].strip()
    tree_plot(root_tree(tree, metadata), metadata, segment, ROOT / "figures" / f"tree_{segment}")
    write_json(work / "analysis.json", {"segment": segment, "backbone_tips": len(metadata),
               "focal_backbone_tips": len(focal), "primary_alignment_nt": len(next(iter(aligned.values()))),
               "strict_alignment_nt": len(next(iter(read_fasta(work / "strict.fasta").values()))),
               "best_fit_model": model, "reference": ref, "reference_cds_length": len(reference),
               "partial_windows": prepared["partial_windows"], "sensitivity": sensitivity, "temporal": temporal,
               "sha256": {"alignment": sha256(work / "backbone.fasta"), "tree": sha256(work / "primary.treefile"), "auspice": sha256(output)}})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "infer", "export"])
    parser.add_argument("--segment", required=True, choices=SEGMENTS)
    args = parser.parse_args()
    {"prepare": prepare, "infer": infer, "export": export}[args.stage](args.segment)


if __name__ == "__main__":
    main()

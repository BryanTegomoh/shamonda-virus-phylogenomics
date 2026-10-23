"""Validate snapshot integrity, tree/sequence consistency, output schemas and privacy."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from Bio.Phylo._io import read as read_tree

from shav.common import COMMUNITY_PREFIX, ROOT, SEGMENTS, read_fasta, read_json, read_tsv, run, sha256, write_json
from shav.ingest import run as verify_snapshot


def mutation_consistency(dataset: dict[str, Any], reference: str, observed: dict[str, str]) -> None:
    """Independently replay tree mutations, including root mutations, in CDS coordinates."""
    def visit(node: dict[str, Any], inherited: str) -> None:
        sequence = list(inherited)
        for mutation in node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []):
            before, position, after = mutation[0], int(mutation[1:-1]) - 1, mutation[-1]
            if not 0 <= position < len(sequence) or sequence[position] != before:
                raise ValueError(f"Inconsistent tree mutation {mutation} at {node['name']}")
            sequence[position] = after
        state = "".join(sequence)
        if node["name"] in observed:
            for actual, expected in zip(state, observed[node["name"]], strict=True):
                if expected in "ACGT" and actual != expected:
                    raise ValueError(f"Tree/alignment disagreement at {node['name']}")
        for child in node.get("children", []):
            visit(child, state)
    visit(dataset["tree"], reference)


def interactive_divergence_consistency(dataset: dict[str, Any]) -> None:
    """Check cumulative inferred nuc event counts independently of ML branch length."""
    def visit(node: dict[str, Any], parent: int, root: bool = False) -> None:
        expected = parent if root else parent + sum(
            event[0] in "ACGT" and event[-1] in "ACGT" and event[0] != event[-1]
            for event in node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []))
        if node.get("node_attrs", {}).get("div") != expected:
            raise ValueError(f"Interactive mutation-count divergence differs at {node['name']}")
        for child in node.get("children", []):
            visit(child, expected)
    visit(dataset["tree"], 0, True)


def placement_tree_agreement(display: Path, placement: Path) -> dict[str, Any]:
    if sha256(display) != sha256(placement):
        raise ValueError("Placement tree differs from validated display tree")
    data: dict[str, Any] = read_json(placement)
    return data


def repository_privacy() -> None:
    paths = ["case-inputs/probe.fasta", "case-results/probe.json", "private/probe.tsv"]
    result = subprocess.run(["git", "check-ignore", "--no-index", *paths], cwd=ROOT, capture_output=True, text=True, check=True)
    if set(result.stdout.splitlines()) != set(paths):
        raise ValueError("Private data directories are not all ignored")
    tracked = subprocess.run(["git", "ls-files", "--", "case-inputs", "case-results", "private"], cwd=ROOT, capture_output=True, text=True, check=True)
    if tracked.stdout.strip():
        raise ValueError("Private data are tracked by Git")


def main() -> None:
    verify_snapshot()
    repository_privacy()
    receipts: dict[str, Any] = {}
    for segment in SEGMENTS:
        work = ROOT / "results" / segment
        tree = read_tree(work / "tree.nwk", "newick")
        tips = {str(t.name) for t in tree.get_terminals()}
        metadata = {r["accession"] for r in read_tsv(work / "metadata.tsv")}
        aligned = read_fasta(work / "export.aligned.fasta")
        if tips != metadata or tips != set(aligned):
            raise ValueError(f"Tree/alignment/metadata identifiers differ: {segment}")
        data = placement_tree_agreement(ROOT / "auspice" / f"shamonda_{segment}.json", ROOT / "nextclade" / segment / "tree.json")
        if sha256(ROOT / "auspice" / f"{COMMUNITY_PREFIX}_{segment}.json") != sha256(ROOT / "auspice" / f"shamonda_{segment}.json"):
            raise ValueError(f"Community export differs from validated display tree: {segment}")
        reference = next(iter(read_fasta(ROOT / "nextclade" / segment / "reference.fasta").values()))
        if data.get("root_sequence", {}).get("nuc") != reference:
            raise ValueError(f"Dataset reference differs from tree reference: {segment}")
        mutation_consistency(data, reference, aligned)
        interactive_divergence_consistency(data)
        run(["augur", "validate", "export-v2", str(ROOT / "nextclade" / segment / "tree.json")], log=work / "schema-validation.log")
        prepared, analysis = read_json(work / "prepared.json"), read_json(work / "analysis.json")
        if prepared["input_sha256"]["genbank"] != sha256(ROOT / "data/public/records.gb"):
            raise ValueError(f"Stale GenBank analysis: {segment}")
        if prepared["input_sha256"]["metadata"] != sha256(ROOT / "data/public/metadata.tsv"):
            raise ValueError(f"Stale metadata analysis: {segment}")
        for name, path in [("alignment", work / "backbone.fasta"), ("tree", work / "primary.treefile"),
                           ("auspice", ROOT / "auspice" / f"shamonda_{segment}.json")]:
            if analysis["sha256"][name] != sha256(path):
                raise ValueError(f"Output checksum differs: {segment}/{name}")
        receipts[segment] = {"tips": len(tips), "schema": "passed", "mutation_replay": "passed", "snapshot": "matched", "placement_tree": "byte-identical and schema-validated", "community_export": "byte-identical to validated display tree", "interactive_divergence": "inferred nucleotide event counts verified"}
    write_json(ROOT / "reports" / "integrity-validation.json", receipts)
    print("Snapshot, three schemas, mutation replay, source/output hashes and private-file exclusions passed")


if __name__ == "__main__":
    main()

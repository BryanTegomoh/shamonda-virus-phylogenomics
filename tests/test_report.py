from io import StringIO

import pytest
from Bio.Phylo._io import read as read_tree

from shav.common import read_fasta
from shav.report import european_summary
from shav.validate import mutation_consistency


def test_europe_exclusive_monophyly():
    metadata = {a: {"group": "SHAV", "country": country, "bioproject": ""}
                for a, country in [("fr", "France"), ("de", "Germany"), ("jp", "Japan"), ("za", "South Africa")]}
    together = read_tree(StringIO("((fr:1,de:1)'90/99':1,(jp:1,za:1):1);"), "newick")
    summary = european_summary(together, metadata)
    assert summary["europe_exclusive_monophyly"] is True
    assert len(summary["maximal_supported_europe_only_groups"]) == 1
    apart = read_tree(StringIO("((fr:1,jp:1)'90/99':1,(de:1,za:1):1);"), "newick")
    summary = european_summary(apart, metadata)
    assert summary["europe_exclusive_monophyly"] is False
    assert summary["maximal_supported_europe_only_groups"] == []


def test_context_tip_breaks_exclusive_cluster():
    metadata = {a: {"group": group, "country": country, "bioproject": ""} for a, group, country in
                [("fr", "SHAV", "France"), ("de", "SHAV", "Germany"), ("context", "Sathuperi", ""), ("out", "SHAV", "Japan")]}
    tree = read_tree(StringIO("((fr:1,de:1,context:1):1,out:1);"), "newick")
    assert european_summary(tree, metadata)["europe_exclusive_monophyly"] is False


def test_mutation_replay_includes_root_and_detects_wrong_ancestor():
    dataset = {"tree": {"name": "root", "branch_attrs": {"mutations": {"nuc": ["A1C"]}},
                        "children": [{"name": "tip", "branch_attrs": {"mutations": {"nuc": ["A2G"]}}}]}}
    mutation_consistency(dataset, "AAA", {"tip": "CGA"})
    with pytest.raises(ValueError, match="Inconsistent"):
        mutation_consistency(dataset, "CAA", {"tip": "CGA"})


def test_bare_fasta_header_fails_clearly(tmp_path):
    path = tmp_path / "bad.fasta"; path.write_text(">\nACGT\n")
    with pytest.raises(ValueError, match="Empty FASTA"):
        read_fasta(path)

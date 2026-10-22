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


def test_placement_copy_mismatch_fails(tmp_path):
    from shav.validate import placement_tree_agreement
    display, placement = tmp_path / 'display.json', tmp_path / 'placement.json'
    display.write_text('{"tree": {"name": "root"}}')
    placement.write_bytes(display.read_bytes())
    assert placement_tree_agreement(display, placement)['tree']['name'] == 'root'
    placement.write_text('{"tree": {"name": "stale"}}')
    with pytest.raises(ValueError, match='Placement tree differs'):
        placement_tree_agreement(display, placement)


def test_human_conclusion_tracks_hosts_and_search_results():
    from shav.report import human_conclusion
    manifest = {'human_checks': {'human': {'uids': ['123'], 'checked_at': '2026-11-01'}}}
    text = human_conclusion({'explicit_human_records': ['AB123456.1'], 'focal_missing_host': 2}, manifest)
    assert 'AB123456.1' in text and '1 results' in text
    assert 'No explicitly human-labelled' not in text
    assert 'Source inspection on 2026-10-22' in text
    text = human_conclusion({'explicit_human_records': [], 'focal_missing_host': 0}, manifest)
    assert 'No explicitly human-labelled' in text


def test_support_conclusion_tracks_refreshed_thresholds():
    from shav.report import support_conclusion
    group = {'mrca_sh_alrt': 100, 'mrca_ufboot': 100}
    summaries = [{'segment': 'L', 'europe': {'project_labels': {'other_European_tip_set': dict(group), 'EU2_German_project': dict(group)}}}]
    assert 'none.' in support_conclusion(summaries)
    summaries[0]['europe']['project_labels']['EU2_German_project']['mrca_ufboot'] = 94
    assert 'L EU2' in support_conclusion(summaries)
    assert 'earlier European group' not in support_conclusion(summaries)


def test_divergence_validation_rejects_per_site_values():
    from shav.validate import interactive_divergence_consistency
    tip = {'name': 'tip', 'node_attrs': {'div': 1}, 'branch_attrs': {'mutations': {'nuc': ['A1C']}}}
    data = {'tree': {'name': 'root', 'node_attrs': {'div': 0}, 'children': [tip]}}
    interactive_divergence_consistency(data)
    tip['node_attrs']['div'] = .01
    with pytest.raises(ValueError, match='mutation-count divergence'):
        interactive_divergence_consistency(data)

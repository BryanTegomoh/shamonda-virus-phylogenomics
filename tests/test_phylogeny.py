from io import StringIO

import pytest
from Bio.Phylo._io import read as read_tree
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from shav.phylogeny import CodingSequence, backtranslate, decimal_date, distance, extract_cds, infer_unannotated_cds, mask_codons, split_set, support


def test_cds_uses_annotated_reverse_strand():
    rec = SeqRecord(Seq("TTATTTCAT"), id="reverse")
    rec.features = [SeqFeature(SimpleLocation(0, 9, strand=-1), type="CDS", qualifiers={"product": ["nucleocapsid protein"]})]
    result = extract_cds(rec, "S")
    assert (result.nucleotide, result.protein, result.strand) == ("ATGAAA", "MK", -1)


def test_internal_stop_raises():
    rec = SeqRecord(Seq("ATGTAGAAATAA"), id="stop")
    rec.features = [SeqFeature(SimpleLocation(0, 12), type="CDS")]
    with pytest.raises(ValueError, match="Internal stop"):
        extract_cds(rec, "L")


def test_missing_annotation_raises():
    with pytest.raises(ValueError, match="No annotated"):
        extract_cds(SeqRecord(Seq("ATGAAA"), id="empty"), "L")


def test_partial_codon_start():
    rec = SeqRecord(Seq("AATGAAATAA"), id="partial")
    rec.features = [SeqFeature(SimpleLocation(0, 10), type="CDS", qualifiers={"codon_start": ["2"]})]
    assert extract_cds(rec, "L").protein == "MK"


def test_backtranslation_preserves_codon_and_gap():
    cds = CodingSequence("a", "ATGAAA", "MK", "", "", 1, 1)
    assert backtranslate({"a": "M-K"}, {"a": cds}) == {"a": "ATG---AAA"}
    with pytest.raises(ValueError, match="changed sequence"):
        backtranslate({"a": "M-R"}, {"a": cds})


def test_masks_entire_codons_and_preserves_threshold_boundary():
    a = "AAA" * 60 + "ANC"
    masked, columns = mask_codons({"a": a, "b": "AAA" * 61}, 0)
    assert len(columns) == 180 and len(masked["a"]) == 180
    relaxed, _ = mask_codons({"a": a, "b": "AAA" * 61}, .5)
    assert len(relaxed["a"]) == 183
    with pytest.raises(ValueError, match="Insufficient"):
        mask_codons({"a": "NNN" * 60}, 0)


def test_distance_does_not_treat_missing_as_matches():
    assert distance("ACGN--", "ATGAAA") == (1, 3, 1 / 3)
    assert distance("NNN", "---") == (0, 0, None)
    with pytest.raises(ValueError, match="equal length"):
        distance("A", "AAA")


def test_support_order_and_both_thresholds():
    tree = read_tree(StringIO("((a:1,b:1)'79/100':1,(c:1,d:1)'90/96':1);"), "newick")
    assert support(tree.get_nonterminals()[1]) == (79, 100)
    assert len(split_set(tree, strong_only=True)) == 1


def test_unrooted_splits_are_invariant_to_root():
    tree = read_tree(StringIO("((a:1,b:2):1,(c:1,d:2):1,e:1);"), "newick")
    before = split_set(tree)
    tree.root_with_outgroup("a")
    assert split_set(tree) == before


def test_dates_keep_missing_and_uncertainty_explicit():
    assert decimal_date("2002") == 2002.5
    assert 2026 < decimal_date("2026-07") < 2027
    assert decimal_date("2026-02-31") is None
    assert decimal_date("unknown") is None


def test_unannotated_public_reference_recovery_on_both_strands():
    from Bio import SeqIO
    from shav.common import ROOT
    record = next(r for r in SeqIO.parse(ROOT / "data/public/records.gb", "genbank") if r.id == "NC_018464.1")
    template = extract_cds(record, "S")
    for strand in (1, -1):
        seq = Seq(template.nucleotide)
        query = SeqRecord(seq if strand == 1 else seq.reverse_complement(), id="query")
        observed = infer_unannotated_cds(query, template)
        assert observed.nucleotide == template.nucleotide
        assert observed.strand == strand
    with pytest.raises(ValueError, match="No credible"):
        infer_unannotated_cds(SeqRecord(Seq("N" * 300), id="bad"), template)


@pytest.mark.parametrize('midpoint', [False, True])
def test_rooting_preserves_support_by_split(midpoint):
    from shav.phylogeny import edge_supports, root_tree
    tree = read_tree(StringIO("(((a:1,b:1)'84/79':1,c:1)'100/95':2,(d:1,e:1)'99/100':1,o:8);"), 'newick')
    rows = {str(t.name): {'group': 'Bunyamwera' if t.name == 'o' else 'SHAV'} for t in tree.get_terminals()}
    assert edge_supports(root_tree(tree, rows, midpoint=midpoint)) == edge_supports(tree)
    assert tree.get_nonterminals()[1].name == '100/95'


def test_duplicate_root_edge_and_conflict():
    from shav.phylogeny import edge_supports, root_tree
    tree = read_tree(StringIO("((a:1,b:1)'84/79':1,(c:1,d:1)'84/79':1);"), 'newick')
    rows = {str(t.name): {'group': 'SHAV'} for t in tree.get_terminals()}
    assert len(edge_supports(tree)) == 1
    assert edge_supports(root_tree(tree, rows)) == edge_supports(tree)
    tree.get_nonterminals()[2].name = '100/100'
    with pytest.raises(ValueError, match='Conflicting support'):
        root_tree(tree, rows)


def test_pruning_collapsed_path_preserves_all_support_provenance():
    from shav.phylogeny import prune_with_support
    tree = read_tree(StringIO("(((a:1,b:1)'70/100':1,x:1)'100/100':1,c:1,d:1);"), 'newick')
    focused, provenance = prune_with_support(tree, {'a', 'b', 'c', 'd'})
    assert support(focused.common_ancestor(['a', 'b'])) == (70, 100)
    assert {r['source_sh_alrt'] for r in provenance} == {'70', '100'}
    assert len(provenance) == 2
    with pytest.raises(ValueError, match='four existing'):
        prune_with_support(tree, {'a', 'b', 'c', 'absent'})
    tree.common_ancestor(['a', 'b']).name = None
    focused, _ = prune_with_support(tree, {'a', 'b', 'c', 'd'})
    assert support(focused.common_ancestor(['a', 'b'])) == (0, 0)


def test_interactive_counts_reversions_excludes_ambiguity_and_root_differences():
    from shav.phylogeny import mutation_count_divergence
    from shav.validate import mutation_consistency
    leaf = {'name': 'tip', 'branch_attrs': {'mutations': {'nuc': ['C1A', 'A2N']}}}
    dataset = {'meta': {}, 'tree': {'name': 'root', 'branch_attrs': {'mutations': {'nuc': ['A1C']}}, 'children': [leaf]}}
    mutation_count_divergence(dataset)
    assert dataset['tree']['node_attrs']['div'] == 0
    assert leaf['node_attrs']['div'] == 1
    mutation_consistency(dataset, 'AAA', {'tip': 'ANA'})

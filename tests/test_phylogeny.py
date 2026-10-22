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

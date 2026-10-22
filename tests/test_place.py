from pathlib import Path

import pytest
from Bio.Seq import Seq

from shav.common import read_json, write_fasta, write_tsv
from shav.place import execute, inputs, place_one


def test_sample_segment_identifier_and_override(tmp_path: Path) -> None:
    write_fasta(tmp_path / 'q.fa', [('x|S', 'ACGT'), ('other', 'ACGT')])
    write_tsv(tmp_path / 'q.tsv', [{'sample': 'x', 'segment': 'S', 'sequence_id': ''},
                                  {'sample': 'x', 'segment': 'M', 'sequence_id': 'other'}])
    assert [r[1] for r in inputs(tmp_path / 'q.fa', tmp_path / 'q.tsv')] == ['x|S', 'other']


@pytest.mark.parametrize('rows', [
    [{'sample': 'x', 'segment': 'Q'}],
    [{'sample': 'x', 'segment': 'S'}, {'sample': 'x', 'segment': 'S'}],
    [{'sample': 'missing', 'segment': 'S'}],
])
def test_invalid_mapping(tmp_path: Path, rows: list[dict[str, str]]) -> None:
    write_fasta(tmp_path / 'q.fa', [('x', 'ACGT')])
    write_tsv(tmp_path / 'q.tsv', rows)
    with pytest.raises(ValueError):
        inputs(tmp_path / 'q.fa', tmp_path / 'q.tsv')


def test_batch_reports_record_errors_and_missing_segments(tmp_path: Path) -> None:
    write_fasta(tmp_path / 'q.fa', [('x|S', 'Z'), ('x|M', 'Z')])
    write_tsv(tmp_path / 'q.tsv', [{'sample': 'x', 'segment': 'S'}, {'sample': 'x', 'segment': 'M'}])
    assert not execute(tmp_path / 'q.fa', tmp_path / 'q.tsv', tmp_path / 'out')
    report = read_json(tmp_path / 'out/report.json')
    assert len(report['errors']) == 2
    assert report['missing_segments'] == {'x': ['L']}
    with pytest.raises(ValueError, match='already exists'):
        execute(tmp_path / 'q.fa', tmp_path / 'q.tsv', tmp_path / 'out')


def test_low_coverage_withheld(tmp_path: Path) -> None:
    reference = 'ATG' + 'GCTGATAACTGTGGTGAA' * 20
    write_fasta(tmp_path / 'nextclade/S/reference.fasta', [('ref', reference)])
    for name in ('pathogen.json', 'tree.json'):
        (tmp_path / 'nextclade/S' / name).write_text('{}')
    result = place_one({'sample': 'x', 'segment': 'S'}, 'x', str(Seq(reference[:240]).reverse_complement()), tmp_path / 'query', tmp_path)
    assert result['status'] == 'rejected_qc'
    assert result['strand'] == -1
    assert result['assignment'] is None

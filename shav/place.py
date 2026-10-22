"""Private, coding-aware query comparison and local phylogenetic placement."""
from __future__ import annotations

import argparse
import datetime
import subprocess
from pathlib import Path
from typing import Any

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from shav.common import ROOT, SEGMENTS, Row, read_fasta, read_json, read_tsv, run, sha256, write_fasta, write_json
from shav.phylogeny import CodingSequence, backtranslate, build_tree, distance, infer_unannotated_cds, mask_codons


def inputs(sequences: Path, metadata: Path) -> list[tuple[Row, str, str]]:
    fasta = read_fasta(sequences)
    rows = read_tsv(metadata)
    seen: set[tuple[str, str]] = set()
    used: set[str] = set()
    result: list[tuple[Row, str, str]] = []
    for row in rows:
        sample, segment = row.get('sample', ''), row.get('segment', '')
        if not sample or segment not in SEGMENTS:
            raise ValueError('Metadata requires nonempty sample and segment L, M or S')
        key = (sample, segment)
        if key in seen:
            raise ValueError(f'Duplicate sample, segment: {key}')
        seen.add(key)
        identifier = row.get('sequence_id') or (f'{sample}|{segment}' if f'{sample}|{segment}' in fasta else sample)
        if identifier not in fasta or identifier in used:
            raise ValueError(f'Missing or reused FASTA identifier: {identifier}')
        used.add(identifier)
        result.append((row, identifier, fasta[identifier]))
    if set(fasta) != used:
        raise ValueError('FASTA records without matching metadata')
    if not rows:
        raise ValueError('No metadata records')
    return result


def rebuild(segment: str, query: CodingSequence, directory: Path, root: Path = ROOT) -> Path:
    """Realign proteins, backtranslate and independently infer all coding sites."""
    coding: dict[str, CodingSequence] = {}
    public = read_fasta(root / 'results' / segment / 'cds.fasta')
    backbone = set(read_fasta(root / 'results' / segment / 'backbone.fasta'))
    for accession in sorted(backbone):
        nt = public[accession]
        coding[accession] = CodingSequence(accession, nt, str(Seq(nt).translate()), '', '', 1, 1)
    coding[query.accession] = query
    write_fasta(directory / 'rebuild.proteins.fasta', ((a, c.protein) for a, c in coding.items()))
    run(['mafft', '--auto', '--thread', '2', str(directory / 'rebuild.proteins.fasta')],
        stdout=directory / 'rebuild.proteins.aligned.fasta', log=directory / 'rebuild.mafft.log')
    aligned = backtranslate(read_fasta(directory / 'rebuild.proteins.aligned.fasta'), coding)
    masked, _ = mask_codons(aligned, .5)
    write_fasta(directory / 'rebuild.coding.fasta', sorted(masked.items()))
    return build_tree(directory / 'rebuild.coding.fasta', directory / 'rebuild',
                      {'bootstrap': 1000, 'shalrt': 1000, 'seed': 159150, 'threads': 2})


def place_one(row: Row, identifier: str, sequence: str, directory: Path,
              root: Path = ROOT, force_rebuild: bool = False) -> dict[str, Any]:
    segment = row['segment']
    directory.mkdir(parents=True, exist_ok=True)
    if set(sequence) - set('ACGTRYSWKMBDHVN'):
        raise ValueError('Invalid nucleotide alphabet, use ungapped IUPAC DNA')
    dataset = root / 'nextclade' / segment
    ref_id, reference = next(iter(read_fasta(dataset / 'reference.fasta').items()))
    template = CodingSequence(ref_id, reference, str(Seq(reference).translate()), '', '', 1, 1)
    # Internal IDs prevent public-accession collisions and unsafe filename interpolation.
    query = infer_unannotated_cds(SeqRecord(Seq(sequence), id='QUERY'), template)
    coverage = len(query.nucleotide) / len(reference)
    ambiguous = sum(b not in 'ACGT' for b in query.nucleotide) / len(query.nucleotide)
    result: dict[str, Any] = {'sample': row['sample'], 'segment': segment, 'sequence_id': identifier,
        'status': 'unresolved', 'metadata': row, 'cds_fraction': coverage, 'ambiguity_fraction': ambiguous,
        'strand': query.strand, 'cds_location': query.location, 'coordinate_system': 'primary CDS',
        'reference': ref_id, 'assignment': None, 'warnings': [],
        'dataset_sha256': {name: sha256(dataset / name) for name in ('reference.fasta', 'tree.json', 'pathogen.json')}}
    if coverage < .9 or ambiguous > .05:
        result['warnings'] = ['Insufficient CDS coverage (<90%) or excessive ambiguity (>5%); placement withheld']
        result['status'] = 'rejected_qc'
        return result
    write_fasta(directory / 'query.cds.fasta', [('QUERY', query.nucleotide)])
    run(['nextclade', 'run', '--input-dataset', str(dataset), '--jobs', '2',
         '--output-json', str(directory / 'nextclade.json'), '--output-fasta', str(directory / 'aligned.fasta'),
         '--output-tree', str(directory / 'placement.auspice.json'), str(directory / 'query.cds.fasta')],
        log=directory / 'nextclade.log')
    data = read_json(directory / 'nextclade.json')
    if data.get('errors') or not data.get('results'):
        raise ValueError(f'Nextclade query failure: {data.get("errors")}')
    next_result = data['results'][0]
    aligned = read_fasta(directory / 'aligned.fasta')['QUERY']
    public = read_fasta(root / 'results' / segment / 'export.aligned.fasta')
    public_metadata = {r['accession']: r for r in read_tsv(root / 'results' / segment / 'metadata.tsv')}
    comparisons: list[dict[str, Any]] = []
    for accession, nt in public.items():
        differences, callable_sites, fraction = distance(aligned, nt)
        comparisons.append({'accession': accession, 'differences': differences, 'callable_sites': callable_sites,
                            'p_distance': fraction, 'context': public_metadata[accession]})
    eligible = [c for c in comparisons if c['p_distance'] is not None and c['callable_sites'] >= .9 * len(reference)]
    if not eligible:
        result['warnings'] = ['No public comparison has 90% callable CDS overlap']
        return result
    minimum = min(c['p_distance'] for c in eligible)
    closest = [c for c in eligible if abs(c['p_distance'] - minimum) < 1e-12]
    result.update({'closest_public': closest, 'comparisons': comparisons, 'nextclade': next_result,
                   'status': 'compared', 'interpretation': 'Genetic similarity and local tree placement only. No transmission, virulence, source or definitive clade assignment.'})
    ambiguous_placement = len(closest) > 1 or minimum > .05 or next_result.get('qc', {}).get('overallStatus') == 'bad'
    if ambiguous_placement or force_rebuild:
        result['warnings'] = ['Ambiguous or divergent relationship; full coding-sequence rebuild required. Tied nearest sequences remain tied.']
        tree = rebuild(segment, query, directory, root)
        result['rebuild_tree'] = str(tree)
        result['rebuild_sha256'] = sha256(tree)
        result['status'] = 'rebuilt_ambiguous' if ambiguous_placement else 'rebuilt'
    return result


def execute(sequences: Path, metadata: Path, output: Path, root: Path = ROOT,
            force_rebuild: bool = False) -> bool:
    if output.exists():
        raise ValueError(f'Output already exists: {output}; choose a fresh directory')
    output.mkdir(parents=True)
    report: dict[str, Any] = {'status': 'failed', 'results': [], 'errors': [], 'missing_segments': {},
                             'privacy': 'Local case output. Do not commit or upload.',
                             'inputs_sha256': {'sequences': sha256(sequences), 'metadata': sha256(metadata)}}
    try:
        records = inputs(sequences, metadata)
        samples = {r['sample'] for r, _, _ in records}
        report['missing_segments'] = {sample: sorted(set(SEGMENTS) - {r['segment'] for r, _, _ in records if r['sample'] == sample}) for sample in sorted(samples)}
        for number, (row, identifier, sequence) in enumerate(records):
            try:
                result = place_one(row, identifier, sequence, output / f'query_{number:04d}', root, force_rebuild)
                report['results'].append(result)
                if result['status'] in ('unresolved', 'rejected_qc'):
                    report['errors'].append({'sequence_id': identifier, 'error': result['warnings']})
            except (ValueError, OSError, subprocess.CalledProcessError) as exc:
                report['errors'].append({'sequence_id': identifier, 'error': str(exc)})
    except (ValueError, OSError) as exc:
        report['errors'].append({'error': str(exc)})
    report['status'] = 'complete' if not report['errors'] else 'failed'
    write_json(output / 'report.json', report)
    return not report['errors']


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequences', required=True, type=Path)
    parser.add_argument('--metadata', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=ROOT / 'case-results' / datetime.datetime.now(datetime.UTC).strftime('%Y%m%dT%H%M%S%fZ'))
    parser.add_argument('--rebuild', action='store_true')
    args = parser.parse_args()
    try:
        success = execute(args.sequences, args.metadata, args.output, force_rebuild=args.rebuild)
    except (ValueError, OSError) as exc:
        parser.exit(1, f'{exc}\n')
    parser.exit(0 if success else 1, f'Report: {args.output / "report.json"}\n')


if __name__ == '__main__':
    main()

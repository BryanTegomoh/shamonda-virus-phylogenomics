"""Public-only held-out-isolate validation. Run with the project runtime."""
from __future__ import annotations

import shutil
from typing import Any

from Bio.Phylo._io import read as read_tree, write as write_tree

from shav.common import ROOT, read_fasta, read_json, read_tsv, run, write_fasta, write_json, write_tsv
from shav.phylogeny import build_tree, mutation_count_divergence, root_tree
from shav.place import place_one
from shav.validate import interactive_divergence_consistency, mutation_consistency

BASE = ROOT / 'results/placement-validation'
SETTINGS = {'bootstrap': 1000, 'shalrt': 1000, 'seed': 159150, 'threads': 2}


def validate(segment: str) -> dict[str, Any]:
    directory = BASE / segment
    directory.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(ROOT / 'results' / segment / 'metadata.tsv')
    metadata = {r['accession']: r for r in rows}
    candidates = [r for r in rows if r['group'] == 'SHAV' and r['isolate'] and not r['accession'].startswith('NC_')]
    held = candidates[-1]
    aliases = {r['accession'] for r in rows if r['isolate'] == held['isolate'] and r['isolate']}
    alignment = read_fasta(ROOT / 'results' / segment / 'backbone.fasta')
    training = {a: nt for a, nt in alignment.items() if a not in aliases}
    write_fasta(directory / 'training.fasta', sorted(training.items()))
    write_fasta(directory / 'full.fasta', sorted(alignment.items()))
    train_tree = build_tree(directory / 'training.fasta', directory / 'training', SETTINGS)
    full_tree = build_tree(directory / 'full.fasta', directory / 'full', SETTINGS)
    rooted = root_tree(read_tree(train_tree, 'newick'), metadata)
    actual_training = {str(t.name) for t in rooted.get_terminals()}
    if actual_training != set(training) or aliases & actual_training:
        raise ValueError('Held-out alias leakage or missing training tips')
    for number, node in enumerate(rooted.get_nonterminals()):
        node.name = f'NODE_{number:05d}'
        node.confidence = None
        node.comment = None
    write_tree(rooted, directory / 'training.nwk', 'newick')
    sandbox = directory / 'sandbox'
    dataset = sandbox / 'nextclade' / segment
    dataset.mkdir(parents=True, exist_ok=True)
    for name in ('reference.fasta', 'pathogen.json', 'genome_annotation.gff3'):
        shutil.copyfile(ROOT / 'nextclade' / segment / name, dataset / name)
    exported = read_fasta(ROOT / 'results' / segment / 'export.aligned.fasta')
    write_fasta(directory / 'training.export.fasta', ((a, exported[a]) for a in training))
    write_tsv(directory / 'metadata.tsv', [metadata[a] for a in training])
    run(['augur', 'ancestral', '--tree', str(directory / 'training.nwk'), '--alignment', str(directory / 'training.export.fasta'),
         '--root-sequence', str(dataset / 'reference.fasta'), '--keep-ambiguous', '--keep-overhangs', '--seed', '159150',
         '--output-node-data', str(directory / 'ancestral.json')], log=directory / 'ancestral.log')
    run(['augur', 'export', 'v2', '--tree', str(directory / 'training.nwk'), '--metadata', str(directory / 'metadata.tsv'),
         '--metadata-id-columns', 'accession', '--node-data', str(directory / 'ancestral.json'), '--auspice-config',
         str(ROOT / 'results' / segment / 'auspice_config.json'), '--include-root-sequence-inline', '--output', str(dataset / 'tree.json')], log=directory / 'export.log')
    tree_data = read_json(dataset / 'tree.json')
    tree_data['meta']['genome_annotations'] = read_json(ROOT / 'nextclade' / segment / 'tree.json')['meta']['genome_annotations']
    mutation_count_divergence(tree_data)
    mutation_consistency(tree_data, next(iter(read_fasta(dataset / 'reference.fasta').values())), {a: exported[a] for a in training})
    interactive_divergence_consistency(tree_data)
    write_json(dataset / 'tree.json', tree_data)
    run(['augur', 'validate', 'export-v2', str(dataset / 'tree.json')], log=directory / 'schema-validation.log')
    work = sandbox / 'results' / segment
    work.mkdir(parents=True, exist_ok=True)
    write_fasta(work / 'export.aligned.fasta', ((a, exported[a]) for a in training))
    write_fasta(work / 'backbone.fasta', sorted(training.items()))
    cds = read_fasta(ROOT / 'results' / segment / 'cds.fasta')
    write_fasta(work / 'cds.fasta', ((a, cds[a]) for a in training))
    write_tsv(work / 'metadata.tsv', [metadata[a] for a in training])
    write_fasta(directory / 'heldout.fasta', [('validation', cds[held['accession']])])
    write_tsv(directory / 'heldout.tsv', [{'sample': 'validation', 'segment': segment}])
    result = place_one({'sample': 'validation', 'segment': segment}, 'validation', cds[held['accession']], directory / 'query', sandbox, public_accession=held['accession'])
    full = read_tree(full_tree, 'newick')
    target = next(t for t in full.get_terminals() if t.name == held['accession'])
    distances = {str(t.name): float(full.distance(target, t)) for t in full.get_terminals() if t.name not in aliases}
    minimum = min(distances.values())
    full_nearest = sorted(a for a, d in distances.items() if abs(d - minimum) < 1e-7)
    comparison_nearest = [c['accession'] for c in result.get('closest_public', [])]
    attachment = result['nextclade']['nearestNodeName']
    attached_node = next(n for n in rooted.find_clades() if n.name == attachment)
    attachment_distances = {str(t.name): float(rooted.distance(attached_node, t)) for t in rooted.get_terminals()}
    attachment_minimum = min(attachment_distances.values())
    placement_nearest = sorted(a for a, d in attachment_distances.items() if abs(d - attachment_minimum) < 1e-7)
    output: dict[str, Any] = {'segment': segment, 'heldout': held, 'removed_explicit_isolate_aliases': sorted(aliases),
        'training_tree_rebuilt': True, 'full_tree_independently_rebuilt': True,
        'training_tip_count': len(actual_training), 'full_tip_count': len(alignment), 'excluded_aliases_present_training_tree': [], 'full_tree_nearest': full_nearest,
        'nextclade_attachment': attachment, 'nextclade_attachment_nearest': placement_nearest,
        'placement_full_nearest_overlap': bool(set(placement_nearest) & set(full_nearest)),
        'tie_compatible_broad_relationship': bool(set(placement_nearest) & set(comparison_nearest)) and bool(set(full_nearest) & set(comparison_nearest)),
        'callable_distance_nearest': comparison_nearest, 'nearest_sets_overlap': bool(set(full_nearest) & set(comparison_nearest)),
        'placement': result, 'interpretation': 'A broad relationship check, not exact-node placement accuracy. Identical sequences and weak branches can make node identity non-identifiable.'}
    write_json(directory / 'metrics.json', output)
    return output


if __name__ == '__main__':
    summaries = [validate(s) for s in ('S', 'M', 'L')]
    write_json(BASE / 'summary.json', [{
        'segment': d['segment'], 'heldout_accession': d['heldout']['accession'],
        'metrics': f"{d['segment']}/metrics.json",
        'training_tip_count': d['training_tip_count'], 'full_tip_count': d['full_tip_count'],
        'placement_full_nearest_overlap': d['placement_full_nearest_overlap'],
        'tie_compatible_broad_relationship': d['tie_compatible_broad_relationship'],
        'placement_status': d['placement']['status'],
    } for d in summaries])

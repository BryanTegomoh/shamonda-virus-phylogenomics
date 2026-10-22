# Local query placement

Run `make place SEQUENCES=/absolute/path/case.fasta METADATA=/absolute/path/case.tsv`.
The TSV requires `sample` and `segment` (L, M, S). Optional columns are
`sequence_id`, `date`, `host`, and `country`. For multiple segments from one sample,
use FASTA IDs `sample|segment`, or give explicit unique `sequence_id` values.
A single segment may use its sample name. Metadata and FASTA must match exactly.

Output defaults to a fresh ignored `case-results/` directory. Do not upload or
commit case data. `python -m shav.place` also accepts `--output` and `--rebuild`.
Existing output directories are rejected to prevent accidental overwrites.

The reference datasets represent primary CDS coordinates: N (S), M polyprotein,
and RdRp (L). Whole deposited segments, coding sequences, and reverse complements
are accepted through six-frame homolog recovery. Recovery requires at least 70%
protein identity and 50 matched residues. Routine placement additionally requires
90% of reference CDS length and at most 5% ambiguous bases. Short amplicons remain
unresolved, with a nonzero batch exit. No missing segment is inferred.

Nextclade aligns each eligible CDS and places it on its segment-specific tree.
The JSON report includes per-record QC, exact mismatches and callable nucleotide
denominators for every public comparison, all tied nearest accessions, public
context, dataset hashes, missing segments and uncertainties. Ambiguous bases and
gaps do not enter callable-site denominators. A comparison needs at least 90%
callable reference-CDS overlap. Ties, distances above 5%, or bad Nextclade QC trigger
a protein MAFFT alignment, codon backtranslation and IQ-TREE rebuild with model
selection, 1,000 ultrafast bootstraps, 1,000 SH-aLRT replicates and seed 159150.
The 5% trigger is an operational conservative threshold, not a validated species
or lineage boundary. Rebuild completion does not resolve a biological tie.

These are genetic comparisons, with local tree artifacts, not definitive clade,
transmission, virulence, origin or clinical predictions. Sparse public sampling
limits interpretation. Errors are collected per record; a partly failed batch
returns nonzero and retains its report and logs.

Public validation is reproducible with
`PYTHONPATH=. bash scripts/runtime.sh python results/placement-validation/validate.py`.
For each segment it removes every matching explicit isolate label, independently
rebuilds training and full trees, reconstructs training ancestral sequences, and
runs real Nextclade placement. It compares broad nearest relationships and reports
ties rather than claiming exact-node accuracy. See the recorded segment metrics.

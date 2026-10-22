# Public-only placement validation

All FASTA and TSV files here derive from public GenBank records already in this
repository. `smoke.fasta` is the public S reference NC_018464.1 relabelled `case`;
it is not a patient or clinical sample. Its real Nextclade run and automatic full
rebuild completed in ignored `case-results/placement-smoke/`.

Run the script with `PYTHONPATH=. bash scripts/runtime.sh python
results/placement-validation/validate.py`. Each segment uses one explicit public
isolate and excludes all of that isolate's labelled accession aliases. Training
trees are inferred from retained sequences. Existing inference is reused only when
the alignment hash and complete inference settings match its signature. Full trees are inferred independently
in a separate prefix, with seed 159150, 1,000 SH-aLRT and 1,000 ultrafast bootstrap
replicates. Ancestral inference and Nextclade use the rebuilt training tree.

`metrics.json` records exact callable comparisons, Nextclade attachment and the
relationship to the independently inferred full tree. Overlap of nearest sets is
a descriptive broad relationship check, not a validated exact-node accuracy rate.
Ties and poor resolution remain explicit. The script's query rebuilds are full
protein realignment, codon backtranslation and nucleotide inference.

Observed held-out checks (one isolate per segment):

| Segment | Held out | Training/full tips | Closest callable-distance ties | Attachment-nearest agrees with full-tree nearest |
|---|---|---|---|---|
| S | PZ860854.1 | 43/44 | 7 | No |
| M | PZ860853.1 | 36/37 | 1 | Yes |
| L | PZ860852.1 | 30/31 | 4 | No |

Every excluded explicit-isolate alias is absent from the rebuilt training tree.
For S and L, both the attachment-nearest and full-tree nearest are members of the
same callable-distance tie set. This supports a broad tied relationship, not an
exact placement claim. Both queries triggered and completed full coding rebuilds
and retained null definitive assignments. M has a unique nearest public sequence,
OZ502472.1, with agreement in these two trees; one example is not an accuracy
estimate. None of these experiments measures clinical performance.

The 2026-10-22 final-review rerun regenerated ancestral/export data and repeated real
Nextclade acceptance against unchanged signature-verified inference. Interactive
count semantics follow [the placement methods](../../docs/placement.md). Every
training export passed schema, mutation replay and cumulative-event checks.

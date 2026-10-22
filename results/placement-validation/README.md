# Public-only placement validation

All FASTA and TSV files here derive from public GenBank records already in this
repository. `smoke.fasta` is the public S reference NC_018464.1 relabelled `case`;
it is not a patient or clinical sample. Its real Nextclade run and automatic full
rebuild completed in ignored `case-results/placement-smoke/`.

Run the script with `PYTHONPATH=. bash scripts/runtime.sh python
results/placement-validation/validate.py`. Each segment uses one explicit public
isolate and excludes all of that isolate's labelled accession aliases. Training
trees are rebuilt from retained sequences. Full trees are inferred independently
in a separate prefix, with seed 159150, 1,000 SH-aLRT and 1,000 ultrafast bootstrap
replicates. Ancestral inference and Nextclade use the rebuilt training tree.

`metrics.json` records exact callable comparisons, Nextclade attachment and the
relationship to the independently inferred full tree. Overlap of nearest sets is
a descriptive broad relationship check, not a validated exact-node accuracy rate.
Ties and poor resolution remain explicit. The script's query rebuilds are full
protein realignment, codon backtranslation and nucleotide inference.

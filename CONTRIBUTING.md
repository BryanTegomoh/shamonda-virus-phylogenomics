# Contributing

For a sequence or metadata correction, identify the accession.version and segment,
link the deposited record or primary publication, and explain the discrepancy.
Preserve original source values and record the reason for a normalized value.
Do not infer isolate identity from accession order or geography alone.

For a bug report, include the command, tool versions and a small public or synthetic
input that reproduces the failure. Do not share clinical sequences or identifying
metadata in issues, pull requests or commits. The public fixtures in
`results/placement-validation/` are deposited reference sequences, not patient data.

Before submitting changes, run `make validate`. Changes affecting curation or
phylogenetic inference also require a rebuilt public analysis and updated findings.
Keep source retrieval dates distinct from rebuild dates. Describe what changed
and report the checks actually executed.

Use short, imperative commit subjects in sentence case, for example
`Preserve collection-date precision in metadata`. Keep each commit focused on one
change. Use the body to explain the reason when the subject is insufficient.

For scientific interpretation, distinguish source observations from hypotheses.
Changes to cluster, reassortment or temporal conclusions need the relevant support
values and sensitivity results. A nearest sequence is not evidence of a geographic
source or a direct transmission event.

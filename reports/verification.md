# Baseline verification, 2026-10-22

The following checks were run on implementation `0b187f4`. Tool versions are
recorded in [runtime.json](runtime.json).

Executed acceptance checks:

- `make validate`: 70 tests passed, strict mypy passed seven modules, all three
  actual placement-tree schemas passed, with reference agreement, mutation replay,
  inferred-substitution counts, source/output checksums and private-file exclusion.
- Strict mypy also passed the installer and public held-out validation driver.
- The final `make build CORES=4` completed all ten workflow steps using valid
  inference signatures; `make report` and integrity validation passed afterward.
- Actual Auspice 2.66.0 rendered L, M and S. The full views and a SHAV-only S zoom
  displayed the documented inferred-substitution scale. Static figures were
  inspected for labels, support annotations and substitutions/site axes.
- The official pinned Nextclade 3.25.1 download passed size and executable-version
  checks. An idempotent installation check passed separately.

The [scientific validation receipt](final-review-validation.json) records invariant
split supports under both roots and unchanged inference trees. The
[held-out experiment](../results/placement-validation/README.md) records the actual
placement outcomes and their limited scope. These checks do not establish clinical
performance, transmission chains, phenotype prediction or exact placement accuracy.

## Public interactive views

The README's Nextstrain links use the documented public-URL sharing service to
load the committed GitHub JSON files. All three rendered in the public viewer:
21 of 31 tips for L, 27 of 37 for M and 34 of 44 for S with the SHAV filter active.
The views were colored by country and zoomed to the selected tips.

The dataset descriptions now link to the public repository. This metadata-only
change preserves the complete tree objects and reference sequences; the matching
placement datasets and recorded export checksums were updated together.
`make validate` passed after this change: 70 tests, strict typing, dataset schemas
and snapshot integrity checks.

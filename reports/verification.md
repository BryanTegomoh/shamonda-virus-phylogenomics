# Baseline verification, 2026-10-22

The implementation through `0b187f4` passed independent specification and code
review after correcting rerooted support assignment, interactive divergence units,
specimen metadata, placement-tree integrity and refresh-dependent conclusions.
No Critical or Important findings remained in the scoped final review.

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

# Public interactive datasets

`shamonda_L.json`, `shamonda_M.json` and `shamonda_S.json` are the canonical display
exports. Their existing public-URL links remain valid.

Nextstrain Community requires dataset filenames to begin with the repository
name. The `community` rule in `phylogenetic/Snakefile` creates byte-identical
`shamonda-virus-phylogenomics_*.json` publication copies. `make build` refreshes
these copies; `make validate` fails if they differ from the canonical exports.
Do not edit the publication copies separately.

The guided interpretation is maintained in `narratives/`. Data selection, tree
inference and findings are documented in the repository's main README and report.

# Shamonda virus phylogenomics

Sequence curation and segment-specific phylogenies for Shamonda virus (SHAV),
with a reproducible workflow for comparing new consensus sequences against
available public data.

[Findings](reports/findings.md) · [Sequence inventory](data/public/metadata.tsv) ·
[Methods and limitations](reports/findings.md#methods-and-reproducibility) ·
[New-sequence placement](docs/placement.md)

## Findings at a glance

**Baseline: 22 October 2026.** The inventory reconciles 97 SHAV nucleotide records
into 94 canonical segment records. These are accession counts, not 94 independent
genomes. Another 30 reference segments provide related-virus context. The raw-read
catalog contains 20 runs; reads were not downloaded or assembled.

- **European sequences occupy two genetic backgrounds.** They do not form one
  Europe-exclusive cluster in any segment. The pattern is consistent with separate
  introductions, but does not identify the source country, entry route or arrival date.
- **Support differs across segments.** The earlier European group is strongly
  supported in L and M; EU2 is strongly supported in L and S. EU2's M branch and
  the earlier group's S branch fall below the prespecified joint support thresholds.
- **The human case cannot be assigned to either group.** No explicitly
  human-labelled sequence was identified in the snapshot. Fifty focal records lack
  host metadata, so this does not establish global absence of a deposited human sequence.

The [full report](reports/findings.md) gives branch-support values, segment
comparisons and the evidence gaps. Clinical associations and countermeasures are
covered in a separate [source-linked evidence table](docs/evidence.tsv). Observed
sequence differences are not validated virulence or transmissibility markers.

## Phylogenetic trees

Select a tree for the full-resolution figure. These views retain SHAV tips from
the context-inclusive maximum-likelihood analysis. Branch lengths are substitutions
per nucleotide site. Dots mark conservative support annotations meeting both
SH-aLRT ≥80 and ultrafast bootstrap ≥95; their provenance is documented in the report.

| L: polymerase | M: polyprotein | S: nucleocapsid N |
|:---:|:---:|:---:|
| [![L segment phylogeny](figures/tree_L_SHAV.png)](figures/tree_L_SHAV.png) | [![M segment phylogeny](figures/tree_M_SHAV.png)](figures/tree_M_SHAV.png) | [![S segment phylogeny](figures/tree_S_SHAV.png)](figures/tree_S_SHAV.png) |
| [Full context PDF](figures/tree_L.pdf) | [Full context PDF](figures/tree_M.pdf) | [Full context PDF](figures/tree_S.pdf) |

S describes the N coding region, not the entire overlapping N/NSs sequence.
Segments are analyzed independently. No concatenated tree is used.

Interactive Auspice datasets are included in [auspice/](auspice/). After installing
the tools below, run `make view` and select L, M or S in the local viewer.
Its mutation axis counts inferred A/C/G/T substitutions; the static trees above
retain model-based substitutions/site. [Coordinate and scale conventions](docs/placement.md)
explain the difference.

## Data and provenance

| File | Contents |
|---|---|
| [metadata.tsv](data/public/metadata.tsv) | Every retrieved accession, original and normalized metadata, inclusion decisions and exclusions |
| [sequences.fasta](data/public/sequences.fasta) | Canonical SHAV and context sequences |
| [records.gb](data/public/records.gb) | Deposited GenBank records with annotations |
| [manifest.json](data/public/manifest.json) | Source inventories, exact queries, retrieval dates and checksums |
| [read_catalog.tsv](data/public/read_catalog.tsv) | Linked SRA/INSDC experiments and runs |
| [Sequence acknowledgments](reports/sequence_acknowledgments.tsv) | Accession-level authors and publications |

The [ingestion notes](ingest/README.md) describe reconciliation of NCBI
Nucleotide, Datasets and publication-linked accessions. Unknown metadata remain
blank. The catalog uses NCBI SRA and shared INSDC accessions; it is not an
independent reconciliation of ENA-specific metadata.

## Install and run

Prerequisites: conda or mamba, uv, curl and Git. The baseline was tested on Apple
Silicon with Python 3.12 and the versions in [environment.yml](environment.yml).
[Recorded runtime versions](reports/runtime.json) identify the tools actually used.

```sh
git clone https://github.com/BryanTegomoh/shamonda-virus-phylogenomics.git
cd shamonda-virus-phylogenomics
conda env create -f environment.yml
conda activate shamonda-phylogenomics-tools
make setup
bash scripts/runtime.sh python scripts/install-nextclade.py
```

### Explore the included results

```sh
make ingest  # Verify the pinned public-data snapshot
make view    # Serve the included trees locally on port 4000
```

The committed figures and reports can be read without running an analysis.
The viewer uses local datasets; no public hosting or upload is configured.

### Reproduce the analysis

```sh
make build CORES=4
make report
make validate
```

The workflow uses protein-guided coding alignments and IQ-TREE model selection
with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates, with fixed seeds.
Complete-case alignments and documented partial-sequence windows provide
sensitivity analyses. [Methods](reports/findings.md#methods-and-reproducibility)
describe the inclusion rules and rooting choices.

To retrieve a new public-data snapshot, run `make ingest REFRESH=1`, then rebuild,
regenerate the report and validate. This replaces the baseline snapshot. Preserve
the previous results in Git before refreshing. Original responses are cached
locally; [cache and replay options](ingest/README.md) are documented separately.

### Compare a new sequence

```sh
make place SEQUENCES=case-inputs/case.fasta METADATA=case-inputs/case.tsv
```

The metadata TSV requires `sample` and `segment`. See the [input specification](docs/placement.md)
for FASTA identifiers, optional metadata and QC thresholds. Outputs include a case
report, callable-site distances and an interactive tree. Ambiguous or divergent
queries trigger a full tree rebuild and retain unresolved alternatives.

Keep clinical inputs in ignored `case-inputs/`; reports default to ignored
`case-results/`. Do not put case sequences or identifying metadata in GitHub issues
or pull requests. The local datasets leave Nextclade QC modules unconfigured;
separate coverage, ambiguity and overlap checks control eligibility.

## Validation and limits

[Verification](reports/verification.md) records the executed tests and integrity
checks. [Held-out placement validation](results/placement-validation/README.md)
uses one public isolate per segment: M agreed on the nearest reference, while L
and S remained tied and required rebuilds. This is not an accuracy estimate.

Public sampling is sparse and uneven. Clustering does not establish a transmission
chain, and these analyses do not resolve introduction dates or calibrated
evolutionary rates. The workflow is for research and surveillance comparisons;
it has not been clinically validated.

## Contributing and citation

Accession corrections, metadata improvements and reproducible bug reports are
welcome. Read [CONTRIBUTING.md](CONTRIBUTING.md) for the required checks and commit
conventions. [CITATION.cff](CITATION.cff) provides the repository citation.
Credit the original sequence authors and retain accession versions when reusing
the data.

The repository layout follows [Nextstrain's pathogen repository guide](https://github.com/nextstrain/pathogen-repo-guide).

## License

Original code, project-authored documentation and figures are available under the
[MIT License](LICENSE). Deposited sequence records, source metadata and other
third-party material retain their source terms and attribution requirements.

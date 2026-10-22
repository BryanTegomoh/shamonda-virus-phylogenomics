# Shamonda virus phylogenomics

Private research workflow for public SHAV sequence curation, independent L/M/S
phylogenies and local placement of new consensus sequences.

Start with the [measured findings](reports/findings.md), [L](figures/tree_L_SHAV.png),
[M](figures/tree_M_SHAV.png), and [S](figures/tree_S_SHAV.png) trees. The report
answers the European-clustering and human-sequence questions and documents
limitations. [Clinical and countermeasure evidence](docs/evidence.tsv) is curated
separately from observed sequence variation.

## Run

Prerequisites: conda or mamba, uv, curl and Git.

Tested on Apple Silicon with Python 3.12.12 for helpers and the scientific tools
pinned in [environment.yml](environment.yml). Create that conda environment and
activate it, then install the Python lock and Nextclade:

```sh
conda env create -f environment.yml
conda activate shamonda-phylogenomics-tools
make setup
bash scripts/runtime.sh python scripts/install-nextclade.py
make ingest
make build CORES=6
make report
make validate
make view
```

`make ingest` verifies the delivered snapshot without changing it. To obtain a new
snapshot use `make ingest REFRESH=1`. Exact queries, annotations, accession versions,
source hashes, exclusions and retrieval times live in the
[manifest](data/public/manifest.json). Raw responses are cached locally. The default
ingest CLI also supports `--cache-dir` for storage outside cloud-synced folders.
Reads are cataloged, not downloaded or assembled.

The runtime script uses a project `.venv` or an existing external cache venv and
scientific executables on PATH. On the original Mac it can reuse the existing
hantavirus tool environment. No packages are added to that environment. Runtime
versions used for delivery are recorded in [runtime.json](reports/runtime.json).

## New case

```sh
make place SEQUENCES=/path/to/case.fasta METADATA=/path/to/case.tsv
```

See [input format, QC and interpretation](docs/placement.md). This accepts deposited
consensus sequences, coding sequences or sufficiently complete segments on either
strand. New-case outputs default to ignored `case-results/`. Exact placement is
not established for ambiguous sequences; the workflow reports ties and rebuilds
the tree. [Held-out validation](results/placement-validation/summary.json) records
actual per-segment outcomes rather than a general accuracy claim.

`make view` opens a local Auspice server on port 4000. Datasets remain private.
No Nextstrain Community upload or public hosting is configured.

## Layout and provenance

- `ingest/`, `data/public/`: source inventories, annotated records and metadata.
- `phylogenetic/`, `results/`, `auspice/`: coding alignments, supported ML trees and interactive exports.
- `nextclade/`: local primary-CDS placement datasets.
- `reports/`, `figures/`: measured findings and reusable outputs.

Organization follows [Nextstrain's pathogen repository guide](https://github.com/nextstrain/pathogen-repo-guide)
and the author's [hantavirus workflow](https://github.com/BryanTegomoh/hantavirus-nextstrain).
The implementation is specific to SHAV; unrelated template automation is not included.
Original sequence authors and submitters are acknowledged in the
[accession-level source table](reports/sequence_acknowledgments.tsv) and GenBank records.

The original requested path `/Users/bryan/Documents/ChatGPT/Shamonda virus` is a
symlink to `/Users/bryan/Projects/shamonda-virus-phylogenomics`. The working files
were moved outside iCloud after sync conflicts replaced active source files.

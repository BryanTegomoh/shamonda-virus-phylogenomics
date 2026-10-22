# Shamonda virus phylogenomics implementation

## Global constraints

Private repository `BryanTegomoh/shamonda-virus-phylogenomics`; canonical checkout
`/Users/bryan/Documents/ChatGPT/Shamonda virus`. Public consensus analysis and raw-read
catalog only. Independent L/M/S trees, no concatenation. Strict typed Python, real
failure-path tests, deterministic analyses, source-backed claims. No unsupported
virulence/transmission scores. Private case inputs and results never enter Git.
The user explicitly authorized implementation, validation, commit, and private push.

## Task 1: public data ingestion

Own `shav/ingest.py`, `tests/test_ingest.py`, `ingest/`, `config/ingest.json`, and
`data/public/`. Implement the union of NCBI Datasets taxon 159150, Nucleotide taxon,
name-search and publication accessions, plus a curated Simbu reference panel. Fetch
all public SHAV nucleotide records including partial sequences. Reconcile duplicated
GenBank/RefSeq records, preserve provenance and annotation. Catalog SRA/ENA links
without reads. Track exact retrieval queries, versions, hashes, original/normalized
metadata and inclusion decisions. Include linked publications, biosamples, bioprojects.
Default output `data/public/metadata.tsv` (one row per accession), `sequences.fasta`,
`records.gb` with complete annotated GenBank records, `manifest.json`, and
`read_catalog.tsv`. Metadata must have accession, segment L/M/S, isolate, country,
host, date (ISO with original precision), date_raw, date_precision, group (SHAV or
reference virus name), taxid, canonical_accession, inclusion, reason, biosample,
bioproject, sra, passage, publication. Preserve blank unknowns. No guessed isolate
links; lab passage is not patient/animal evidence. `inclusion=include` marks canonical
records; all others have explicit reasons. Sequence IDs are accession.version.
The wider reference panel is annotated near-complete references for Schmallenberg,
Sathuperi, Douglas, Peaton, Sango, Shuni, Aino, Akabane, Simbu, plus Bunyamwera as an
outgroup. Prefer curated RefSeq; resolve names/taxonomy from live NCBI and retain
source identity labels separately from modern species names. Compare taxon vs
keyword/Datasets inventories. Check explicitly labelled human records, preserving
unknown hosts as unknown. Use pinned retrieval snapshot with refresh opt-in. Rate
limit, retry only transient failures, atomic cache writes, verify complete retrieval,
fail if expected accessions are missing. Tests must cover pagination, failures,
deduplication, dates, segment conflicts, orientation and malformed metadata.
Can use existing `shav/common.py` but do not change it; coordinate needed changes.
Provide `python -m shav.ingest [--refresh]` and save results of a full successful run.

## Task 2: clinical and epidemiological evidence

Own `docs/evidence.tsv` and `docs/evidence-notes.md`. Independently verify primary
sources and curate one evidence table covering clinical associations, transmission,
SHAV gene functions, human infection announced 2026-10-22, countermeasure availability,
cross-protection, European variant labels and sequence accession mapping. Columns:
topic, statement, evidence_level, applicability, source_title, source_url,
publication_date, checked_date, limitations. Separate SHAV direct evidence, related
virus evidence and computational hypotheses. No unsupported phenotype prediction.
Verify source URLs and connected author/journal/DOI triples for papers. Explain
unresolved European nomenclature and any primary source inconsistencies. No email,
external messages, handbook edits or biological experimentation.

## Task 3: phylogeny and reports

Controller owns `shav/phylogeny.py`, `shav/report.py`, `shav/validate.py`, tests for
these, workflow wiring, runtime, README, reports and figures. Extract annotated
primary CDS (L polymerase, M polyprotein, S N), translate and protein-align with MAFFT,
backtranslate. Include near-complete >=90% homologous CDS and partial analyses over
explicit common regions. ModelFinder + 1000 UFBoot + 1000 SH-aLRT, fixed seed. Compare
alignment masks and roots, supported clusters and pairwise callable-site distances.
Reassortment candidates only for matched isolates with supported incompatible
splits; uncertainty is explicit. Temporal signal and randomization precede time
estimation, and sparse metadata can preclude dates/rates. Export validated Auspice
trees and static publication figures. Run full public analysis and report exact
counts/exclusions and limitations.

## Task 4: placement and acceptance

Controller or a fresh implementer owns `shav/place.py` and placement tests, local
Nextclade segment datasets and held-out placement validation. Interface:
`make place SEQUENCES=case.fasta METADATA=case.tsv` (sample, segment required; date,
host, country optional). QC, reverse orientation, callable-site distances, missing
segments, ambiguity and divergent queries explicit. Reference and tree roots must
agree. Validate local Nextclade placements against full tree rebuilds on held-out
references. Unresolved cases require a rebuild and must not produce definitive
assignment. Accept one, two or three segments, never invent missing ones. Outputs
under ignored `case-results/`. Final checks: tests, mypy strict, full workflows,
Snakemake DAG, Auspice schema, actual held-out runs, privacy, commits and private push.

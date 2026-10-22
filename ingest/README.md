# Public retrieval snapshot

Run `python -m shav.ingest` to verify the pinned output hashes. Explicit
`--refresh` retrieves a new snapshot and replaces the manifest. Missing outputs,
changed configuration, malformed records and incomplete requests fail loudly.

`data/public/records.gb` retains every retrieved nucleotide record and its full
annotation. `metadata.tsv` contains one row per accession.version, including
excluded keyword-context records. `sequences.fasta` contains only canonical,
explicitly assigned focal and reference segments. Unknown metadata stays blank.
No records are linked into a putative isolate from their accession order.

Sources are the union of NCBI Nucleotide focal taxonomy, broad name search,
publication accessions in configuration, and NCBI Datasets. The manifest records
exact queries, source inventories, source differences, retrieval time, Datasets
version and SHA-256 digests. Raw E-utilities responses are cached atomically under
`.cache/ingest/` by URL SHA-256. `datasets-summary.jsonl` preserves the separate
Datasets source response. Cached responses can be removed only after preserving
the snapshot if an exact replay is needed.

The panel uses three annotated reference segments for each requested historical
virus label. RefSeq is preferred. Douglas virus has no RefSeq returned under its
source label, so the panel uses HE795090.1, HE795091.1 and HE795092.1 (isolate 93-6).
Source labels and current annotated organism names remain separate. Selecting a
reference for each segment is not evidence that segment genealogies agree.

Exact duplicate nucleotide sequences (including reverse complements) collapse
only when explicit isolate/strain, segment and virus group annotations match.
Unlabelled records remain separate. As-deposited orientation is retained; strand
normalization belongs to the downstream annotated-CDS alignment workflow.

SRA/INSDC runs are catalogued without downloading reads. A nucleotide record gets
run links only through matching explicit BioSample identifiers. SRA run accessions
are shared across INSDC archives including ENA; this snapshot queries NCBI SRA,
not a separately reconciled ENA-specific metadata inventory.

The 2026-10-22 retrieval contains no explicitly human-labelled nucleotide source
host. This does not establish absence of human infection or absence of deposited
human sequences with missing metadata. The IHUMI_Shav000001 records
PZ860852.1/PZ860853.1/PZ860854.1 explicitly annotate Bos taurus. Cell-culture passage
metadata is stored separately and is never treated as a host observation.

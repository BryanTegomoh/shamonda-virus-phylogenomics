# Public retrieval snapshot

Run `python -m shav.ingest` to verify the pinned output hashes. Explicit
`--refresh` retrieves a new snapshot and replaces the manifest. Use `--rebuild`
to regenerate outputs after a code change from the preserved source cache,
fetching only missing requests. `--cache-dir PATH` selects an external cache when
cloud storage has evicted local responses. Such a cache fetches fresh bytes for
missing requests rather than reconstructing responses. Cache entries record
retrieval and access times. Missing outputs,
changed configuration, malformed records and incomplete requests fail loudly.

`data/public/records.gb` retains every retrieved nucleotide record and its full
annotation. `metadata.tsv` contains one row per accession.version, including
excluded keyword-context records. `sequences.fasta` contains only canonical,
annotated focal and reference segments. Missing source segment qualifiers can be resolved only from explicit orthobunyavirus nucleocapsid/NSs CDS products, with assignment provenance retained. Unknown metadata stays blank.
No records are linked into a putative isolate from their accession order.

Sources are the union of NCBI Nucleotide focal taxonomy, broad name search,
publication accessions in configuration, and NCBI Datasets. The manifest records
exact UID-to-accession mappings, ingestion source hash, queries, source inventories, source differences, retrieval time, Datasets
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

SRA/INSDC runs are catalogued without downloading reads. Discovery searches
the focal name/taxon plus explicit BioSample, BioProject and focal nucleotide
identifiers. Requested SRA UIDs are resolved to experiments and samples before
package retrieval; missing or substituted packages fail. Candidates found only
by identifiers must match a verified focal BioSample or BioProject; unmatched
hits are retained with exclusion reasons in the manifest. A nucleotide record gets
run links only through matching explicit BioSample identifiers. SRA run accessions
are shared across INSDC archives including ENA; this snapshot queries NCBI SRA,
not a separately reconciled ENA-specific metadata inventory.

The 2026-10-22 retrieval contains no explicitly human-labelled nucleotide source
host. Nuccore EInfo field discovery and taxon/name searches combined with human
keywords are preserved as timestamped sensitivity checks. Host classification
uses only explicit GenBank source qualifiers. This does not establish absence of human infection or absence of deposited
human sequences with missing metadata. The IHUMI_Shav000001 records
PZ860852.1/PZ860853.1/PZ860854.1 explicitly annotate Bos taurus. Cell-culture passage
metadata is stored separately and is never treated as a host observation.

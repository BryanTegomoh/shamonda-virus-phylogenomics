"""Reproducible public nucleotide inventory. Never infer hosts or isolate links."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from shav.common import ROOT, sha256, write_fasta, write_json, write_tsv

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
FIELDS = "accession segment isolate country host date date_raw date_precision group taxid canonical_accession inclusion reason biosample bioproject sra passage publication length organism orientation provenance geo_loc_raw location assignment_provenance".split()


class Client:
    def __init__(self, cache: Path, refresh: bool = False) -> None:
        self.cache = cache
        self.refresh = refresh
        self.requests: list[dict[str, str]] = []
        self.last = 0.0

    def get(self, endpoint: str, **params: str) -> bytes:
        url = BASE + endpoint + "?" + urllib.parse.urlencode(params)
        key = hashlib.sha256(url.encode()).hexdigest()
        path = self.cache / (key + ".raw")
        self.cache.mkdir(parents=True, exist_ok=True)
        if self.refresh or not path.exists():
            for attempt in range(5):
                time.sleep(max(0, 0.36 - (time.monotonic() - self.last)))
                try:
                    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "SHAV-public-inventory/0.1"}), timeout=30) as response:
                        body = response.read()
                    self.last = time.monotonic()
                    if not body:
                        raise ValueError(f"Empty response: {url}")
                    temp = path.with_suffix(".tmp")
                    temp.write_bytes(body)
                    temp.replace(path)
                    break
                except urllib.error.HTTPError as exc:
                    self.last = time.monotonic()
                    if exc.code != 429 and exc.code < 500:
                        raise
                    if attempt == 4:
                        raise
                    time.sleep(2 ** (attempt + 1))
                except (urllib.error.URLError, TimeoutError):
                    self.last = time.monotonic()
                    if attempt == 4:
                        raise
                    time.sleep(2 ** (attempt + 1))
        body = path.read_bytes()
        self.requests.append({"url": url, "cache": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(body).hexdigest()})
        return body

    def search(self, db: str, term: str, page_size: int = 100) -> list[str]:
        ids: list[str] = []
        total: int | None = None
        while total is None or len(ids) < total:
            result = json.loads(self.get("esearch.fcgi", db=db, term=term, retmode="json", retstart=str(len(ids)), retmax=str(page_size)))
            data = result["esearchresult"]
            if "errorlist" in data or "ERROR" in result:
                raise ValueError(f"Search failed: {result}")
            count = int(data["count"])
            if total is not None and count != total:
                raise ValueError("Inventory changed during pagination; refresh again")
            total = count
            page = data["idlist"]
            if not isinstance(page, list) or (not page and len(ids) < total):
                raise ValueError("Truncated search")
            ids.extend(page)
        if len(ids) != total or len(set(ids)) != total:
            raise ValueError("Duplicate or incomplete search pagination")
        return ids

    def fetch(self, ids: list[str]) -> list[SeqRecord]:
        records: list[SeqRecord] = []
        for start in range(0, len(ids), 40):
            batch = ids[start:start + 40]
            body = self.get("efetch.fcgi", db="nuccore", id=",".join(batch), rettype="gbwithparts", retmode="text")
            parsed: list[SeqRecord] = list(SeqIO.parse(io.StringIO(body.decode()), "genbank"))
            if len(parsed) != len(batch):
                raise ValueError(f"Incomplete GenBank retrieval: expected {len(batch)}, got {len(parsed)}")
            records.extend(parsed)
        return records


def normalize_date(raw: str) -> tuple[str, str]:
    if not raw or raw.lower() in {"unknown", "not collected", "not provided", "missing"}:
        return "", ""
    if re.fullmatch(r"\d{4}", raw):
        return raw, "year"
    for fmt, precision, output in [("%b-%Y", "month", "%Y-%m"), ("%d-%b-%Y", "day", "%Y-%m-%d"), ("%Y-%m-%d", "day", "%Y-%m-%d"), ("%Y-%m", "month", "%Y-%m")]:
        try:
            return datetime.strptime(raw, fmt).strftime(output), precision
        except ValueError:
            pass
    return "", "unparsed"


def segment(record: SeqRecord) -> str:
    candidates: set[str] = set()
    for feature in record.features:
        if feature.type == "source":
            for value in feature.qualifiers.get("segment", []):
                match = re.fullmatch(r"(?:segment\s*)?([LMS])(?:\s*(?:RNA|segment))?", value, re.I)
                if match:
                    candidates.add(match.group(1).upper())
    for value in re.findall(r"\b([LMS])\s+(?:RNA\s+)?segment\b|\bsegment\s+([LMS])\b", record.description, re.I):
        candidates.update(v.upper() for v in value if v)
    # Orthobunyavirus N/nucleocapsid and NSs are encoded on S. Only explicit CDS
    # products supply this fallback, never accession order or isolate name.
    is_orthobunyavirus = "Shamonda" in record.description or "Orthobunyavirus" in cast(list[str], record.annotations.get("taxonomy", []))
    if is_orthobunyavirus:
        for feature in record.features:
            if feature.type == "CDS" and any(re.search(r"nucleocapsid|nucleoprotein|^NSs$", product, re.I) for product in feature.qualifiers.get("product", [])):
                candidates.add("S")
            if feature.type == "CDS" and any(re.search(r"RNA-dependent RNA polymerase", product, re.I) for product in feature.qualifiers.get("product", [])):
                candidates.add("L")
            if feature.type == "CDS" and any(re.search(r"glycoprotein", product, re.I) for product in feature.qualifiers.get("product", [])):
                candidates.add("M")
    if len(candidates) > 1:
        raise ValueError(f"Conflicting segment annotations: {record.id}: {candidates}")
    return next(iter(candidates), "")


def metadata(record: SeqRecord, group: str, provenance: list[str]) -> dict[str, str]:
    if record.id is None or record.seq is None:
        raise ValueError("Missing record identifier or sequence")
    if not re.fullmatch(r"[A-Z][A-Z_0-9]*\.\d+", record.id) or not len(record.seq):
        raise ValueError(f"Malformed accession or empty sequence: {record.id}")
    sources = [f for f in record.features if f.type == "source"]
    if len(sources) != 1:
        raise ValueError(f"Expected one source feature: {record.id}")
    qualifiers = sources[0].qualifiers
    def q(name: str) -> str:
        return "; ".join(qualifiers.get(name, []))
    raw_date = q("collection_date")
    date, precision = normalize_date(raw_date)
    links = list(record.dbxrefs) + qualifiers.get("db_xref", [])
    def link(prefix: str) -> str:
        return ";".join(x.split(":", 1)[1] for x in links if x.lower().startswith(prefix.lower() + ":"))
    refs = cast(list[Any], record.annotations.get("references", []))  # Biopython heterogeneous annotations boundary.
    publications = ["PMID:" + ref.pubmed_id if getattr(ref, "pubmed_id", "") else ref.title for ref in refs]
    row = dict.fromkeys(FIELDS, "")
    row.update(accession=record.id, segment=segment(record), isolate=q("isolate") or q("strain"), country=q("geo_loc_name") or q("country"), host=q("host"), date=date, date_raw=raw_date, date_precision=precision, group=group, taxid=link("taxon"), canonical_accession=record.id, inclusion="include", reason="canonical public nucleotide record", biosample=link("BioSample"), bioproject=link("BioProject") or link("Project"), sra=link("SRA"), passage=q("lab_host") or q("passage_history"), publication="; ".join(publications), length=str(len(record.seq)), organism=str(record.annotations.get("organism", "")), orientation="as deposited", provenance=";".join(sorted(provenance)))
    geo = row["country"]
    row["geo_loc_raw"] = geo
    row["country"] = geo.partition(":")[0].strip()
    row["location"] = geo.partition(":")[2].strip()
    has_explicit_segment = any(feature.qualifiers.get("segment") for feature in sources) or bool(re.search(r"\b[LMS]\s+(?:RNA\s+)?segment\b|\bsegment\s+[LMS]\b", record.description, re.I))
    row["assignment_provenance"] = "explicit source/title segment annotation" if has_explicit_segment else ("annotated orthobunyavirus segment-specific CDS product" if row["segment"] else "")
    return row


def reconcile(records: list[SeqRecord], rows: list[dict[str, str]]) -> None:
    """Collapse exact duplicate sequences only within explicit same-isolate annotations.

    RefSeq source accession links take precedence. Reverse complements collapse only
    when annotation establishes the same isolate, segment, and virus group.
    """
    by_id = {record.id: record for record in records}
    ordered = sorted(rows, key=lambda r: (not r["accession"].startswith("NC_"), r["accession"]))
    canonical: list[dict[str, str]] = []
    for row in ordered:
        current_seq = by_id[row["accession"]].seq
        if current_seq is None:
            raise ValueError("Missing sequence during reconciliation")
        sequence = str(current_seq).upper()
        match: dict[str, str] | None = None
        for previous in canonical:
            if not row["isolate"] or any(row[k] != previous[k] for k in ("isolate", "segment", "group")):
                continue
            prev_sequence = str(by_id[previous["accession"]].seq).upper()
            if sequence == prev_sequence or str(current_seq.reverse_complement()).upper() == prev_sequence:  # type: ignore[no-untyped-call]  # Biopython boundary.
                match = previous
                break
        if match:
            row.update(canonical_accession=match["accession"], inclusion="exclude", reason="exact sequence duplicate with matching explicit isolate, segment and group")
        else:
            canonical.append(row)


def run(refresh: bool = False) -> None:
    config_path = ROOT / "config/ingest.json"
    config = json.loads(config_path.read_text())
    output = ROOT / "data/public"
    manifest_path = output / "manifest.json"
    if manifest_path.exists() and not refresh:
        manifest = json.loads(manifest_path.read_text())
        for name, digest in manifest["outputs"].items():
            if sha256(output / name) != digest:
                raise ValueError(f"Pinned snapshot changed: {name}; restore or use --refresh")
        if manifest["config_sha256"] != sha256(config_path):
            raise ValueError("Config changed; explicit --refresh required")
        print("Verified pinned public snapshot")
        return
    client = Client(ROOT / ".cache/ingest", refresh)
    inventories: dict[str, list[str]] = {}
    queries = {"taxon": "txid159150[Organism:exp]", "keyword": '"Shamonda"[All Fields] OR "SHAV"[Title]', "publication": " OR ".join(f'{accession}[Accession]' for accession in config["publication_accessions"])}
    provenance: dict[str, list[str]] = {}
    for label, query in queries.items():
        inventories[label] = client.search("nuccore", query)
        for uid in inventories[label]:
            provenance.setdefault(uid, []).append(label)
    datasets_path = ROOT / "ingest/datasets-summary.jsonl"
    datasets_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([config["datasets_executable"], "summary", "virus", "genome", "taxon", config["taxid"], "--as-json-lines"], check=True, capture_output=True, text=True)
    datasets_path.write_text(result.stdout)
    dataset_accessions = [json.loads(line)["accession"] for line in result.stdout.splitlines() if line.strip()]
    inventories["datasets"] = dataset_accessions
    # Query by accession as Datasets and Entrez expose different inventory units.
    dataset_records = client.fetch(dataset_accessions)
    focal_records = client.fetch(sorted(provenance))
    records_by_id: dict[str, SeqRecord] = {str(r.id): r for r in focal_records + dataset_records}
    origins = {str(r.id): provenance[uid] for uid, r in zip(sorted(provenance), focal_records, strict=True)}
    for r in dataset_records:
        origins.setdefault(str(r.id), []).append("datasets")
    expected = set(config["publication_accessions"]) | set(dataset_accessions)
    if expected - records_by_id.keys():
        raise ValueError(f"Missing expected accessions: {expected - records_by_id.keys()}")
    groups = {str(r.id): "SHAV" if any("taxon:" + config["taxid"] in f.qualifiers.get("db_xref", []) for f in r.features if f.type == "source") else "keyword context" for r in focal_records + dataset_records}
    reference_choices: dict[str, list[str]] = {}
    for name in config["references"]:
        query = f'"{name}"[Organism] AND srcdb_refseq[PROP]'
        ids = client.search("nuccore", query)
        if not ids:
            # Some historically named viruses have no RefSeq under their source name.
            query = f'"{name}"[Title]'
            ids = client.search("nuccore", query)
        candidates = client.fetch(ids)
        chosen: list[SeqRecord] = []
        for seg in ("L", "M", "S"):
            options = [r for r in candidates if segment(r) == seg and (name not in config.get("reference_accessions", {}) or str(r.id) in config["reference_accessions"][name])]
            if not options:
                raise ValueError(f"Missing reference segment {seg}: {name}")
            chosen.append(max(options, key=lambda r: (len(r), str(r.id).startswith("NC_"), str(r.id))))
        reference_choices[name] = [str(r.id) for r in chosen]
        for r in chosen:
            records_by_id[str(r.id)] = r
            groups[str(r.id)] = name
            origins.setdefault(str(r.id), []).append("reference:" + query)
    records = [records_by_id[k] for k in sorted(records_by_id)]
    rows = [metadata(r, groups[str(r.id)], origins[str(r.id)]) for r in records]
    reconcile(records, rows)
    for row in rows:
        if row["group"] == "keyword context":
            row.update(inclusion="exclude", reason="keyword hit outside exact focal taxonomy; retained for audit only")
        elif not row["segment"]:
            row.update(inclusion="exclude", reason="segment not explicitly annotated; retained without inferred assignment")
    sra_ids = client.search("sra", 'txid159150[Organism:exp] OR "Shamonda"[All Fields]')
    read_rows: list[dict[str, str]] = []
    for start in range(0, len(sra_ids), 20):
        xml = client.get("efetch.fcgi", db="sra", id=",".join(sra_ids[start:start + 20]), retmode="xml")
        root = ET.fromstring(xml)
        packages = root.findall(".//EXPERIMENT_PACKAGE")
        if not packages:
            raise ValueError("SRA retrieval returned no experiment packages")
        for package in packages:
            sample = package.find(".//SAMPLE")
            experiment = package.find(".//EXPERIMENT")
            for run_element in package.findall(".//RUN"):
                read_rows.append({"run": run_element.get("accession", ""), "experiment": experiment.get("accession", "") if experiment is not None else "", "sample": sample.get("accession", "") if sample is not None else "", "biosample": ";".join(e.text or "" for e in package.findall(".//EXTERNAL_ID") if e.get("namespace") == "BioSample"), "bioproject": ";".join(e.text or "" for e in package.findall(".//EXTERNAL_ID") if e.get("namespace") == "BioProject"), "source": "NCBI SRA linked INSDC catalog", "reads_downloaded": "false"})
    for row in rows:
        sample_ids = set(row["biosample"].split(";")) - {""}
        matched_runs = sorted({read["run"] for read in read_rows if sample_ids & set(read["biosample"].split(";"))})
        if matched_runs:
            row["sra"] = ";".join(sorted(set(filter(None, row["sra"].split(";"))) | set(matched_runs)))
    write_tsv(output / "metadata.tsv", rows, FIELDS)
    write_fasta(output / "sequences.fasta", [(str(r.id), str(r.seq)) for r in records if next(row for row in rows if row["accession"] == str(r.id))["inclusion"] == "include"])
    SeqIO.write(records, output / "records.gb", "genbank")
    write_tsv(output / "read_catalog.tsv", read_rows, ["run", "experiment", "sample", "biosample", "bioproject", "source", "reads_downloaded"])
    version = subprocess.run([config["datasets_executable"], "version"], check=True, capture_output=True, text=True).stdout.strip()
    write_json(output / "manifest.json", {"retrieved_at": datetime.now(timezone.utc).isoformat(), "config_sha256": sha256(config_path), "queries": queries, "inventories": inventories, "reference_choices": reference_choices, "inventory_comparison": {"accessions_by_source": {label: sorted(accession for accession, labels in origins.items() if label in labels) for label in queries}, "datasets_only": sorted(set(dataset_accessions) - {accession for accession, labels in origins.items() if "taxon" in labels}), "keyword_outside_taxon": sorted(accession for accession, labels in origins.items() if "keyword" in labels and "taxon" not in labels)}, "sra_ids": sra_ids, "requests": client.requests, "datasets_version": version, "datasets_summary_sha256": sha256(datasets_path), "outputs": {name: sha256(output / name) for name in ["metadata.tsv", "sequences.fasta", "records.gb", "read_catalog.tsv"]}, "counts": {"records": len(rows), "included": sum(r["inclusion"] == "include" for r in rows), "explicit_human_records": [r["accession"] for r in rows if re.search(r"human|Homo sapiens", r["host"], re.I)], "unknown_host_records": sum(not r["host"] for r in rows)}, "interpretation": "Host is verbatim source host only. Lab host never establishes patient/animal origin. Different accessions are not linked into isolates without explicit source annotations. Reference source labels are retained separately from current organism names. Dates without recognized precision remain blank and raw values persist."})
    print(f"Saved {len(rows)} nucleotide records, {len(read_rows)} public read catalog runs")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    run(parser.parse_args().refresh)


if __name__ == "__main__":
    main()

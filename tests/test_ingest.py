from __future__ import annotations

import io
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from shav.ingest import Client, metadata, normalize_date, reconcile, segment, discover_reads, parse_sra_packages, select_linked_reads


def record(accession: str = "AB123456.1", sequence: str = "ACGTAA", seg: str = "L") -> SeqRecord:
    result = SeqRecord(Seq(sequence), id=accession, description=f"Shamonda virus {seg} segment")
    result.features = [SeqFeature(SimpleLocation(0, len(sequence)), type="source", qualifiers={"segment": [seg], "isolate": ["named-isolate"], "db_xref": ["taxon:159150"]})]
    return result


@pytest.mark.parametrize(("raw", "expected"), [("2020", ("2020", "year")), ("Jan-2020", ("2020-01", "month")), ("02-Jan-2020", ("2020-01-02", "day")), ("2020-02-30", ("", "unparsed")), ("unknown", ("", "")), ("2020/2021", ("", "unparsed"))])
def test_dates(raw: str, expected: tuple[str, str]) -> None:
    assert normalize_date(raw) == expected


def test_segment_conflict() -> None:
    value = record()
    value.features[0].qualifiers["segment"] = ["S"]
    with pytest.raises(ValueError, match="Conflicting"):
        segment(value)


def test_metadata_preserves_unknown_host() -> None:
    value = record()
    value.features[0].qualifiers["lab_host"] = ["human cell line"]
    row = metadata(value, "SHAV", ["taxon"])
    assert row["host"] == ""
    assert row["lab_host"] == "human cell line"
    assert row["passage"] == ""
    assert row["date"] == ""


def test_malformed_metadata() -> None:
    value = record()
    value.features = []
    with pytest.raises(ValueError, match="one source"):
        metadata(value, "SHAV", [])
    with pytest.raises(ValueError, match="Malformed"):
        metadata(record("not_an_accession"), "SHAV", [])


def test_reverse_complement_deduplication_requires_isolate() -> None:
    values = [record("NC_123456.1"), record("AB123456.1", "TTACGT"), record("AB123457.1", "ACGTAA")]
    values[2].features[0].qualifiers.pop("isolate")
    rows = [metadata(value, "SHAV", []) for value in values]
    reconcile(values, rows)
    assert rows[0]["inclusion"] == "include"
    assert rows[1]["canonical_accession"] == "NC_123456.1"
    assert rows[1]["inclusion"] == "exclude"
    assert rows[2]["inclusion"] == "include"


def test_pagination(tmp_path: Path) -> None:
    client = Client(tmp_path)
    pages = [{"esearchresult": {"count": "3", "idlist": ["1", "2"]}}, {"esearchresult": {"count": "3", "idlist": ["3"]}}]
    with patch.object(client, "get", side_effect=[json.dumps(p).encode() for p in pages]) as get:
        assert client.search("nuccore", "query", 2) == ["1", "2", "3"]
        assert get.call_args_list[1].kwargs["retstart"] == "2"


@pytest.mark.parametrize("pages", [[{"count": "2", "idlist": ["1"]}, {"count": "2", "idlist": []}], [{"count": "2", "idlist": ["1", "1"]}], [{"count": "3", "idlist": ["1"]}, {"count": "2", "idlist": ["2"]}]])
def test_incomplete_pagination_fails(tmp_path: Path, pages: list[dict[str, object]]) -> None:
    client = Client(tmp_path)
    with patch.object(client, "get", side_effect=[json.dumps({"esearchresult": p}).encode() for p in pages]):
        with pytest.raises(ValueError):
            client.search("nuccore", "query")


def test_fetch_truncation_fails(tmp_path: Path) -> None:
    client = Client(tmp_path)
    with patch.object(client, "summary", return_value={"1": {"accessionversion": "AB123456.1"}}), patch.object(client, "get", return_value=b""):
        with pytest.raises(ValueError, match="Incomplete"):
            client.fetch(["1"])


def test_http_permanent_failure_not_retried(tmp_path: Path) -> None:
    import urllib.error
    client = Client(tmp_path)
    error = urllib.error.HTTPError("https://example.invalid", 400, "bad input", {}, None)  # type: ignore[arg-type]  # Headers unused by test.
    with patch("shav.ingest.urllib.request.urlopen", side_effect=error) as fetch:
        with pytest.raises(urllib.error.HTTPError):
            client.get("esearch.fcgi", db="nuccore")
        assert fetch.call_count == 1


def test_transient_failure_retries_then_cache(tmp_path: Path) -> None:
    import urllib.error
    from unittest.mock import MagicMock
    from shav.common import ROOT
    # Place cache under project root to exercise manifest-relative paths.
    cache = ROOT / ".cache" / ("test-client-" + tmp_path.name)
    client = Client(cache)
    response = MagicMock()
    response.__enter__.return_value.read.return_value = b'{"ok":true}'
    error = urllib.error.HTTPError("https://example.invalid", 503, "unavailable", {}, None)  # type: ignore[arg-type]  # Headers unused by test.
    try:
        with patch("shav.ingest.urllib.request.urlopen", side_effect=[error, response]) as fetch, patch("shav.ingest.time.sleep"):
            assert client.get("test", db="nuccore") == b'{"ok":true}'
            assert fetch.call_count == 2
            assert client.get("test", db="nuccore") == b'{"ok":true}'
            assert fetch.call_count == 2
        assert len(client.requests[0]["sha256"]) == 64
    finally:
        for file in cache.glob("*"):
            file.unlink()
        if cache.exists():
            cache.rmdir()


def test_source_geography_normalization() -> None:
    value = record()
    value.features[0].qualifiers["geo_loc_name"] = ["Japan: Kagoshima"]
    row = metadata(value, "SHAV", [])
    assert row["country"] == "Japan"
    assert row["location"] == "Kagoshima"
    assert row["geo_loc_raw"] == "Japan: Kagoshima"


def test_cds_product_assigns_s_and_rejects_conflict() -> None:
    value = record()
    value.description = "Shamonda virus nucleocapsid gene"
    value.features[0].qualifiers.pop("segment")
    value.features.append(SeqFeature(type="CDS", qualifiers={"product": ["nucleocapsid protein"]}))
    row = metadata(value, "SHAV", [])
    assert row["segment"] == "S"
    assert "CDS product" in row["assignment_provenance"]
    value.features[0].qualifiers["segment"] = ["L"]
    with pytest.raises(ValueError, match="Conflicting"):
        segment(value)


def test_cds_polymerase_assigns_l() -> None:
    value = record()
    value.description = "Shamonda virus partial polymerase gene"
    value.features[0].qualifiers.pop("segment")
    value.features.append(SeqFeature(type="CDS", qualifiers={"product": ["RNA-dependent RNA polymerase"]}))
    assert segment(value) == "L"


def genbank_bytes(values: list[SeqRecord]) -> bytes:
    output = io.StringIO()
    for value in values:
        value.annotations["molecule_type"] = "RNA"
    SeqIO.write(values, output, "genbank")
    return output.getvalue().encode()


@pytest.mark.parametrize("actual", [["AB123457.1"], ["AB123456.2"], ["AB123456.1", "AB123456.1"]])
def test_fetch_rejects_substituted_accessions(tmp_path: Path, actual: list[str]) -> None:
    client = Client(tmp_path)
    with patch.object(client, "get", return_value=genbank_bytes([record(identity) for identity in actual])):
        with pytest.raises(ValueError, match="substituted"):
            client.fetch(["AB123456.1"])
    assert client.identities == {}


def test_fetch_uid_mapping_preserves_reordered_provenance(tmp_path: Path) -> None:
    client = Client(tmp_path)
    summary = {"1": {"accessionversion": "AB123456.1"}, "2": {"accessionversion": "AB123457.1"}}
    with patch.object(client, "summary", return_value=summary), patch.object(client, "get", return_value=genbank_bytes([record("AB123457.1"), record()])):
        values = client.fetch(["1", "2"])
    origins = {client.identities[uid]: labels for uid, labels in {"1": ["taxon"], "2": ["keyword"]}.items()}
    assert [metadata(value, "SHAV", origins[str(value.id)])["provenance"] for value in values] == ["keyword", "taxon"]


def test_summary_rejects_wrong_uid_set(tmp_path: Path) -> None:
    client = Client(tmp_path)
    with patch.object(client, "get", return_value=b'{"result":{"uids":["2"],"2":{}}}'):
        with pytest.raises(ValueError, match="identity mapping"):
            client.summary("nuccore", ["1"])


def sra_package(experiment: str = "SRX1", sample: str = "SRS1", runs: str = '<RUN accession="SRR1"/>') -> bytes:
    return f'<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><EXPERIMENT accession="{experiment}"><IDENTIFIERS><PRIMARY_ID>{experiment}</PRIMARY_ID></IDENTIFIERS></EXPERIMENT><SAMPLE accession="{sample}"><IDENTIFIERS><PRIMARY_ID>{sample}</PRIMARY_ID><EXTERNAL_ID namespace="BioSample">SAMN1</EXTERNAL_ID></IDENTIFIERS></SAMPLE><RUN_SET>{runs}</RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>'.encode()


@pytest.mark.parametrize("body", [b'<EXPERIMENT_PACKAGE_SET/>', sra_package("SRX2"), sra_package(sample="SRS2"), sra_package(runs=""), sra_package(runs='<RUN accession="bad"/>'), sra_package(runs='<RUN accession="SRR1"/><RUN accession="SRR1"/>')])
def test_sra_incomplete_or_substituted_packages_fail(body: bytes) -> None:
    expected = {"SRX1": {"runs": "SRR1", "uid": "1", "sample": "SRS1", "biosample": "SAMN1", "bioproject": ""}}
    with pytest.raises(ValueError):
        parse_sra_packages(body, expected)


def test_sra_verified_package() -> None:
    expected = {"SRX1": {"runs": "SRR1", "uid": "1", "sample": "SRS1", "biosample": "SAMN1", "bioproject": ""}}
    assert parse_sra_packages(sra_package(), expected)[0]["run"] == "SRR1"


def test_discover_reads_explicit_link_finds_taxon_missed_run(tmp_path: Path) -> None:
    client = Client(tmp_path)
    row = metadata(record(), "SHAV", ["taxon"])
    row.update(biosample="SAMN1", bioproject="PRJNA1")
    def search(db: str, query: str) -> list[str]:
        assert db == "sra"
        return ["42"] if query == '\"SAMN1\"[All Fields]' else []
    context = metadata(record("AB999999.1"), "keyword context", [])
    context.update(biosample="SAMN_CONTEXT", bioproject="PRJNA_CONTEXT")
    with patch.object(client, "search", side_effect=search):
        ids, inventory = discover_reads(client, [row, context])
    assert ids == ["42"]
    assert inventory["taxon_name"]["count"] == 0
    assert inventory["biosample:SAMN1"]["uids"] == ["42"]
    assert "bioproject:PRJNA1" in inventory
    assert "biosample:SAMN_CONTEXT" not in inventory
    assert "bioproject:PRJNA_CONTEXT" not in inventory
    assert "AB123456.1" in inventory["nucleotide_accessions:0"]["query"]


def test_sra_missing_run_within_present_experiment_fails() -> None:
    expected = {"SRX1": {"runs": "SRR1;SRR2", "uid": "1", "sample": "SRS1", "biosample": "SAMN1", "bioproject": ""}}
    with pytest.raises(ValueError, match="Incomplete or substituted SRA runs"):
        parse_sra_packages(sra_package(), expected)


def test_pinned_snapshot_rejects_changed_ingestion_code(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from shav import ingest
    from shav.common import sha256
    (tmp_path / "config").mkdir()
    (tmp_path / "data/public").mkdir(parents=True)
    config = tmp_path / "config/ingest.json"
    config.write_text("{}")
    (tmp_path / "data/public/manifest.json").write_text(json.dumps({"outputs": {}, "config_sha256": sha256(config), "ingest_source_sha256": "changed"}))
    monkeypatch.setattr(ingest, "ROOT", tmp_path)
    with pytest.raises(ValueError, match="Ingestion code changed"):
        ingest.run()


def test_identifier_search_false_positive_is_audited_and_excluded() -> None:
    false_hit = {"experiment": "ERX17003618", "biosample": "SAMEA123171705", "bioproject": "PRJEB121521", "run": "ERR17611905"}
    verified = {"experiment": "SRX1", "biosample": "SAMN1", "bioproject": "PRJNA1", "run": "SRR1"}
    substring_hit = {"experiment": "SRX2", "biosample": "SAMN10", "bioproject": "PRJEB1231740", "run": "SRR2"}
    included, excluded = select_linked_reads([false_hit, substring_hit, verified], set(), {"SAMN1"}, {"PRJEB123174"})
    assert included == [verified]
    assert excluded[0]["run"] == "ERR17611905"
    assert "lacks matching" in excluded[0]["reason"]
    assert excluded[1]["run"] == "SRR2"


def test_specimen_lab_host_and_passage_remain_separate() -> None:
    value = record()
    value.features[0].qualifiers.update(isolation_source=['blood', 'serum'], lab_host=['Vero', 'BHK-21'], passage_history=['P1', 'P2'])
    row = metadata(value, 'SHAV', [])
    assert row['isolation_source'] == 'blood; serum'
    assert row['lab_host'] == 'Vero; BHK-21'
    assert row['passage'] == 'P1; P2'
    assert row['host'] == ''
    blank = metadata(record(), 'SHAV', [])
    assert all(blank[key] == '' for key in ['isolation_source', 'lab_host', 'passage', 'host'])

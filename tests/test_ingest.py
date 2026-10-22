from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from shav.ingest import Client, metadata, normalize_date, reconcile, segment


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
    assert row["passage"] == "human cell line"
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
    with patch.object(client, "get", return_value=b""):
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

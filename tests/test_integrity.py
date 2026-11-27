import hashlib
import json
from unittest.mock import Mock

import pytest
import requests

from cosmic_anomaly_lab.literature import deduplicate
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest
from cosmic_anomaly_lab.registry import Registry


def response(content):
    r = Mock()
    r.headers = {"Content-Length": str(len(content))}
    r.url = "https://example.org/catalogue"
    r.iter_content.return_value = iter([content])
    r.__enter__ = Mock(return_value=r)
    r.__exit__ = Mock(return_value=False)
    return r


def test_download_budget_checksum_and_cache(tmp_path):
    cache = DownloadCache(tmp_path, budget_bytes=10)
    cache.session.get = Mock(return_value=response(b"abc"))
    path = cache.fetch(
        "https://example.org/catalogue",
        "x",
        survey="test",
        version="1",
        kind="synthetic_test",
        checksum="md5:" + hashlib.md5(b"abc").hexdigest(),
    )
    assert path.read_bytes() == b"abc"
    cache.fetch(
        "https://example.org/catalogue", "x", survey="test", version="1", kind="synthetic_test"
    )
    assert cache.session.get.call_count == 1
    records = json.loads(cache.manifest.read_text())
    assert records[0]["bytes"] == 3
    path.write_bytes(b"def")
    with pytest.raises(ValueError, match="Cache checksum"):
        cache.fetch(
            "https://example.org/catalogue", "x", survey="test", version="1", kind="synthetic_test"
        )
    cache.session.get = Mock(return_value=response(b"123456789"))
    with pytest.raises(ValueError, match="budget"):
        cache.fetch("https://example.org/y", "y", survey="test", version="1", kind="synthetic_test")
    assert not (cache.raw / "y").exists()


def test_post_download_reuses_the_recorded_query(tmp_path):
    cache = DownloadCache(tmp_path, budget_bytes=10)
    cache.session.post = Mock(return_value=response(b"[1]"))
    cache.session.get = Mock()
    url = "https://example.org/catalogue"
    metadata = dict(survey="test", version="1", kind="synthetic_test")
    request_a = {"search": [["band", "i"]]}
    request_b = {"search": [["band", "z"]]}
    first = cache.fetch(url, "i.json", json_payload=request_a, **metadata)
    assert first.read_bytes() == b"[1]"
    assert cache.fetch(url, "i.json", json_payload=request_a, **metadata) == first
    cache.session.post.assert_called_once_with(
        url, json=request_a, stream=True, timeout=(20, 90)
    )
    cache.session.post.return_value = response(b"[2]")
    second = cache.fetch(url, "z.json", json_payload=request_b, **metadata)
    assert second.read_bytes() == b"[2]"
    records = json.loads(cache.manifest.read_text())
    assert [record["request_json"] for record in records] == [request_a, request_b]
    cache.session.get.assert_not_called()


def test_failed_checksum_leaves_no_committed_file(tmp_path):
    cache = DownloadCache(tmp_path)
    cache.session.get = Mock(return_value=response(b"abc"))
    with pytest.raises(ValueError, match="Published checksum"):
        cache.fetch(
            "https://example.org/x",
            "x",
            survey="t",
            version="1",
            kind="test",
            checksum="sha256:wrong",
        )
    assert not (cache.raw / "x").exists() and not (cache.raw / "x.part").exists()


def test_form_post_records_and_reuses_its_request_body(tmp_path):
    cache = DownloadCache(tmp_path, budget_bytes=100)
    cache.session.post = Mock(side_effect=[response(b"first"), response(b"second"), response(b"third")])
    url = "https://example.org/position"
    metadata = dict(survey="test", version="1", kind="synthetic_test")
    first = {"ra": "53.19", "dec": "-27.81"}
    second = {"ra": "53.20", "dec": "-27.81"}
    path = cache.fetch(url, "first", form_payload=first, **metadata)
    assert cache.fetch(url, "first", form_payload=first, **metadata) == path
    cache.session.post.assert_called_once_with(url, data=first, stream=True, timeout=(20, 90))
    assert cache.fetch(url, "second", form_payload=second, **metadata).read_bytes() == b"second"
    assert cache.fetch(url, "json", json_payload=first, **metadata).read_bytes() == b"third"
    assert [r.get("request_form") for r in cache.records] == [first, second, None]
    assert [r.get("request_json") for r in cache.records] == [None, None, first]
    with pytest.raises(ValueError, match="only one POST"):
        cache.fetch(url, "invalid", form_payload=first, json_payload=first, **metadata)


def test_prepared_get_uses_the_same_budget_and_provenance(tmp_path):
    cache = DownloadCache(tmp_path, budget_bytes=4)
    cache.session.get = Mock(side_effect=AssertionError("Should reuse the prepared GET"))
    url = "https://example.org/catalogue"
    fetched = response(b"abc")
    fetched.history = []
    fetched.request = requests.Request("GET", url).prepare()
    metadata = dict(survey="test", version="1", kind="synthetic_test")
    path = cache.fetch(url, "file", response=fetched, **metadata)
    assert path.read_bytes() == b"abc"
    assert cache.records[0]["sha256"] == digest(path)
    with pytest.raises(ValueError, match="does not match"):
        cache.fetch("https://example.org/other", "wrong", response=fetched, **metadata)
    too_large = response(b"123")
    too_large.history = []
    too_large.request = fetched.request
    with pytest.raises(ValueError, match="budget"):
        cache.fetch(url, "larger", response=too_large, **metadata)
    assert not (cache.raw / "larger").exists()


def test_download_path_traversal_refused(tmp_path):
    cache = DownloadCache(tmp_path)
    with pytest.raises(ValueError, match="inside"):
        cache.fetch("https://example.org/x", "../../escape", survey="t", version="1", kind="test")


def test_registry_multiple_observations_and_idempotency(tmp_path):
    r = Registry(tmp_path / "r.sqlite")
    sid = r.add_source("survey", "1", 359.9, 10, redshift=None)
    assert r.add_source("survey", "1", 359.9, 10, redshift=None) == sid
    r.add_observation("a", sid, "survey", "field", "https://example.org/a")
    r.add_observation("b", sid, "survey", "field", "https://example.org/b")
    r.add_evidence(sid, "run", "INSUFFICIENT_EVIDENCE", {"observed": True})
    assert len(r.search("survey")) == 1
    assert r.db.execute("select count(*) from observations").fetchone()[0] == 2
    with pytest.raises(ValueError):
        r.add_evidence(sid, "run", "DISCOVERY", {})
    r.close()


def test_deduplicate_identifiers_not_titles():
    out = deduplicate(
        [
            {"title": "A", "doi": "https://doi.org/10.1/ABC", "arxiv": "2401.00001v2"},
            {"title": "B", "doi": "10.1/abc"},
            {"title": "A", "arxiv": "2402.00001"},
        ]
    )
    assert len(out) == 2


def test_canonical_text_hash_survives_git_newlines(tmp_path):
    a, b = tmp_path / "windows", tmp_path / "unix"
    a.write_bytes(b"one\r\ntwo\r\n")
    b.write_bytes(b"one\ntwo\n")
    assert text_digest(a) == text_digest(b)
    assert digest(a) != digest(b)

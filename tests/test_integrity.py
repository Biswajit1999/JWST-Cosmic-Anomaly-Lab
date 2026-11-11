import hashlib
import json
from unittest.mock import Mock

import pytest

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

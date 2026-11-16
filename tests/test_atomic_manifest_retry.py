import json
from pathlib import Path

import pytest

from cosmic_anomaly_lab.provenance import write_json


def test_transient_lock_preserves_atomic_json(tmp_path, monkeypatch):
    path = tmp_path / "manifest.json"
    path.write_text('{"old": true}')
    replace = Path.replace
    calls = []

    def flaky(source, destination):
        calls.append(1)
        if len(calls) < 3:
            assert json.loads(path.read_text()) == {"old": True}
            raise PermissionError("simulated temporary file lock")
        return replace(source, destination)

    monkeypatch.setattr(Path, "replace", flaky)
    monkeypatch.setattr("cosmic_anomaly_lab.provenance.time.sleep", lambda _: None)
    write_json(path, {"new": True})
    assert len(calls) == 3
    assert json.loads(path.read_text()) == {"new": True}


def test_persistent_lock_fails_without_losing_previous_manifest(tmp_path, monkeypatch):
    path = tmp_path / "manifest.json"
    path.write_text('{"old": true}')

    def locked(*args):
        raise PermissionError("locked")

    monkeypatch.setattr(Path, "replace", locked)
    monkeypatch.setattr("cosmic_anomaly_lab.provenance.time.sleep", lambda _: None)
    with pytest.raises(PermissionError):
        write_json(path, {"new": True})
    assert json.loads(path.read_text()) == {"old": True}

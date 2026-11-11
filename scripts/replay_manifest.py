"""Restore the recorded product bytes, refusing unannounced upstream changes."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cosmic_anomaly_lab.provenance import DownloadCache

root = Path(__file__).resolve().parents[1]
records = json.loads((root / "data/manifests/downloads.json").read_text())
latest = {r["path"]: r for r in records}
cache = DownloadCache(root)
for r in latest.values():
    name = Path(r["path"]).relative_to("data/raw").as_posix()
    cache.fetch(
        r["url"],
        name,
        survey=r["survey"],
        version=r["version"],
        kind=r["kind"],
        checksum="sha256:" + r["sha256"],
    )
    print(name, flush=True)

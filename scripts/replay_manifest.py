"""Restore the recorded product bytes, refusing unannounced upstream changes."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cosmic_anomaly_lab.provenance import DownloadCache

root = Path(__file__).resolve().parents[1]
records = json.loads((root / "data/manifests/downloads.json").read_text())
latest = {r["path"]: r for r in records}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--budget-mib",
    type=int,
    default=4096,
    help="Finite raw-data cap; expanded search requires about 2 GiB",
)
parser.add_argument(
    "--balmer-comparison",
    action="store_true",
    help="Restore only the published DJA/DR4 inputs for the four-target Balmer comparison.",
)
args = parser.parse_args()
if args.balmer_comparison:
    dja = json.loads((root / "results/multiwavelength/dja_products.json").read_text())
    dr4 = json.loads((root / "results/jades_balmer_census/spectral_products.json").read_text())
    target_ids = {record["dr4_id"] for record in dja}
    inputs = dja + [record for record in dr4 if record["dr4_id"] in target_ids
                    and record["status"] == "DOWNLOADED"]
    selected = {}
    for item in inputs:
        record = latest[item["path"]]
        if record["sha256"] != item["sha256"]:
            raise ValueError(f"Conflicting checksum receipts for {item['path']}")
        selected[item["path"]] = record
    latest = selected
cache = DownloadCache(root, budget_bytes=args.budget_mib * 1024**2)
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

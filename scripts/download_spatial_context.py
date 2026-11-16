"""Reproduce image and nearby-spectrum retrieval from stored source URLs."""

import json
from pathlib import Path

from cosmic_anomaly_lab.provenance import DownloadCache, digest

root = Path(__file__).resolve().parents[1]
cache = DownloadCache(root, 4096*1024**2)
products = json.loads((root / "results/spatial_context/image_products.json").read_text())["products"]
products.append(json.loads((root / "results/spatial_context/neighbor_product.json").read_text()))
for p in products:
    if p.get("status", "DOWNLOADED") != "DOWNLOADED":
        continue
    path = cache.fetch(p["url"], p["path"].removeprefix("data/raw/"), survey="DJA",
                       version="2026-11-16-snapshot", kind="image_or_neighbor_spectrum")
    assert digest(path) == p["sha256"]

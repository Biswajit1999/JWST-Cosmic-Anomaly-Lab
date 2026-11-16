"""Provenance and reporting checks for spatial and exoplanet research outputs."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
exo = root / "results/exoplanet_pilot"
spatial = root / "results/spatial_context"
products = json.loads((exo / "products.json").read_text())["products"]
products += json.loads((spatial / "image_products.json").read_text())["products"]
products += [json.loads((p).read_text()) for p in [exo / "pixel_product.json", exo / "atmosphere_product.json", spatial / "neighbor_product.json"]]
verified = 0
for p in products:
    if p.get("status", "DOWNLOADED") != "DOWNLOADED":
        continue
    assert digest(root / p["path"]) == p["sha256"]
    if "published_checksum" in p:
        algorithm, checksum = p["published_checksum"].split(":")
        assert digest(root / p["path"], algorithm) == checksum
    verified += 1
for output, script in [(exo / "summary.json", "analyse_exoplanet_pilot.py"),
                        (exo / "injection_summary.json", "check_transit_recovery.py"),
                        (exo / "control_audit.json", "audit_transit_controls.py"),
                        (exo / "pixel_summary.json", "image_kepler_control.py"),
                        (exo / "atmosphere_summary.json", "analyse_atmosphere_benchmark.py"),
                        (spatial / "summary.json", "analyse_spatial_context.py"),
                        (spatial / "neighbor_summary.json", "analyse_neighbor_spectra.py")]:
    record = json.loads(output.read_text())
    assert record["script_sha256"] == text_digest(root / "scripts" / script)
    if "module_sha256" in record:
        assert record["module_sha256"] == text_digest(root / "src/cosmic_anomaly_lab/transits.py")
    inputs = record.get("input_sha256", {})
    if isinstance(inputs, dict):
        for path, expected in inputs.items():
            assert digest(root / path) == expected
screen = pd.read_csv(exo / "screening.csv")
metadata = pd.read_csv(exo / "input_metadata.csv")
assert len(screen) == 13 and len(metadata) == 24
assert int(metadata.retained.sum()) == 80650
assert (screen.status == "INSUFFICIENT_TWO_QUARTER_DATA").sum() == 1
assert not screen.loc[screen.role == "search", "screen_pass"].any()
known = screen[screen.kic == 11446443].iloc[0]
assert known.status == "KNOWN_SIGNAL_RECOVERY" and not known.screen_pass
assert json.loads(known.matched_koi)[0]["disposition"] == "CONFIRMED"
injections = pd.read_csv(exo / "injection_recovery.csv")
assert injections.period_recovered.all() and len(injections) == 3
assert injections.kind.str.startswith("SYNTHETIC").all()
geometry = pd.read_csv(spatial / "slit_geometry.csv")
assert geometry.loc[geometry.slit_target == geometry.tested_target, "centroid_inside"].all()
assert not geometry.loc[geometry.slit_target != geometry.tested_target, "centroid_inside"].any()
neighbor = pd.read_csv(spatial / "catalogue_neighbors.csv")
assert neighbor.loc[neighbor.neighbor.str.endswith("_27359"), "small_dz_velocity_proxy_kms"].isna().all()
spectra = pd.read_csv(exo / "wasp39_transmission.csv")
assert len(spectra) == 344 and np.isfinite(spectra).all().all()
assert not json.loads((exo / "wasp39_properties.json").read_text())["new_abundances_retrieved"]
raw_bytes = sum(p.stat().st_size for p in (root / "data/raw").rglob("*") if p.is_file())
assert raw_bytes < 4096*1024**2
artifacts = []
for folder in [exo, spatial, root / "figures/exoplanet_pilot", root / "figures/spatial_context"]:
    artifacts.extend(p for p in folder.iterdir() if p.is_file() and p.name != "verification.json")
artifacts += [root / "reports/spatial_context.md", root / "reports/exoplanet_pilot.md"]
write_json(exo / "verification.json", dict(executed_utc=utcnow(), checks="PASS", raw_products_verified=verified,
    raw_storage_bytes=raw_bytes, raw_storage_budget_bytes=4096*1024**2,
    artifact_sha256={p.relative_to(root).as_posix(): text_digest(p) if p.suffix in [".csv", ".json", ".md"] else digest(p)
                     for p in artifacts}, verifier_sha256=text_digest(Path(__file__))))
print(f"PASS: {verified} raw products, experiment scripts, all reported counts and missing cases, storage budget {raw_bytes} bytes")

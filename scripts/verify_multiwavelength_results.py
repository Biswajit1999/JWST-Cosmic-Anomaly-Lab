"""Audit product fingerprints, fixed-sample counts and missing-data outcomes."""

import json
from pathlib import Path

import pandas as pd

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
out = root / "results/multiwavelength"
dja = json.loads((out / "dja_summary.json").read_text())
context = json.loads((out / "context_summary.json").read_text())
for name, summary in [("check_balmer_dja_reductions.py", dja), ("analyse_multiwavelength_context.py", context)]:
    assert text_digest(root / "scripts" / name) == summary["script_sha256"]
    for path, expected in summary["input_sha256"].items():
        fingerprint = (digest(root / path) if name == "analyse_multiwavelength_context.py"
                       or path.endswith(".gz") else text_digest(root / path))
        assert fingerprint == expected, path
products = json.loads((out / "dja_products.json").read_text())
for product in products:
    assert digest(root / product["path"]) == product["sha256"]
    assert product["separation_arcsec"] < .2
assert len(products) == 6
assert (dja["core_attempts"], dja["measured_cores"]) == (72, 60)
assert (context["photometric_band_rows"], context["valid_signed_band_rows"]) == (88, 74)
assert context["negative_valid_flux_rows"] == 17
associations = pd.read_csv(out / "association_summary.csv")
assert len(associations) == 12 and associations.xray_candidates.sum() == 0
assert associations.photometric_candidates.sum() == 8
ratio = pd.read_csv(out / "dja_core_ratios.csv")
assert ratio.loc[~ratio.eligible, "ratio"].isna().all()
assert not ratio[ratio.eligible][["ratio", "error"]].isna().any().any()
used = sum(p.stat().st_size for p in (root / "data/raw").rglob("*") if p.is_file())
assert used < 4096 * 1024**2
write_json(out / "verification.json", dict(executed_utc=utcnow(), status="PASS",
    raw_storage_bytes=used, raw_storage_budget_bytes=4096*1024**2,
    checks=["script and input fingerprints", "six product fingerprints and positions",
            "fixed-sample measurement counts", "signed negative photometry preserved",
            "ineligible ratios remain missing", "no inferred X-ray upper limits",
            "raw cache within budget"]))
print("PASS", used, "raw bytes; no physical classification inferred")

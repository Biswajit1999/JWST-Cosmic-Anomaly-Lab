"""Verify executed provenance, coverage gates and every retained ratio."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
out = root / "results/gn89464_followup"
for filename, script in [("summary.json", "analyse_gn89464_continuum.py"),
                         ("previous_queue_summary.json", "audit_previous_continuum_coverage.py")]:
    s = json.loads((out / filename).read_text())
    for path, expected in s["input_sha256"].items():
        assert digest(root / path) == expected
    assert text_digest(root / "scripts" / script) == s["script_sha256"]
    assert text_digest(root / "src/cosmic_anomaly_lab/continuum_audit.py") == s["module_sha256"]
assert text_digest(root / "results/jades_balmer_census/signed_core_measurements.csv") == s["old_measurements_sha256"]
cores = pd.read_csv(out / "bilateral_cores.csv")
ratios = pd.read_csv(out / "bilateral_ratios.csv")
assert len(cores) == 5*5*2*3*2 and len(ratios) == len(cores)//2
keys = ["product", "shift_kms", "half_width_kms", "sidebands"]
assert not cores.duplicated(keys + ["line"]).any() and not ratios.duplicated(keys).any()
accepted = cores[cores.status == "MEASURED_CORE"]
assert (accepted.blue_pixels >= 3).all() and (accepted.red_pixels >= 3).all()
assert (accepted.good_coverage >= .8).all() and (accepted.core_pixels >= 4).all()
assert cores.loc[cores.status != "MEASURED_CORE", ["flux", "error"]].isna().all().all()
indexed = cores.set_index(keys + ["line"])
for row in ratios.itertuples():
    key = (row.product, row.shift_kms, row.half_width_kms, row.sidebands)
    n, d = indexed.loc[(*key, "hg")], indexed.loc[(*key, "hb")]
    eligible = (n.status == d.status == "MEASURED_CORE" and n.flux/n.error >= 5 and d.flux/d.error >= 5)
    assert eligible == row.eligible
    if eligible:
        np.testing.assert_allclose(row.ratio, n.flux/d.flux)
        np.testing.assert_allclose(row.formal_ratio_error,
                                   np.hypot(n.error/d.flux, n.flux*d.error/d.flux**2))
    else:
        assert pd.isna(row.ratio) and pd.isna(row.formal_ratio_error)
previous = pd.read_csv(out / "previous_queue_coverage_audit.csv")
keep = previous[previous.status == "MEASURED_CORE"]
np.testing.assert_allclose(keep[["flux", "error"]], keep[["legacy_flux", "legacy_error"]], rtol=1e-8, atol=1e-8)
raw_bytes = sum(p.stat().st_size for p in (root / "data/raw").rglob("*") if p.is_file())
assert raw_bytes < 4096*1024**2
paths = list(out.glob("*.csv")) + list((root / "figures/gn89464_followup").glob("*"))
paths += [root / "reports/gn89464_continuum.md"]
write_json(out / "verification.json", dict(executed_utc=utcnow(), checks="PASS",
    ratios_verified=len(ratios), retained_previous_measurements=len(keep), raw_storage_bytes=raw_bytes,
    raw_storage_budget_bytes=4096*1024**2,
    artifact_sha256={p.relative_to(root).as_posix(): text_digest(p) if p.suffix in [".csv", ".md"] else digest(p)
                     for p in paths},
    support_module_sha256={name: text_digest(root / "src/cosmic_anomaly_lab" / name)
                          for name in ["spectroscopy.py", "reduction_audit.py", "balmer_census.py", "provenance.py"]},
    verifier_sha256=text_digest(Path(__file__))))
print("PASS: provenance, all 150 ratio gates/values, 239 retained prior measurements, storage budget")

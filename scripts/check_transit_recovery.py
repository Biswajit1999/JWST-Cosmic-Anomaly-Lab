"""Explicit synthetic injections into real data; never count as sky detections."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.transits import flatten_segments, phase_distance, scan_transits

root = Path(__file__).resolve().parents[1]
out = root / "results/exoplanet_pilot"
kid = json.loads((out / "selection.json").read_text())["search_targets"][0]
products = json.loads((out / "products.json").read_text())["products"]
records = []
for period in [1.7, 5.3, 12.7]:
    arrays = []
    for p in products:
        if p["kic"] != kid or p["status"] != "DOWNLOADED":
            continue
        with fits.open(root / p["path"]) as h:
            d = h[1].data
            transit = abs(phase_distance(d["TIME"], period, 150.3)) < .05
            flux = d["PDCSAP_FLUX"] * (1-.001*transit)
            ix, y, error = flatten_segments(d["TIME"], flux, d["PDCSAP_FLUX_ERR"], d["SAP_QUALITY"])
            arrays.append((np.array(d["TIME"][ix]), y, error))
    t, y, e = [np.concatenate([a[j] for a in arrays]) for j in range(3)]
    result, _, _, _ = scan_transits(t, y, e)
    records.append(dict(kic=kid, kind="SYNTHETIC_INJECTION_NOT_OBSERVED_PLANET", injected_period=period,
        injected_depth=.001, injected_duration=.1, recovered_period=result["period_days"],
        recovered_depth=result["depth"], depth_snr=result["depth_snr"],
        period_recovered=bool(abs(result["period_days"]/period-1)<.01)))
    pd.DataFrame(records).to_csv(out / "injection_recovery.csv", index=False)
    print(records[-1], flush=True)
write_json(out / "injection_summary.json", dict(executed_utc=utcnow(), injections=len(records),
    recovered=sum(r["period_recovered"] for r in records), script_sha256=text_digest(Path(__file__)),
    module_sha256=text_digest(root / "src/cosmic_anomaly_lab/transits.py"),
    full_completeness_or_false_alarm_calibration=False))

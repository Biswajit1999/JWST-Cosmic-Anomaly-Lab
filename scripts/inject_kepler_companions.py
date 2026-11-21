"""Labelled transit injections into real Kepler cadences for one fixed host."""

import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from astropy.io import fits
from search_kepler_companions import estimate_candidates, events_at, meets, out, root, split_curves

from cosmic_anomaly_lab.companion_search import prepare_companion_quarter
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths
from cosmic_anomaly_lab.transits import phase_distance


def experiment(task):
    star, products, template = task
    curves = []
    period, epoch, duration, amplitude = (template[k] for k in ["period_days", "epoch_bkjd", "duration_days", "depth"])
    for p in products:
        file = root / p["path"]
        assert digest(file) == p["sha256"]
        with fits.open(file) as h:
            assert h[0].header["KEPLERID"] == star["kic"] and h[1].header["BJDREFI"] == 2454833
            assert h[1].header["TIMESYS"] == "TDB"
            q, cadence = int(h[0].header["QUARTER"]), float(h[1].header["TIMEDEL"])
            d = h[1].data
            x = phase_distance(d["TIME"], period, epoch)
            overlap = np.maximum(0, np.minimum(x+cadence/2, duration/2)-np.maximum(x-cadence/2, -duration/2))
            flux = d["PDCSAP_FLUX"]*(1-amplitude*overlap/cadence)
            t, y, err, _, _ = prepare_companion_quarter(d["TIME"], flux, d["PDCSAP_FLUX_ERR"],
                                                       d["SAP_QUALITY"], star["masks"], cadence)
        curves.append(dict(time=t, flux=y, error=err, cadence=cadence, quarter=q))
    halves = split_curves(sorted(curves, key=lambda c: c["quarter"]))
    assert all(halves.values())
    checks = []
    for fit in estimate_candidates(halves["early"], 1):
        event_lists = {label: events_at(group, fit) for label, group in halves.items()}
        summaries = {label: combine_event_depths(events) for label, events in event_lists.items()}
        mids = np.array([r["midpoint"] for half in event_lists.values() for r in half if r["status"] == "MEASURED"])
        offset = abs(phase_distance(mids, period, epoch))
        timing = bool(len(mids) and (offset <= duration/2).mean() >= .8)
        same_signal = bool(abs(fit["period_days"]/period - 1) < .005 and timing)
        checks.append(dict(**fit, period_and_timing_match=same_signal,
            later_depth=summaries["later"]["depth"], later_contrast=summaries["later"]["contrast"],
            replicated=all(meets(s) for s in summaries.values())))
    return dict(kic=star["kic"], kind="SYNTHETIC_INJECTION", injected=template,
                recovered=any(t["period_and_timing_match"] and t["replicated"] for t in checks), trials=checks)


if __name__ == "__main__":
    roster = json.loads((out / "selection.json").read_text())["targets"]
    star = next(s for s in roster if s["role"] == "search_host")
    products = next(t["products"] for t in json.loads((out / "products.json").read_text())["targets"]
                    if t["kic"] == star["kic"])
    injections = [dict(period_days=p, epoch_bkjd=e, duration_days=d, depth=a)
                  for p, e, d in [(8.137, 151.394, .14), (24.691, 153.753, .24)]
                  for a in [.0002, .0005]]
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(experiment, [(star, products, trial) for trial in injections]))
    write_json(out / "injections.json", dict(executed_utc=utcnow(), experiments=results,
        script_sha256=text_digest(__file__), input_sha256={p["path"]: p["sha256"] for p in products},
        scope="Two periods and two depths in one star. This is not a completeness measurement."))
    summary = pd.DataFrame([dict(**r["injected"], recovered=r["recovered"]) for r in results])
    summary.to_csv(out / "injections.csv", index=False)
    print(summary.to_string(index=False))

"""Search masked Kepler light curves, preserving disjoint replication data."""

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from matplotlib.figure import Figure

from cosmic_anomaly_lab.companion_search import (
    duration_resolved_period_grid,
    prepare_companion_quarter,
)
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths, measure_events
from cosmic_anomaly_lab.transits import phase_distance
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_companions"
config_path = root / "configs/kepler_companions.yaml"
config = yaml.safe_load(config_path.read_text())


def tasks(wait):
    selection = json.loads((out / "selection.json").read_text())
    assert selection["config_sha256"] == text_digest(config_path)
    for i, star in enumerate(selection["targets"]):
        deadline = time.monotonic() + (300 if wait else 0)
        while True:
            ready = json.loads((out / "products.json").read_text())["targets"]
            if len(ready) > i:
                assert ready[i]["kic"] == star["kic"]
                yield star, ready[i]["products"]
                break
            if time.monotonic() >= deadline:
                raise ValueError(f"Missing inventory entry for {star['kic']}")
            time.sleep(2)


def load_curve(product, star, mask_multiplier):
    path = root / product["path"]
    assert digest(path) == product["sha256"]
    with fits.open(path) as h:
        header, timing, data = h[0].header, h[1].header, h[1].data
        assert header["KEPLERID"] == star["kic"]
        assert timing["BJDREFI"] == 2454833 and timing.get("BJDREFF", 0) == 0
        assert timing["TIMESYS"] == "TDB" and timing["TIMEPIXR"] == .5
        cadence = float(timing["TIMEDEL"])
        assert .019 < cadence < .023
        t, f, e, masked, valid = prepare_companion_quarter(
            data["TIME"], data["PDCSAP_FLUX"], data["PDCSAP_FLUX_ERR"], data["SAP_QUALITY"],
            star["masks"], cadence, mask_multiplier=mask_multiplier
        )
    return dict(quarter=int(header["QUARTER"]), cadence=cadence, time=t, flux=f, error=e,
                quality_cadences=valid, masked_cadences=masked)


def split_curves(curves):
    return {label: [c for c in curves if lo <= c["quarter"] <= hi and len(c["time"]) > 100]
            for label, (lo, hi) in [("early", config["discovery_quarters"]),
                                    ("later", config["replication_quarters"])]}


def estimate_candidates(curves, sign):
    t, y, error = (np.concatenate([c[key] for c in curves]) for key in ["time", "flux", "error"])
    flux = 1 + sign*(y-1)
    diff = np.diff(flux)
    scatter = 1.4826*np.median(abs(diff-np.median(diff)))/np.sqrt(2)
    model = BoxLeastSquares(t, flux, np.maximum(error, scatter))
    lo, hi = config["trial_period_days"]
    durations = config["trial_durations_days"]
    span = np.ptp(t)
    grid = duration_resolved_period_grid(lo, hi, span, min(durations))
    result = model.power(grid, durations, oversample=5, objective="snr")
    eligible = np.ones(len(grid), bool)
    for rank in range(1, config["trial_peaks_per_sign"]+1):
        if not eligible.any():
            break
        idx = int(np.argmax(np.where(eligible, result.power, -np.inf)))
        local = np.geomspace(grid[max(0, idx-3)], grid[min(len(grid)-1, idx+3)], 91)
        fit = model.power(local, durations, objective="snr", oversample=20)
        i = int(np.argmax(fit.power))
        period, epoch, duration = (float(fit[k][i]) for k in ["period", "transit_time", "duration"])
        yield dict(rank=rank, period_days=period, epoch_bkjd=epoch, duration_days=duration,
                   bls_contrast=float(fit.power[i]), trial_periods=len(grid))
        # Remove neighbouring periods within two fitted transit durations of drift over training.
        drift = abs(1/grid - 1/period)*span*period
        eligible &= drift > 2*duration


def events_at(curves, fit, sign=1, phase=0):
    return [dict(quarter=c["quarter"], **row) for c in curves for row in measure_events(
        c["time"], 1+sign*(c["flux"]-1), c["error"], fit["period_days"],
        fit["epoch_bkjd"]+phase*fit["period_days"], fit["duration_days"], c["cadence"])]


def meets(summary):
    return (summary["events"] >= config["minimum_events_per_half"]
            and summary["contrast"] >= config["depth_contrast_threshold"])


def plot(kic, fit, halves, label):
    fig = Figure(figsize=(10, 4.8), layout="constrained")
    for ax, (period_set, curves) in zip(fig.subplots(1, 2), halves.items(), strict=True):
        phase, rel_flux = [], []
        for c in curves:
            offset = phase_distance(c["time"], fit["period_days"], fit["epoch_bkjd"])
            selected = abs(offset) <= 3*fit["duration_days"]
            phase.extend(offset[selected]*24)
            rel_flux.extend((c["flux"][selected]-1)*1e6)
        ax.plot(phase, rel_flux, ".", ms=2, color="#2c5072", alpha=.25)
        axes = np.linspace(-3*fit["duration_days"]*24, 3*fit["duration_days"]*24, 25)
        points = pd.cut(phase, axes)
        binned = pd.DataFrame(dict(offset=phase, flux=rel_flux)).groupby(points, observed=True).median()
        ax.plot(binned.offset, binned.flux, "o", ms=3, color="#b45c27", label="Bin medians")
        ax.axvspan(-12*fit["duration_days"], 12*fit["duration_days"], color="#db8c32", alpha=.15)
        ax.set(title=period_set.capitalize()+" quarters", xlabel="Hours from predicted transit",
               ylabel="Observed relative flux (ppm)")
        ax.legend(fontsize=8)
    fig.suptitle(f"KIC {kic} · observed Kepler photometry · {fit['period_days']:.6f} days")
    save_figure(fig, root / f"figures/kepler_companions/kic{kic}_{label}")


def run_target(task):
    star, products = task
    kic = star["kic"]
    curves = sorted((load_curve(p, star, 1) for p in products), key=lambda c: c["quarter"])
    assert len({c["quarter"] for c in curves}) == len(curves)
    halves = split_curves(curves)
    spans = [v[-1]["time"][-1]-v[0]["time"][0] if v else 0 for v in halves.values()]
    data = dict(kic=kic, role=star["role"], files=len(curves), baseline_early_days=spans[0],
                baseline_later_days=spans[1], retained_cadences=sum(len(c["time"]) for c in curves),
                quality_cadences=sum(c["quality_cadences"] for c in curves),
                masked_cadences=sum(c["masked_cadences"] for c in curves))
    if min(spans) < 100:
        return dict(data, status="INSUFFICIENT_BASELINE"), [], []
    assert halves["early"][-1]["time"][-1] < halves["later"][0]["time"][0]
    fits_rows, event_rows = [], []
    wide_mask = None
    for sign, label in [(1, "observed"), (-1, "inverted")]:
        candidates = list(estimate_candidates(halves["early"], sign))
        for fit in candidates:
            record = dict(kic=kic, role=star["role"], sign=label, **fit)
            stats = {}
            for split, subset in halves.items():
                events = events_at(subset, fit, sign)
                stats[split] = combine_event_depths(events)
                event_rows.extend(dict(kic=kic, sign=label, rank=fit["rank"], split=split,
                                       variant="primary", **ev) for ev in events)
            passes = all(meets(s) for s in stats.values())
            record.update(replicated=passes,
                          **{f"{split}_{key}": v for split, s in stats.items() for key, v in s.items()})
            if passes and sign == 1:
                if wide_mask is None:
                    wide_mask = split_curves(sorted((load_curve(p, star, 1.5) for p in products),
                                                   key=lambda c: c["quarter"]))
                widened = {label: combine_event_depths(events_at(half, fit)) for label, half in wide_mask.items()}
                secondary = combine_event_depths(events_at(halves["later"], fit, phase=.5))
                record.update(wide_mask_replicated=all(meets(s) for s in widened.values()),
                              wide_mask_later_contrast=widened["later"]["contrast"],
                              half_phase_later_contrast=secondary["contrast"])
                plot(kic, fit, halves, f"rank{fit['rank']}")
            fits_rows.append(record)
    print(kic, star["role"], data["masked_cadences"],
          [(round(r["period_days"], 5), r["replicated"]) for r in fits_rows if r["sign"] == "observed"],
          flush=True)
    return dict(data, status="EXECUTED"), fits_rows, event_rows


def collect(results):
    inventory, fits, events = [], [], []
    for target, trial, measurements in results:
        inventory.append(target)
        fits.extend(trial)
        events.extend(measurements)
    pd.DataFrame(inventory).to_csv(out / "coverage.csv", index=False)
    pd.DataFrame(fits).to_csv(out / "candidates.csv", index=False)
    pd.DataFrame(events).to_csv(out / "events.csv", index=False)
    write_json(out / "summary.json", dict(executed_utc=utcnow(), hosts=len(inventory), peaks=len(fits),
        script_sha256=text_digest(__file__), config_sha256=text_digest(config_path),
        preparation_sha256=text_digest(root / "src/cosmic_anomaly_lab/companion_search.py"),
        products_sha256=text_digest(out / "products.json"), selection_sha256=text_digest(out / "selection.json"),
        heldout_ephemeris_refitted=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--wait-for-downloads", action="store_true")
    args = parser.parse_args()
    if args.workers == 1:
        collect(map(run_target, tasks(args.wait_for_downloads)))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            collect(executor.map(run_target, tasks(args.wait_for_downloads)))

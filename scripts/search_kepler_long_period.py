"""Fit early-quarter transit ephemerides and measure them in untouched quarters."""

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

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths, measure_events
from cosmic_anomaly_lab.transits import flatten_segments, phase_distance
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_long_period"
config_path = root / "configs/kepler_long_period.yaml"
config = yaml.safe_load(config_path.read_text())
roster = json.loads((out / "selection.json").read_text())["roster"]


def completed_targets(wait_for_downloads):
    for index, target in enumerate(roster):
        deadline = time.monotonic() + (360 if wait_for_downloads else 0)
        while True:
            inventory = json.loads((out / "products.json").read_text())["targets"]
            if len(inventory) > index:
                ready = inventory[index]
                assert (ready["kic"], ready["role"]) == (target["kic"], target["role"])
                yield ready
                break
            if time.monotonic() >= deadline:
                raise ValueError(f"Download inventory is missing KIC {target['kic']}")
            time.sleep(2)


def load_curve(record, kic):
    path = root / record["path"]
    assert digest(path) == record["sha256"]
    with fits.open(path) as product:
        header, timing, data = product[0].header, product[1].header, product[1].data
        assert header["KEPLERID"] == kic
        assert timing["BJDREFI"] == 2454833 and timing.get("BJDREFF", 0) == 0
        assert timing["TIMESYS"] == "TDB" and timing["TIMEPIXR"] == .5
        cadence = float(timing["TIMEDEL"])
        assert .019 < cadence < .023
        idx, flux, err = flatten_segments(data["TIME"], data["PDCSAP_FLUX"],
                                         data["PDCSAP_FLUX_ERR"], data["SAP_QUALITY"],
                                         window_days=4)
        time = np.asarray(data["TIME"][idx], float)
        return dict(quarter=int(header["QUARTER"]), cadence=cadence, time=time,
                    flux=flux, error=err, file=record)


def trial_grid(model, durations):
    grid = model.autoperiod(durations, minimum_period=config["period_days"][0],
                           maximum_period=config["period_days"][1], frequency_factor=3)
    return grid[(grid >= config["period_days"][0]) & (grid <= config["period_days"][1])]


def search(curves, sign=1):
    time = np.concatenate([c["time"] for c in curves])
    flux = np.concatenate([c["flux"] for c in curves])
    errors = []
    for curve in curves:
        differences = np.diff(curve["flux"])
        floor = 1.4826 * np.median(abs(differences - np.median(differences))) / np.sqrt(2)
        errors.append(np.maximum(curve["error"], floor))
    effective_error = np.concatenate(errors)
    model = BoxLeastSquares(time, 1 + sign*(flux - 1), effective_error)
    durations = config["durations_days"]
    span = np.ptp(time)
    periods = trial_grid(model, durations)
    periodogram = model.power(periods, durations, objective="snr", oversample=5)
    idx = int(np.argmax(periodogram.power))
    frequency = 1 / periodogram.period[idx]
    step = 3 * min(durations) / span**2
    local = 1 / np.linspace(frequency - 2*step, frequency + 2*step, 101)
    lo, hi = config["period_days"]
    local = local[(local >= lo) & (local <= hi)]
    fit = model.power(local, durations, objective="snr", oversample=20)
    ix = int(np.argmax(fit.power))
    return dict(period_days=float(fit.period[ix]), epoch_bkjd=float(fit.transit_time[ix]),
                duration_days=float(fit.duration[ix]), search_contrast=float(fit.power[ix]),
                grid_periods=len(periodogram.period))


def events_at(curves, fit, sign=1, phase=0):
    period, epoch, duration = (fit[k] for k in ["period_days", "epoch_bkjd", "duration_days"])
    return [dict(quarter=curve["quarter"], **ev) for curve in curves
            for ev in measure_events(curve["time"], 1 + sign*(curve["flux"] - 1), curve["error"],
                                     period, epoch + phase*period, duration, curve["cadence"])]


def gate(summary):
    return summary["events"] >= config["minimum_events_per_half"] and (
        summary["contrast"] >= config["minimum_contrast_per_half"]
    )


def plot_candidate(kic, role, fit, curves, event_rows):
    fig = Figure(figsize=(11, 6.5), layout="constrained")
    panels = fig.subplots(2, 2)
    for i, (split, values) in enumerate(curves.items()):
        ax, event_ax = panels[i]
        for curve in values:
            phase = phase_distance(curve["time"], fit["period_days"], fit["epoch_bkjd"])
            use = abs(phase) <= 3 * fit["duration_days"]
            ax.plot(phase[use]*24, (curve["flux"][use]-1)*1e6, ".", ms=3, alpha=.38)
        ax.axvspan(-12*fit["duration_days"], 12*fit["duration_days"], color="#d58a38", alpha=.12)
        ax.set(title=f"{split} · folded observed light curve",
               xlabel="Hours from predicted transit", ylabel="Relative flux (ppm)")
        valid = [r for r in event_rows[split] if r["status"] == "MEASURED"]
        event_ax.errorbar([r["midpoint"] for r in valid], [r["depth"]*1e6 for r in valid],
                          yerr=[r["error"]*1e6 for r in valid], fmt="o", ms=3, color="#254d68")
        event_ax.axhline(0, color="#777", lw=.7)
        event_ax.set(title=f"{split} · signed event depths", xlabel="Event midpoint (BKJD)", ylabel="Depth (ppm)")
    fig.suptitle(f"KIC {kic} · {role.replace('_', ' ')} · {fit['period_days']:.5f} days")
    save_figure(fig, root / f"figures/kepler_long_period/kic{kic}")


def search_target(target):
    fits_rows, event_rows, snapshot = [], [], {}
    kic, role = target["kic"], target["role"]
    curves = sorted([load_curve(p, kic) for p in target["products"]], key=lambda curve: curve["quarter"])
    quarters = [c["quarter"] for c in curves]
    assert len(set(quarters)) == len(quarters)
    halves = {label: [c for c in curves if lo <= c["quarter"] <= hi and len(c["time"]) > 100]
              for label, (lo, hi) in [("discovery", config["discovery_quarters"]),
                                      ("replication", config["replication_quarters"])]}
    snapshot.update({p["path"]: p["sha256"] for p in target["products"]})
    span = {key: float(values[-1]["time"][-1] - values[0]["time"][0]) if values else 0
            for key, values in halves.items()}
    record = dict(kic=kic, role=role, downloaded_quarters=len(curves),
                  discovery_quarters=len(halves["discovery"]), replication_quarters=len(halves["replication"]),
                  discovery_days=span["discovery"], replication_days=span["replication"],
                  retained_cadences=sum(len(c["time"]) for half in halves.values() for c in half))
    if min(span.values()) < 100:
        fits_rows.append(dict(**record, status="INSUFFICIENT_BASELINE", primary_replicated=False))
        return fits_rows, event_rows, snapshot
    assert halves["discovery"][-1]["time"][-1] < halves["replication"][0]["time"][0]
    for sign, name in [(1, "observed"), (-1, "inverted_control")]:
        fit = search(halves["discovery"], sign)
        rows = {label: events_at(half, fit, sign) for label, half in halves.items()}
        summaries = {label: combine_event_depths(events) for label, events in rows.items()}
        replicated = all(gate(s) for s in summaries.values())
        ancillary = {}
        if sign == 1:
            for phase, label in [(.5, "half_phase"), (.25, "quarter_phase"), (.75, "three_quarter_phase")]:
                trial = events_at(halves["replication"], fit, phase=phase)
                ancillary.update({label + "_" + k: v for k, v in combine_event_depths(trial).items()})
                event_rows.extend(dict(kic=kic, sign=name, split="replication", phase=phase, **row)
                                  for row in trial)
        fits_rows.append(dict(**record, sign=name, status="EXECUTED", **fit,
                              primary_replicated=replicated, **ancillary,
                              **{label + "_" + k: v for label, s in summaries.items() for k, v in s.items()}))
        event_rows.extend(dict(kic=kic, sign=name, split=label, phase=0., **row)
                          for label, group in rows.items() for row in group)
        if sign == 1 and (replicated or role == "confirmed_control"):
            plot_candidate(kic, role, fit, halves, rows)
        print(kic, name, f"{fit['period_days']:.5f}", replicated,
              *(summaries[label]["contrast"] for label in halves), flush=True)
    return fits_rows, event_rows, snapshot


def collect(results):
    fits_rows, event_rows, snapshot = [], [], {}
    for rows, events, hashes in results:
        fits_rows.extend(rows)
        event_rows.extend(events)
        snapshot.update(hashes)
    pd.DataFrame(fits_rows).to_csv(out / "search.csv", index=False)
    pd.DataFrame(event_rows).to_csv(out / "events.csv", index=False)
    write_json(out / "search_summary.json", dict(executed_utc=utcnow(), total_targets=len(roster),
        partial_run=False, catalogues_used_to_fit=False, replication_parameters_frozen=True,
        input_sha256=snapshot, script_sha256=text_digest(Path(__file__)),
        config_sha256=text_digest(config_path),
        transit_followup_sha256=text_digest(root / "src/cosmic_anomaly_lab/transit_followup.py"),
        detrending_sha256=text_digest(root / "src/cosmic_anomaly_lab/transits.py")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--wait-for-downloads", action="store_true",
                        help="Read each completed target while the single cache writer downloads the rest")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.workers == 1:
        collect(map(search_target, completed_targets(args.wait_for_downloads)))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            collect(executor.map(search_target, completed_targets(args.wait_for_downloads)))

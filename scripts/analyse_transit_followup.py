"""Catalogue-directed Kepler audit and locked-ephemeris TESS replication."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths, measure_events
from cosmic_anomaly_lab.transits import flatten_segments, phase_distance
from cosmic_anomaly_lab.visual_evidence import save_figure

parser = argparse.ArgumentParser()
parser.add_argument("mission", choices=["kepler", "tess"])
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
out = root / "results" / (args.mission + "_followup")
products = json.loads((out / "products.json").read_text())
input_hashes = {}


def checked_path(record):
    path = root / record["path"]
    if path not in input_hashes:
        actual = digest(path)
        assert record["sha256"] == actual
        input_hashes[path] = actual
    return path


def lightcurve(record, column, detrend=False):
    with fits.open(checked_path(record)) as hdus:
        header, timing = hdus[0].header, hdus[1].header
        data = hdus[1].data
        time, flux, error = [np.asarray(data[k], dtype=float) for k in ["TIME", column, column + "_ERR"]]
        if args.mission == "tess":
            assert header["TICID"] == record["tic"] and header["SECTOR"] == record["sector"]
            quality = data["QUALITY"]
            segment = header["SECTOR"]
            reference = 2457000
        else:
            assert header["KEPLERID"] == 11446424
            quality = data["SAP_QUALITY"]
            segment = header["QUARTER"]
            reference = 2454833
        assert timing["BJDREFI"] == reference and timing.get("BJDREFF", 0) == 0
        assert timing["TIMESYS"] == "TDB" and timing.get("TIMEPIXR", .5) == .5
        cadence = float(timing["TIMEDEL"])
        assert 0 < cadence < .03
    if detrend:
        idx, flux, error = flatten_segments(time, flux, error, quality)
        time = time[idx]
    else:
        good = (quality == 0) & np.isfinite(time + flux + error) & (flux > 0) & (error > 0)
        time, flux, error = time[good], flux[good], error[good]
        scale = np.median(flux)
        flux, error = flux / scale, error / scale
    assert len(time) > 100 and np.all(np.diff(time) > 0)
    return time, flux, error, cadence, int(segment)


def at_ephemeris(curve, period, epoch, duration):
    time, flux, error, cadence, segment = curve
    return [dict(segment=segment, **row) for row in
            measure_events(time, flux, error, period, epoch, duration, cadence)]


if args.mission == "kepler":
    catalogue_path = root / "data/raw/exoplanets/koi_cumulative.json"
    old = json.loads((root / "results/exoplanet_pilot/summary.json").read_text())
    assert digest(catalogue_path) == old["koi_sha256"]
    catalogue = pd.read_json(catalogue_path)
    row = catalogue[catalogue.kepoi_name == "K08223.01"].iloc[0]
    assert row.kepid == 11446424
    period, epoch, duration = float(row.koi_period), float(row.koi_time0bk), float(row.koi_duration / 24)
    all_events, summaries, nulls, omissions = [], {}, [], []
    fig = Figure(figsize=(12, 7), layout="constrained")
    axes = fig.subplots(2, 2)
    for j, column in enumerate(["PDCSAP_FLUX", "SAP_FLUX"]):
        curves = [lightcurve(p, column) for p in products["products"]]
        events = [dict(reduction=column, **e) for curve in curves
                  for e in at_ephemeris(curve, period, epoch, duration)]
        all_events.extend(events)
        later = [e for e in events if e["segment"] > 3]
        summaries[column] = dict(all=combine_event_depths(events), q4_and_later=combine_event_depths(later),
                                q2_q3=combine_event_depths([e for e in events if e["segment"] <= 3]),
                                quarters=sorted(x[-1] for x in curves), cadences=sum(len(x[0]) for x in curves))
        for quarter in sorted({e["segment"] for e in later if e["status"] == "MEASURED"}):
            omissions.append(dict(reduction=column, omitted_quarter=quarter,
                                  **combine_event_depths([e for e in later if e["segment"] != quarter])))
        for fraction in np.linspace(.02, .98, 193):
            shifted = [e for curve in curves if curve[-1] > 3 for e in
                       at_ephemeris(curve, period, epoch + period * fraction, duration)]
            nulls.append(dict(reduction=column, phase=float(fraction), **combine_event_depths(shifted)))
        good = [e for e in events if e["status"] == "MEASURED"]
        ax = axes[0, j]
        for subset, color, name in [(False, "#66789a", "Kepler Q2–Q3"), (True, "#c2702f", "Kepler Q4 onward")]:
            points = [e for e in good if (e["segment"] > 3) == subset]
            ax.errorbar([e["midpoint"] for e in points], [e["depth"]*1e6 for e in points],
                        yerr=[e["error"]*1e6 for e in points], fmt="o", ms=3, color=color, label=name)
        ax.axhline(0, color="gray", lw=.7)
        ax.set(title=column, xlabel="Predicted event midpoint (BKJD)", ylabel="Signed fitted depth (ppm)")
        ax.legend(fontsize=8)
        info = summaries[column]["q4_and_later"]
        ax = axes[1, j]
        comparable = [e["contrast"] for e in nulls if e["reduction"] == column and e["events"] >= info["events"]]
        ax.hist(comparable, bins=np.linspace(-6, 8, 29), color="#617b9a")
        if info["contrast"] is not None:
            ax.axvline(info["contrast"], color="#b55f25", label="Catalogue phase")
        ax.set(xlabel="Conditional depth / error", ylabel="Phase offsets", title="Q4 onward, fixed-period phase comparison")
        ax.legend(fontsize=8)
        info["phase_offsets_with_at_least_as_many_events"] = len(comparable)
        info["phase_offsets_at_least_as_deep_in_contrast"] = (
            sum(c >= info["contrast"] for c in comparable) if info["contrast"] is not None else None)
    fig.suptitle("K08223.01 · catalogue-directed check of actual Kepler observations")
    save_figure(fig, root / "figures/kepler_followup/k08223_event_audit")
    pd.DataFrame(all_events).to_csv(out / "events.csv", index=False)
    pd.DataFrame(nulls).to_csv(out / "phase_offsets.csv", index=False)
    pd.DataFrame(omissions).to_csv(out / "quarter_omissions.csv", index=False)
    result = dict(candidate="K08223.01", catalogue_disposition=row.koi_disposition,
                  period_days=period, epoch_bkjd=epoch, duration_days=duration,
                  catalogue_sha256=old["koi_sha256"], reductions=summaries,
                  independent_of_catalogue=False, planet_validated=False)
else:
    toi, planets = (json.loads(checked_path(products["ancillary"][k]).read_text()) for k in ["toi", "planets"])
    rows, events = [], []
    for sic in sorted({p["tic"] for p in products["products"]}):
        first, second = sorted([p for p in products["products"] if p["tic"] == sic], key=lambda p: p["sector"])
        train, hold = [lightcurve(p, "PDCSAP_FLUX", detrend=True) for p in [first, second]]
        time, flux, error = train[:3]
        scatter = 1.4826 * np.median(abs(flux - np.median(flux)))
        durations = [.05, .08, .12, .18]
        bls = BoxLeastSquares(time, flux, np.maximum(error, scatter))
        grid = bls.autopower(durations, minimum_period=.8, maximum_period=8, frequency_factor=2,
                             objective="snr")
        index = int(np.argmax(grid.power))
        period, epoch, duration = (float(grid[k][index]) for k in ["period", "transit_time", "duration"])
        refine = bls.power(np.linspace(period - 2 * period**2 * .05 / np.ptp(time)**2,
                                       period + 2 * period**2 * .05 / np.ptp(time)**2, 121),
                           durations, objective="snr", oversample=30)
        ix = int(np.argmax(refine.power))
        period, epoch, duration = (float(refine[k][ix]) for k in ["period", "transit_time", "duration"])
        a, b = [at_ephemeris(curve, period, epoch, duration) for curve in [train, hold]]
        events.extend(dict(tic=sic, **ev) for ev in a + b)
        ra, rb = [combine_event_depths(group) for group in [a, b]]
        on_star = [x for x in toi if int(x["tid"]) == sic]
        hosts = [x for x in planets if x["tic_id"] == f"TIC {sic}"]
        passed = all(x["events"] >= 3 and x["contrast"] >= 5 for x in [ra, rb])
        rows.append(dict(tic=sic, role=first["role"], period_days=period, epoch_btjd=epoch,
                         duration_days=duration, training_bls_snr=float(refine.power[ix]),
                         training_points=len(train[0]), holdout_points=len(hold[0]),
                         **{"s2_" + k: v for k, v in ra.items()}, **{"s3_" + k: v for k, v in rb.items()},
                         preliminary_replication=passed, toi_entries=json.dumps(on_star),
                         confirmed_on_star=json.dumps(hosts)))
        print(sic, period, ra["contrast"], rb["contrast"], passed, flush=True)
        if passed or first["role"] == "control":
            fig = Figure(figsize=(10, 4.5), layout="constrained")
            for ax, curve in zip(fig.subplots(1, 2), [train, hold], strict=True):
                time, flux, _, _, sector = curve
                phase = phase_distance(time, period, epoch)
                use = abs(phase) <= 3 * duration
                ax.plot(phase[use] * 24, (flux[use] - 1)*1e6, ".", ms=2, alpha=.35)
                ax.axvspan(-duration*12, duration*12, color="orange", alpha=.12)
                ax.set(title=f"Sector {sector}" + (" · search" if sector == 2 else " · fixed ephemeris"),
                       xlabel="Hours from predicted midpoint", ylabel="Flux relative to baseline (ppm)")
            fig.suptitle(f"Observed TESS light curves · TIC {sic} · {period:.6f} d")
            save_figure(fig, root / f"figures/tess_followup/tic{sic}")
    pd.DataFrame(rows).to_csv(out / "replication.csv", index=False)
    pd.DataFrame(events).to_csv(out / "events.csv", index=False)
    result = dict(targets=len(rows), total_cadences=sum(r["training_points"] + r["holdout_points"] for r in rows),
                  search_replications=sum(r["preliminary_replication"] for r in rows if r["role"] == "search"),
                  control_replication=next(r["preliminary_replication"] for r in rows if r["role"] == "control"),
                  holdout_searched=False, planet_validated=False)
write_json(out / "summary.json", dict(executed_utc=utcnow(), **result,
    script_sha256=text_digest(Path(__file__)),
    event_model_sha256=text_digest(root / "src/cosmic_anomaly_lab/transit_followup.py"),
    input_sha256={p.relative_to(root).as_posix(): checksum for p, checksum in sorted(input_hashes.items())}))
print(json.dumps(result, indent=2), flush=True)

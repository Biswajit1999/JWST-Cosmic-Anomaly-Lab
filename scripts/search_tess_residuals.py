"""Search for an additional repeating dip after accounting for stellar variability."""

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

root = Path(__file__).resolve().parents[1]
old = root / "results/tess_followup"
out = root / "results/tess_residual_search"
out.mkdir(parents=True, exist_ok=True)
products = json.loads((old / "products.json").read_text())["products"]
replications = pd.read_csv(old / "replication.csv").query("preliminary_replication and role == 'search'")
known_models = pd.read_csv(old / "heldout_models.csv")
rows, events, inputs = [], [], {}
fig = Figure(figsize=(10, 8), layout="constrained")
axes = fig.subplots(len(replications), 2, squeeze=False)
for i, star in enumerate(replications.itertuples()):
    options = known_models[(known_models.tic == star.tic) & known_models.model.str.startswith("harmonics")]
    # Use training RMS alone to choose the nuisance frequency; remove the same
    # frequency basis in both sectors, allowing coefficients to differ.
    model = options.loc[options.training_rms_ppm.idxmin(), "model"]
    variable_period = star.period_days * (2 if model == "harmonics_2P" else 1)
    curves, originals = [], []
    for product in sorted((p for p in products if p["tic"] == star.tic), key=lambda p: p["sector"]):
        path = root / product["path"]
        assert digest(path) == product["sha256"]
        inputs[product["path"]] = product["sha256"]
        with fits.open(path) as hdus:
            header, timing, data = hdus[0].header, hdus[1].header, hdus[1].data
            assert header["TICID"] == star.tic and header["SECTOR"] == product["sector"]
            assert timing["BJDREFI"] == 2457000 and timing["TIMESYS"] == "TDB"
            idx, flux, errors = flatten_segments(data["TIME"], data["PDCSAP_FLUX"],
                                                 data["PDCSAP_FLUX_ERR"], data["QUALITY"])
            time = np.asarray(data["TIME"][idx], float)
            cadence = float(timing["TIMEDEL"])
        phase = 2*np.pi * (time - star.epoch_btjd) / variable_period
        design = np.column_stack([np.ones(len(time)), np.sin(phase), np.cos(phase),
                                  np.sin(2*phase), np.cos(2*phase)])
        coefficients = np.linalg.lstsq(design, flux, rcond=None)[0]
        residual = flux - design @ coefficients
        floor = 1.4826 * np.median(abs(residual - np.median(residual)))
        curves.append((time, residual + 1, np.maximum(errors, floor), cadence, product["sector"]))
        originals.append(float(np.std(flux)))
    for sign, label in [(1, "observed_residual"), (-1, "inverted_residual")]:
        time, flux, error, _, _ = curves[0]
        bls = BoxLeastSquares(time, 1 + sign*(flux - 1), error)
        durations = [.05, .08, .12, .18]
        period_grid = bls.autoperiod(durations, minimum_period=.8, maximum_period=8, frequency_factor=2)
        aliases = {a/b for a in range(1, 5) for b in range(1, 5)}
        use = (period_grid >= .8) & (period_grid <= 8)
        for harmonic in aliases:
            use &= abs(1/period_grid - harmonic/variable_period) > 1/np.ptp(time)
        if not use.any():
            rows.append(dict(tic=star.tic, sign=label, status="NO_UNMASKED_PERIODS", replicated=False))
            continue
        fit = bls.power(period_grid[use], durations, objective="snr")
        ix = int(np.argmax(fit.power))
        period, epoch, duration = (float(fit[k][ix]) for k in ["period", "transit_time", "duration"])
        summaries = []
        for j, (t, y, e, cadence, sector) in enumerate(curves):
            group = list(measure_events(t, 1 + sign*(y-1), e, period, epoch, duration, cadence))
            events.extend(dict(tic=star.tic, sign=label, sector=sector, **ev) for ev in group)
            summaries.append(combine_event_depths(group))
            if sign == 1:
                ax = axes[i, j]
                v = phase_distance(t, period, epoch)
                inside = abs(v) <= 3*duration
                ax.plot(v[inside]*24, (y[inside]-1)*1e6, ".", ms=2, color="#315676", alpha=.38)
                ax.axvspan(-12*duration, 12*duration, color="#db9443", alpha=.15)
                ax.set(title=f"TIC {star.tic} · Sector {sector} · {period:.4f} d",
                       xlabel="Hours from predicted midpoint", ylabel="Residual flux (ppm)")
        passed = all(s["events"] >= 3 and s["contrast"] >= 5 for s in summaries)
        rows.append(dict(tic=star.tic, sign=label, status="EXECUTED", replicated=passed,
                         nuisance_period_days=variable_period, period_days=period, epoch_btjd=epoch,
                         duration_days=duration, tested_periods=int(use.sum()), bls_contrast=float(fit.power[ix]),
                         **{f"s{j+2}_raw_rms_ppm": rms*1e6 for j, rms in enumerate(originals)},
                         **{f"s{j+2}_residual_rms_ppm": float(np.std(c[1])*1e6) for j, c in enumerate(curves)},
                         **{f"s{j+2}_{key}": value for j, summary in enumerate(summaries)
                            for key, value in summary.items()}))
fig.suptitle("Observed TESS residuals after fitting the dominant variability\nSector 2 dip search, same candidate timing in Sector 3", fontsize=12)
save_figure(fig, root / "figures/tess_residual_search/phase_comparison")
frame = pd.DataFrame(rows)
frame.to_csv(out / "search.csv", index=False)
pd.DataFrame(events).to_csv(out / "events.csv", index=False)
write_json(out / "summary.json", dict(executed_utc=utcnow(), input_sha256=inputs,
    script_sha256=text_digest(Path(__file__)), primary_sample_selection="Two prior TESS variable lookalikes",
    holdout_transit_timing_refitted=False, fitted_nuisance_in_each_sector=True,
    new_planets_validated=0, candidates_tested=int(len(replications))))
print(frame[["tic", "sign", "nuisance_period_days", "period_days",
             "s2_contrast", "s3_contrast", "replicated"]].to_string(index=False))

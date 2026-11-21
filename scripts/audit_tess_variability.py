"""Compare periodic box and harmonic predictions on the same held-out samples."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transits import flatten_segments, phase_distance
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/tess_followup"
products = json.loads((out / "products.json").read_text())["products"]
ranking = pd.read_csv(out / "replication.csv")
records = []
for row in ranking[ranking.preliminary_replication].itertuples():
    curves = []
    for p in sorted([p for p in products if p["tic"] == row.tic], key=lambda p: p["sector"]):
        path = root / p["path"]
        assert digest(path) == p["sha256"]
        with fits.open(path) as h:
            data = h[1].data
            idx, flux, _ = flatten_segments(data["TIME"], data["PDCSAP_FLUX"],
                                             data["PDCSAP_FLUX_ERR"], data["QUALITY"])
            time = np.asarray(data["TIME"][idx], float)
            curves.append((time, flux))
    period, epoch, duration = row.period_days, row.epoch_btjd, row.duration_days

    def design(time, model):
        phase = (time - epoch) / period
        unit = np.ones(len(time))
        if model == "box":
            return np.column_stack([unit, abs(phase_distance(time, period, epoch)) <= duration/2])
        base = phase / 2 if model == "harmonics_2P" else phase
        return np.column_stack([unit, np.sin(2*np.pi*base), np.cos(2*np.pi*base),
                                np.sin(4*np.pi*base), np.cos(4*np.pi*base)])

    train, holdout = curves
    fig = Figure(figsize=(10, 4.8), layout="constrained")
    ax = fig.subplots()
    phase = (holdout[0] - epoch) / period % 2
    ax.plot(phase, (holdout[1] - 1)*1e6, ".", ms=2, alpha=.22, color="#334d6b",
            label="Sector 3 observed")
    grid = epoch + period*np.linspace(0, 2, 1801)
    for model, color in zip(["box", "harmonics_P", "harmonics_2P"],
                            ["#c47234", "#328568", "#7d65a6"], strict=True):
        x = design(train[0], model)
        beta = np.linalg.lstsq(x, train[1], rcond=None)[0]
        prediction = design(holdout[0], model) @ beta
        residual = holdout[1] - prediction
        offset = float(np.mean(residual))
        records.append(dict(tic=row.tic, role=row.role, model=model, parameters=len(beta),
                            training_rms_ppm=float(np.std(train[1] - x @ beta)*1e6),
                            holdout_rms_ppm=float(np.std(residual)*1e6), holdout_offset_ppm=offset*1e6))
        ax.plot((grid - epoch) / period, (design(grid, model) @ beta + offset - 1)*1e6,
                color=color, lw=1.5, label=model.replace("_", " "))
    ax.set(title=f"TIC {row.tic} · predictions from Sector 2, timing and shapes fixed",
           xlabel="Cycles of the fitted box period (two cycles shown)",
           ylabel="Relative flux (ppm)", xlim=(0, 2))
    ax.legend(fontsize=9)
    save_figure(fig, root / f"figures/tess_followup/tic{row.tic}_model_check")
data = pd.DataFrame(records)
data.to_csv(out / "heldout_models.csv", index=False)
write_json(out / "model_audit.json", dict(executed_utc=utcnow(), script_sha256=text_digest(Path(__file__)),
    searched_heldout_period=False, fitted_heldout_offset_only=True,
    comparison="Unweighted residual RMS on same held-out cadences; template diagnostics, not class inference",
    ranking_sha256=text_digest(out / "replication.csv")))
print(data.to_string(index=False))

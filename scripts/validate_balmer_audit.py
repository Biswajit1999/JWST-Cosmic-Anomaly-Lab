"""Finite synthetic controls, including a deliberately incorrect instrumental LSF."""

import json
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_absorption import BalmerModel, fit_balmer
from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.validation import binomial_interval
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/balmer_audit"
sources = json.loads((out / "inputs.json").read_text())
fits = json.loads((out / "fits.json").read_text())
template = next(
    f
    for f in fits
    if f["label"] == "DJA G395H" and f["family"] == "absorbed" and f["resolution_factor"] == 1
)
rng = np.random.default_rng(20261115)
trials, summaries = [], []
for label, true_resolution, fitted_resolution in [
    ("DJA G395H", 2700, 2700),
    ("DJA G395M", 1000, 1000),
    ("DJA G395M", 1500, 1000),
]:
    s = next(s for s in sources if s["label"] == label)
    v, error, good = (
        np.array(s["velocity_kms"]),
        np.array(s["error_uJy"], float),
        np.array(s["good"]),
    )
    for truth in ["co_centered", "absorbed"]:
        p = (
            template["parameter_vector"][:11]
            if truth == "co_centered"
            else template["parameter_vector"]
        )
        noiseless = BalmerModel(v, true_resolution, truth).predict(p)
        records = []
        for trial in range(10):
            flux = noiseless + rng.normal(size=len(v)) * np.where(good, error, 0)
            pair = [
                fit_balmer(v, flux, error, good, fitted_resolution, family, starts=2)
                for family in ["co_centered", "absorbed"]
            ]
            row = dict(
                spectrum=label,
                truth=truth,
                true_resolution=true_resolution,
                fitted_resolution=fitted_resolution,
                trial=trial,
                delta_bic=pair[0]["bic"] - pair[1]["bic"],
                absorbed_bound_parameters=pair[1]["parameters_near_bounds"],
            )
            trials.append(row)
            records.append(row)
        count = sum(r["delta_bic"] > 10 for r in records)
        summaries.append(
            dict(
                spectrum=label,
                truth=truth,
                true_resolution=true_resolution,
                fitted_resolution=fitted_resolution,
                trials=len(records),
                delta_bic_over_10=count,
                binomial_interval_95=binomial_interval(count, len(records)),
            )
        )
        write_json(
            out / "synthetic_validation.json",
            dict(
                data_kind="SIMULATED",
                seed=20261115,
                executed_utc=utcnow(),
                completed_scenarios=len(summaries),
                summaries=summaries,
                trials=trials,
                template_parameters=template["parameter_vector"],
                limitation="Source-inspired local template, independent Gaussian pixel errors, only ten trials per cell; not a survey false-positive rate or population completeness",
                code_sha256=text_digest(root / "src/cosmic_anomaly_lab/balmer_absorption.py"),
                script_sha256=text_digest(Path(__file__)),
            ),
        )
        print(
            label, truth, true_resolution, fitted_resolution, count, "/", len(records), flush=True
        )
fig = Figure(figsize=(11, 5), layout="constrained")
ax = fig.subplots()
for i, summary in enumerate(summaries):
    values = [
        r["delta_bic"]
        for r in trials
        if all(
            r[k] == summary[k]
            for k in ["spectrum", "truth", "true_resolution", "fitted_resolution"]
        )
    ]
    ax.scatter(
        np.full(len(values), i) + np.linspace(-0.12, 0.12, len(values)),
        values,
        color="#bb5932" if summary["truth"] == "absorbed" else "#377a93",
        s=22,
    )
ax.axhline(10, color="gray", ls="--", label="Descriptive ΔBIC = 10 reference")
ax.set_xticks(
    range(len(summaries)),
    [
        f"{s['truth'].replace('_', ' ')}\nR true {s['true_resolution']}\nR fit {s['fitted_resolution']}"
        for s in summaries
    ],
    fontsize=8,
)
ax.set(
    ylabel="BIC(emission only) − BIC(absorbed)",
    title="SIMULATED controls · ten trials per cell · not survey detection calibration",
)
ax.legend(fontsize=8)
ax.grid(axis="y", alpha=0.15)
save_figure(fig, root / "figures/balmer_synthetic_controls")

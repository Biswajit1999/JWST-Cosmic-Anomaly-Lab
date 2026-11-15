"""Run bounded local emission/absorption comparisons for a known-source benchmark."""

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_absorption import BalmerModel, fit_balmer
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
config = yaml.safe_load((root / "configs/balmer_audit.yaml").read_text())
out = root / "results/balmer_audit"
out.mkdir(parents=True, exist_ok=True)
inputs = []
for grating in ["G395H", "G395M"]:
    path = root / f"data/raw/spectra/jades-gds06-v4_{grating.lower()}-f290lp_1286_159717.spec.fits"
    inputs.append((f"DJA {grating}", grating, path, read_spectrum(path)))
dr4 = (
    root
    / "data/raw/independent/hlsp_jades_jwst_nirspec_goods-s-mediumjwst-00159717_f290lp-g395m_v1.0_x1d.fits"
)
for aperture in [3, 5]:
    for nods in [3, 2]:
        spec = read_dr4_spectrum(dr4, aperture, nods)
        reference = inputs[0][3].metadata
        separation = (
            SkyCoord(spec.metadata["RA"] * u.deg, spec.metadata["DEC"] * u.deg)
            .separation(SkyCoord(reference["SRCRA"] * u.deg, reference["SRCDEC"] * u.deg))
            .arcsec
        )
        assert separation < 0.01 and int(spec.metadata["PROPOSID"]) == 1286
        assert int(spec.metadata["TARGET"]) == 159717
        inputs.append((f"DR4 G395M {aperture}px/{nods}-nod", "G395M", dr4, spec))

rows, full, source_records = [], [], []
fig = Figure(figsize=(14, 13), layout="constrained")
axes = fig.subplots(3, 2).flat
colors = {"co_centered": "#4682a0", "absorbed": "#b74a27", "split_emission": "#7c6197"}
for (label, grating, path, spec), ax in zip(inputs, axes, strict=True):
    velocity = (
        spec.wavelength.value / (config["halpha_vacuum_um"] * (1 + config["redshift_reference"]))
        - 1
    ) * 299792.458
    window = abs(velocity) <= config["velocity_limit_kms"]
    v, flux, error, good = (
        velocity[window],
        spec.flux.value[window],
        spec.error.value[window],
        spec.good[window],
    )
    source_records.append(
        dict(
            label=label,
            path=str(path.relative_to(root)).replace("\\", "/"),
            sha256=digest(path),
            metadata=spec.metadata,
            velocity_kms=v.tolist(),
            flux_uJy=[float(f) if np.isfinite(f) else None for f in flux],
            error_uJy=[float(e) if np.isfinite(e) else None for e in error],
            good=good.tolist(),
        )
    )
    ax.errorbar(
        v[good], flux[good], yerr=error[good], fmt=".", ms=3, color="black", alpha=0.65, lw=0.6
    )
    for factor in config["resolution_factors"]:
        resolution = config["resolution"][grating] * factor
        for family in config["families"]:
            fit = fit_balmer(
                v,
                flux,
                error,
                good,
                resolution,
                family,
                starts=config["fitting_starts"],
                step=config["grid_step_kms"],
            )
            record = dict(label=label, resolution_factor=factor, **fit)
            full.append(record)
            rows.append(
                {
                    k: v
                    for k, v in record.items()
                    if k not in ["model", "fit_mask", "parameter_vector", "attempts", "parameters"]
                }
                | {
                    "tau": fit["parameters"].get("tau"),
                    "abs_sigma": fit["parameters"].get("abs_sigma"),
                    "abs_velocity": fit["parameters"].get("abs_velocity"),
                }
            )
            if factor == 1:
                ax.plot(
                    v, fit["model"], lw=1.2, color=colors[family], label=family.replace("_", " ")
                )
            print(
                label,
                factor,
                family,
                round(fit["chi2"], 2),
                fit["parameters_near_bounds"],
                flush=True,
            )
        # Save progress after each completed set; never confuse partial run with summary.
        write_json(out / "fits.json", full)
    ax.set(
        xlim=(-2200, 2200),
        xlabel="Velocity at z=5.07781 (km/s)",
        ylabel="Flux density (µJy)",
        title=label,
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.15)
fig.suptitle(
    "OBSERVED GS-159717 · local Hα emission / absorption model comparison\n"
    "Nominal Gaussian LSF; fixed covering fraction; targeted known-source benchmark",
    fontsize=13,
)
save_figure(fig, root / "figures/balmer_model_comparison")
pd.DataFrame(rows).to_csv(out / "fit_comparisons.csv", index=False)
write_json(out / "inputs.json", source_records)

# Inspect standardized residuals separately, rather than reporting only a ranking.
rf = Figure(figsize=(13, 8), layout="constrained")
rax = rf.subplots(2, 2)
for i, source in enumerate(source_records[:2]):
    v, flux, error = map(
        np.array, [source["velocity_kms"], source["flux_uJy"], source["error_uJy"]]
    )
    good = np.array(source["good"])
    for fit in full:
        if fit["label"] == source["label"] and fit["resolution_factor"] == 1:
            for j, limits in enumerate([(-4500, 4500), (-650, 650)]):
                rax[i, j].plot(
                    v[good],
                    ((flux - np.array(fit["model"])) / error)[good],
                    ".-",
                    ms=2,
                    lw=0.7,
                    color=colors[fit["family"]],
                    label=fit["family"].replace("_", " "),
                )
                rax[i, j].set(
                    xlim=limits,
                    xlabel="Velocity (km/s)",
                    ylabel="(data − model) / released error",
                    title=source["label"],
                )
                rax[i, j].axhline(0, color="gray", lw=0.6)
                rax[i, j].legend(fontsize=8)
rf.suptitle("OBSERVED residuals · diagonal statistical errors only; no class inference")
save_figure(rf, root / "figures/balmer_model_residuals")

# Numerical quadrature check on all nominal best fits, without refitting.
checks = []
for fit in full:
    if fit["resolution_factor"] != 1:
        continue
    source = next(s for s in source_records if s["label"] == fit["label"])
    fine = BalmerModel(source["velocity_kms"], fit["resolution"], fit["family"], step=2.5).predict(
        fit["parameter_vector"]
    )
    good = np.array(source["good"])
    error = np.array(source["error_uJy"], float)
    checks.append(
        dict(
            label=fit["label"],
            family=fit["family"],
            maximum_difference_over_error=float(
                np.max(abs(fine[good] - np.array(fit["model"])[good]) / error[good])
            ),
        )
    )
write_json(out / "numerical_convergence.json", checks)
write_json(
    out / "summary.json",
    dict(
        executed_utc=utcnow(),
        research_snapshot=config["research_snapshot"],
        source_id=config["source_id"],
        unique_source_count=1,
        extraction_variants=len(inputs),
        fits=len(full),
        physical_classification="NOT_ESTABLISHED",
        design="Post-selection known-source local-model benchmark",
        input_sha256={
            str(p.relative_to(root)).replace("\\", "/"): text_digest(p)
            for p in [
                root / "configs/balmer_audit.yaml",
                Path(__file__),
                root / "src/cosmic_anomaly_lab/balmer_absorption.py",
            ]
        },
    ),
)
print("Completed", len(full), "fits", flush=True)

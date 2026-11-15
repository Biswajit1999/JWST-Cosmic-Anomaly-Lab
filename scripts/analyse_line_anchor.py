"""Measure [OIII] independently of Halpha and retain the instrument/gas degeneracy."""

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.line_anchor import OIII, fit_oiii, resolution_lower_bound
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/line_anchor"
config = yaml.safe_load((root / "configs/line_anchor.yaml").read_text())
inputs = []
for grating, filt in [("G395H", "f290lp"), ("G235M", "f170lp"), ("G395M", "f290lp")]:
    path = root / f"data/raw/spectra/jades-gds06-v4_{grating.lower()}-{filt}_1286_159717.spec.fits"
    inputs.append((f"DJA {grating}", grating, path, read_spectrum(path)))
for grating, filt in [("G235M", "f170lp"), ("G395M", "f290lp")]:
    path = (
        root
        / f"data/raw/independent/hlsp_jades_jwst_nirspec_goods-s-mediumjwst-00159717_{filt}-{grating.lower()}_v1.0_x1d.fits"
    )
    for aperture in [3, 5]:
        for nods in [3, 2]:
            spec = read_dr4_spectrum(path, aperture, nods)
            ref = inputs[0][3].metadata
            separation = (
                SkyCoord(spec.metadata["RA"] * u.deg, spec.metadata["DEC"] * u.deg)
                .separation(SkyCoord(ref["SRCRA"] * u.deg, ref["SRCDEC"] * u.deg))
                .arcsec
            )
            assert separation < 0.01 and int(spec.metadata["TARGET"]) == 159717
            assert int(spec.metadata["PROPOSID"]) == 1286
            inputs.append((f"DR4 {grating} {aperture}px/{nods}-nod", grating, path, spec))
full, rows, source_records = [], [], []
fig = Figure(figsize=(13, 11), layout="constrained")
axes = iter(fig.subplots(3, 2).flat)
profile_fig = Figure(figsize=(12, 5), layout="constrained")
left, right = profile_fig.subplots(1, 2)
for label, grating, path, spec in inputs:
    source_records.append(
        dict(
            label=label,
            path=path.relative_to(root).as_posix(),
            sha256=digest(path),
            metadata=spec.metadata,
        )
    )
    ref_path = root / f"data/raw/instrument/{grating.lower()}_disp.fits"
    with fits.open(ref_path) as h:
        w, r = np.asarray(h[1].data["WAVELENGTH"]), np.asarray(h[1].data["R"])
        nominal = float(np.interp(OIII[1] * (1 + config["redshift"]), w, r))
        halpha_in_range = w.min() <= 0.656461 * (1 + config["redshift"]) <= w.max()
        nominal_halpha = (
            float(np.interp(0.656461 * (1 + config["redshift"]), w, r)) if halpha_in_range else None
        )
    for window in config["half_windows_kms"]:
        fit = fit_oiii(spec, config["redshift"], window)
        record = dict(
            label=label,
            grating=grating,
            reference_R_at_oiii=nominal,
            reference_R_at_halpha=nominal_halpha,
            reference_sha256=digest(ref_path),
            **fit,
        )
        row = {
            k: v
            for k, v in record.items()
            if k
            not in [
                "rest_wave_um",
                "flux_uJy",
                "error_uJy",
                "mask",
                "model_uJy",
                "profile_sigma_kms",
                "profile_delta_chi2",
                "parameter_vector",
            ]
        }
        if fit["status"] == "CONDITIONAL_GAUSSIAN_WIDTH":
            bound = resolution_lower_bound(fit["sigma_delta3p84_upper"])
            row["conditional_R_lower_at_oiii"] = bound
            row["conditional_q_lower"] = bound / nominal if bound else None
            # Transfer only within G395H; never extrapolate G235M to Halpha/G395M.
            row["conditional_same_grating_R_lower_at_halpha"] = (
                bound / nominal * nominal_halpha if grating == "G395H" and bound else None
            )
            if window == 900:
                ax = next(axes)
                rest, flux, error = map(
                    np.array, [fit["rest_wave_um"], fit["flux_uJy"], fit["error_uJy"]]
                )
                mask = np.array(fit["mask"])
                velocity = (rest / OIII[1] - 1) * 299792.458
                ax.errorbar(
                    velocity[mask],
                    flux[mask],
                    yerr=error[mask],
                    fmt=".",
                    color="#303030",
                    ms=3,
                    lw=0.6,
                )
                ax.plot(velocity, fit["model_uJy"], color="#ac552e", lw=1)
                ax.set(
                    title=label,
                    xlabel="Velocity relative to [OIII] 5008 (km/s)",
                    ylabel="Flux density (µJy)",
                )
                ax.grid(alpha=0.15)
                panel = left if grating == "G395H" else right
                panel.plot(fit["profile_sigma_kms"], fit["profile_delta_chi2"], label=label)
                panel.set(
                    xlim=(max(2, fit["observed_sigma_kms"] - 30), fit["observed_sigma_kms"] + 30),
                    ylim=(0, 16),
                    xlabel="Observed Gaussian σ (km/s)",
                    ylabel="Profile Δχ²",
                )
        full.append(record)
        rows.append(row)
        print(
            label,
            window,
            fit["status"],
            fit.get("observed_sigma_kms"),
            fit.get("reduced_chi2"),
            flush=True,
        )
        write_json(out / "fits.json", full)
for panel in [left, right]:
    panel.axhline(3.84, color="gray", ls="--", lw=0.8, label="Formal Δχ²=3.84")
    panel.legend(fontsize=8)
    panel.grid(alpha=0.15)
left.set_title("DJA high-resolution control")
right.set_title("Medium-resolution extractions")
fig.suptitle(
    "OBSERVED GS-159717 · [OIII] widths measured without Hα\nFree doublet amplitudes; common observed width; released diagonal errors",
    fontsize=12,
)
profile_fig.suptitle(
    "Conditional observed-width profiles · not a unique instrument calibration", fontsize=12
)
save_figure(fig, root / "figures/oiii_line_anchor")
save_figure(profile_fig, root / "figures/oiii_width_profiles")
pd.DataFrame(rows).to_csv(out / "width_constraints.csv", index=False)
write_json(out / "inputs.json", source_records)
write_json(
    out / "summary.json",
    dict(
        executed_utc=utcnow(),
        research_snapshot=config["research_snapshot"],
        source_id=config["source_id"],
        extraction_variants=len(inputs),
        attempted_windows=len(rows),
        fitted_windows=sum(r["status"] == "CONDITIONAL_GAUSSIAN_WIDTH" for r in rows),
        unique_LSF_calibrated=False,
        physical_classification="NOT_ESTABLISHED",
        input_sha256={
            p.relative_to(root).as_posix(): text_digest(p)
            for p in [
                Path(__file__),
                root / "configs/line_anchor.yaml",
                root / "src/cosmic_anomaly_lab/line_anchor.py",
                root / "src/cosmic_anomaly_lab/line_profiles.py",
            ]
        },
    ),
)
print("Completed", flush=True)

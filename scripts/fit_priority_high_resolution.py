"""Conditional narrow/blended-broad comparisons for coverage-qualified Halpha spectra."""

from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.line_profiles import fit_complex
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/high_resolution_coverage"
coverage = pd.read_csv(out / "line_coverage.csv")
targets = coverage[(coverage.line == "Hα") & (coverage.status == "COVERAGE_QUALIFIED")]
rows, details = [], []


def serial(value):
    if isinstance(value, np.ndarray):
        return [serial(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, list):
        return [serial(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


for target in targets.itertuples():
    path = root / "data/raw/spectra" / target.file
    spec = read_spectrum(path)
    rest = spec.wavelength.to_value(u.AA) / (1 + target.fixed_redshift)
    velocity = (rest / 6564.61 - 1) * 299792.458
    window = abs(velocity) <= 4500
    conversion = (np.ones(len(rest)) * u.uJy).to_value(
        u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(spec.wavelength)
    ) / 1e-20
    reference = root / f"data/raw/instrument/{target.grating.lower()}_disp.fits"
    with fits.open(reference) as h:
        nominal = float(
            np.interp(target.center_observed_um, h[1].data["WAVELENGTH"], h[1].data["R"])
        )
    fig = Figure(figsize=(11, 7), layout="constrained")
    ax, rx = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    mask = window & spec.good
    ax.errorbar(
        velocity[mask],
        (spec.flux.value * conversion)[mask],
        yerr=(spec.error.value * conversion)[mask],
        fmt=".",
        ms=3,
        color="black",
        lw=0.6,
        label="Observed",
    )
    for factor in [0.8, 1, 1.5]:
        pair = []
        for broad in [False, True]:
            model = fit_complex(
                rest[window],
                (spec.flux.value * conversion)[window],
                (spec.error.value * conversion)[window],
                spec.good[window],
                main_line=6564.61,
                components=[[(6564.61, 1)], [(6549.86, 1 / 2.96), (6585.27, 1)]],
                resolution=nominal * factor,
                family="gaussian",
                broad=broad,
            )
            details.append(
                serial(
                    dict(
                        file=target.file,
                        source_id=target.source_id,
                        reference_sha256=digest(reference),
                        raw_sha256=digest(path),
                        resolution_factor=factor,
                        **model,
                    )
                )
            )
            pair.append(model)
            if factor == 1:
                color = "#ab562c" if broad else "#327c94"
                ax.plot(
                    velocity[window],
                    model["model"],
                    color=color,
                    label="Added broad Hα" if broad else "Narrow Hα + [NII]",
                )
                rx.plot(velocity[window], model["residual"], ".-", color=color, ms=2, lw=0.7)
        rows.append(
            dict(
                source_id=target.source_id,
                file=target.file,
                grating=target.grating,
                resolution_factor=factor,
                resolution=nominal * factor,
                n=pair[0]["n"],
                delta_bic_narrow_minus_broad=pair[0]["bic"] - pair[1]["bic"],
                narrow_reduced_chi2=pair[0]["reduced_chi2"],
                broad_reduced_chi2=pair[1]["reduced_chi2"],
                broad_sigma_kms=pair[1]["broad_sigma_kms"],
                narrow_bounds=pair[0]["parameter_indices_near_bound"],
                broad_bounds=pair[1]["parameter_indices_near_bound"],
            )
        )
    ax.set(ylabel="Fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)")
    ax.legend(fontsize=9)
    rx.set(
        xlabel="Velocity at fixed catalogue z (km/s)", ylabel="Residual / error", xlim=(-4500, 4500)
    )
    rx.axhline(0, color="gray", lw=0.6)
    for panel in [ax, rx]:
        panel.grid(alpha=0.15)
    fig.suptitle(
        f"OBSERVED {target.source_id} · conditional local Hα comparison\n{target.file}\nNominal reference resolution; no absorption model or physical-class assignment",
        fontsize=10,
    )
    save_figure(fig, root / "figures/high_resolution_fits" / target.file.removesuffix(".spec.fits"))
    write_json(out / "halpha_fits.json", details)
    print(target.file, round(rows[-2]["delta_bic_narrow_minus_broad"], 2), flush=True)
pd.DataFrame(rows).to_csv(out / "halpha_comparisons.csv", index=False)
write_json(
    out / "fit_summary.json",
    dict(
        executed_utc=utcnow(),
        research_snapshot="2026-11-15",
        products=len(targets),
        source_records=int(targets.source_id.nunique()),
        fits=len(details),
        physical_classification="NOT_ESTABLISHED",
        design="Coverage-qualified local narrow Halpha+[NII] versus added Gaussian broad Halpha; three LSF scales; no significance calibration",
        input_sha256={
            p.relative_to(root).as_posix(): text_digest(p)
            for p in [
                Path(__file__),
                root / "src/cosmic_anomaly_lab/line_profiles.py",
                out / "line_coverage.csv",
            ]
        },
    ),
)

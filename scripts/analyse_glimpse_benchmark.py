"""Conditional line-shape comparisons in the authors' published GLIMPSE spectrum."""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.line_profiles import fit_complex
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
path = root / "data/raw/benchmarks/stacked_1_2_17775_CORR.spec.fits"
spec = read_spectrum(path, wavelength_range=(2.68, 5.6))
redshift = 3.50102
wave = spec.wavelength.to_value(u.AA) / (1 + redshift)
conversion = (np.ones(len(wave)) * u.uJy).to_value(
    u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(spec.wavelength)
) / 1e-20
flux, error = spec.flux.value * conversion, spec.error.value * conversion
with fits.open(root / "data/raw/instrument/g395m_disp.fits") as hdus:
    res_wave, res_r = hdus[1].data["WAVELENGTH"], hdus[1].data["R"]
windows = {
    "Halpha": dict(
        bounds=(6420, 6800),
        main_line=6564.61,
        components=[
            [(6564.61, 1)],
            [(6549.86, 1 / 3), (6585.27, 1)],
            [(6679.995, 1)],
            [(6718.29, 1)],
            [(6732.67, 1)],
        ],
    ),
    "Paschen_delta": dict(
        bounds=(9850, 10250),
        main_line=10052.1,
        components=[[(10052.1, 1)], [(9997, 1)], [(10131, 1)], [(10173, 1)]],
    ),
}
records, plot_fits, quality_audit = [], {}, {}
for name, settings in windows.items():
    selected = (wave >= settings["bounds"][0]) & (wave <= settings["bounds"][1])
    strict_good = spec.good[selected]
    finite_good = np.isfinite(flux[selected]) & np.isfinite(error[selected]) & (error[selected] > 0)
    provisional = strict_good.sum() < 25
    quality_audit[name] = dict(
        strict_valid_pixels=int(strict_good.sum()),
        finite_flux_error_pixels=int(finite_good.sum()),
        provisional_mask_override=bool(provisional),
        policy="PROVISIONAL finite flux/error only; zero NPIX unresolved"
        if provisional
        else "Strict released spectral quality mask",
    )
    for lsf_scale in [1, 1.3, 1.7]:
        resolution = float(
            np.interp(settings["main_line"] * (1 + redshift) / 1e4, res_wave, res_r) * lsf_scale
        )
        models = {}
        for family in ["gaussian", "exponential"]:
            fit = fit_complex(
                wave[selected],
                flux[selected],
                error[selected],
                finite_good if provisional else strict_good,
                main_line=settings["main_line"],
                components=settings["components"],
                resolution=resolution,
                family=family,
            )
            models[family] = fit
            if lsf_scale == 1.7:
                plot_fits[(name, family)] = fit
            summary = {k: v for k, v in fit.items() if not isinstance(v, np.ndarray)}
            records.append(
                dict(window=name, lsf_scale=lsf_scale, provisional=bool(provisional), **summary)
            )
        gaussian, exponential = models["gaussian"], models["exponential"]
        print(
            name,
            lsf_scale,
            "deltaBIC G-E",
            gaussian["bic"] - exponential["bic"],
            "reduced chi2 G/E",
            gaussian["reduced_chi2"],
            exponential["reduced_chi2"],
            flush=True,
        )
fig = Figure(figsize=(13, 9), layout="constrained")
grid = fig.add_gridspec(3, 2, height_ratios=[2, 1, 1])
for col, name in enumerate(windows):
    ax = fig.add_subplot(grid[0, col])
    g, e = plot_fits[(name, "gaussian")], plot_fits[(name, "exponential")]
    mask = g["mask"]
    x = g["wavelength"]
    ax.plot(x, np.where(mask, g["observed_flux"], np.nan), color="#20252d", lw=1, label="Observed")
    ax.fill_between(
        x,
        np.where(mask, g["observed_flux"] - g["observed_error"], np.nan),
        np.where(mask, g["observed_flux"] + g["observed_error"], np.nan),
        color="#dadfe4",
    )
    for family, color in [("gaussian", "#c77a32"), ("exponential", "#397aa5")]:
        fit = plot_fits[(name, family)]
        ax.plot(x, fit["model"], color=color, lw=1.4, label=family + " broad component")
    prefix = "PROVISIONAL · " if quality_audit[name]["provisional_mask_override"] else ""
    ax.set(
        title=f"{prefix}{name.replace('_', ' ')} · ΔBIC(G−E) = {g['bic'] - e['bic']:.1f}",
        ylabel="fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)",
    )
    ax.legend(fontsize=8)
    for i, (family, color) in enumerate([("gaussian", "#c77a32"), ("exponential", "#397aa5")], 1):
        ax = fig.add_subplot(grid[i, col])
        fit = plot_fits[(name, family)]
        ax.plot(x, fit["residual"], color=color, lw=1)
        ax.axhline(0, color="black", lw=0.6)
        ax.axhspan(-2, 2, color="#eceff2")
        ax.set(
            ylabel="(data−model) / σ",
            xlabel="Rest wavelength (Å)" if i == 2 else "",
            title=f"{family} residuals · reduced χ² = {fit['reduced_chi2']:.2f}",
        )
fig.suptitle(
    "GLIMPSE-17775 · OBSERVED G395M · CONDITIONAL PROFILE COMPARISON\n"
    "Published candidate, not a new discovery. Empirical shapes are not competing complete physical models.\n"
    "Hα: finite flux/error only despite zero NPIX; unresolved metadata conflict. Paδ: strict mask.\n"
    "Shown: resolving power ×1.7; shared narrow-line nuisance components; no absorption component.",
    fontsize=11,
)
save_figure(fig, root / "figures/glimpse_profile_comparison")
comparisons = []
for name in windows:
    for scale in [1, 1.3, 1.7]:
        g, e = [
            next(
                r
                for r in records
                if r["window"] == name and r["lsf_scale"] == scale and r["family"] == f
            )
            for f in ["gaussian", "exponential"]
        ]
        comparisons.append(
            dict(
                window=name,
                lsf_scale=scale,
                delta_bic_gaussian_minus_exponential=g["bic"] - e["bic"],
                delta_bic_with_doubled_errors=(g["chi2"] - e["chi2"]) / 4
                + (g["k"] - e["k"]) * np.log(g["n"]),
            )
        )
write_json(
    root / "results/glimpse_profile_comparison.json",
    dict(
        executed_utc=utcnow(),
        data_sha256=digest(path),
        redshift=redshift,
        code_sha256={
            "script": text_digest(Path(__file__)),
            "profiles": text_digest(root / "src/cosmic_anomaly_lab/line_profiles.py"),
        },
        benchmark_reference="https://arxiv.org/abs/2511.07515",
        data_release="https://doi.org/10.5281/zenodo.17575391",
        method="Pixel-integrated Gaussian versus Gaussian convolved with symmetric exponential (normal-Laplace) broad profile; Gaussian LSF; tied narrow nuisance lines; diagonal errors; multi-start bounded least squares",
        settings=windows,
        quality_audit=quality_audit,
        fits=records,
        comparisons=comparisons,
        limitations=[
            "Published selected benchmark, not a blind discovery",
            "Same published reduction, not independent re-reduction",
            "Empirical shape comparison, not complete scattering or radiative-transfer calculation",
            "Halpha NPIX is zero: finite-flux/error fit is explicitly provisional pending metadata resolution",
            "No absorption component; nuisance lines and continuum can affect wings",
            "No lens-corrected luminosity, mass or accretion rate inferred",
            "No stellar/AGN/composite model evidence computed",
        ],
        physical_classification="NOT_ESTABLISHED_BY_THIS_REANALYSIS",
    ),
)

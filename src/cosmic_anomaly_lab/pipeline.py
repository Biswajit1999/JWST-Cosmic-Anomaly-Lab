"""Executed pilot: observed diagnostics, held-out ranking, simulations and figures."""

import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval

from .anomaly_detection import ErrorAwareDistance, calibrated_pvalues, line_ratio_features
from .crossmatching import candidates
from .imaging import measure_cutout
from .provenance import text_digest, utcnow, write_json
from .registry import Registry
from .spectroscopy import continuum_break, fit_halpha_blend, read_spectrum
from .validation import simulate

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 220,
        "figure.constrained_layout.use": True,
    }
)
COLORS = {"CEERS": "#2563a6", "JADES": "#ba5227", "RUBIES": "#407c57", "UNCOVER": "#7954a1"}


def savefig(root, name, figure):
    for extension in ["png", "pdf"]:
        figure.savefig(root / "figures" / f"{name}.{extension}", bbox_inches="tight")
    plt.close(figure)


def analyse_spectra(root, run_id):
    config = yaml.safe_load((root / "configs/pilot.yaml").read_text())
    pilot = pd.read_csv(root / "data/catalogues/pilot_selection.csv")
    reference = fits.getdata(root / "data/raw/instrument/prism_disp.fits", 1)
    registry = Registry(root / "data/processed/registry.sqlite")
    measurements, fits_rows, failures = [], [], []
    nrows = max(1, (len(pilot) + 3) // 4)
    fig, axes = plt.subplots(nrows, 4, figsize=(15, 2.8 * nrows), squeeze=False)
    for ax, row in zip(axes.flat, pilot.itertuples()):
        try:
            spec = read_spectrum(
                root / "data/raw/spectra" / row.file,
                wavelength_range=config["spectroscopy"]["observed_wavelength_um"],
            )
            spec.as_specutils()  # Verify interoperability with the supported spectral API.
            brk = continuum_break(
                spec,
                row.z_best,
                blue=config["continuum_proxy"]["rest_blue_um"],
                red=config["continuum_proxy"]["rest_red_um"],
            )
            line = fit_halpha_blend(spec, row.z_best, reference["WAVELENGTH"], reference["R"])
            variants = {}
            for scale in [1.0, 1.6]:
                variants[str(scale)] = fit_halpha_blend(
                    spec, row.z_best, reference["WAVELENGTH"], reference["R"], lsf_scale=scale
                )
            inflated = fit_halpha_blend(
                spec, row.z_best, reference["WAVELENGTH"], reference["R"], error_scale=1.5
            )
            best = min(line["fits"], key=lambda f: f["bic"]) if line["fits"] else None
            measurement = dict(
                source_id=row.source_id,
                file=row.file,
                survey=row.survey,
                field=row.field,
                sample=row.sample,
                z_spec=row.z_best,
                good_pixels=int(spec.good.sum()),
                median_pixel_sn=float(
                    np.median(spec.flux.value[spec.good] / spec.error.value[spec.good])
                ),
                break_ratio=brk["ratio"],
                break_ratio_error=brk["ratio_error"],
                break_status=brk["status"],
                blend_fit_status=line["status"],
                preferred_empirical_profile=best["model"] if best else None,
            )
            measurements.append(measurement)
            for fit in line["fits"]:
                fits_rows.append(dict(source_id=row.source_id, file=row.file, **fit))
            evidence = dict(
                data_kind="OBSERVED",
                measurement=measurement,
                continuum_break=brk,
                empirical_blend=line,
                lsf_sensitivity=variants,
                error_inflation_1p5=inflated,
                literature=row.classification_reference
                if isinstance(row.classification_reference, str)
                else None,
                status="INSUFFICIENT_EVIDENCE",
                interpretation="No physical family comparison performed",
                systematics={
                    k: "NOT_TESTED"
                    for k in [
                        "foreground",
                        "alternate_redshift",
                        "blending",
                        "lensing",
                        "independent_reduction",
                        "calibration",
                        "model_library_completeness",
                    ]
                },
            )
            registry.add_evidence(row.source_id, run_id, "INSUFFICIENT_EVIDENCE", evidence)
            write_json(root / "reports/sources" / f"{row.source_id}.json", evidence)
            source_report = root / "reports/sources" / f"{row.source_id}.md"
            source_report.write_text(
                f"# {row.source_id}\n\nObserved {row.survey} / {row.field}; DJA v4.4 reduction.\n\n"
                f"- Product: [{row.file}]({row.url})\n- Archive redshift: {row.z_best:.6f}; quality grade {row.grade:g}. No propagated redshift uncertainty.\n"
                f"- Sample {row.sample}: {row.selection}\n- Usable spectral pixels: {spec.good.sum()}\n"
                f"- Continuum ratio: {brk['ratio']} +/- {brk['ratio_error']} (dimensionless; conditional statistical error).\n"
                f"- Local blend model: {measurement['preferred_empirical_profile']}. This is an empirical profile, not a physical classification.\n"
                "- Assessment: **INSUFFICIENT_EVIDENCE** for deciding among AGN, stars, galaxy composites and reprocessing models.\n\n"
                "Review the paired JSON for measured fluxes (erg/s/cm2), equivalent widths (angstrom), uncertainties, and LSF sensitivity. "
                "Remaining checks: 2D extraction, calibration, foreground, alternate redshifts, blends, lensing, independent reduction, and full physical model library. "
                "Follow-up: higher-resolution hydrogen/metal-line profiles and PSF-matched imaging. This is a qualitative proposal, not an exposure forecast.\n",
                encoding="utf-8",
            )
            good = spec.good
            ax.plot(
                spec.wavelength.value,
                np.where(good, spec.flux.value, np.nan),
                color=COLORS[row.survey],
                lw=0.7,
            )
            ax.fill_between(
                spec.wavelength.value,
                np.where(good, (spec.flux - spec.error).value, np.nan),
                np.where(good, (spec.flux + spec.error).value, np.nan),
                color=COLORS[row.survey],
                alpha=0.15,
                lw=0,
            )
            ax.set(
                title=f"{row.survey} {row.srcid} | z={row.z_best:.2f} | {row.sample}",
                xlabel="Observed wavelength (um)",
                ylabel="Flux density (uJy)",
            )
            ax.set_xlim(0.7, 5.2)
            ax.tick_params(labelsize=8)
        except Exception as exc:
            failures.append(dict(file=row.file, error=str(exc)))
            ax.set_title(f"Failed: {row.file}")
    registry.close()
    for ax in axes.flat[len(pilot) :]:
        ax.set_visible(False)
    fig.suptitle(
        f"OBSERVED: {len(pilot)}-source stratified JWST PRISM pilot; shaded regions are pixel errors",
        fontsize=15,
    )
    savefig(root, "observed_spectra", fig)
    pd.DataFrame(measurements).to_csv(root / "results/spectral_measurements.csv", index=False)
    pd.DataFrame(fits_rows).to_csv(root / "results/empirical_blend_fits.csv", index=False)
    write_json(root / "results/spectral_failures.json", failures)
    return pd.DataFrame(measurements), pd.DataFrame(fits_rows)


def analyse_images(root):
    pilot = pd.read_csv(root / "data/catalogues/pilot_selection.csv").set_index("source_id")
    records, failures = [], []
    files = sorted((root / "data/raw/images").glob("*.fits"))
    fig, axes = plt.subplots(len(files), 3, figsize=(10, 3.2 * len(files)), squeeze=False)
    for i, path in enumerate(files):
        row = pilot.loc[path.stem]
        try:
            measurements, panels = measure_cutout(path, row.ra, row.dec)
            for measurement in measurements:
                records.append(dict(source_id=path.stem, survey=row.survey, **measurement))
            for j, (band, data, wcs) in enumerate(panels):
                scale = measurements[j]["pixel_scale_arcsec"]
                ny, nx = data.shape
                x0, y0 = wcs.world_to_pixel_values(row.ra, row.dec)
                norm = ImageNormalize(
                    data, interval=PercentileInterval(99.5), stretch=AsinhStretch()
                )
                axes[i, j].imshow(
                    data,
                    origin="lower",
                    norm=norm,
                    cmap="magma",
                    extent=[
                        (-0.5 - x0) * scale,
                        (nx - 0.5 - x0) * scale,
                        (-0.5 - y0) * scale,
                        (ny - 0.5 - y0) * scale,
                    ],
                )
                axes[i, j].set(
                    title=f"{row.survey} {row.srcid} | {band}",
                    xlabel="Offset (arcsec)",
                    ylabel="Offset (arcsec)",
                )
                axes[i, j].add_patch(plt.Circle((0, 0), 0.3, fill=False, color="white", lw=0.8))
        except Exception as exc:
            failures.append(dict(source_id=path.stem, error=str(exc)))
    fig.suptitle("OBSERVED: DJA NIRCam cutouts; independent display stretches, PSFs unmatched")
    savefig(root, "observed_cutouts", fig)
    pd.DataFrame(records).to_csv(root / "results/aperture_diagnostics.csv", index=False)
    write_json(root / "results/imaging_failures.json", failures)
    return records


def compare_catalogues(root):
    sources = pd.read_csv(root / "data/catalogues/sources.csv.gz")
    comparisons = []
    for survey, filename, zcol, quality_col in [
        ("JADES", "jades_dr4_reference.csv.gz", "z_Spec", "z_Spec_flag"),
        ("UNCOVER", "uncover_dr4_1_reference.csv", "z_spec", "flag_zspec_qual"),
    ]:
        left = sources[sources.survey.str.contains(survey)].reset_index(drop=True)
        right = pd.read_csv(root / "data/catalogues" / filename)
        match = candidates(left, right, 0.3)
        match["survey"] = survey
        match["source_id"] = left.iloc[match.left_index].source_id.to_numpy()
        match["dja_z"] = left.iloc[match.left_index].z_spec.to_numpy()
        match["dja_quality"] = left.iloc[match.left_index].z_quality.to_numpy()
        match["reference_z"] = right.iloc[match.right_index][zcol].to_numpy()
        match["reference_quality"] = right.iloc[match.right_index][quality_col].to_numpy()
        match["delta_z_over_1pz"] = (match.dja_z - match.reference_z) / (1 + match.reference_z)
        if survey == "UNCOVER":
            match["background_advisory"] = right.iloc[
                match.right_index
            ].flag_potential_local_background_issue.to_numpy()
        comparisons.append(match)
    matches = pd.concat(comparisons, ignore_index=True)
    matches.to_csv(root / "results/catalogue_counterpart_candidates.csv.gz", index=False)
    robust_ref = ((matches.survey == "JADES") & matches.reference_quality.isin(["A", "B"])) | (
        (matches.survey == "UNCOVER")
        & pd.to_numeric(matches.reference_quality, errors="coerce").ge(2)
    )
    valid = matches[
        (~matches.ambiguous)
        & (matches.dja_quality == 3)
        & robust_ref
        & matches.dja_z.ge(0)
        & matches.reference_z.ge(0)
    ]
    fig, ax = plt.subplots(figsize=(7, 5))
    for survey, group in valid.groupby("survey"):
        ax.scatter(
            group.reference_z,
            group.delta_z_over_1pz,
            s=8,
            alpha=0.4,
            label=f"{survey}: n={len(group)}",
            color=COLORS[survey],
        )
    ax.axhline(0, color="black", lw=0.7)
    ax.set(
        xlabel="Reference spectroscopic redshift",
        ylabel="(z_DJA - z_reference) / (1 + z_reference)",
        title="OBSERVED: positional candidates with robust archive grades",
    )
    ax.legend()
    savefig(root, "catalogue_redshift_comparison", fig)
    summary = dict(
        candidate_pairs=len(matches),
        unambiguous_robust_pairs=len(valid),
        abs_normalized_difference_over_0p02=int((valid.delta_z_over_1pz.abs() > 0.02).sum()),
        median_abs_normalized_difference=float(valid.delta_z_over_1pz.abs().median()),
        caveat="Positional candidates only; quality definitions differ. No automatic source merge or redshift adjudication.",
    )
    write_json(root / "results/redshift_comparison.json", summary)
    return summary


def rank_archive(root):
    frame = pd.read_csv(root / "data/catalogues/observations.csv.gz", low_memory=False)
    # The archive changes from blended to separate [OIII] columns with resolution.
    # Combine BOTH components; never compare 5007 alone to a doublet sum.
    missing = frame.line_oiii.isna()
    frame["oiii_representation"] = "archive_doublet_sum"
    frame.loc[missing, "line_oiii"] = (
        frame.loc[missing, "line_oiii_5007"] + frame.loc[missing, "line_oiii_4959"]
    )
    frame.loc[missing, "line_oiii_err"] = np.hypot(
        frame.loc[missing, "line_oiii_5007_err"], frame.loc[missing, "line_oiii_4959_err"]
    )
    frame.loc[missing, "oiii_representation"] = "sum_of_components_covariance_unavailable"
    frame = frame[
        (frame.grating == "PRISM")
        & (frame.grade == 3)
        & frame.z_best.between(2, 6)
        & (((frame.zline - frame.z_best) / (1 + frame.z_best)).abs() < 0.005)
    ].copy()
    frame = frame.sort_values(["sn50", "file"], ascending=[False, True]).drop_duplicates(
        "source_id"
    )
    mask, x, cov = line_ratio_features(frame)
    frame = frame.loc[mask].reset_index(drop=True)
    frame["feature_oiii_hb_dex"], frame["feature_oiii_oii_dex"] = x[:, 0], x[:, 1]
    frame["z_bin"] = np.where(frame.z_best < 4, "2-4", "4-6")
    frame["hash"] = frame.source_id.map(
        lambda x: int(hashlib.sha256(x.encode()).hexdigest()[:8], 16)
    )
    rankings, folds = [], []
    for zbin in ["2-4", "4-6"]:
        for field in sorted(frame.field.unique()):
            test = (frame.z_bin == zbin) & (frame.field == field)
            train_pool = (frame.z_bin == zbin) & (frame.field != field)
            # Independent calibration split, common IDs cannot leak between fields.
            training = train_pool & (frame["hash"] % 3 != 0)
            calibration = train_pool & (frame["hash"] % 3 == 0)
            counts = dict(
                z_bin=zbin,
                held_out_field=field,
                training=int(training.sum()),
                calibration=int(calibration.sum()),
                test=int(test.sum()),
            )
            if min(training.sum(), calibration.sum()) < 30 or test.sum() == 0:
                folds.append(counts | {"status": "INSUFFICIENT_SAMPLE"})
                continue
            model = ErrorAwareDistance().fit(x[training], cov[training])
            cs = model.score(x[calibration], cov[calibration])
            score = model.score(x[test], cov[test])
            p = calibrated_pvalues(cs, score)
            out = frame.loc[
                test,
                [
                    "source_id",
                    "file",
                    "survey",
                    "field",
                    "z_best",
                    "url",
                    "feature_oiii_hb_dex",
                    "feature_oiii_oii_dex",
                ],
            ].copy()
            out["distance"] = score
            out["calibration_p"] = p
            out["calibration_n"] = int(calibration.sum())
            out["assessment"] = "INSUFFICIENT_EVIDENCE"
            out["interpretation"] = "UNVETTED_STATISTICAL_OUTLIER_NOT_PHYSICAL_DISCOVERY"
            rankings.append(out)
            folds.append(
                counts
                | {
                    "status": "EXECUTED",
                    "fraction_p_le_0p01": float(np.mean(p <= 0.01)),
                    "null_FPR_measured": False,
                }
            )
    ranked = pd.concat(rankings, ignore_index=True) if rankings else pd.DataFrame()
    if len(ranked):
        ranked = ranked.sort_values(["calibration_p", "distance"], ascending=[True, False])
        from scipy.stats import false_discovery_control

        ranked["bh_adjusted_p_descriptive"] = false_discovery_control(
            ranked.calibration_p.to_numpy()
        )
        ranked.to_csv(root / "results/archive_ranking.csv", index=False)
    write_json(root / "results/anomaly_folds.json", folds)
    frame.drop(columns="hash").to_csv(root / "results/anomaly_feature_sample.csv.gz", index=False)
    fig, ax = plt.subplots(figsize=(7, 5))
    for survey, group in frame.groupby("survey"):
        ax.scatter(
            group.feature_oiii_hb_dex,
            group.feature_oiii_oii_dex,
            s=10,
            alpha=0.45,
            label=f"{survey} (n={len(group)})",
            color=COLORS[survey],
        )
    ax.set(
        xlabel="log10([OIII] / Hbeta)",
        ylabel="log10([OIII] / [OII])",
        title="OBSERVED archive line fits: z=2-6; all three lines S/N >= 5",
    )
    ax.legend()
    savefig(root, "archive_line_ratio_population", fig)
    return dict(
        feature_sources=len(frame),
        ranked_sources=len(ranked),
        executed_folds=sum(f["status"] == "EXECUTED" for f in folds),
        nominal_p_le_0p01_unvetted=int((ranked.calibration_p <= 0.01).sum()) if len(ranked) else 0,
        bh_adjusted_p_le_0p05=int((ranked.bh_adjusted_p_descriptive <= 0.05).sum())
        if len(ranked)
        else 0,
        interpretation="Calibration assumes exchangeability; field shift and shared fits limit p-value interpretation",
    )


def run_simulations(root):
    metrics, pulls, control_p, injection_p, wave, flux = simulate()
    write_json(root / "results/simulation_validation.json", metrics)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(pulls, bins=25, density=True, color="#2563a6", alpha=0.7)
    axes[0].set(
        xlabel="(Recovered line flux - truth) / error",
        ylabel="Probability density",
        title="SIMULATED: 500 line recoveries",
    )
    axes[1].hist(
        control_p, bins=np.linspace(0, 1, 21), alpha=0.65, density=True, label="Null controls"
    )
    axes[1].hist(
        injection_p,
        bins=np.linspace(0, 1, 21),
        alpha=0.65,
        density=True,
        label="Injected 1.5 dex shift",
    )
    axes[1].set(
        xlabel="Held-out calibration p",
        ylabel="Probability density",
        title="SIMULATED: anomaly injection / recovery",
    )
    axes[1].legend()
    savefig(root, "simulation_validation", fig)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.loglog(wave.value * 5, flux.value)
    ax.set(
        xlabel="Observed wavelength at z=4 (um)",
        ylabel="Flux density (uJy)",
        title="SIMULATED: 5000 K, 1000 solar-radius blackbody; no IGM or dust",
    )
    savefig(root, "analytic_blackbody_baseline", fig)
    return metrics


def run(root):
    root = Path(root).resolve()
    for sub in ["results", "figures", "reports/sources"]:
        (root / sub).mkdir(parents=True, exist_ok=True)
    tracked_inputs = list((root / "src/cosmic_anomaly_lab").glob("*.py")) + list(
        (root / "configs").glob("*.yaml")
    )
    tracked_inputs += [
        root / "data/catalogues/pilot_selection.csv",
        root / "data/manifests/downloads.json",
    ]
    input_hashes = {p.relative_to(root).as_posix(): text_digest(p) for p in sorted(tracked_inputs)}
    dependencies = {
        name: version(name)
        for name in [
            "numpy",
            "scipy",
            "astropy",
            "astroquery",
            "specutils",
            "pandas",
            "scikit-learn",
            "matplotlib",
        ]
    }
    signature = hashlib.sha256(
        json.dumps({"inputs": input_hashes, "dependencies": dependencies}, sort_keys=True).encode()
    ).hexdigest()
    run_id = "pilot-v0.1-" + signature[:10]
    measurements, fits_frame = analyse_spectra(root, run_id)
    print(f"Spectra analysed: {len(measurements)}", flush=True)
    images = analyse_images(root)
    comparisons = compare_catalogues(root)
    ranking = rank_archive(root)
    simulations = run_simulations(root)
    summary = dict(
        run_id=run_id,
        input_text_sha256=input_hashes,
        input_hash_convention="UTF-8 with universal-newline normalization; download checksums remain raw bytes",
        dependency_versions=dependencies,
        executed_utc=utcnow(),
        spectra_analysed=len(measurements),
        break_measurements=int(measurements.break_ratio.notna().sum()),
        blend_fit_sources=int(fits_frame.source_id.nunique()) if len(fits_frame) else 0,
        aperture_measurements=len(images),
        redshift_comparison=comparisons,
        anomaly=ranking,
        simulations=simulations,
        figures=sorted(p.name for p in (root / "figures").glob("*")),
    )
    write_json(root / "results/run_summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)
    return summary

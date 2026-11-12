"""Fixed-sample observational screening, with null outcomes and failures retained."""

import gzip
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from matplotlib.figure import Figure

from .provenance import digest, text_digest, utcnow, write_json
from .spectroscopy import continuum_break, fit_halpha_blend, read_spectrum
from .visual_evidence import save_figure

matplotlib.use("Agg")


def triage_flags(record, config):
    """Return screening predicates, not posterior probabilities or classes."""

    def exceeds(key, threshold):
        value = record.get(key)
        return value is not None and np.isfinite(value) and value > threshold

    continuum = (
        exceeds("continuum_ratio", config["continuum_ratio_minimum"])
        and exceeds("continuum_lower_3sigma", config["continuum_ratio_lower_3sigma_minimum"])
        and exceeds(
            "alternate_continuum_lower_3sigma", config["continuum_ratio_lower_3sigma_minimum"]
        )
    )
    broad = exceeds(
        "minimum_broad_delta_bic", config["broad_over_narrow_delta_BIC_minimum"]
    ) and exceeds("minimum_broad_sn", config["broad_blend_SN_minimum"])
    return dict(
        continuum_screen=bool(continuum),
        broad_blend_screen=bool(broad),
        joint_followup_priority=bool(continuum and broad),
    )


def screen_source(path, row, resolution_wave, resolution_r, config):
    spec = read_spectrum(path)
    brk = continuum_break(spec, row.z_best)
    alternate = continuum_break(
        spec,
        row.z_best,
        blue=config["alternate_continuum_blue_um"],
        red=config["alternate_continuum_red_um"],
    )
    record = dict(
        source_id=row.source_id,
        file=row.file,
        survey=row.survey,
        field=row.field,
        z_best=float(row.z_best),
        archive_sn50=float(row.sn50),
        status="ANALYSED",
        usable_pixels=int(spec.good.sum()),
        raw_sha256=digest(path),
        continuum_status=brk["status"],
        continuum_ratio=brk["ratio"],
        continuum_ratio_error=brk["ratio_error"],
        continuum_lower_3sigma=None
        if brk["ratio"] is None
        else brk["ratio"] - 3 * brk["ratio_error"],
        alternate_continuum_lower_3sigma=None
        if alternate["ratio"] is None
        else alternate["ratio"] - 3 * alternate["ratio_error"],
        minimum_broad_delta_bic=None,
        minimum_broad_sn=None,
        nominal_broad_delta_bic=None,
        nominal_broad_reduced_chi2=None,
        line_status="INSUFFICIENT_COVERAGE",
        lsf_preference_changes=False,
    )
    fits_by_lsf = {}
    for scale in config["lsf_scales"]:
        fit = fit_halpha_blend(spec, row.z_best, resolution_wave, resolution_r, lsf_scale=scale)
        fits_by_lsf[str(scale)] = fit
    available = [f["fits"] for f in fits_by_lsf.values() if f["fits"]]
    if len(available) == len(config["lsf_scales"]):
        differences = [f[1]["bic"] - f[2]["bic"] for f in available]
        record.update(
            minimum_broad_delta_bic=float(min(differences)),
            minimum_broad_sn=float(min(f[2]["sn"] for f in available)),
            nominal_broad_delta_bic=float(differences[1]),
            nominal_broad_reduced_chi2=float(available[1][2]["chi2"] / available[1][2]["dof"]),
            line_status="EMPIRICAL_BLEND_DIAGNOSTIC",
            lsf_preference_changes=bool(min(differences) < 0 < max(differences)),
        )
    record.update(triage_flags(record, config))
    return record, dict(continuum=brk, alternate_continuum=alternate, fits_by_lsf=fits_by_lsf)


def run_search(root, limit=None):
    root = Path(root).resolve()
    config = yaml.safe_load((root / "configs/expanded_search.yaml").read_text())
    rows = pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz")
    denominator = len(rows)
    if limit is not None:
        rows = rows.head(limit)
    destination = root / ("results/expanded" if limit is None else "data/processed/expanded_smoke")
    destination.mkdir(parents=True, exist_ok=True)
    with fits.open(root / "data/raw/instrument/prism_disp.fits") as hdus:
        resolution_wave = np.array(hdus[1].data["WAVELENGTH"])
        resolution_r = np.array(hdus[1].data["R"])
    measurements = []
    detail_stream = gzip.open(destination / "source_fits.jsonl.gz", "wt", encoding="utf-8")
    for i, row in enumerate(rows.itertuples(), 1):
        try:
            record, detail = screen_source(
                root / "data/raw/spectra" / row.file,
                row,
                resolution_wave,
                resolution_r,
                config["triage"],
            )
            detail_stream.write(json.dumps(dict(**record, **detail), allow_nan=False) + "\n")
        except (OSError, ValueError, KeyError, np.linalg.LinAlgError) as exc:
            record = dict(
                source_id=row.source_id,
                file=row.file,
                survey=row.survey,
                field=row.field,
                z_best=float(row.z_best),
                status="FAILED",
                failure_reason=str(exc),
            )
        measurements.append(record)
        if i % 100 == 0:
            print(f"Analysed {i}/{len(rows)}", flush=True)
    detail_stream.close()
    frame = pd.DataFrame(measurements)
    frame.to_csv(destination / "measurements.csv", index=False)
    good = frame[frame.status == "ANALYSED"].copy()
    priority = good[good.joint_followup_priority.fillna(False)].sort_values(
        ["continuum_lower_3sigma", "minimum_broad_delta_bic"], ascending=False
    )
    priority.to_csv(destination / "followup_priority.csv", index=False)
    summary = dict(
        executed_utc=utcnow(),
        eligible_denominator=denominator,
        attempted=len(frame),
        analysed=len(good),
        fraction_of_eligible_analysed=len(good) / denominator,
        failed=int((frame.status == "FAILED").sum()),
        continuum_measured=int(good.continuum_ratio.notna().sum()),
        line_blends_measured=int(good.minimum_broad_delta_bic.notna().sum()),
        continuum_screen_pass=int(good.continuum_screen.sum()),
        broad_blend_screen_pass=int(good.broad_blend_screen.sum()),
        joint_followup_priorities=len(priority),
        lsf_preference_changes=int(good.lsf_preference_changes.sum()),
        joint_retained_with_doubled_line_errors=int(
            (
                (priority.minimum_broad_delta_bic / 4 > 10) & (priority.minimum_broad_sn / 2 > 5)
            ).sum()
        ),
        by_survey=good.groupby("survey")
        .agg(analysed=("source_id", "size"), joint_priorities=("joint_followup_priority", "sum"))
        .reset_index()
        .to_dict("records"),
        input_sha256={
            "selection": digest(root / "data/catalogues/expanded_selection.csv.gz"),
            "config": text_digest(root / "configs/expanded_search.yaml"),
            "analysis": text_digest(Path(__file__)),
            "spectroscopy": text_digest(root / "src/cosmic_anomaly_lab/spectroscopy.py"),
        },
        physical_classification="NOT_ESTABLISHED",
        interpretation="Threshold-selected follow-up queue, not detections or discovery probabilities",
        limitations=[
            "Targeted high-SN successful-redshift subset",
            "Common DJA reduction",
            "Unresolved blends, fixed archive z, conditional diagonal errors",
            "No physical stellar/AGN/envelope evidence comparison",
            "No calibrated multiple-testing significance",
        ],
    )
    write_json(destination / "summary.json", summary)
    if limit is None:
        plot_search(root, good, priority)
    return summary


def plot_search(root, frame, priority):
    fig = Figure(figsize=(12, 5.8), layout="constrained")
    axes = fig.subplots(1, 2)
    colors = dict(CEERS="#286b9b", JADES="#4c8d75", RUBIES="#b46e3e", UNCOVER="#845f9c")
    for survey, group in frame.groupby("survey"):
        mask = group.continuum_ratio.between(0.1, 30)
        axes[0].scatter(
            group.loc[mask, "z_best"],
            group.loc[mask, "continuum_ratio"],
            s=8,
            alpha=0.45,
            color=colors[survey],
            label=survey,
        )
        mask = group.continuum_ratio.between(0.1, 30) & group.minimum_broad_delta_bic.notna()
        axes[1].scatter(
            group.loc[mask, "continuum_ratio"],
            group.loc[mask, "minimum_broad_delta_bic"],
            s=8,
            alpha=0.35,
            color=colors[survey],
        )
    axes[0].set(xlabel="Archive redshift", ylabel="Red / blue continuum ratio", yscale="log")
    axes[0].legend(fontsize=9)
    axes[0].axhline(2, color="black", ls="--", lw=0.8)
    axes[1].scatter(
        priority.continuum_ratio,
        priority.minimum_broad_delta_bic,
        s=45,
        facecolors="none",
        edgecolors="#151515",
        label="Joint follow-up screen",
    )
    axes[1].set(
        xlabel="Red / blue continuum ratio",
        ylabel="Minimum ΔBIC (narrow − broad), across LSF tests",
        xscale="log",
        yscale="symlog",
    )
    axes[1].axhline(10, color="black", ls="--", lw=0.8)
    axes[1].axvline(2, color="black", ls="--", lw=0.8)
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=0.15)
    fig.suptitle(
        f"OBSERVED · {len(frame):,} spectra screened | {len(priority)} follow-up priorities\n"
        "Descriptive thresholds; neither a physical classification nor a calibrated discovery significance",
        fontsize=12,
    )
    save_figure(fig, root / "figures/expanded_search")

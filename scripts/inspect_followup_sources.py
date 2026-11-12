"""Images and higher-resolution checks for a bounded, explicitly selected queue."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.line_profiles import fit_complex
from cosmic_anomaly_lab.provenance import DownloadCache, digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import evidence_sheet, read_images, save_figure

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--resume", action="store_true", help="Reuse completed source inspections")
args = parser.parse_args()
measured = pd.read_csv(root / "results/expanded/measurements.csv")
priority = pd.read_csv(root / "results/expanded/followup_priority.csv")
selection = pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz")
observations = pd.read_csv(root / "data/catalogues/observations.csv.gz")
chosen = priority.head(12).source_id.tolist()
# Also inspect comparison sources of similar redshift and S/N, without claiming
# they are physically normal or matched in every confounding variable.
controls, pairs = [], []
for sid in chosen[:4]:
    candidate = measured[measured.source_id == sid].iloc[0]
    pool = measured[
        (measured.status == "ANALYSED")
        & ~measured.source_id.isin(chosen + controls)
        & (measured.survey == candidate.survey)
        & ~measured.joint_followup_priority.fillna(False).astype(bool)
        & ((measured.z_best - candidate.z_best).abs() / (1 + candidate.z_best) < 0.1)
        & measured.archive_sn50.between(candidate.archive_sn50 / 2, candidate.archive_sn50 * 2)
    ].copy()
    if len(pool):
        distance = (
            abs(pool.z_best - candidate.z_best) / (1 + candidate.z_best)
            + abs(np.log(pool.archive_sn50 / candidate.archive_sn50)) * 0.2
        )
        control = pool.loc[distance.idxmin()]
        controls.append(control.source_id)
        pairs.append(
            dict(
                priority_source_id=sid,
                comparison_source_id=control.source_id,
                normalized_redshift_difference=float(
                    abs(control.z_best - candidate.z_best) / (1 + candidate.z_best)
                ),
                archive_SN_ratio=float(control.archive_sn50 / candidate.archive_sn50),
            )
        )
benchmarks = ["CAL-ebc2a6a997b3535d"]  # The Cliff, RUBIES-UDS-154183; de Graaff et al. 2025
queue = selection[selection.source_id.isin(chosen + controls + benchmarks)].copy()
queue["inspection_role"] = np.where(
    queue.source_id.isin(chosen), "JOINT_SCREEN_PRIORITY", "COMPARISON"
)
queue.loc[queue.source_id.isin(benchmarks), "inspection_role"] = "PUBLISHED_BENCHMARK"
queue.to_csv(root / "results/expanded/inspection_selection.csv", index=False)
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
audit = []
old_path = root / "results/expanded/inspection_audit.json"
previous = (
    {r["source_id"]: r for r in json.loads(old_path.read_text())["sources"]}
    if args.resume and old_path.exists()
    else {}
)
for row in queue.itertuples():
    if row.source_id in previous:
        audit.append(previous[row.source_id])
        continue
    record = dict(
        source_id=row.source_id,
        survey=row.survey,
        source_number=int(row.srcid),
        role=row.inspection_role,
        redshift=float(row.z_best),
    )
    try:
        url = f"https://grizli-cutout.herokuapp.com/thumb?ra={row.ra}&dec={row.dec}&size=3&filters=f150w-clear,f277w-clear,f444w-clear&output=fits_weight"
        image_path = cache.fetch(
            url,
            "images/" + row.source_id + ".fits",
            survey=row.survey,
            version="DJA-cutout-live",
            kind="observed_image_cutout",
        )
        images, rgb, detail = read_images(image_path, row.ra, row.dec)
        record.update(
            image_status="RGB_AVAILABLE" if rgb is not None else "PARTIAL_COVERAGE",
            image_sha256=digest(image_path),
            image_diagnostics=detail,
        )
    except Exception as exc:
        images, rgb = {}, None
        record.update(image_status="UNAVAILABLE", image_error=str(exc))
    record.update(evidence_sheet(root, row, images, rgb, "figures/expanded_sources"))
    # Prefer one medium-resolution observation with the archive Halpha region
    # bracketed. This remains an upstream source association, not an independent
    # positional merge or an independent reduction.
    ha = 0.656461 * (1 + row.z_best)
    medium = observations[
        (observations.source_id == row.source_id)
        & observations.grating.isin(["G140M", "G235M", "G395M"])
        & (observations.grade == 3)
        & (observations.wmin < ha * 0.975)
        & (observations.wmax > ha * 1.036)
    ].sort_values("sn50", ascending=False)
    record["medium_resolution_status"] = "NO_ELIGIBLE_MEDIUM_RESOLUTION_COUNTERPART"
    if len(medium):
        m = medium.iloc[0]
        try:
            path = cache.fetch(
                m.url,
                "spectra/" + m.file,
                survey=m.survey,
                version="DJA-v4.4",
                kind="observed_spectrum",
            )
            instrument = cache.fetch(
                f"https://raw.githubusercontent.com/gbrammer/msaexp/main/msaexp/data/jwst_nirspec_{m.grating.lower()}_disp.fits",
                "instrument/" + m.grating.lower() + "_disp.fits",
                survey="instrument-reference",
                version="msaexp-main-snapshot",
                kind="resolving_power",
            )
            spec = read_spectrum(path)
            rest = spec.wavelength.to_value(u.AA) / (1 + row.z_best)
            window = (rest >= 6420) & (rest <= 6800)
            if np.count_nonzero(spec.good & window) < 25:
                raise ValueError("Insufficient medium-resolution Halpha coverage")
            conversion = (np.ones(len(rest)) * u.uJy).to_value(
                u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(spec.wavelength)
            ) / 1e-20
            with fits.open(instrument) as hdus:
                resolution = float(
                    np.interp(ha, hdus[1].data["WAVELENGTH"], hdus[1].data["R"]) * 1.3
                )
            components = [
                [(6564.61, 1)],
                [(6549.86, 1 / 3), (6585.27, 1)],
                [(6679.995, 1)],
                [(6718.29, 1)],
                [(6732.67, 1)],
            ]
            models = {}
            for name, broad in [("narrow_only", False), ("gaussian", True), ("exponential", True)]:
                models[name] = fit_complex(
                    rest[window],
                    (spec.flux.value * conversion)[window],
                    (spec.error.value * conversion)[window],
                    spec.good[window],
                    main_line=6564.61,
                    components=components,
                    resolution=resolution,
                    family="gaussian" if name == "narrow_only" else name,
                    broad=broad,
                )
            record.update(
                medium_resolution_status="CONDITIONAL_PROFILE_FITS",
                medium_file=m.file,
                medium_sha256=digest(path),
                medium_grating=m.grating,
                medium_fits={
                    name: {k: v for k, v in f.items() if not isinstance(v, np.ndarray)}
                    for name, f in models.items()
                },
                medium_delta_bic_narrow_minus_gaussian=models["narrow_only"]["bic"]
                - models["gaussian"]["bic"],
                medium_delta_bic_gaussian_minus_exponential=models["gaussian"]["bic"]
                - models["exponential"]["bic"],
            )
            fig = Figure(figsize=(10, 7), layout="constrained")
            axes = fig.subplots(2, 1, gridspec_kw={"height_ratios": [2, 1]})
            g = models["gaussian"]
            x = g["wavelength"]
            axes[0].errorbar(
                x[g["mask"]],
                g["observed_flux"][g["mask"]],
                g["observed_error"][g["mask"]],
                fmt=".",
                color="#454545",
                ms=3,
                label="Observed",
            )
            for name, color in [
                ("narrow_only", "#8a8a8a"),
                ("gaussian", "#bf7c37"),
                ("exponential", "#397aa5"),
            ]:
                axes[0].plot(x, models[name]["model"], color=color, label=name.replace("_", " "))
                axes[1].plot(x, models[name]["residual"], color=color, lw=0.8)
            axes[0].legend(fontsize=9)
            axes[0].set(ylabel="fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)")
            axes[1].axhspan(-2, 2, color="#edf0f2")
            axes[1].axhline(0, color="black", lw=0.6)
            axes[1].set(xlabel="Rest wavelength (Å)", ylabel="(data−model) / σ")
            fig.suptitle(
                f"OBSERVED · {row.survey} {row.srcid} · {m.grating} · {row.inspection_role}\n"
                "Conditional shapes, nominal LSF ×1.3; no absorption component; no physical classification",
                fontsize=11,
            )
            save_figure(fig, root / "figures/expanded_sources" / (row.source_id + "_medium"))
        except Exception as exc:
            record.update(medium_resolution_status="UNAVAILABLE_OR_FAILED", medium_reason=str(exc))
    audit.append(record)
    write_json(
        root / "results/expanded/inspection_audit.json",
        dict(
            executed_utc=utcnow(),
            sources=audit,
            comparison_pairs=pairs,
            limitation="Bounded inspection of up to twelve highest continuum-lower-bound joint priorities and up to four approximate comparison sources; not validation of the whole screening queue",
        ),
    )
    print(row.source_id, record["image_status"], record["medium_resolution_status"], flush=True)
write_json(
    old_path,
    dict(
        executed_utc=utcnow(),
        sources=audit,
        comparison_pairs=pairs,
        limitation="All eleven joint priorities, three approximate comparison sources and The Cliff as an explicitly added published benchmark; not a physical classification",
    ),
)

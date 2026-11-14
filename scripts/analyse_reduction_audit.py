"""Compare fixed-peer shapes across DR4 extractions and inspect UV spatial data."""

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib.colors import AsinhNorm
from matplotlib.figure import Figure

from cosmic_anomaly_lab.open_search import band_measurements, normalized_shape, shape_rarity
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import (
    collapse_spatial,
    read_dr4_spectrum,
    signed_spatial_moment,
)
from cosmic_anomaly_lab.spectroscopy import pixel_edges, read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/reduction_audit"
config = yaml.safe_load((root / "configs/reduction_audit.yaml").read_text())
downloads = json.loads((out / "downloads.json").read_text())["targets"]
queue = pd.read_csv(root / "results/open_search/shape_followup.csv").set_index("source_id")
with gzip.open(root / "results/open_search/shape_details.jsonl.gz", "rt", encoding="utf-8") as f:
    details = {(r["source_id"], r["mode"]): r for r in map(json.loads, f)}
with gzip.open(
    root / "results/open_search/band_measurements.jsonl.gz", "rt", encoding="utf-8"
) as f:
    bands = {(r["source_id"], r["mode"]): r for r in map(json.loads, f)}
rng = np.random.default_rng(20261114)
rows, full = [], []
fig = Figure(figsize=(14, 13), layout="constrained")
axes = fig.subplots(4, 2).flat
colors = {(3, 3): "#207a8b", (5, 3): "#a45c34", (3, 2): "#835aa3", (5, 2): "#67934c"}
for product, ax in zip(downloads, axes, strict=True):
    source = product["source_id"]
    target = queue.loc[source]
    original = details[(source, target["mode"])]
    peer_shapes = normalized_shape(
        np.array([bands[(sid, target["mode"])]["flux_uJy"] for sid in original["peer_ids"]])
    )
    edges = np.array(original["edges_um"])
    centers = np.sqrt(edges[:-1] * edges[1:])
    ax.plot(centers, original["normalized_flux"], "k.-", label="DJA selected product", lw=1.5)
    med, spread = np.array(original["peer_median"]), np.array(original["peer_scatter"])
    ax.fill_between(
        centers, med - spread, med + spread, color="#dce2e5", label="Original peer spread"
    )
    product_fig = Figure(figsize=(12, 7), layout="constrained")
    bx, cx = product_fig.subplots(2, 1)
    dja = read_spectrum(root / "data/raw/spectra" / product["dja_file"])
    bx.plot(
        dja.wavelength.value,
        np.where(dja.good, dja.flux.value, np.nan),
        color="#202020",
        lw=0.9,
        label="DJA selected product",
    )
    cx.plot(centers, original["normalized_flux"], "k.-", label="DJA selected product")
    for aperture in [3, 5]:
        for nods in [3, 2]:
            record = dict(
                source_id=source,
                dr4_id=product["dr4_id"],
                mode=target["mode"],
                aperture_pixels=aperture,
                nods=nods,
                status="UNAVAILABLE",
            )
            try:
                if product["status"] != "DOWNLOADED":
                    raise ValueError("No downloaded product")
                path = root / product["path"]
                spec = read_dr4_spectrum(path, aperture, nods)
                header_separation = (
                    SkyCoord(spec.metadata["RA"] * u.deg, spec.metadata["DEC"] * u.deg)
                    .separation(
                        SkyCoord(dja.metadata["SRCRA"] * u.deg, dja.metadata["SRCDEC"] * u.deg)
                    )
                    .arcsec
                )
                assert header_separation <= config["coordinate_tolerance_arcsec"]
                assert int(spec.metadata["PROPOSID"]) == product["pid"]
                assert int(spec.metadata["TARGET"]) == int(product["dr4_id"].rsplit("_", 1)[1])
                measured = band_measurements(spec, target.z_best, edges)
                record.update(
                    status=measured["status"],
                    raw_sha256=digest(path),
                    header_separation_arcsec=float(header_separation),
                    minimum_band_coverage=min(measured["coverage"]),
                    minimum_band_pixel_count=min(measured["pixel_count"]),
                )
                color = colors[(aperture, nods)]
                label = f"DR4 {aperture}px / {nods}-nod"
                bx.plot(
                    spec.wavelength.value,
                    np.where(spec.good, spec.flux.value, np.nan),
                    color=color,
                    lw=0.7,
                    alpha=0.8,
                    label=label,
                )
                if measured["status"] == "MEASURED":
                    median_band_sn = float(
                        np.median(
                            np.array(measured["flux_uJy"])
                            / np.sqrt(np.diag(measured["covariance_uJy2"]))
                        )
                    )
                    fit = shape_rarity(
                        measured["flux_uJy"], measured["covariance_uJy2"], peer_shapes, rng
                    )
                    record.update(
                        score=fit["score"],
                        score_q10=fit["perturbation_q10"],
                        dja_score=original["score"],
                        median_band_sn=median_band_sn,
                        meets_original_band_sn_floor=median_band_sn >= 3,
                        normalized_shape_RMS_difference=float(
                            np.sqrt(
                                np.mean(
                                    (np.array(fit["normalized_flux"]) - original["normalized_flux"])
                                    ** 2
                                )
                            )
                        ),
                    )
                    for panel in [ax, cx]:
                        panel.plot(
                            centers,
                            fit["normalized_flux"],
                            ".-",
                            color=color,
                            lw=1,
                            ms=3,
                            label=label,
                        )
                    full.append(dict(**record, bands=measured, metadata=spec.metadata, rarity=fit))
            except (OSError, ValueError, KeyError, AssertionError) as exc:
                record.update(status="UNAVAILABLE", reason=str(exc))
            rows.append(record)
    for panel in [ax, cx]:
        panel.set(xlabel="Rest wavelength (µm; fixed DJA z)", ylabel="fν / RMS(fν)")
        panel.grid(alpha=0.15)
    ax.set_title(f"{product['dr4_id']}\nDJA z={target.z_best:.3f} · {target['mode']}", fontsize=10)
    if product == downloads[0]:
        ax.legend(fontsize=7)
    bx.set(xlabel="Observed wavelength (µm)", ylabel="Released flux density (µJy)", xlim=(0.7, 5.2))
    bx.legend(fontsize=8, ncol=2)
    bx.grid(alpha=0.15)
    product_fig.suptitle(
        f"{product['dr4_id']} · OBSERVED reduction / aperture comparison\n"
        "Shared photons and systematic errors: curves are not independent detections",
        fontsize=12,
    )
    save_figure(product_fig, root / "figures/reduction_sources" / product["dr4_id"])
    print("Compared", product["dr4_id"], flush=True)
fig.suptitle(
    "OBSERVED cross-reduction shape check · fixed original peers, no population refitting\n"
    "JADES DR4 3/5-pixel and 2/3-nod extractions; no new-class classification",
    fontsize=13,
)
save_figure(fig, root / "figures/reduction_shape_comparison")
pd.DataFrame(rows).to_csv(out / "extraction_comparisons.csv", index=False)
write_json(out / "shape_details.json", full)

# A predeclared, signed spatial audit of the two DJA PRISM products of GS-13905.
spatial = []
uv_files = [
    "gds-deep-v4_prism-clear_1210_13905.spec.fits",
    "jades-gds-wide-v4_prism-clear_1180_13905.spec.fits",
]
sf = Figure(figsize=(13, 9), layout="constrained")
sa = sf.subplots(3, 2)
exposure_sets = []
image_fig = Figure(figsize=(13, 10), layout="constrained")
image_axes = image_fig.subplots(3, 1)
observed_band_rows = []


def record_observed_bands(spec, label):
    for band, bounds in config["spatial_check"]["observed_bands_um"].items():
        measured = band_measurements(spec, 0, bounds)
        row = dict(product=label, band=band, observed_um=bounds, **measured)
        if measured["status"] == "MEASURED":
            row["statistical_error_uJy"] = float(np.sqrt(measured["covariance_uJy2"][0][0]))
        observed_band_rows.append(row)


for col, filename in enumerate(uv_files):
    path = root / "data/raw/spectra" / filename
    with fits.open(path) as h:
        wave = np.asarray(h["SPEC1D"].data["wave"], float)
        science = np.asarray(h["SCI"].data, float)
        weight = np.asarray(h["WHT"].data, float)
        y = np.arange(science.shape[0]) - h["SCI"].header["YTRACE"]
        shown = np.where(np.isfinite(science) & (weight > 0), science, np.nan)
        mesh = image_axes[col].pcolormesh(
            pixel_edges(wave),
            pixel_edges(y),
            shown,
            shading="flat",
            cmap="RdBu_r",
            norm=AsinhNorm(linear_width=0.002, vmin=-0.08, vmax=0.08),
            rasterized=True,
        )
        image_axes[col].set(
            xlim=(0.7, 2.7),
            ylim=(-8, 8),
            ylabel="Trace offset (pixels)",
            title="DJA " + ("deep / PID 1210" if col == 0 else "wide / PID 1180"),
        )
        for boundary in [-3, 3]:
            image_axes[col].axhline(boundary, ls="--", lw=0.6, color="black")
        image_fig.colorbar(
            mesh, ax=image_axes[col], label="SCI (µJy; common asinh stretch)", extend="both"
        )
        dja_spec = read_spectrum(path)
        label = "DJA " + ("deep" if col == 0 else "wide")
        record_observed_bands(dja_spec, label)
        image_axes[2].plot(
            dja_spec.wavelength.value,
            np.where(dja_spec.good, dja_spec.flux.value, np.nan),
            lw=0.9,
            label=label,
        )
        slit = h["SLITS"].data
        names = sorted(set(str(x) for x in slit["filename"]))
        exposure_sets.append(set(names))
        entry = dict(
            file=filename,
            raw_sha256=digest(path),
            slit_files=names,
            mjd_start_range=[float(np.min(slit["start_time"])), float(np.max(slit["start_time"]))],
            pa_v3_range_deg=[float(np.min(slit["pa_v3"])), float(np.max(slit["pa_v3"]))],
            source_xpos_range=[
                float(np.min(slit["source_xpos"])),
                float(np.max(slit["source_xpos"])),
            ],
            source_ypos_range=[
                float(np.min(slit["source_ypos"])),
                float(np.max(slit["source_ypos"])),
            ],
            spatial_bands=[],
        )
        for line, (name, (lo, hi)) in enumerate(
            config["spatial_check"]["observed_bands_um"].items()
        ):
            profile, variance = collapse_spatial(wave, science, weight, lo, hi)
            moments = [
                signed_spatial_moment(y, profile, variance, half_width=radius)
                for radius in config["spatial_check"]["aperture_half_width_pixels"]
            ]
            entry["spatial_bands"].append(
                dict(
                    name=name,
                    observed_um=[lo, hi],
                    moments=moments,
                    y_pixels=y.tolist(),
                    profile_uJy=[float(v) if np.isfinite(v) else None for v in profile],
                    variance_uJy2=[float(v) if np.isfinite(v) else None for v in variance],
                )
            )
            panel = sa[line, col]
            panel.errorbar(
                y, profile, yerr=np.sqrt(variance), fmt="o-", ms=3, lw=0.8, color="#286f88"
            )
            panel.axvspan(-3, 3, color="#e4ccb6", alpha=0.3)
            panel.axhline(0, color="grey", lw=0.6)
            visible = (abs(y) <= 8) & np.isfinite(profile + variance)
            lower = min(0, float(np.min((profile - np.sqrt(variance))[visible])))
            upper = max(0, float(np.max((profile + np.sqrt(variance))[visible])))
            margin = max((upper - lower) * 0.08, 1e-6)
            panel.set(
                xlim=(-8, 8),
                ylim=(lower - margin, upper + margin),
                xlabel="Row offset from released trace (pixels)",
                ylabel="Mean SCI (µJy)",
                title=f"{'Deep / PID 1210' if col == 0 else 'Wide / PID 1180'} · {name.replace('_', ' ')} {lo}–{hi} µm",
            )
            panel.grid(alpha=0.15)
        spatial.append(entry)
sf.suptitle(
    "GS-13905 · signed spatial profiles in the released DJA 2D products\n"
    "Negative nod/background pixels retained; shaded aperture ±3 pixels; errors omit resampling covariance",
    fontsize=13,
)
save_figure(sf, root / "figures/uv_spatial_audit")
for product in downloads:
    if product["source_id"] != config["spatial_check"]["source_id"]:
        continue
    for aperture in [3, 5]:
        for nods in [3, 2]:
            spec = read_dr4_spectrum(root / product["path"], aperture, nods)
            label = f"{product['dr4_id']} {aperture}px/{nods}-nod"
            record_observed_bands(spec, label)
            if aperture == 3 and nods == 3:
                image_axes[2].plot(
                    spec.wavelength.value,
                    np.where(spec.good, spec.flux.value, np.nan),
                    lw=1,
                    label=label,
                )
image_axes[2].set(xlim=(0.7, 2.7), ylim=(-0.04, 0.35), ylabel="Released 1D flux density (µJy)")
image_axes[2].legend(fontsize=8)
for panel in image_axes:
    panel.axvline(
        0.121567 * (1 + queue.loc[config["spatial_check"]["source_id"], "z_best"]),
        ls=":",
        color="#855700",
        lw=1,
    )
    panel.set_xlabel("Observed wavelength (µm)")
image_fig.suptitle(
    "GS-13905 · OBSERVED rectified spectral images and released extractions\n"
    "Dotted line: Lyα at catalogue z; 2D rows are not sky-registered between visits",
    fontsize=12,
)
save_figure(image_fig, root / "figures/uv_reduction_images")
write_json(out / "uv_observed_bands.json", observed_band_rows)
write_json(
    out / "spatial_audit.json",
    dict(
        executed_utc=utcnow(),
        source_id=config["spatial_check"]["source_id"],
        products=spatial,
        shared_slit_filenames=sorted(exposure_sets[0] & exposure_sets[1]),
        interpretation="Different listed input exposures, slit placements and roll angles; "
        "not independent reductions, and not a controlled variability experiment",
    ),
)
summary = dict(
    executed_utc=utcnow(),
    associated_products=len(downloads),
    associated_sources=len(set(x["source_id"] for x in downloads)),
    extraction_variants=len(rows),
    scorable_variants=sum("score" in r for r in rows),
    variants_meeting_original_band_sn_floor=sum(
        r.get("meets_original_band_sn_floor", False) for r in rows
    ),
    original_sample_refit=False,
    physical_classification="NOT_ESTABLISHED",
    input_sha256={
        str(p.relative_to(root)).replace("\\", "/"): text_digest(p)
        for p in [
            root / "configs/reduction_audit.yaml",
            root / "src/cosmic_anomaly_lab/reduction_audit.py",
            Path(__file__),
            out / "downloads.json",
            root / "src/cosmic_anomaly_lab/open_search.py",
            root / "results/open_search/shape_followup.csv",
            root / "scripts/download_reduction_audit.py",
        ]
    },
    binary_input_sha256={
        str(p.relative_to(root)).replace("\\", "/"): digest(p)
        for p in [
            root / "results/open_search/shape_details.jsonl.gz",
            root / "results/open_search/band_measurements.jsonl.gz",
        ]
    },
)
write_json(out / "summary.json", summary)
print(summary)

"""Observed image/evidence atlas. Display transforms do not classify sources."""

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.colors import AsinhNorm
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
from matplotlib.ticker import FuncFormatter

from .imaging import aperture_measurement
from .provenance import digest, text_digest, utcnow, write_json
from .spectroscopy import continuum_break, pixel_edges, read_spectrum

matplotlib.use("Agg")

BANDS = ("F150W-CLEAR", "F277W-CLEAR", "F444W-CLEAR")


def check_registration(reference, other, shape, tolerance_pixels=0.05):
    """Test all pixel centres, including distortion, before composing channels."""
    yy, xx = np.indices(shape)
    ra, dec = reference.pixel_to_world_values(xx, yy)
    ox, oy = other.world_to_pixel_values(ra, dec)
    displacement = np.hypot(ox - xx, oy - yy)
    if not np.all(np.isfinite(displacement)) or displacement.max() > tolerance_pixels:
        raise ValueError("Band WCS grids differ; reprojection required before RGB")
    return float(displacement.max())


def colour_composite(channels, valid, scale, q=8.0):
    """One intensity stretch for calibrated R,G,B; preserve positive colour ratios.

    Clip negative flux ONLY for display. Saturation divides all channels by the
    same pixel maximum. Missing pixels are grey, never an invented detection.
    """
    channels = np.asarray(channels, dtype=float)
    if channels.ndim != 3 or channels.shape[-1] != 3 or valid.shape != channels.shape[:2]:
        raise ValueError("Expected H x W x 3 channels and H x W mask")
    if not np.isfinite(scale) or scale <= 0 or not np.isfinite(q) or q <= 0:
        raise ValueError("Invalid RGB stretch")
    positive = np.maximum(np.where(valid[..., None], channels, 0), 0)
    if not np.all(np.isfinite(positive)):
        raise ValueError("Nonfinite valid RGB input")
    intensity = positive.mean(axis=-1)
    factor = np.divide(
        np.arcsinh(q * intensity / scale) / np.arcsinh(q),
        intensity,
        out=np.zeros_like(intensity),
        where=intensity > 0,
    )
    rgb = positive * factor[..., None]
    rgb /= np.maximum(1, rgb.max(axis=-1))[..., None]
    rgb[~valid] = 0.22
    return rgb


def read_images(path, ra, dec):
    """Calibrate from PHOTFNU and mask nonfinite/nonpositive-weight pixels."""
    images, details, apertures = {}, {}, []
    with fits.open(path) as hdus:
        hdus.verify("exception")
        for hdu in hdus:
            if hdu.header.get("EXTVER") != "SCI":
                continue
            band = hdu.header["FILTER"]
            if band not in BANDS:
                continue
            if band in images:
                raise ValueError(f"Duplicate science band {band}")
            weight_hdus = [
                h
                for h in hdus
                if h.header.get("EXTVER") == "WHT" and h.header.get("FILTER") == band
            ]
            if len(weight_hdus) != 1:
                raise ValueError(f"Missing/ambiguous weight for {band}")
            data, weight = np.array(hdu.data, float), np.array(weight_hdus[0].data, float)
            conversion = float(hdu.header["PHOTFNU"]) * 1e6
            if (
                data.ndim != 2
                or data.shape != weight.shape
                or not np.isfinite(conversion)
                or conversion <= 0
            ):
                raise ValueError("Invalid image shape or PHOTFNU")
            wcs = WCS(hdu.header).celestial
            valid = np.isfinite(data) & np.isfinite(weight) & (weight > 0)
            x0, y0 = wcs.world_to_pixel_values(ra, dec)
            scales = proj_plane_pixel_scales(wcs) * 3600
            yy, xx = np.indices(data.shape)
            radius = np.hypot((xx - x0) * scales[0], (yy - y0) * scales[1])
            sky = sigma_clip(data[valid & (radius > 1.5) & (radius < 2.5)], sigma=3).compressed()
            if len(sky) < 20:
                reason = "Insufficient valid pixels in display background annulus"
                details[band] = dict(
                    status="UNAVAILABLE", reason=reason, valid_fraction=float(valid.mean())
                )
                apertures.append(dict(filter=band, status="UNAVAILABLE", reason=reason))
                continue
            background = float(np.median(sky) * conversion)
            calibrated = np.where(valid, data * conversion - background, np.nan)
            images[band] = dict(data=calibrated, valid=valid, wcs=wcs, radius=radius)
            details[band] = dict(
                photfnu_Jy_per_unit=conversion / 1e6,
                display_background_uJy_pixel=background,
                valid_fraction=float(valid.mean()),
                pixel_scale_arcsec=scales.tolist(),
                shape=list(data.shape),
            )
            try:
                result = aperture_measurement(
                    data, weight, wcs, ra, dec, jy_per_unit=conversion / 1e6
                )
                apertures.append(dict(filter=band, status="MEASURED", **result))
            except ValueError as exc:
                apertures.append(dict(filter=band, status="UNAVAILABLE", reason=str(exc)))
    if not images:
        raise ValueError("No requested SCI bands")
    reference = next(iter(images.values()))
    for band, im in images.items():
        if im["data"].shape != reference["data"].shape:
            raise ValueError("Band shapes differ; reprojection required")
        details[band]["registration_max_offset_pixels"] = check_registration(
            reference["wcs"], im["wcs"], im["data"].shape
        )
    rgb, scale = None, None
    if all(b in images for b in BANDS):
        channels = np.stack([images[b]["data"] for b in BANDS[::-1]], axis=-1)
        valid = np.logical_and.reduce([images[b]["valid"] for b in BANDS])
        intensity = np.maximum(channels, 0).mean(axis=-1)
        central = valid & (reference["radius"] < 1) & (intensity > 0)
        if central.sum() >= 20:
            scale = float(np.percentile(intensity[central], 99.7))
            rgb = colour_composite(channels, valid, scale)
    return images, rgb, dict(bands=details, rgb_scale_uJy_pixel=scale, apertures=apertures)


def mark_image(ax, im, ra, dec):
    wcs, data = im["wcs"], im["data"]
    x, y = wcs.world_to_pixel_values(ra, dec)
    scales = proj_plane_pixel_scales(wcs) * 3600
    ax.add_patch(
        Ellipse((x, y), 0.8 / scales[0], 0.8 / scales[1], fill=False, edgecolor="#66e0ed", lw=0.8)
    )
    ax.plot([8, 8 + 1 / scales[0]], [8, 8], color="white", lw=2)
    ax.text(8, 11, "1″", color="white", fontsize=9)
    ax.set(xlim=(-0.5, data.shape[1] - 0.5), ylim=(-0.5, data.shape[0] - 0.5), xticks=[], yticks=[])


def image_panel(ax, images, rgb, row, band=None):
    ax.set_box_aspect(1)
    if band is None and rgb is not None:
        ax.imshow(rgb, origin="lower", interpolation="nearest")
        mark_image(ax, next(iter(images.values())), row.ra, row.dec)
        return
    if band in images:
        im = images[band]
        values = im["data"][im["valid"]]
        high = max(float(np.percentile(values, 99.7)), 1e-8)
        cmap = matplotlib.colormaps["gray"].with_extremes(bad="#383838")
        ax.imshow(
            im["data"],
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            norm=AsinhNorm(linear_width=high / 12, vmin=-high / 20, vmax=high),
        )
        mark_image(ax, im, row.ra, row.dec)
    else:
        ax.set_facecolor("#e9edf0")
        ax.text(
            0.5,
            0.5,
            "RGB unavailable\nincomplete usable bands"
            if band is None
            else "No usable band coverage",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10,
        )
        ax.set(xticks=[], yticks=[])


def save_figure(fig, base):
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=170, facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white")
    fig.clear()


def evidence_sheet(root, row, images, rgb, output_directory="figures/sources"):
    path = root / "data/raw/spectra" / row.file
    spec = read_spectrum(path)
    wave, flux, error = spec.wavelength.value, spec.flux.value, spec.error.value
    continuum = continuum_break(spec, row.z_best)
    fig = Figure(figsize=(13, 10.4), layout="constrained")
    grid = fig.add_gridspec(3, 4, height_ratios=[1.15, 1, 0.65])
    for i, band in enumerate([None, *BANDS]):
        ax = fig.add_subplot(grid[0, i])
        image_panel(ax, images, rgb, row, band)
        ax.set_title("NIRCam · false colour" if band is None else band.split("-")[0], fontsize=11)
    ax = fig.add_subplot(grid[1, :3])
    goodflux = np.where(spec.good, flux, np.nan)
    ax.fill_between(
        wave, goodflux - error, goodflux + error, color="#ccdce5", label="±1σ statistical"
    )
    ax.plot(wave, goodflux, lw=0.9, color="#173d57", label="Released 1D spectrum")
    ax.axhline(0, color="#b6b6b6", lw=0.6)
    ax.set(
        xlim=(0.7, 5.2),
        xlabel="Observed wavelength (µm)",
        ylabel="Flux density (µJy)",
        title="NIRSpec / PRISM · data and uncertainty",
    )
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.15)
    ax = fig.add_subplot(grid[1, 3])
    rest = spec.rest_wavelength(row.z_best).value
    mask = spec.good & (rest > 0.30) & (rest < 0.46)
    ax.axvspan(0.34, 0.36, color="#b8dce8", alpha=0.65)
    ax.axvspan(0.40, 0.42, color="#ecc1b6", alpha=0.65)
    if mask.sum() >= 3:
        ax.plot(rest, np.where(mask, flux, np.nan), lw=0.9, color="#173d57")
        ax.fill_between(
            rest,
            np.where(mask, flux - error, np.nan),
            np.where(mask, flux + error, np.nan),
            color="#ccdce5",
        )
    ratio = continuum.get("ratio")
    label = (
        f"Red / blue = {ratio:.2f} ± {continuum['ratio_error']:.2f}"
        if ratio is not None
        else continuum["status"].replace("_", " ")
    )
    ax.set(
        xlim=(0.30, 0.46),
        xlabel="Rest wavelength (µm)",
        ylabel="Observed fν (µJy)",
        title="Continuum windows",
    )
    ax.text(0.04, 0.96, label, transform=ax.transAxes, fontsize=8, va="top")
    ax.grid(alpha=0.15)
    ax = fig.add_subplot(grid[2, :])
    with fits.open(path) as hdus:
        sci = hdus["SCI"]
        data = np.asarray(sci.data, float) * u.Unit(sci.header["BUNIT"]).to(u.uJy)
        weight = np.asarray(hdus["WHT"].data, float)
        if data.shape != weight.shape or data.shape[1] != wave.size:
            raise ValueError("2D spectrum and wavelength grid disagree")
        valid = np.isfinite(data) & np.isfinite(weight) & (weight > 0)
        valid &= ((wave >= 0.7) & (wave <= 5.2))[None, :]
        trace = float(sci.header["YTRACE"])
        values = data[valid]
        high = max(float(np.percentile(values, 99.5)), 1e-8)
        mesh = ax.pcolormesh(
            pixel_edges(wave),
            np.arange(data.shape[0] + 1) - 0.5 - trace,
            np.where(valid, data, np.nan),
            cmap="magma",
            norm=AsinhNorm(linear_width=high / 12, vmin=-high / 12, vmax=high),
            rasterized=True,
            shading="flat",
        )
        fig.colorbar(
            mesh,
            ax=ax,
            label="Released SCI (µJy)",
            fraction=0.025,
            pad=0.015,
            ticks=[0, high / 10, high],
            format=FuncFormatter(lambda x, pos: f"{x:.3g}"),
        )
    ax.set(
        xlim=(0.7, 5.2),
        xlabel="Observed wavelength (µm)",
        ylabel="Pixels from trace",
        title="Observed 2D spectral image · resampled wavelength grid; display stretch only",
    )
    fig.suptitle(
        f"OBSERVED EVIDENCE  |  {row.survey} · {row.srcid}  |  z = {row.z_best:.4f}\n"
        f"{row.source_id}   •   RA {row.ra:.6f}°  Dec {row.dec:.6f}°\n"
        "RGB = F444W / F277W / F150W. Circle: catalogue position (0.4″ radius). Native image orientation.\n"
        "PSFs unmatched; mono panels stretched independently. No physical classification established.",
        fontsize=11,
    )
    base = root / output_directory / row.source_id
    save_figure(fig, base)
    return dict(
        continuum=continuum,
        spectrum_sha256=digest(path),
        evidence_png=base.relative_to(root).as_posix() + ".png",
    )


def build_atlas(root):
    root = Path(root).resolve()
    selection = root / "data/catalogues/pilot_selection.csv"
    rows = pd.read_csv(selection)
    records, prepared = [], []
    for row in rows.itertuples():
        path = root / "data/raw/images" / f"{row.source_id}.fits"
        try:
            images, rgb, detail = read_images(path, row.ra, row.dec)
            status = "RGB_AVAILABLE" if rgb is not None else "PARTIAL_BAND_COVERAGE"
        except (OSError, ValueError, KeyError) as exc:
            images, rgb, detail = {}, None, dict(image_error=str(exc))
            status = "IMAGE_UNAVAILABLE"
        sheet = evidence_sheet(root, row, images, rgb)
        records.append(
            dict(
                source_id=row.source_id,
                survey=row.survey,
                source_number=int(row.srcid),
                redshift=float(row.z_best),
                image_status=status,
                image_sha256=digest(path) if path.exists() else None,
                **detail,
                **sheet,
            )
        )
        prepared.append((row, images, rgb))
        print(f"{row.source_id}: {status}", flush=True)
    fig = Figure(figsize=(13, 16), layout="constrained")
    axes = fig.subplots(5, 4)
    for ax, (row, images, rgb) in zip(axes.flat, prepared, strict=True):
        fallback = "F444W-CLEAR" if rgb is None else None
        image_panel(ax, images, rgb, row, fallback)
        ax.set_title(
            f"{row.survey} · {row.srcid}  |  z {row.z_best:.2f}\n"
            + ("F444W ONLY — RGB unavailable" if rgb is None else row.source_id),
            fontsize=8,
        )
    fig.suptitle(
        "OBSERVED JWST SOURCE ATLAS\n"
        "R / G / B = F444W / F277W / F150W  •  1″ scale bars  •  catalogue positions circled\n"
        "Per-source intensity stretch; brightness is not comparable between panels. PSFs unmatched.",
        fontsize=12,
    )
    save_figure(fig, root / "figures/observed_rgb_atlas")
    selected = [
        x
        for x in prepared
        if x[0].source_id
        in ["CAL-3539b1d648b05582", "CAL-3c82037c6ff15eda", "CAL-cf2762adfeb35e39"]
    ]
    fig = Figure(figsize=(13, 4.9), layout="constrained")
    for ax, (row, images, rgb) in zip(fig.subplots(1, 3), selected, strict=True):
        image_panel(ax, images, rgb, row)
        label = "Published LRD/AGN benchmark" if row.sample == "A" else "Comparison source"
        ax.set_title(f"{row.survey} · {row.srcid}  |  z {row.z_best:.2f}\n{label}", fontsize=11)
    fig.suptitle(
        "ACTUAL NIRCAM IMAGES · FALSE COLOUR\n"
        "F444W → red   F277W → green   F150W → blue  |  1″ scale bars\n"
        "Colour and apparent compactness do not establish a black-hole star. PSFs unmatched.",
        fontsize=12,
    )
    save_figure(fig, root / "figures/observed_image_comparison")
    manifest = dict(
        created_utc=utcnow(),
        selection_sha256=text_digest(selection),
        renderer_sha256=text_digest(Path(__file__)),
        observed_sources=len(records),
        rgb_sources=sum(r["image_status"] == "RGB_AVAILABLE" for r in records),
        measured_apertures=sum(
            a["status"] == "MEASURED" for r in records for a in r.get("apertures", [])
        ),
        method="PHOTFNU calibration; clipped median sky at 1.5–2.5 arcsec; common RGB intensity asinh Q=8; "
        "scale=99.7 percentile of positive mean RGB intensity inside 1 arcsec; "
        "negative pixels clipped for RGB display only; grey=invalid; native WCS grid checked at all pixels",
        limitations=[
            "No PSF matching or deblending",
            "Per-source display scaling, not comparable brightness",
            "2D image is a visual diagnostic, not an independent extraction or S/N significance",
            "Aperture statistics omit drizzle covariance and aperture corrections",
            "Physical model comparison not performed by this renderer",
        ],
        sources=records,
    )
    write_json(root / "results/visual_atlas_run.json", manifest)
    gallery = [
        "# Observed image and spectrum gallery",
        "",
        f"{len(records)} observed source evidence sheets; {manifest['rgb_sources']} three-band RGB cutouts. "
        "These are telescope data, not illustrations or simulated sources.",
        "",
        "![Observed source comparison](../figures/observed_image_comparison.png)",
        "",
        "Open an individual sheet below for the three filter images, the 1D spectrum with errors, "
        "the continuum windows, and the actual 2D spectral image.",
        "",
        "**No source here is established as a black-hole star by this analysis.** "
        "Read the [physical-test status](black_hole_star_status.md) and "
        "[rendering methods and limitations](../docs/visual_evidence.md).",
        "",
        "[Full atlas PNG](../figures/observed_rgb_atlas.png) · "
        "[Atlas PDF](../figures/observed_rgb_atlas.pdf) · "
        "[Run manifest / measurements](../results/visual_atlas_run.json)",
        "",
        "| Source | Archive z | Image coverage | Evidence |",
        "| --- | --- | --- | --- |",
    ]
    for record in records:
        sid = record["source_id"]
        gallery.append(
            f"| {record['survey']} {record['source_number']} | {record['redshift']:.4f} | "
            f"{record['image_status'].replace('_', ' ')} | "
            f"[PNG](../figures/sources/{sid}.png) · [PDF](../figures/sources/{sid}.pdf) · "
            f"[Report](sources/{sid}.md) |"
        )
        source_report = root / "reports/sources" / f"{sid}.md"
        if source_report.exists():
            marker = "<!-- observed-atlas -->"
            original = source_report.read_text(encoding="utf-8").split(marker)[0].rstrip()
            source_report.write_text(
                original + f"\n\n{marker}\n"
                f"[Observed image + spectrum sheet](../../figures/sources/{sid}.png) · "
                f"[PDF](../../figures/sources/{sid}.pdf) · [Gallery](../visual_gallery.md)\n",
                encoding="utf-8",
            )
    gallery += ["", "## Individual evidence sheets", ""]
    for record in sorted(records, key=lambda r: r["source_id"] != "CAL-cf2762adfeb35e39"):
        sid = record["source_id"]
        gallery += [
            f"### {record['survey']} {record['source_number']}",
            "",
            f"![Observed evidence for {sid}](../figures/sources/{sid}.png)",
            "",
        ]
    (root / "reports/visual_gallery.md").write_text("\n".join(gallery), encoding="utf-8")
    return manifest

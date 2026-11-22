"""Measure the selected GOODS-S candidate at fixed sky coordinates."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from scipy.ndimage import maximum_filter

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jwst_epoch_search"
FIGURES = ROOT / "figures/jwst_epoch_search"
RADII = (.12, .18, .24)


def read_receipt(name):
    return json.loads((OUT / name).read_text())


def verified_path(product):
    path = ROOT / product["path"]
    if digest(path) != product["sha256"]:
        raise ValueError(f"Archive checksum mismatch: {path.name}")
    return path


def pixel_scale(wcs):
    return float(np.mean(proj_plane_pixel_scales(wcs)) * 3600)


def aperture(image, x, y, scale, radius, annulus):
    half = int(np.ceil(annulus[1] / scale)) + 1
    ix, iy = round(x), round(y)
    if min(ix, iy, image.shape[1]-1-ix, image.shape[0]-1-iy) < half:
        return np.nan
    stamp = image[iy-half:iy+half+1, ix-half:ix+half+1]
    yy, xx = np.ogrid[iy-half:iy+half+1, ix-half:ix+half+1]
    distance = np.hypot(xx-x, yy-y) * scale
    source = distance <= radius + 1e-8
    sky = (distance >= annulus[0]) & (distance <= annulus[1])
    if not np.isfinite(stamp[source]).all() or np.isfinite(stamp[sky]).mean() < .8:
        return np.nan
    return float(np.sum(stamp[source] - np.nanmedian(stamp[sky])))


def photometry(image, x, y, scale, radii, annulus):
    step = round(.56 / scale)
    offsets = np.arange(-int(6 / scale), int(6 / scale)+1, step)
    controls = [(dx, dy) for dy in offsets for dx in offsets
                if 1.5 <= np.hypot(dx, dy) * scale <= 6]
    result = []
    for radius in radii:
        blank = np.array([aperture(image, x+dx, y+dy, scale, radius, annulus) for dx, dy in controls])
        blank = blank[np.isfinite(blank)]
        flux = aperture(image, x, y, scale, radius, annulus)
        result.append({
            "radius_arcsec": radius,
            "flux_njy": flux if np.isfinite(flux) else None,
            "local_mad_sigma_njy": float(1.4826*np.median(np.abs(blank-np.median(blank)))) if len(blank) >= 40 else None,
            "control_median_njy": float(np.median(blank)) if len(blank) >= 40 else None,
            "control_count": len(blank),
        })
    return result


def load_mosaic(record, sky, half=1050):
    with fits.open(verified_path(record["sci"]), memmap=False) as h:
        header = h[0].header
        wcs = WCS(header).celestial
        x, y = map(float, wcs.world_to_pixel_values(*sky))
        x0, y0 = round(x)-half, round(y)-half
        if min(x0, y0) < 0 or x0+2*half >= header["NAXIS1"] or y0+2*half >= header["NAXIS2"]:
            raise ValueError("Candidate is too close to mosaic array boundary for reference cutout")
        region = np.s_[y0:y0+2*half+1, x0:x0+2*half+1]
        image = np.array(h[0].data[region], dtype=np.float32) * (float(header["PHOTFNU"]) * 1e9)
    with fits.open(verified_path(record["wht"]), memmap=False) as h:
        weight = h[0].data[region]
        image[~np.isfinite(weight) | (weight <= 0)] = np.nan
    return image, wcs, x-x0, y-y0, x0, y0


def show_stamp(ax, image, x, y, scale, bounds, cmap, radius, half_arcsec=.72):
    half = round(half_arcsec / scale)
    ix, iy = round(x), round(y)
    stamp = image[iy-half:iy+half+1, ix-half:ix+half+1]
    extent = np.array([ix-half-.5-x, ix+half+.5-x, iy-half-.5-y, iy+half+.5-y])*scale
    ax.imshow(stamp, origin="lower", extent=extent, cmap=cmap, vmin=bounds[0], vmax=bounds[1],
              interpolation="nearest")
    ax.add_patch(Circle((0, 0), radius, fill=False, color="#4b6270" if cmap == "RdBu_r" else "#69e2ec",
                        linewidth=1.1, linestyle="--"))
    ax.set(xticks=[], yticks=[])


def format_flux(row):
    flux, sigma = row["flux_njy"], row["local_mad_sigma_njy"]
    return "Insufficient coverage" if flux is None or sigma is None else f"{flux:+.1f} ± {sigma:.1f} nJy"


def flux_change(older, newer):
    fluxes = [m["flux_njy"] for m in (older, newer)]
    scatters = [m["local_mad_sigma_njy"] for m in (older, newer)]
    return {
        "radius_arcsec": older["radius_arcsec"],
        "delta_njy": fluxes[1]-fluxes[0] if None not in fluxes else None,
        "quadrature_local_scatter_njy": float(np.hypot(*scatters)) if None not in scatters else None,
    }


def centroids(image, xy, core, search=2):
    result = []
    half = round(core)+search
    for x, y in xy:
        ix, iy = round(x), round(y)
        if min(ix, iy, image.shape[1]-ix-1, image.shape[0]-iy-1) < half:
            result.append([np.nan, np.nan])
            continue
        yy, xx = np.mgrid[iy-half:iy+half+1, ix-half:ix+half+1]
        stamp = image[yy, xx]
        region = (abs(xx-x) <= search) & (abs(yy-y) <= search) & np.isfinite(stamp)
        if not np.any(region):
            result.append([np.nan, np.nan])
            continue
        peak = np.unravel_index(np.argmax(np.where(region, stamp, -np.inf)), region.shape)
        distance = np.hypot(xx-xx[peak], yy-yy[peak])
        local = (distance <= core) & np.isfinite(stamp)
        sky = (distance > core) & np.isfinite(stamp)
        if sky.sum() < 10 or local.sum() < 5:
            result.append([np.nan, np.nan])
            continue
        weights = np.maximum(stamp[local]-np.median(stamp[sky]), 0)
        result.append([np.average(xx[local], weights=weights), np.average(yy[local], weights=weights)]
                      if weights.sum() > 0 else [np.nan, np.nan])
    return np.array(result)


def reference_sky(mosaic):
    image, wcs, x, y, x0, y0 = mosaic
    yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
    sep = np.hypot(xx-x, yy-y)*pixel_scale(wcs)
    clean = np.where(np.isfinite(image), image, -np.inf)
    peaks = (clean == maximum_filter(clean, size=23)) & (clean > 4) & (sep > 2) & (sep < 15)
    py, px = np.nonzero(peaks)
    selected = np.argsort(clean[py, px])[-32:]
    centers = centroids(image, np.column_stack([px[selected], py[selected]]), 4.5)
    centers = centers[np.isfinite(centers).all(axis=1)]
    return wcs.pixel_to_world_values(centers[:, 0]+x0, centers[:, 1]+y0)


def load_dithers(receipt, refs):
    frames = []
    for record in receipt["products"]:
        with fits.open(verified_path(record), memmap=False) as h:
            header = h["SCI"].header
            if header["BUNIT"] != "MJy/sr":
                raise ValueError("Unexpected detector image units")
            scale_factor = float(header["PIXAR_SR"])*1e15
            wcs = WCS(header).celestial
            dq = h["DQ"].data
            image = h["SCI"].data*scale_factor
        # Jumps are excluded from the astrometric fit. Retain their science
        # values in a separate array so detector-locked structure is visible.
        clean = np.where(dq == 0, image, np.nan)
        predicted = np.column_stack(wcs.world_to_pixel_values(*refs))
        shifts = centroids(clean, predicted, 2.8)-predicted
        shifts = shifts[np.isfinite(shifts).all(axis=1)]
        if len(shifts):
            shifts = shifts[np.linalg.norm(shifts-np.median(shifts, axis=0), axis=1) < .85]
        if len(shifts) < 5:
            raise ValueError(f"Too few local reference matches for {record['file']}: {len(shifts)}")
        translation = np.median(shifts, axis=0)
        frames.append({
            "record": record, "wcs": wcs, "scale": pixel_scale(wcs), "image": image,
            "dq": dq, "clean": clean, "translation": translation, "reference_count": len(shifts),
            "residual": float(np.sqrt(np.mean(np.sum((shifts-translation)**2, axis=1)))),
        })
    return frames


def native_peak(image, x, y, scale):
    half = int(np.ceil(.36 / scale)) + 1
    ix, iy = round(x), round(y)
    if min(ix, iy, image.shape[1]-1-ix, image.shape[0]-1-iy) < half:
        return None
    yy, xx = np.mgrid[iy-half:iy+half+1, ix-half:ix+half+1]
    stamp = image[yy, xx]
    window = (np.hypot(xx-x, yy-y)*scale <= .36) & np.isfinite(stamp)
    if not window.any():
        return None
    index = np.unravel_index(np.argmax(np.where(window, stamp, -np.inf)), stamp.shape)
    return int(xx[index]), int(yy[index])


def audit_detector_target(frames, sky):
    rows = []
    for frame in frames:
        wcs, scale, translation = (frame[k] for k in ("wcs", "scale", "translation"))
        record, image, dq, clean = (frame[k] for k in ("record", "image", "dq", "clean"))
        x, y = np.array(wcs.world_to_pixel_values(*sky)) + translation
        measurements = photometry(clean, x, y, scale, (.16, .22), (.32, .52))
        # A native-pixel extremum is a diagnostic, not a source detection.
        peak = native_peak(image, x, y, scale)
        px, py = peak if peak else (None, None)
        flux_peak = aperture(image, x, y, scale, .16, (.32, .52))
        rows.append({
            "file": record["file"], "mjd": record["expstart_mjd"], "usable_reference_count": frame["reference_count"],
            "offset_pixels_xy": translation.tolist(), "forced_center_pixels_xy": [float(x), float(y)],
            "reference_residual_rms_pixels": frame["residual"], "unmasked_peak_pixels_xy": peak,
            "unmasked_peak_dq": int(dq[py, px]) if peak else None,
            "unmasked_peak_distance_arcsec": float(np.hypot(px-x, py-y)*scale) if peak else None,
            "apertures_dq_zero_only": measurements,
            "diagnostic_flux_including_flagged_pixels_njy": flux_peak if np.isfinite(flux_peak) else None,
        })
    return rows


def plot_detectors(groups):
    figure = Figure(figsize=(8.4, 8.9), layout="constrained")
    axes = figure.subplots(len(groups), 3, squeeze=False)
    for ax_row, (band, frames, measurements) in zip(axes, groups, strict=True):
        centers = np.array([m["forced_center_pixels_xy"] for m in measurements])
        ix, iy = np.rint(np.mean(centers, axis=0)).astype(int)
        half = 7
        for index, (ax, frame, values) in enumerate(zip(ax_row, frames, measurements, strict=True), 1):
            image = frame["image"][iy-half:iy+half+1, ix-half:ix+half+1]
            dq = frame["dq"][iy-half:iy+half+1, ix-half:ix+half+1]
            extent = (ix-half-.5, ix+half+.5, iy-half-.5, iy+half+.5)
            floor = np.nanmedian(image)
            ax.imshow(image-floor, origin="lower", extent=extent, cmap="gray", vmin=-2.8, vmax=8.0,
                      interpolation="nearest")
            x, y = values["forced_center_pixels_xy"]
            ax.add_patch(Circle((x, y), .16/frame["scale"], fill=False, edgecolor="#36ddeb",
                                linestyle="--", linewidth=1.35))
            bad_y, bad_x = np.nonzero((dq & 4) != 0)
            ax.scatter(ix-half+bad_x, iy-half+bad_y, marker=".", color="#ffbd55", s=8, alpha=.7)
            peak = values["unmasked_peak_pixels_xy"]
            if peak:
                ax.scatter(*peak, marker="+", color="#ff637d", s=76, linewidths=1.45)
            ax.set(xticks=[], yticks=[], title=f"{band} · exposure {index}")
    figure.suptitle("GOODS-S brightening on native NIRCam pixels\n"
                    "Native pixels · cyan = fixed sky aperture · pink = brightest pixel\n"
                    "Amber = jump flag · same detector coordinates within each row", fontsize=12)
    return figure


def summarize_detector(rows):
    centers = np.array([r["forced_center_pixels_xy"] for r in rows])
    motion = float(np.max(np.linalg.norm(centers[:, None]-centers, axis=-1)))
    peaks = [r["unmasked_peak_pixels_xy"] for r in rows]
    spread = None
    if all(p is not None for p in peaks):
        pixels = np.array(peaks)
        spread = float(np.max(np.linalg.norm(pixels[:, None]-pixels, axis=-1)))
    return {
        "sky_track_span_pixels": motion, "native_peak_span_pixels": spread,
        "jump_flagged_peaks": sum(bool((r["unmasked_peak_dq"] or 0) & 4) for r in rows),
        "exposures": rows,
    }


def fading_shortlist():
    table = pd.read_csv(OUT / "matched_peaks.csv")
    selected = table[(table.sign == -1) & table.compactness_gate]
    for band in ("f356w", "f444w"):
        old, new = selected[f"{band}_old_njy"], selected[f"{band}_new_njy"]
        selected = selected[(old > 0) & (old-new >= .25*old)]
    return selected.sort_values("minimum_absolute_score", ascending=False)


def measure_fade(row, older_frames, stamps, products, axes):
    sky = row.ra, row.dec
    result = {"candidate": row.candidate, "ra": row.ra, "dec": row.dec,
              "matched_score": row.minimum_absolute_score,
              "detector_audit": summarize_detector(audit_detector_target(older_frames, sky)),
              "mosaics": []}
    stamp_data = []
    for band, year in (("F444W", 2023), ("F444W", 2025), ("F444W", 2026),
                       ("F356W", 2025), ("F356W", 2026)):
        image, wcs, _, _, x0, y0 = stamps[(band, year)]
        x, y = map(float, wcs.world_to_pixel_values(*sky))
        if min(x-x0, y-y0, image.shape[1]-1-x+x0, image.shape[0]-1-y+y0) < 6.5/pixel_scale(wcs):
            image, wcs, _, _, x0, y0 = load_mosaic(products[(band, year)], sky, half=175)
            x, y = map(float, wcs.world_to_pixel_values(*sky))
        x, y = x-x0, y-y0
        scale = pixel_scale(wcs)
        meas = photometry(image, x, y, scale, RADII, (.28, .44))
        cx, cy = centroids(image, [(x, y)], 4.5, search=6)[0]
        result["mosaics"].append({
            "band": band, "year": year, "apertures": meas,
            "local_centroid_offset_arcsec_xy": [(cx-x)*scale, (cy-y)*scale] if np.isfinite(cx+cy) else None,
            "recentered_024arcsec_flux_njy": (
                m if np.isfinite(cx+cy) and np.isfinite(m := aperture(image, cx, cy, scale, .24, (.28, .44))) else None
            ),
        })
        if band == "F444W":
            stamp_data.append((year, image, x, y, scale, meas[2]))
    local = np.concatenate([i[round(y)-5:round(y)+6, round(x)-5:round(x)+6].ravel()
                            for _, i, x, y, _, _ in stamp_data])
    limits = (np.nanpercentile(local, 3), np.nanpercentile(local, 99))
    for ax, (year, im, x, y, scale, aperture_row) in zip(axes, stamp_data, strict=True):
        show_stamp(ax, im, x, y, scale, limits, "gray", .24)
        ax.set_title(f"{year} · {format_flux(aperture_row)}", fontsize=9)
    axes[0].set_ylabel(f"{row.ra:.6f}°\n{row.dec:+.6f}°", fontsize=9)
    return result


def main():
    primary = read_receipt("products.json")
    followup = read_receipt("target_followup_products.json")
    sky = followup["ra"], followup["dec"]
    records = primary["images"] + followup["images"]
    stamps, measured, indexed = {}, {}, {}
    for record in records:
        band, year = record["band"].split("-")[0], int(record["dateobs"][:4])
        key = band, year
        indexed[key] = record
        mosaic = load_mosaic(record, sky)
        image, wcs, x, y, _, _ = mosaic
        scale = pixel_scale(wcs)
        stamps[key] = mosaic
        measured[key] = {
            "band": band, "year": year, "utc_start": record["dateobs"], "association": record["assoc_name"],
            "apertures": photometry(image, x, y, scale, RADII, (.28, .44)),
        }
        print(band, year, format_flux(measured[key]["apertures"][0]), flush=True)

    figure = Figure(figsize=(9.3, 9.4), layout="constrained")
    axes = figure.subplots(3, 3)
    changes = []
    for row, band in enumerate(("F444W", "F356W", "F277W")):
        old, new = (stamps[(band, year)] for year in (2025, 2026))
        older, newer = (measured[(band, year)]["apertures"] for year in (2025, 2026))
        changes.append({
            "band": band,
            "apertures": [flux_change(o, n) for o, n in zip(older, newer, strict=True)],
        })
        scale = pixel_scale(new[1])
        if np.hypot(old[2]-new[2], old[3]-new[3]) > .01 or abs(pixel_scale(old[1])-scale) > 1e-6:
            raise ValueError(f"{band}: mosaics differ in registration near target")
        stretch = .7 if band == "F277W" else 1.0
        show_stamp(axes[row, 0], old[0], old[2], old[3], scale, (-.18, stretch), "gray", RADII[0])
        show_stamp(axes[row, 1], new[0], new[2], new[3], scale, (-.18, stretch), "gray", RADII[0])
        show_stamp(axes[row, 2], new[0]-old[0], new[2], new[3], scale, (-stretch, stretch), "RdBu_r", RADII[0])
        axes[row, 0].set_title(f"{band} · October 2025\n{format_flux(older[0])}", fontsize=10)
        axes[row, 1].set_title(f"{band} · October 2026\n{format_flux(newer[0])}", fontsize=10)
        change = changes[-1]["apertures"][0]
        change_label = format_flux({
            "flux_njy": change["delta_njy"], "local_mad_sigma_njy": change["quadrature_local_scatter_njy"],
        })
        axes[row, 2].set_title(f"{band} · 2026 minus 2025\n{change_label}", fontsize=10)
    axes[0, 0].plot([-.60, -.20], [-.60, -.60], color="#69e2ec", lw=2)
    axes[0, 0].text(-.40, -.53, "0.4″", ha="center", color="#69e2ec", fontsize=8)
    figure.suptitle("GOODS-S brightening in the JWST mosaics\n"
                    f"RA {sky[0]:.6f}°, Dec {sky[1]:.6f}° · fixed 0.12″ aperture", fontsize=14)
    save_figure(figure, FIGURES / "target_evidence")

    dither_receipt = read_receipt("dither_products.json")
    if (dither_receipt["target_ra"], dither_receipt["target_dec"]) != sky:
        raise ValueError("Dither and mosaic receipt positions differ")
    groups = dither_receipt["groups"]
    native_groups, detector_audits = [], {}
    refs = reference_sky(stamps[("F444W", 2026)])
    for band in ("F444W", "F356W", "F277W"):
        frames = load_dithers(groups[band.lower()], refs)
        rows = audit_detector_target(frames, sky)
        native_groups.append((band, frames, rows))
        detector_audits[band] = summarize_detector(rows)
        print(band, "sky/peak separation spans", detector_audits[band]["sky_track_span_pixels"],
              detector_audits[band]["native_peak_span_pixels"], flush=True)
    save_figure(plot_detectors(native_groups), FIGURES / "target_dithers")
    actual_peaks = [tuple(r["unmasked_peak_pixels_xy"]) for group in detector_audits.values()
                    for r in group["exposures"] if r["unmasked_peak_pixels_xy"] is not None]
    detector_locked = (len(actual_peaks) == 9 and len(set(actual_peaks)) == 1
                       and all(g["jump_flagged_peaks"] == 3 and g["sky_track_span_pixels"] > 3
                               for g in detector_audits.values()))

    older_refs = reference_sky(stamps[("F444W", 2025)])
    older_frames = load_dithers(groups["f444w_2025"], older_refs)
    fading = fading_shortlist()
    fades = []
    if len(fading):
        fade_figure = Figure(figsize=(8.3, 2.52*len(fading)+.45), layout="constrained")
        fade_axes = fade_figure.subplots(len(fading), 3, squeeze=False)
        fades = [measure_fade(row, older_frames, stamps, indexed, axes)
                 for row, axes in zip(fading.itertuples(), fade_axes, strict=True)]
        fade_figure.suptitle("F444W imaging of the four fading positions\n"
                            "Common intensity scale per row · fixed 0.24″ apertures", fontsize=13)
        save_figure(fade_figure, FIGURES / "fading_followup")
    write_json(OUT / "candidate_photometry.json", {
        "executed_utc": utcnow(), "candidate": followup["candidate"], "ra": sky[0], "dec": sky[1],
        "brightening_disposition": "DETECTOR_FIXED_ARTIFACT" if detector_locked else "UNRESOLVED",
        "errors": "Local aperture scatter estimated from MAD; no trials correction or systematic error model",
        "mosaic_sky_annulus_arcsec": [.28, .44], "dither_sky_annulus_arcsec": [.32, .52],
        "mosaics": list(measured.values()), "changes": changes, "detector_checks_2026": detector_audits,
        "fading_shortlist_2025": fades,
        "input_receipt_sha256": {name: text_digest(OUT / name) for name in
                                 ("products.json", "matched_peaks.csv", "target_followup_products.json",
                                  "dither_products.json")},
        "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

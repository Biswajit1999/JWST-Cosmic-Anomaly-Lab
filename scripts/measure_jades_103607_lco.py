"""Measure the optical transient position on public Las Cumbres BANZAI images."""

import json
from collections import defaultdict
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from download_jades_103607_lco_eso import BUDGET
from download_jades_103607_lco_eso import TAG as INVENTORY
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_lco_r_history"
PS1 = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
PS1_COLS = (
    "objID", "raMean", "decMean", "nDetections", "gMeanPSFMag", "rMeanPSFMag",
    "rMeanPSFMagErr", "iMeanPSFMag", "rMeanKronMag",
)
COLOR_PIVOT = 0.7
GRID = np.arange(-72, 73, 12)
OFFSETS = [(east, north) for north in GRID for east in GRID if np.hypot(east, north) >= 12]
DISPLAY_SCALE = 0.45
DISPLAY_EDGES = (np.arange(-78, 79) * DISPLAY_SCALE)
DISPLAY_X, DISPLAY_Y = np.meshgrid(DISPLAY_EDGES, DISPLAY_EDGES)
SUBPIXELS = np.linspace(-0.4, 0.4, 5)


def mad(values):
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def reference_stars(cache, coord):
    url = PS1 + "?" + urlencode({
        "ra": f"{coord.ra.deg:.12f}", "dec": f"{coord.dec.deg:.12f}",
        "radius": 0.34, "nDetections.gte": 3, "rMeanPSFMag.lte": 21,
        "columns": f"[{','.join(PS1_COLS)}]", "pagesize": 5000,
    })
    path = cache.fetch(
        url, f"metadata/{TAG}_ps1_034deg.csv", survey="STScI MAST / Pan-STARRS1",
        version="PS1 DR2 mean photometry", kind="reference_star_catalogue",
    )
    table = pd.read_csv(path)
    if not len(table) or len(table) >= 5000 or not set(PS1_COLS) <= set(table.columns):
        raise ValueError("Missing or truncated Pan-STARRS reference catalogue")
    color = table["gMeanPSFMag"] - table["iMeanPSFMag"]
    use = (
        table["rMeanPSFMag"].between(14.5, 19.2)
        & table["rMeanPSFMagErr"].between(1e-6, 0.06)
        & table["gMeanPSFMag"].between(10, 24) & table["iMeanPSFMag"].between(10, 24)
        & (table["rMeanPSFMag"] - table["rMeanKronMag"]).abs().lt(0.14)
        & color.between(0, 2.4)
    )
    stars = table[use].copy()
    stars["color"] = color[use] - COLOR_PIVOT
    positions = SkyCoord(stars["raMean"], stars["decMean"], unit="deg")
    return stars, positions, dict(receipt=receipt(path, url), rows=len(table), selected_stars=len(stars))


def aperture(science, mask, error, center, radius):
    x, y = center
    half = int(np.ceil(2.8 * radius)) + 2
    ix, iy = round(x), round(y)
    ny, nx = science.shape
    if not (half <= ix < nx - half and half <= iy < ny - half):
        return None
    ys = slice(iy - half, iy + half + 1)
    xs = slice(ix - half, ix + half + 1)
    pixels = science[ys, xs]
    dy, dx = np.indices(pixels.shape, dtype=float)
    dx += ix - half - x
    dy += iy - half - y
    rr = np.hypot(dx, dy)
    core = sum(
        (dx + sx)**2 + (dy + sy)**2 <= radius**2 for sx in SUBPIXELS for sy in SUBPIXELS
    ) / SUBPIXELS.size**2
    valid = (mask[ys, xs] == 0) & np.isfinite(pixels)
    ring = (rr >= 1.8 * radius) & (rr <= 2.8 * radius)
    sky = ring & valid
    if np.any((core > 0) & ~valid) or sky.sum() < max(8, 0.75 * ring.sum()):
        return None
    count = float(np.dot(core[core > 0], pixels[core > 0] - np.median(pixels[sky])))
    formal = None
    if error is not None:
        sigma = error[ys, xs]
        used = (core > 0) | sky
        if np.isfinite(sigma[used]).all() and np.all(sigma[used] > 0):
            source_variance = np.sum((core[core > 0] * sigma[core > 0])**2)
            sky_variance = np.pi / 2 * np.sum(sigma[sky]**2) * (core.sum() / sky.sum())**2
            formal = float(np.sqrt(source_variance + sky_variance))
    return count, formal


def calibrate(science, mask, catalog, positions, wcs, radius, exposure, stars, star_coords):
    nearest, separation, _ = positions.match_to_catalog_sky(star_coords)
    sx, sy = wcs.world_to_pixel(positions)
    measurements, seen = [], set()
    for i in np.argsort(separation):
        k = nearest[i]
        if catalog["FLAG"][i] != 0 or separation[i].arcsec > 0.9 or k in seen:
            continue
        measured = aperture(science, mask, None, (sx[i], sy[i]), radius)
        if measured is None or measured[0] <= 0:
            continue
        if not np.isfinite(catalog["FLUXERR"][i]) or catalog["FLUXERR"][i] <= 0:
            continue
        if catalog["FLUX"][i] / catalog["FLUXERR"][i] < 15:
            continue
        reference = stars.iloc[k]
        zp = reference["rMeanPSFMag"] + 2.5 * np.log10(measured[0] / exposure)
        measurements.append((int(stars["objID"].iloc[k]), reference["color"], zp, separation[i].arcsec))
        seen.add(k)
    if len(measurements) < 8:
        return None
    ids, colors, zero, offset = zip(*measurements, strict=True)
    colors, zero = np.asarray(colors), np.asarray(zero)
    keep = abs(zero - np.median(zero)) <= max(0.12, 3.5 * mad(zero))
    design = np.column_stack((np.ones(len(colors)), colors))
    for _ in range(2):
        if np.count_nonzero(keep) < 8:
            return None
        intercept, slope = np.linalg.lstsq(design[keep], zero[keep], rcond=None)[0]
        residual = zero - design @ [intercept, slope]
        keep &= abs(residual) <= max(0.05, 3.5 * mad(residual[keep]))
    if np.count_nonzero(keep) < 8 or np.ptp(colors[keep]) < 0.4:
        return None
    intercept, slope = np.linalg.lstsq(design[keep], zero[keep], rcond=None)[0]
    scatter = mad((zero - design @ [intercept, slope])[keep])
    if scatter > 0.2:
        return None
    return dict(
        aperture_ab_zeropoint_at_pivot=float(intercept), color_coefficient=float(slope),
        zeropoint_residual_scatter_mag=scatter, matched_stars=int(len(zero)),
        used_stars=int(keep.sum()), reference_objids_used=[str(ids[i]) for i in np.flatnonzero(keep)],
        median_astrometry_separation_arcsec=float(np.median(np.asarray(offset)[keep])),
    )


def blank_fluxes(science, mask, error, wcs, sky, catalog_coords, radius):
    control_x, control_y = wcs.world_to_pixel(sky)
    sx, sy = wcs.world_to_pixel(catalog_coords)
    coords = np.column_stack([sx, sy])
    detected = coords[np.all(np.isfinite(coords), axis=1)]
    nearest, _ = cKDTree(detected).query(np.column_stack([control_x, control_y]))
    fluxes = np.full(len(sky), np.nan)
    for index in np.flatnonzero(nearest > 2.8 * radius):
        measure = aperture(science, mask, error, (control_x[index], control_y[index]), radius)
        if measure is not None:
            fluxes[index] = measure[0]
    return fluxes


def image_stamp(science, mask, wcs, display, scale, factor):
    x, y = wcs.world_to_pixel(display)
    x0 = max(0, int(np.floor(x.min())) - 2)
    y0 = max(0, int(np.floor(y.min())) - 2)
    local = science[y0:int(np.ceil(y.max())) + 3, x0:int(np.ceil(x.max())) + 3].astype(float)
    local[mask[y0:y0+local.shape[0], x0:x0+local.shape[1]] != 0] = np.nan
    local -= np.nanmedian(local)
    sampled = map_coordinates(local, (y-y0, x-x0), order=1, cval=np.nan, prefilter=False)
    return sampled * factor * (DISPLAY_SCALE / scale)**2


def process(row, coord, control, display, stars, star_coords):
    if row["status"] != "downloaded":
        return row, None, None
    meta = row["archive"]
    path = ROOT / meta["path"]
    if digest(path) != meta["sha256"]:
        raise ValueError(f"LCO image checksum mismatch for frame {row['id']}")
    result = dict(
        frame_id=row["id"], date_obs_utc=row["observation_date"], site=row["site_id"],
        instrument=row["instrument_id"], filter=row["primary_optical_element"],
        title=row["target_name"], proposal=row["proposal_id"], archive=meta,
    )
    with fits.open(path, memmap=True) as file:
        header = file["SCI"].header
        science, mask = file["SCI"].data, file["BPM"].data
        error = file["ERR"].data if "ERR" in file else None
        catalog = file["CAT"].data
        exposure, seeing = float(header["EXPTIME"]), float(header.get("L1FWHM", 0))
        date_offset = (Time(header["DATE-OBS"]) - Time(row["observation_date"])).sec
        if abs(date_offset) > 0.1 or header["FILTER"].strip() != result["filter"]:
            raise ValueError(f"LCO FITS exposure identity differs from the archive: {row['id']}")
        result.update(
            date_obs_utc=header["DATE-OBS"], exposure_seconds=exposure,
            fwhm_arcsec=seeing, pipeline=header.get("PIPEVER"), has_error_plane=error is not None,
            wcs_error=header.get("WCSERR"), banzai_zeropoint=header.get("L1ZP"),
            moon_illuminated_fraction=header.get("MOONFRAC"),
        )
        if result["filter"] != "rp" or result["wcs_error"] != 0 or not 0.6 < seeing < 7:
            return dict(**result, status="unsuitable band, image quality or astrometry"), None, None
        if science.shape != mask.shape or (error is not None and science.shape != error.shape):
            raise ValueError("Misaligned BANZAI pixel arrays")
        wcs = WCS(header).celestial
        scale = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600)
        x, y = map(float, wcs.world_to_pixel(coord))
        radius = seeing / scale
        target = aperture(science, mask, error, (x, y), radius)
        result.update(
            pixel_scale_arcsec=scale, target_xy_pixels=[float(x), float(y)],
            radius_arcsec=seeing, science_shape=list(science.shape),
        )
        if target is None or not len(catalog):
            return dict(**result, status="invalid native aperture or extraction catalogue"), None, None
        finite = np.isfinite(catalog["RA"]) & np.isfinite(catalog["DEC"])
        catalog = catalog[finite]
        positions = SkyCoord(catalog["RA"], catalog["DEC"], unit="deg")
        controls = blank_fluxes(science, mask, error, wcs, control, positions, radius)
        use = np.isfinite(controls)
        if np.count_nonzero(use) < 30 or mad(controls[use]) <= 0:
            return dict(**result, status="insufficient usable blank apertures"), None, None
        offset = np.median(controls[use])
        scatter_counts = mad(controls[use])
        comparison = positions.separation(coord).arcsec
        result.update(
            target_aperture_counts=target[0],
            local_control_center_counts=float(offset),
            signed_flux_counts=float(target[0] - offset),
            pipeline_formal_aperture_error_counts=target[1],
            local_aperture_error_counts=scatter_counts, controls_used=int(use.sum()),
            closest_catalogue_extraction_arcsec=float(comparison.min()),
        )
        calibration = calibrate(science, mask, catalog, positions, wcs, radius, exposure, stars, star_coords)
        if calibration is None:
            return dict(**result, status="insufficient reference calibration"), None, None
        factor = 3631e9 * 10**(-0.4 * calibration["aperture_ab_zeropoint_at_pivot"]) / exposure
        scatter = scatter_counts * factor
        raw_flux = (target[0] - offset) * factor
        image = image_stamp(science, mask, wcs, display, scale, factor)
        result.update(
            status="measured", calibration=calibration, signed_flux_njy=float(raw_flux),
            pipeline_formal_aperture_error_njy=None if target[1] is None else target[1] * factor,
            local_aperture_error_njy=scatter,
        )
    return result, (controls - offset) * factor, image


def five_scatter(control, target):
    center = float(np.median(control))
    scatter = mad(control)
    upper = 5 * scatter
    amplitude = upper - float(np.quantile(control - center, 0.05, method="lower"))
    return dict(
        signed_flux_njy=float(target - center), local_aperture_error_njy=scatter,
        five_scatter_point_source_depth_njy=upper, five_scatter_point_source_depth_ab=-2.5 * np.log10(upper / 3631e9),
        fixed_aperture_recovery_95pct_5scatter_njy=amplitude,
        fixed_aperture_recovery_95pct_5scatter_ab=-2.5 * np.log10(amplitude / 3631e9),
        control_count=len(control),
        controls_at_or_above_target=int(np.count_nonzero(control >= target)),
    )


def combine(group):
    frames, controls, images = zip(*group, strict=True)
    weights = 1 / np.array([frame["local_aperture_error_njy"] for frame in frames])**2
    common = np.all(np.isfinite(controls), axis=0)
    if np.count_nonzero(common) < 30:
        raise ValueError("Insufficient common blank apertures for a measured LCO epoch")
    control = np.average(np.array(controls)[:, common], weights=weights, axis=0)
    signal = np.average([frame["signed_flux_njy"] for frame in frames], weights=weights)
    measure = five_scatter(control, signal)
    image = np.asarray(images)
    valid = np.isfinite(image)
    weight = weights[:, None, None] * valid
    stack = np.divide(
        np.nansum(image * weight, axis=0), weight.sum(axis=0),
        out=np.full_like(image[0], np.nan), where=weight.sum(axis=0) > 0,
    )
    result = dict(
        first_utc=frames[0]["date_obs_utc"], last_utc=frames[-1]["date_obs_utc"],
        site=frames[0]["site"], instrument=frames[0]["instrument"],
        summed_exposure_seconds=sum(row["exposure_seconds"] for row in frames),
        aperture_radii_arcsec=[row["radius_arcsec"] for row in frames],
        frame_ids=[row["frame_id"] for row in frames], **measure,
    )
    if len(frames) > 2:
        frame_flux = np.array([frame["signed_flux_njy"] for frame in frames])
        loo = (weights @ frame_flux - weights * frame_flux) / (weights.sum() - weights)
        result["leave_one_exposure_out_flux_range_njy"] = [float(loo.min()), float(loo.max())]
    return result, stack


def plot(epochs, frame_rows):
    fig = Figure(figsize=(11, 7.8), layout="constrained")
    layout = fig.add_gridspec(2, 3, height_ratios=[1.15, 1])
    strongest = min(epochs, key=lambda pair: pair[0]["five_scatter_point_source_depth_njy"])
    current = epochs[-1]
    control = next((pair for pair in reversed(epochs[:-1]) if pair[0]["summed_exposure_seconds"] < 30), epochs[0])
    for i, (row, image) in enumerate((strongest, control, current)):
        ax = fig.add_subplot(layout[0, i])
        low, high = np.nanpercentile(image, [4, 98.8])
        edge = DISPLAY_EDGES[-1] + DISPLAY_SCALE / 2
        ax.imshow(image, origin="lower", cmap="gray", extent=(-edge, edge, -edge, edge), vmin=low, vmax=high)
        radius = float(np.median(row["aperture_radii_arcsec"]))
        ax.add_patch(Circle((0, 0), radius, edgecolor="#e58150", facecolor="none", lw=1.4))
        ax.set(
            title=f"{row['first_utc'][:10]} · {row['summed_exposure_seconds']:.0f} s",
            xlim=(20, -20), ylim=(-20, 20), xlabel="East offset (arcsec)",
        )
        if i == 0:
            ax.set_ylabel("North offset (arcsec)")
    axes = fig.add_subplot(layout[1, :])
    ordered = [row for row, _ in epochs]
    x = np.arange(len(ordered))
    values = [row["signed_flux_njy"] / 1e3 for row in ordered]
    scatter = [row["local_aperture_error_njy"] / 1e3 for row in ordered]
    depth = [row["five_scatter_point_source_depth_njy"] / 1e3 for row in ordered]
    axes.errorbar(x, values, yerr=scatter, fmt="o", capsize=3, color="#1c6b7a", label="Signed flux ± local scatter")
    axes.scatter(x, depth, marker="_", s=150, color="#ba642f", label="Five times local scatter")
    axes.axhline(0, color="#606060", lw=0.8)
    axes.set_yscale("symlog", linthresh=3)
    axes.set(
        ylabel="Calibrated r′ flux (µJy)", xlabel="Exposure date (UTC)",
        title=f"{len(frame_rows)} usable exposures grouped by UTC date, instrument, and site",
    )
    axes.set_xticks(x, [row["first_utc"][:10] for row in ordered], rotation=24)
    axes.grid(axis="y", alpha=0.25)
    axes.legend(loc="best", frameon=False)
    fig.suptitle("Las Cumbres r′ images at the JADES 103607 transient position", fontsize=16)
    save_figure(fig, FIGURES / TAG)


def main():
    inv = OUT / f"{INVENTORY}.json"
    data = json.loads(inv.read_text(encoding="utf-8"))
    coord = SkyCoord(*data["position_icrs_deg"], unit="deg")
    offsets = np.asarray(OFFSETS)
    controls = coord.spherical_offsets_by(offsets[:, 0] * u.arcsec, offsets[:, 1] * u.arcsec)
    display = coord.spherical_offsets_by(DISPLAY_X * u.arcsec, DISPLAY_Y * u.arcsec)
    stars, star_coords, source = reference_stars(DownloadCache(ROOT, BUDGET), coord)
    groups = defaultdict(list)
    rows = []
    for product in data["lco_frames"]:
        output = process(product, coord, controls, display, stars, star_coords)
        row = output[0]
        rows.append(row)
        if row["status"] == "measured":
            groups[row["date_obs_utc"][:10], row["site"], row["instrument"]].append(output)
        print(f"LCO {product['id']}: {row['status']}", flush=True)
    if not groups:
        raise ValueError("No qualifying Las Cumbres exposures")
    epochs = [combine(group) for _, group in sorted(groups.items())]
    plot(epochs, [row for row in rows if row["status"] == "measured"])
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__), inventory_sha256=text_digest(inv),
        position_icrs_deg=data["position_icrs_deg"], calibration_reference=source,
        ps1_g_minus_i_pivot=COLOR_PIVOT,
        control_offsets_arcsec=offsets.tolist(), frames=rows, epochs=[row for row, _ in epochs],
        photometry_scope="Native fixed circular aperture of one image FWHM radius, "
        "5x5 subpixel area samples, local median annulus at 1.8–2.8 radii, BPM=0. "
        "Each exposure is calibrated against matched PS1 r PSF magnitudes, fitting the "
        "optical color term with a pivot at g-i=0.7; the transient color is unknown. "
        "The epoch value is a weighted average with errors from common blank apertures. "
        "The five-scatter depth is an empirical sensitivity in a single band and time window. "
        "The 95% recovery estimate adds fixed-aperture flux to measured control residuals "
        "at fixed weights. It does not include mismatched source profiles, color corrections "
        "for an unknown transient spectrum, or calibration correlated across frames.",
    ))


if __name__ == "__main__":
    main()

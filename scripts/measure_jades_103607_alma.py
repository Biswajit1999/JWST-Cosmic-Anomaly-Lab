"""Perform fixed-coordinate continuum photometry on public ALMA image cutouts."""

import json
import re
from collections import defaultdict

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from download_jades_103607_alma import TAG as INVENTORY
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.dates import DateFormatter, YearLocator
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_alma_continuum"
CONTEXT = OUT / "module_b_001_redshift_context.json"
GRID = range(-16, 17, 2)
OFFSETS = np.asarray([(x, y) for y in GRID for x in GRID if 3 <= np.hypot(x, y) <= 16])
HALF_STAMP = 6
PB_LIMIT = 0.20
MIN_CONTROLS = 25


def read_plane(item):
    path = ROOT / item["path"]
    if digest(path) != item["sha256"]:
        raise ValueError(f"ALMA product SHA256 mismatch: {path.name}")
    with fits.open(path, memmap=False) as hdus:
        header = hdus[0].header.copy()
        image = np.asarray(hdus[0].data, dtype=float).squeeze()
    if image.ndim != 2:
        raise ValueError(f"ALMA image is not a single sky plane: {path.name}")
    # Some ALMA SODA headers retain legacy zero-padded PC matrix indices.
    for keyword in list(header):
        match = re.fullmatch(r"PC(\d+)_(\d+)", keyword)
        if match is None:
            continue
        canonical = f"PC{int(match[1])}_{int(match[2])}"
        if canonical == keyword:
            continue
        value = header.pop(keyword)
        if canonical in header and not np.isclose(header[canonical], value):
            raise ValueError(f"Conflicting ALMA WCS matrix entries: {path.name}")
        header[canonical] = value
    return image, header, WCS(header, fix=False).celestial


def sample(image, xy):
    return map_coordinates(image, [xy[1], xy[0]], order=1, prefilter=False, cval=np.nan)


def beam_model(image, xy, matrix, beam):
    x, y = xy
    major, minor, angle = beam
    scale = float(np.min(np.linalg.svd(matrix, compute_uv=False)))
    extent = int(np.ceil(1.65 * major / scale)) + 2
    ix, iy = round(x), round(y)
    if not (extent <= ix < image.shape[1] - extent and extent <= iy < image.shape[0] - extent):
        return np.nan
    values = image[iy-extent:iy+extent+1, ix-extent:ix+extent+1]
    yy, xx = np.indices(values.shape, dtype=float)
    xx += ix - extent - x
    yy += iy - extent - y
    east = matrix[0, 0] * xx + matrix[0, 1] * yy
    north = matrix[1, 0] * xx + matrix[1, 1] * yy
    sin, cos = np.sin(angle), np.cos(angle)
    ell2 = ((east * sin + north * cos) / major)**2 + ((east * cos - north * sin) / minor)**2
    support = ell2 <= 1.65**2
    if not np.all(np.isfinite(values[support])):
        return np.nan
    design = np.column_stack(
        (np.exp(-4 * np.log(2) * ell2[support]), np.ones(support.sum()), east[support], north[support])
    )
    return float(np.linalg.lstsq(design, values[support], rcond=None)[0][0])


def calibration(header):
    frequency = []
    for axis in range(1, header["NAXIS"] + 1):
        key = header.get(f"CTYPE{axis}", "")
        value = header[f"CRVAL{axis}"]
        if key == "STOKES" and not np.isclose(value, 1):
            raise ValueError("Expected ALMA Stokes I")
        if key == "FREQ":
            frequency.append((value * u.Unit(header[f"CUNIT{axis}"])).to_value(u.GHz))
    if len(frequency) != 1:
        raise ValueError("Missing ALMA frequency axis")
    beam = (
        header["BMAJ"] * 3600, header["BMIN"] * 3600, np.deg2rad(float(header["BPA"])),
    )
    factor = (1 * u.Unit(header["BUNIT"])).to_value(u.uJy / u.beam)
    return frequency[0], beam, factor


def mad(values):
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


def measure(row, sky, neighbor, controls):
    cutout = row["science_cutout"]
    metadata = {k: row[k] for k in (
        "obs_id", "member_ous_uid", "proposal_id", "target_name", "band_list", "t_min", "t_max",
    )}
    if cutout is None or cutout["status"] != "downloaded":
        return dict(**metadata, status="no public continuum cutout in this pass"), None
    image, header, wcs = read_plane(cutout)
    frequency, beam, factor = calibration(header)
    image *= factor
    matrix = wcs.pixel_scale_matrix * 3600
    coord = SkyCoord([sky, neighbor, *controls])
    x, y = wcs.world_to_pixel(coord)
    if np.any(~np.isfinite(x) | ~np.isfinite(y)):
        raise ValueError("Invalid local coordinate transform in ALMA cutout")
    metadata.update(
        frequency_ghz=frequency, restoring_beam_fwhm_arcsec=list(beam[:2]),
        restoring_beam_pa_deg=float(np.rad2deg(beam[2])),
        fits_date_obs_utc=header.get("DATE-OBS"), science_cutout=cutout,
        wcs_pixel_scale_arcsec=float(np.sqrt(abs(np.linalg.det(matrix)))),
    )
    pb = row["primary_beam_image"]
    pb_at_sky = None
    acceptable = np.ones(len(controls), dtype=bool)
    if pb and pb["status"] == "downloaded":
        pb_map, _, beam_wcs = read_plane(pb)
        bx, by = beam_wcs.world_to_pixel(coord)
        responses = sample(pb_map, (bx, by))
        pb_at_sky = float(responses[0]) if np.isfinite(responses[0]) else None
        metadata.update(primary_beam_product=pb, primary_beam_at_position=pb_at_sky)
        if pb_at_sky is None or pb_at_sky < PB_LIMIT:
            return dict(**metadata, status="target outside qualified primary beam"), None
        acceptable &= (responses[2:] >= max(PB_LIMIT, 0.7 * pb_at_sky)) & (responses[2:] <= 1.35 * pb_at_sky)
    else:
        metadata.update(primary_beam_at_position=None)
    target = beam_model(image, (x[0], y[0]), matrix, beam)
    near = beam_model(image, (x[1], y[1]), matrix, beam)
    if not np.isfinite(target):
        return dict(**metadata, status="missing or masked target beam support"), None
    acceptable &= controls.separation(neighbor).arcsec >= 2 * beam[0]
    acceptable &= controls.separation(sky).arcsec >= 2 * beam[0]
    amplitudes = np.full(len(controls), np.nan)
    for index in np.flatnonzero(acceptable):
        amplitudes[index] = beam_model(image, (x[index+2], y[index+2]), matrix, beam)
    valid = np.isfinite(amplitudes)
    if np.count_nonzero(valid) < MIN_CONTROLS:
        return dict(**metadata, status="insufficient usable local controls", controls=int(valid.sum())), None
    centered = np.median(amplitudes[valid])
    # A near-target subset catches steep noise changes in primary-beam-corrected edge pixels.
    local = valid & (np.hypot(*OFFSETS.T) <= 10)
    uncertainty = mad(amplitudes[valid])
    if np.count_nonzero(local) >= MIN_CONTROLS:
        uncertainty = max(uncertainty, mad(amplitudes[local]))
    if uncertainty <= 0:
        return dict(**metadata, status="zero local noise"), None
    target -= centered
    near -= centered
    metadata.update(
        status="measured", signed_flux_ujy=float(target), local_error_ujy=uncertainty,
        score=float(target / uncertainty), controls=int(valid.sum()),
        neighbor_signed_flux_ujy=float(near) if np.isfinite(near) else None,
        noise_qualification=(
            "measured primary beam and local residuals" if pb_at_sky is not None
            else "local residuals in the primary-beam-corrected image"
        ),
        control_median_ujy=float(centered),
    )
    return metadata, (amplitudes - centered, (image, wcs))


def mjd_iso(mjd):
    return Time(mjd, format="mjd", scale="utc").isot


def combine(group):
    rows, controls = zip(*group, strict=True)
    data = np.asarray(controls)
    weights = 1 / np.array([row["local_error_ujy"] for row in rows])**2
    valid = np.all(np.isfinite(data), axis=0)
    if np.count_nonzero(valid) < MIN_CONTROLS:
        # Keep a sensitivity-qualified result with the deepest single window.
        idx = int(np.argmax(weights))
        selected = [rows[idx]]
        ctrl = data[idx, np.isfinite(data[idx])]
        weights = weights[[idx]]
    else:
        selected = list(rows)
        ctrl = np.average(data[:, valid], axis=0, weights=weights)
    center = float(np.median(ctrl))
    sigma = float(max(mad(ctrl), 1 / np.sqrt(np.sum(weights))))
    flux = float(np.average([row["signed_flux_ujy"] for row in selected], weights=weights) - center)
    neighbor_values = np.array([
        row["neighbor_signed_flux_ujy"] if row["neighbor_signed_flux_ujy"] is not None else np.nan
        for row in selected
    ])
    neighbor = None
    if np.isfinite(neighbor_values).all():
        neighbor = float(np.average(neighbor_values, weights=weights) - center)
    first = selected[0]
    return dict(
        member_ous_uid=first["member_ous_uid"], proposal_id=first["proposal_id"],
        target_name=first["target_name"], band=int(first["band_list"]),
        start_utc=mjd_iso(min(row["t_min"] for row in selected)),
        end_utc=mjd_iso(max(row["t_max"] for row in selected)),
        image_ids=[row["obs_id"] for row in selected],
        candidate_windows=len(rows), control_count=len(ctrl),
        effective_frequency_ghz=float(np.average([row["frequency_ghz"] for row in selected], weights=weights)),
        signed_flux_ujy=flux, local_error_ujy=sigma, score=flux / sigma, five_scatter_depth_ujy=5 * sigma,
        bright_compact_flux_90pct_recovery_ujy=float(5 * sigma - np.quantile(ctrl - center, 0.10, method="lower")),
        controls_at_or_above_target=int(np.count_nonzero(ctrl >= flux + center)),
        largest_positive_control_score=float((np.max(ctrl) - center) / sigma),
        neighbor_signed_flux_ujy=neighbor, neighbor_score=None if neighbor is None else neighbor / sigma,
    )


def plot(epochs, images, sky, neighbor):
    figure = Figure(figsize=(11.4, 7.7), layout="constrained")
    grid = figure.add_gridspec(2, 3, height_ratios=(1.20, 1))
    selected = [
        min((row for row in epochs if row["band"] == band), key=lambda row: row["five_scatter_depth_ujy"])
        for band in sorted({row["band"] for row in epochs})
    ]
    latest = max(epochs, key=lambda row: row["start_utc"])
    if latest not in selected:
        selected.append(latest)
    if len(selected) < 3:
        others = [row for row in epochs if row not in selected]
        selected.extend(sorted(others, key=lambda row: row["local_error_ujy"])[:3-len(selected)])
    offsets = np.linspace(-HALF_STAMP, HALF_STAMP, 181)
    east, north = np.meshgrid(offsets, offsets)
    coords = sky.spherical_offsets_by(east * u.arcsec, north * u.arcsec)
    dx, dy = sky.spherical_offsets_to(neighbor)
    for i, epoch in enumerate(selected[:3]):
        row, (img, wcs) = min(
            [images[image_id] for image_id in epoch["image_ids"]],
            key=lambda v: v[0]["local_error_ujy"],
        )
        x, y = wcs.world_to_pixel(coords)
        panel = sample(img, (x, y)) / row["local_error_ujy"]
        ax = figure.add_subplot(grid[0, i])
        rendered = ax.imshow(
            panel, extent=(-HALF_STAMP, HALF_STAMP, -HALF_STAMP, HALF_STAMP),
            cmap="RdBu_r", origin="lower", vmin=-4, vmax=4,
        )
        ax.scatter(0, 0, marker="+", s=95, color="#303030")
        ax.scatter(dx.arcsec, dy.arcsec, s=45, facecolors="none", edgecolor="#2d8055")
        major, minor = row["restoring_beam_fwhm_arcsec"]
        ax.add_patch(Ellipse((4.5, -4.5), major, minor, angle=90-row["restoring_beam_pa_deg"],
                             fill=False, linewidth=1.25, color="#212121"))
        ax.set(title=f"{epoch['start_utc'][:10]} · {row['frequency_ghz']:.1f} GHz",
               xlim=(HALF_STAMP, -HALF_STAMP), xlabel="East offset (arcsec)")
        if i == 0:
            ax.set_ylabel("North offset (arcsec)")
        figure.colorbar(rendered, ax=ax, shrink=0.76, label="Local scatter units")
    ax = figure.add_subplot(grid[1, :])
    for band, color in ((3, "#236977"), (6, "#a95235")):
        points = [row for row in epochs if row["band"] == band]
        if not points:
            continue
        date = Time([row["start_utc"] for row in points], scale="utc").datetime
        ax.errorbar(
            date, [row["signed_flux_ujy"] for row in points],
            yerr=[row["local_error_ujy"] for row in points], fmt="o",
            lw=0.9, capsize=3, ms=4.5, label=f"Band {band}", color=color,
        )
    ax.axhline(0, color="#666666", lw=0.8)
    ax.xaxis.set_major_locator(YearLocator(1))
    ax.xaxis.set_major_formatter(DateFormatter("%Y"))
    ax.grid(axis="y", alpha=0.23)
    ax.set(ylabel="Compact continuum flux (µJy)", xlabel="Dataset start (UTC)",
           title="Signed beam-fit flux and local uncertainty by observing dataset")
    ax.legend(frameon=False)
    figure.suptitle("ALMA images and continuum measurements at the JADES 103607 transient position", fontsize=15)
    save_figure(figure, FIGURES / TAG)


def main():
    path = OUT / f"{INVENTORY}.json"
    source = json.loads(path.read_text(encoding="utf-8"))
    sky = SkyCoord(*source["position_icrs_deg"], unit="deg")
    other = json.loads(CONTEXT.read_text(encoding="utf-8"))["simbad"]["closest"]
    comparison = next(row for row in other if row["main_id"].endswith("103483"))
    neighbor = SkyCoord(comparison["ra"], comparison["dec"], unit="deg")
    controls = sky.spherical_offsets_by(OFFSETS[:, 0] * u.arcsec, OFFSETS[:, 1] * u.arcsec)
    groups, images, rows = defaultdict(list), {}, []
    for product in source["public_images"]:
        measured, diagnostic = measure(product, sky, neighbor, controls)
        rows.append(measured)
        if diagnostic is not None:
            amplitudes, image = diagnostic
            groups[measured["member_ous_uid"]].append((measured, amplitudes))
            images[measured["obs_id"]] = measured, image
        print(product["obs_id"], measured["status"], flush=True)
    if not groups:
        raise ValueError("No ALMA continuum cutout satisfied the local photometry checks")
    epochs = sorted([combine(group) for group in groups.values()], key=lambda row: row["start_utc"])
    plot(epochs, images, sky, neighbor)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__), inventory_sha256=text_digest(path),
        neighbor_context_sha256=text_digest(CONTEXT), position_icrs_deg=source["position_icrs_deg"],
        neighbor_simbad_position_icrs_deg=[neighbor.ra.deg, neighbor.dec.deg],
        neighbor_simbad_separation_arcsec=float(sky.separation(neighbor).arcsec),
        control_offsets_east_north_arcsec=OFFSETS.tolist(), windows=rows, datasets=epochs,
        method="Fit a fixed-position unit-peak restoring-beam profile plus a planar background on "
        "native primary-beam-corrected image pixels. Restrict measured primary beams to response "
        "of at least 0.20 at the target and compare controls at similar response. For larger beam "
        "maps, the depth comes from the finite, local corrected pixels without a separate beam gate. "
        "Use the robust local blank-position dispersion and require at least 25 controls. "
        "Combine only within each observation member, weighing spectral windows by inverse local "
        "variance and including the measured scatter of combination controls. The neighbor flux "
        "uses a separate fixed-position fit, with the target noise scale for comparison.",
        scope="Point-source sensitivity and signed continuum measurements in sampled historical "
        "datasets. The recovery estimate transfers an added matching-profile flux onto measured "
        "control residuals. Neighbor flux estimates may covary with the target in broad beams. "
        "Absolute calibration and resolved-source sensitivity are not included. Different spectral "
        "windows and beams do not measure a common passband. No spectral lines are inferred.",
    ))


if __name__ == "__main__":
    main()

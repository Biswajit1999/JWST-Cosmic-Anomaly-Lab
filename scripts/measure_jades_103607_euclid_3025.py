"""Test the 6 August 2024 Euclid excess with individual H images and adjacent bands."""

import json
import re
from urllib.parse import quote, urlencode
from xml.etree import ElementTree

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from measure_jades_103607_euclid_history import (
    CENTER,
    GRID,
    RADIUS,
    SCALE,
    SOURCE,
    EuclidSections,
    aperture,
    archive_rows,
    photometry,
    read_cutout,
    reproject_cutout,
)
from read_jwst_i2d_cutout import FITS_BLOCK, end_of_block
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_euclid_3025"
OBS = "3025"
Q1 = "https://nasa-irsa-euclid-q1.s3.us-east-1.amazonaws.com"
SIDE = 2040
SPAN = end_of_block(SIDE**2 * 4)
BANDS = ("Y", "J", "H")


class NativeHSections(EuclidSections):
    """Read only the parts of a Q1 detector exposure needed at the fixed sky position."""

    range_survey = "NASA IRSA Euclid Q1"
    range_cache_prefix = "cross_detector/euclid_q1_stripes"
    # Q1 DPDD defines bit 0 as INVALID. Higher bits can mark otherwise usable pixels.
    rejected_quality_bits = 1

    def __init__(self, cache, entry):
        self.cache = cache
        self.filename = entry["key"].rsplit("/", 1)[-1]
        self.url = f"{Q1}/{quote(entry['key'])}"
        self.etag = entry["etag"]
        self.size = entry["bytes"]
        self.ranges = []

    def target_detector(self, sky):
        overhead = self.size - 48 * SPAN
        if not 0 < overhead < 512 * 1024 or overhead % FITS_BLOCK:
            raise ValueError(f"Unexpected Q1 detector layout in {self.filename}")
        # Each detector group has three equal-sized padded data arrays. Locate the SCI
        # header in the bounded interval where cumulative headers can place it.
        for index in (14, 15, *range(14)):
            start = 3 * index * SPAN
            window = self.read(start, start + overhead).read_bytes()
            for cursor in range(0, overhead, FITS_BLOCK):
                block = window[cursor : cursor + FITS_BLOCK]
                if not block.startswith(b"XTENSION") or not re.search(rb"DET[1-4][1-4]\.SCI", block):
                    continue
                header, data = self.header(start + cursor)
                x, y = WCS(header, fix=False).celestial.world_to_pixel(sky)
                if not (-0.5 <= x < SIDE-0.5 and -0.5 <= y < SIDE-0.5):
                    continue
                detector = header["EXTNAME"].split(".")[0]
                planes = {"SCI": (header, data)}
                for name in ("RMS", "DQ"):
                    planes[name] = self.header(data + SPAN)
                    data = planes[name][1]
                for name, (hdr, _) in planes.items():
                    if (hdr["EXTNAME"], hdr["NAXIS1"], hdr["NAXIS2"], abs(hdr["BITPIX"])) != (
                        f"{detector}.{name}", SIDE, SIDE, 32
                    ):
                        raise ValueError(f"Unexpected Q1 detector HDU {detector}.{name}")
                return detector, planes
        return None, None


def native_inventory(cache):
    query = urlencode({"list-type": 2, "prefix": f"q1/NIR/{OBS}/EUC_NIR_W-CAL-IMAGE_H-{OBS}-"})
    url = f"{Q1}/?{query}"
    path = cache.fetch(
        url,
        f"metadata/{TAG}_q1_listing.xml",
        survey="NASA IRSA Euclid Q1",
        version="S3 ListObjectsV2 H calibrated exposures",
        kind="observed_visit_inventory",
    )
    root = ElementTree.parse(path).getroot()
    if root.findtext("{*}IsTruncated") != "false":
        raise ValueError("Truncated Q1 calibrated-exposure inventory")
    entries = []
    for obj in root.findall("{*}Contents"):
        key = obj.findtext("{*}Key")
        pattern = rf"q1/NIR/{OBS}/EUC_NIR_W-CAL-IMAGE_H-{OBS}-([0-3])_\d{{8}}T\d{{6}}\.\d+Z\.fits"
        match = re.fullmatch(pattern, key)
        if match:
            entries.append(dict(
                dither=int(match[1]),
                key=key,
                etag=obj.findtext("{*}ETag"),
                bytes=int(obj.findtext("{*}Size")),
            ))
    if not entries or len({r["dither"] for r in entries}) != len(entries):
        raise ValueError("Empty or nonunique Q1 H exposure inventory")
    return sorted(entries, key=lambda row: row["dither"]), receipt(path, url)


def native_cutout(cache, row, sky, grid):
    remote = NativeHSections(cache, row)
    primary, _ = remote.header(0)
    if (primary["NAXIS"], primary["FILTER"], int(primary["OBS_ID"]), primary["DITHOBS"]) != (
        0, "H", int(OBS), row["dither"]
    ):
        raise ValueError("Unexpected Q1 primary header")
    detector, planes = remote.target_detector(sky)
    row.update(
        start_date_utc=primary["UTC-OBS"],
        exposure_seconds=primary["EXPTIME"],
        file_ranges=remote.ranges,
        detector=detector,
    )
    if detector is None:
        row["status"] = "outside all 16 detectors"
        return None
    header, _ = planes["SCI"]
    if header["BUNIT"] != "ELECTRON":
        raise ValueError("Unexpected Q1 science unit")
    wcs = WCS(header, fix=False).celestial
    x, y = wcs.world_to_pixel(grid)
    # Use the local pixel-area Jacobian, including the native TPV distortion.
    c = CENTER
    jac = np.array([[p[c, c+1] - p[c, c-1], p[c+1, c] - p[c-1, c]] for p in (x, y)]) / 2
    area_ratio = abs(float(np.linalg.det(jac)))
    njy_per_electron = 3631e9 * 10 ** (-0.4 * header["ZPAB"])
    factor_check = (
        1000 * primary["PHOTFNU"] * primary["PHRELOB"] * primary["PHRELEX"]
        * header["PHRELDT"] / primary["EXPTIME"]
    )
    if not np.isclose(njy_per_electron, factor_check, rtol=0.01, atol=0):
        raise ValueError(f"Q1 calibration conventions disagree: {njy_per_electron}, {factor_check}")
    bounds = (
        max(0, int(np.floor(x.min())) - 1),
        min(SIDE, int(np.ceil(x.max())) + 2),
        max(0, int(np.floor(y.min())) - 1),
        min(SIDE, int(np.ceil(y.max())) + 2),
    )
    image = reproject_cutout(remote, planes, bounds, x, y) * area_ratio * njy_per_electron
    row.update(
        status="target on detector",
        detector_target_xy=[float(x[c, c]), float(y[c, c])],
        detector_edge_pixels=round(float(min(x[c, c], y[c, c], SIDE-1-x[c, c], SIDE-1-y[c, c])), 4),
        cutout_bounds=list(bounds),
        cutout_valid_fraction=float(np.isfinite(image).mean()),
        header_zpab=header["ZPAB"],
        njy_per_electron=njy_per_electron,
        njy_per_electron_from_photfnu=factor_check,
        output_to_detector_pixel_area_ratio=area_ratio,
    )
    print(f"H dither {row['dither']} {detector}, usable fraction {row['cutout_valid_fraction']:.2f}", flush=True)
    return image


def model_reference(images, minimum):
    values = np.asarray(images)
    ref = np.ma.median(np.ma.masked_invalid(values), axis=0).filled(np.nan)
    return np.where(np.isfinite(values).sum(axis=0) >= minimum, ref, np.nan)


def measure_band(cache, grid, band, selection):
    peak, *reference_rows = selection
    reference_images = [image for row in reference_rows if (image := read_cutout(cache, row, grid, band)) is not None]
    if len(reference_images) < 5:
        raise ValueError(f"Too few usable {band} reference visits")
    ref = model_reference(reference_images, len(reference_images) // 2 + 1)
    image = read_cutout(cache, peak, grid, band)
    if image is None:
        raise ValueError(f"No usable {band} data for {OBS}")
    residual = image - ref
    peak["difference"] = photometry(residual)
    print(band, peak["difference"]["signed_aperture_excess_njy"], flush=True)
    return ref, image, residual


def plot(panels):
    figure = Figure(figsize=(13.8, 6.9), layout="constrained")
    figure.suptitle("Euclid 6 August 2024 at the JADES 103607 transient position")
    axes = figure.subplots(2, 4, sharex=True, sharey=True).ravel()
    half = (GRID[-1] + 0.5) * SCALE
    for axis, (label, image, measurement) in zip(axes, panels, strict=False):
        if image is None or not np.isfinite(image).any():
            axis.set_title(label, fontsize=10)
            axis.set_axis_off()
            continue
        centered = image - np.nanmedian(image)
        limit = np.percentile(np.abs(centered[np.isfinite(centered)]), 95)
        render = axis.imshow(
            centered, cmap="RdBu_r", vmin=-limit, vmax=limit,
            extent=(-half, half, -half, half), origin="lower", interpolation="nearest",
        )
        axis.add_patch(Circle((0, 0), RADIUS, edgecolor="black", facecolor="none", lw=1.15))
        result = "Aperture unavailable" if "dither" in label else "Local median removed"
        if measurement is not None:
            excess = measurement["signed_aperture_excess_njy"]
            scatter = measurement["empirical_scatter_njy"]
            result = f"{excess:+.0f} ± {scatter:.0f} nJy"
        axis.set_title(f"{label}\n{result}", fontsize=10)
        axis.set_xticks([-5, 0, 5])
        axis.set_yticks([-5, 0, 5])
        figure.colorbar(render, ax=axis, shrink=0.73, label="nJy per grid pixel")
    for axis in axes[len(panels):]:
        axis.set_axis_off()
    figure.supxlabel("Offset east (arcsec)")
    figure.supylabel("Offset north (arcsec)")
    FIGURES.mkdir(parents=True, exist_ok=True)
    save_figure(figure, FIGURES / TAG)


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    east, north = np.meshgrid((GRID * SCALE / 3600) * sky.ra.unit, (GRID * SCALE / 3600) * sky.ra.unit)
    grid = sky.spherical_offsets_by(east, north)
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    native_rows, q1_receipt = native_inventory(cache)
    native = {row["dither"]: native_cutout(cache, row, sky, grid) for row in native_rows}
    inventories = {}
    for band in BANDS:
        rows, _, dr1_receipt = archive_rows(cache, sky, band)
        inventories[band] = {row["obs_id"]: row for row in rows if row["polygon_contains"]}
    shared = set.intersection(*(set(rows) for rows in inventories.values()))
    target = next(obs_id for obs_id in shared if obs_id.split("_")[0] == OBS)
    dates = {obs_id: inventories["H"][obs_id]["start_date_utc"] for obs_id in shared}
    early = sorted((obs_id for obs_id in shared if dates[obs_id] < dates[target]), key=dates.get)[-4:]
    late = sorted((obs_id for obs_id in shared if dates[obs_id] > dates[target]), key=dates.get)[:4]
    selections = {band: [inventories[band][obs_id] for obs_id in [target, *early, *late]] for band in BANDS}
    panels, reference_by_band = [], {}
    for band in BANDS:
        ref, image, residual = measure_band(cache, grid, band, selections[band])
        reference_by_band[band] = ref
        panels.append((f"DR1 {band} difference", residual, selections[band][0]["difference"]))
        if band == "H":
            panels.insert(0, ("DR1 H target visit", image, None))
            panels.append(("DR1 H reference visits", ref, None))
    for row in native_rows:
        image = native[row["dither"]]
        residual = None if image is None else image - reference_by_band["H"]
        row["difference"] = None if residual is None else photometry(residual)
        if residual is not None and row["difference"] is None:
            row["status"] = (
                "masked target aperture or local annulus"
                if not np.isfinite(aperture(residual, CENTER, CENTER))
                else "insufficient local controls"
            )
        label = f"Q1 H dither {row['dither']} · {row['detector']}"
        panels.append((label, residual, row["difference"]))
        if row["difference"] is not None:
            print(row["dither"], row["detector"], row["difference"]["signed_aperture_excess_njy"], flush=True)
    plot(panels)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            common_euclid_code_sha256=text_digest(ROOT / "scripts/measure_jades_103607_euclid_history.py"),
            position_source_sha256=text_digest(SOURCE),
            position_icrs_deg=position,
            q1_listing=q1_receipt,
            dr1_sia=dr1_receipt,
            aperture_radius_arcsec=RADIUS,
            reference_selection="Up to four nearest earlier and four nearest later covering visits "
            "with all three filters in SIA, requiring a majority of usable reference images per pixel.",
            reference_obs_ids=[*early, *late],
            dr1_rows_by_band=selections,
            q1_native_h_rows=native_rows,
            q1_h_dither_indices_without_public_file=sorted(set(range(4)) - native.keys()),
            q1_dq_rejected_mask=NativeHSections.rejected_quality_bits,
            scope="Observed signed aperture fluxes and robust aperture scatter measured on local "
            "control positions. Native Q1 H exposures use the same independent-visit DR1 H "
            "reference as the same-date stack. Q1 and DR1 use overlapping August 6 exposures "
            "and cannot count as independent observations. No PSF matching or aperture correction.",
        ),
    )


if __name__ == "__main__":
    main()

"""Inventory released VISTA VIDEO images whose published footprints contain the transient."""

import csv
import json
from collections import Counter
from urllib.parse import urlencode

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from requests import RequestException
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

SOURCE = OUT / "module_b_001_native_november_fade.json"
TAG = "module_b_001_vista_video_coverage"
TAP = "https://archive.eso.org/tap_obs/sync"
BANDS = {"8.29E-7": "Z", "9.74E-7": "Y", "1.166E-6": "J", "1.499E-6": "H", "1.992E-6": "Ks"}


SODA = "https://dataportal.eso.org/dataPortal/soda/sync"
H_VEGA_TO_AB = 1.3663
APERTURE_ARCSEC = 1.0


def aperture(image, x, y, scale):
    half = int(np.ceil(3.5 / scale)) + 1
    ix, iy = round(x), round(y)
    if not (half <= ix < image.shape[1] - half and half <= iy < image.shape[0] - half):
        return np.nan
    stamp = image[iy - half : iy + half + 1, ix - half : ix + half + 1]
    yy, xx = np.indices(stamp.shape)
    distance = np.hypot(xx + ix - half - x, yy + iy - half - y) * scale
    source = stamp[distance <= APERTURE_ARCSEC]
    background = stamp[(distance >= 2.0) & (distance <= 3.5)]
    if not np.isfinite(source).all() or not np.isfinite(background).all():
        return np.nan
    return float((source - np.median(background)).sum())


def inspect_cutout(path, sky):
    with fits.open(path) as hdus:
        if hdus[0].header.get("FILTER") != "H":
            raise ValueError("Unexpected pilot cutout filter")
        for hdu in hdus[1:]:
            if hdu.data is None or hdu.data.ndim != 2:
                continue
            wcs = WCS(hdu.header, fix=False).celestial
            x, y = map(float, wcs.world_to_pixel(sky))
            if not np.isfinite([x, y]).all():
                continue
            scale = float(np.sqrt(np.prod(proj_plane_pixel_scales(wcs))) * 3600)
            if hdu.header["PHOTSYS"].strip().lower() != "vega":
                raise ValueError("Unexpected pilot cutout calibration")
            zp = hdu.header["PHOTZP"] + H_VEGA_TO_AB
            image = np.array(hdu.data, dtype=float) * (3631e9 * 10 ** (-0.4 * zp))
            measured = aperture(image, x, y, scale)
            if not np.isfinite(measured):
                continue
            offsets = range(-20, 21, 5)
            controls = np.asarray(
                [
                    aperture(image, x + dx / scale, y + dy / scale, scale)
                    for dx in offsets
                    for dy in offsets
                    if dx or dy
                ]
            )
            controls = controls[np.isfinite(controls)]
            if len(controls) < 40:
                continue
            center = float(np.median(controls))
            scatter = float(1.4826 * np.median(np.abs(controls - center)))
            return dict(
                status="measured",
                chip=hdu.name,
                photzp_ab=zp,
                seeing_arcsec=hdu.header.get("PSF_FWHM"),
                signed_aperture_njy=measured - center,
                local_aperture_scatter_njy=scatter,
                local_control_median_njy=center,
                five_scatter_aperture_njy=5 * scatter,
                control_count=len(controls),
            ), (image, x, y, scale)
    return dict(status="cutout lacks usable local aperture and controls"), None


def pilot(cache, images, sky):
    # Take dates from the first, median and last distinct H-band observing blocks.
    h_blocks = {}
    for row in images:
        if row["band"] == "H":
            h_blocks.setdefault(row["observing_block"], row)
    h_rows = list(h_blocks.values())
    attempts, panels = [], []
    for index in (0, len(h_rows) // 2, len(h_rows) - 1):
        row = h_rows[index]
        dataset = row["product_id"]
        params = dict(
            ID=f"ivo://eso.org/ID?{dataset}", CIRCLE=f"{sky.ra.deg:.10f} {sky.dec.deg:.10f} 0.010"
        )
        url = f"{SODA}?{urlencode(params)}"
        try:
            path = cache.fetch(
                url,
                f"cross_detector/vista_video/{dataset.replace(':', '-')}_cutout.fits",
                survey="ESO VISTA VIDEO",
                version="DR5 H-band pawprint SODA cutout",
                kind="observed_survey_images",
            )
        except RequestException as error:
            attempts.append(dict(**row, status=type(error).__name__, archive_query=url))
            continue
        result, panel = inspect_cutout(path, sky)
        attempts.append(dict(**row, **result, archive=receipt(path, url)))
        if panel:
            panels.append((row["start_utc"][:10], result, panel))
    if panels:
        plot_cutouts(panels)
    return attempts


def plot_cutouts(panels):
    fig = Figure(figsize=(8.2, 3.85), layout="constrained")
    axes = fig.subplots(1, len(panels), squeeze=False)[0]
    for axis, (date, row, (image, x, y, scale)) in zip(axes, panels, strict=True):
        low, high = np.nanpercentile(image, [15, 99.6])
        extent = (
            (-x - 0.5) * scale,
            (image.shape[1] - x - 0.5) * scale,
            (-y - 0.5) * scale,
            (image.shape[0] - y - 0.5) * scale,
        )
        axis.imshow(image, extent=extent, origin="lower", cmap="gray", vmin=low, vmax=high)
        axis.scatter([0], [0], marker="+", c="#eb7348", s=35)
        depth = row["five_scatter_aperture_njy"] / 1000
        axis.set(
            xlim=(-6, 6),
            ylim=(-6, 6),
            xticks=(-5, 0, 5),
            yticks=(-5, 0, 5),
        )
        axis.set_title(f"{date}\n5-scatter depth {depth:.2f} µJy", fontsize=9)
    fig.supxlabel("Detector x offset (arcsec)", fontsize=10)
    fig.supylabel("Detector y offset (arcsec)", fontsize=10)
    fig.suptitle("JADES 103607 position | observed VISTA H-band pawprints")
    save_figure(fig, FIGURES / TAG)


def main():
    ra, dec = json.loads(SOURCE.read_text())["icrs_deg"]
    # Only individual pawprint images, without catalogue duplicates or stacked tiles.
    query = (
        "SELECT obs_id, dp_id, em_min, em_max, t_min, t_max, t_exptime "
        "FROM ivoa.ObsCore WHERE obs_collection='VIDEO' AND dataproduct_subtype='pawprint' "
        f"AND 1=INTERSECTS(s_region, POINT('ICRS',{ra:.10f},{dec:.10f})) "
        "ORDER BY t_min, dp_id"
    )
    url = f"{TAP}?{urlencode(dict(REQUEST='doQuery', LANG='ADQL', MAXREC=10000, FORMAT='csv', QUERY=query))}"
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    path = cache.fetch(
        url,
        f"metadata/{TAG}.csv",
        survey="ESO VISTA VIDEO",
        version="VIDEO DR5 ObsCore",
        kind="observed_image_footprint_metadata",
    )
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not 0 < len(rows) < 10_000 or len({row["dp_id"] for row in rows}) != len(rows):
        raise ValueError("Expected a complete nonempty set of distinct pawprint products")
    # MJD metadata supply an observation interval, not continuous exposure.
    images = []
    for row in rows:
        band = BANDS[row["em_min"].upper()]
        start_mjd, stop_mjd = float(row["t_min"]), float(row["t_max"])
        start_utc, stop_utc = Time([start_mjd, stop_mjd], format="mjd", scale="utc").isot
        images.append(
            dict(
                product_id=row["dp_id"],
                observing_block=row["obs_id"],
                band=band,
                start_mjd=start_mjd,
                stop_mjd=stop_mjd,
                start_utc=start_utc,
                stop_utc=stop_utc,
                exposure_s=float(row["t_exptime"]),
            )
        )
    summary = {}
    for band in BANDS.values():
        band_rows = [image for image in images if image["band"] == band]
        summary[band] = dict(
            pawprints=len(band_rows),
            observing_blocks=len({row["observing_block"] for row in band_rows}),
            first_start_utc=min(row["start_utc"] for row in band_rows),
            last_stop_utc=max(row["stop_utc"] for row in band_rows),
            pawprints_by_year=dict(Counter(row["start_utc"][:4] for row in band_rows)),
        )
    tests = pilot(cache, images, SkyCoord(ra, dec, unit="deg"))
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=[ra, dec],
            position_sha256=text_digest(SOURCE),
            archive=receipt(path, url),
            adql=query,
            distinct_observing_blocks=len({row["observing_block"] for row in images}),
            summary=summary,
            h_band_pilot=tests,
            h_band_pilot_aperture_radius_arcsec=APERTURE_ARCSEC,
            observations=images,
            scope="Published pawprint polygons contain the exact query point. Per-epoch pixel "
            "masks, usable exposure and sensitivity require image-level checks. Overlapping "
            "pawprints can belong to the same observing block. The finite H-band aperture "
            "pilot subtracts a local annulus and uses signed, control-centered photometry. "
            "Its spatial scatter omits aperture correction, per-pixel confidence maps and "
            "between-epoch PSF matching. It measures sensitivity, not variability.",
        ),
    )


if __name__ == "__main__":
    main()

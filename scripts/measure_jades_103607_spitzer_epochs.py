"""Compare independently dated IRAC observations at JADES 103607."""

import json
import re
from pathlib import PurePosixPath
from urllib.parse import urlencode, urlsplit

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from check_jades_103607_chandra_spitzer_history import IRSA as GOODS_IRSA
from check_jades_103607_chandra_spitzer_history import irac_aperture
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from read_jwst_i2d_cutout import MastImageSections
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

STEM = "module_b_001_spitzer_epochs"
SOURCE = OUT / "module_b_001_native_november_fade.json"
SIMPLE = "https://irsa.ipac.caltech.edu/data/SPITZER/SIMPLE/images"
SHA = "https://irsa.ipac.caltech.edu/SIA"
PROGRAM = "GO8_S2_GOODS_S_80217"
WINDOWS = {"2011": "55821 55866", "2012": "56015 56024"}
MIN_COUNT = 2
MAX_COUNT = 3
SCALE = 0.6
GRID = np.arange(-50, 51)
RADIUS = 1.8
SIMPLE_NJY_PER_DN_S = {"1": 3922.0, "2": 4808.0}


class SpitzerSections(MastImageSections):
    def __init__(self, cache, url):
        self.cache = cache
        self.url = url
        self.filename = PurePosixPath(urlsplit(url).path).name
        with cache.session.head(url, timeout=(15, 40)) as response:
            response.raise_for_status()
            self.etag = response.headers["ETag"]
            self.size = int(response.headers["Content-Length"])
        self.ranges = []


def checked_header(source, sky):
    h, offset = source.header(0)
    w = WCS(h, fix=False).celestial
    px, py = map(float, w.world_to_pixel(sky))
    return h, offset, w, (px, py)


def central_exposure(cache, url, sky):
    source = SpitzerSections(cache, url)
    h, offset, wcs, (px, py) = checked_header(source, sky)
    x, y = round(px), round(py)
    if not (0 <= x < h["NAXIS1"] and 0 <= y < h["NAXIS2"]):
        return dict(exposure_map_value=None, ranges=source.ranges), None
    if h["BITPIX"] != -32 or h.get("BZERO", 0) or h.get("BSCALE", 1) != 1:
        raise ValueError("Unexpected exposure image encoding")
    begin = offset + 4 * (y * h["NAXIS1"] + x)
    local = float(np.fromfile(source.read(begin, begin + 4), dtype=">f4")[0])
    return dict(exposure_map_value=local if np.isfinite(local) else None, ranges=source.ranges), (
        source,
        h,
        offset,
        wcs,
    )


def read_pbcd(cache, url, sky, grid, order):
    name = PurePosixPath(urlsplit(url).path).name
    path = cache.fetch(
        url,
        f"cross_detector/{name}",
        survey="Spitzer S-CANDELS",
        version="SHA Post-BCD science or coverage mosaic",
        kind="observed_epoch_fits",
    )
    with fits.open(path, memmap=False) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
        wcs = WCS(header, fix=False).celestial
        px, py = map(float, wcs.world_to_pixel(sky))
        exposure_at_center = None
        if 0 <= round(px) < data.shape[1] and 0 <= round(py) < data.shape[0]:
            value = float(data[round(py), round(px)])
            if np.isfinite(value):
                exposure_at_center = value
        if exposure_at_center is None or (order == 0 and exposure_at_center <= 0):
            return None, header, exposure_at_center, receipt(path, url)
        x, y = wcs.world_to_pixel(grid)
        sample = map_coordinates(
            data, [y, x], order=order, mode="constant", cval=np.nan, prefilter=False
        )
    return sample, header, exposure_at_center, receipt(path, url)


def simple_coverage(cache, sky):
    records = []
    for epoch, dates in ((1, ["2005-08-19", "2005-08-23"]), (2, ["2006-02-06", "2006-02-11"])):
        for ch in (1, 2):
            name = f"SIMPLE_ECDFS_IRAC_ch{ch}_ep{epoch}_exp.fits"
            url = f"{SIMPLE}/epoch{epoch}/{name}"
            record, _ = central_exposure(cache, url, sky)
            records.append(dict(epoch=epoch, channel=ch, dates_utc=dates, url=url, **record))
    return records


def candidates(cache, sky, year):
    args = dict(
        POS=f"CIRCLE {sky.ra.deg:.10f} {sky.dec.deg:.10f} 0.00028",
        COLLECTION="spitzer_sha",
        INSTRUMENT="IRAC",
        CALIB=2,
        BAND="3.1e-6 5.0e-6",
        FORMAT="image/fits",
        TIME=WINDOWS[year],
        RESPONSEFORMAT="CSV",
        MAXREC=10_000,
    )
    url = SHA + "?" + urlencode(args)
    path = cache.fetch(
        url,
        f"metadata/{STEM}_{year}.csv",
        survey="Spitzer Heritage Archive",
        version="IRSA SIA2 position and date query",
        kind="archive_metadata",
    )
    table = pd.read_csv(path)
    keep = table.loc[
        table.proposal_id.eq(PROGRAM)
        & table.access_url.str.contains(r"_(?:maic|mcov)\.fits$", na=False)
    ]
    paths = keep.access_url.str.extract(
        r"^(?P<base>https://irsa\.ipac\.caltech\.edu/.*/r(?P<aor>\d+)/ch(?P<channel>[12])/pbcd)/[^/]+_(?P<kind>maic|mcov)\.fits$"
    )
    keep = keep.join(paths).dropna(subset=["aor"])
    return keep, dict(
        product=receipt(path, url),
        query_rows=len(table),
        paired_product_candidates=len(keep),
        potentially_truncated=len(table) >= args["MAXREC"],
    )


def inspect_year(cache, sky, grid, year):
    table, query = candidates(cache, sky, year)
    output, records = {}, {}
    for channel in ("1", "2"):
        attempted, stack = [], []
        channel_table = table.loc[table.channel.eq(channel)].sort_values(
            ["t_min", "aor", "access_url"]
        )
        for aor, matches in channel_table.groupby("aor", sort=False):
            found = {}
            for kind in ("maic", "mcov"):
                options = matches.loc[matches.kind.eq(kind)].access_url.unique()
                if len(options) == 1:
                    found[kind] = options[0]
            if set(found) != {"maic", "mcov"}:
                continue
            exposure, _, local, coverage_receipt = read_pbcd(
                cache, found["mcov"], sky, grid, order=0
            )
            record = dict(
                aor=aor,
                dates_mjd=[float(matches.t_min.min()), float(matches.t_max.max())],
                product_urls=found,
                coverage=coverage_receipt,
                exposure_map_value=local,
            )
            attempted.append(record)
            if exposure is None:
                record["status"] = "no target coverage"
                continue
            center = len(GRID) // 2
            if not np.isfinite(irac_aperture(exposure, exposure, center, center)):
                record["status"] = "incomplete aperture footprint"
                continue
            image, header, _, science_receipt = read_pbcd(cache, found["maic"], sky, grid, order=1)
            unit = re.sub(r"\s+", "", header["BUNIT"]).upper()
            if unit not in {"MJY/SR", "MJY/STER"}:
                raise ValueError(f"Unexpected Spitzer PBCD units {unit}")
            record["science"] = science_receipt
            if image is None:
                record["status"] = "no valid science pixel at target"
                continue
            factor = (SCALE * u.arcsec).to_value(u.rad) ** 2 * 1e15
            mask = np.isfinite(image) & np.isfinite(exposure) & (exposure > 0)
            im = np.where(mask, image * factor, np.nan)
            if not np.isfinite(irac_aperture(im, exposure, center, center)):
                record["status"] = "science pixels invalid at target"
                continue
            stack.append(im - np.nanmedian(im))
            record["status"] = "included"
            print(year, channel, aor, "covered", flush=True)
            if len(stack) >= MAX_COUNT:
                break
        records[channel] = dict(attempts=attempted, included=len(stack))
        if len(stack) >= MIN_COUNT:
            valid = np.sum(np.isfinite(stack), axis=0) >= MIN_COUNT
            comb = np.ma.median(np.ma.masked_invalid(stack), axis=0).filled(np.nan)
            output[channel] = np.where(valid, comb, np.nan), valid
    return dict(year=year, channels=records, inventory=query), output


def released_grid(cache, url, grid):
    source = SpitzerSections(cache, url)
    header, offset = source.header(0)
    wcs = WCS(header, fix=False).celestial
    native_scale = abs(np.linalg.det(wcs.pixel_scale_matrix)) ** 0.5 * 3600
    if abs(native_scale - SCALE) > 0.001:
        raise ValueError("Unexpected GOODS or SIMPLE mosaic sampling")
    x, y = wcs.world_to_pixel(grid)
    xmin, xmax = int(np.floor(x.min())) - 2, int(np.ceil(x.max())) + 3
    ymin, ymax = int(np.floor(y.min())) - 2, int(np.ceil(y.max())) + 3
    image = source.stripe((header, offset), (xmin, xmax, ymin, ymax))
    sample = map_coordinates(image, [y - ymin, x - xmin], order=1, prefilter=False)
    return sample, header, source.ranges


def released_mosaics(cache, grid):
    output, records = {}, {}
    for channel, epoch in (("1", 1), ("2", 2)):
        output[channel], records[channel] = {}, {}
        urls = {
            "2004": [
                f"{GOODS_IRSA}/s_irac_{channel}_s2_v0.30_{name}.fits" for name in ("sci", "exp")
            ],
            "2005" if epoch == 1 else "2006": [
                f"{SIMPLE}/epoch{epoch}/SIMPLE_ECDFS_IRAC_ch{channel}_ep{epoch}{suffix}.fits"
                for suffix in ("", "_exp")
            ],
        }
        for year, (science, expmap) in urls.items():
            data, h, science_rows = released_grid(cache, science, grid)
            exposure, w, weight_rows = released_grid(cache, expmap, grid)
            if h["BUNIT"].strip() != "DN/sec" or w["BUNIT"].strip() != "sec":
                raise ValueError("Unexpected source mosaic photometric units")
            factor = (
                float(h.get("FLUXCONV", 0)) * 1000
                if year == "2004"
                else SIMPLE_NJY_PER_DN_S[channel]
            )
            if not factor > 0:
                raise ValueError("GOODS science calibration missing")
            valid = np.isfinite(data) & np.isfinite(exposure) & (exposure > 0)
            output[channel][year] = np.where(valid, data * factor, np.nan), valid
            records[channel][year] = dict(
                data_ranges=science_rows + weight_rows,
                calibration_njy_per_dn_second=factor,
                calibration_source="GOODS DR2 FLUXCONV"
                if year == "2004"
                else "Damen et al. 2011 Table 3",
            )
    return output, records


def difference(epochs):
    results, images = {}, {}
    n = len(GRID) // 2
    for channel, observations in epochs.items():
        results[channel] = {}
        early_image, early_coverage = observations["2004"]
        for year, (late_image, late_coverage) in observations.items():
            if year == "2004":
                continue
            diff = late_image - early_image
            coverage = early_coverage & late_coverage
            target = irac_aperture(diff, coverage, n, n)
            controls = np.asarray(
                [
                    irac_aperture(diff, coverage, x, y)
                    for x in range(9, len(GRID) - 9, 10)
                    for y in range(9, len(GRID) - 9, 10)
                    if np.hypot(x - n, y - n) > 10
                ]
            )
            controls = controls[np.isfinite(controls)]
            if len(controls) < 15 or not np.isfinite(target):
                continue
            median = float(np.median(controls))
            sigma = float(1.4826 * np.median(abs(controls - median)))
            results[channel][year] = dict(
                signed_aperture_difference_njy=target,
                control_median_njy=median,
                differential_njy=target - median,
                scatter_njy=sigma,
                baseline_2004_aperture_njy=irac_aperture(early_image, early_coverage, n, n),
                later_aperture_njy=irac_aperture(late_image, late_coverage, n, n),
                controls_njy=controls.tolist(),
            )
            images[channel, year] = early_image, late_image, diff
    return results, images


def plot(results, images):
    fig = Figure(figsize=(9.1, 2.3 * len(images)), layout="constrained")
    axes = fig.subplots(len(images), 3, squeeze=False)
    n = len(GRID) // 2
    r = 15
    edge = (r + 0.5) * SCALE
    for row, ((channel, year), stamps) in zip(axes, images.items(), strict=True):
        m = results[channel][year]
        vmax = np.nanpercentile(
            abs(np.array(stamps[:2])[:, n - r : n + r + 1, n - r : n + r + 1]), 98
        )
        vdiff = np.nanpercentile(abs(stamps[2][n - r : n + r + 1, n - r : n + r + 1]), 98)
        for i, (ax, im, title) in enumerate(
            zip(row, stamps, ("2004 GOODS", year, f"{year} minus 2004"), strict=True)
        ):
            lim = vmax if i < 2 else vdiff
            ax.imshow(
                im[n - r : n + r + 1, n - r : n + r + 1],
                origin="lower",
                cmap="magma" if i < 2 else "RdBu_r",
                vmin=-0.13 * lim if i < 2 else -lim,
                vmax=lim,
                extent=(-edge, edge, -edge, edge),
            )
            ax.add_patch(
                Circle((0, 0), RADIUS, fill=False, color="#5bc6cd" if i < 2 else "#254c58", lw=1.2)
            )
            ax.set(title=title, xticks=[], yticks=[])
        row[0].set_ylabel(f"IRAC {3.6 if channel == '1' else 4.5:.1f} µm")
        row[2].set_xlabel(f"{m['differential_njy']:+.0f} ± {m['scatter_njy']:.0f} nJy")
    fig.suptitle("JADES 103607 · historical Spitzer infrared imaging", fontsize=13)
    save_figure(fig, FIGURES / STEM)


def main():
    coords = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*coords, unit="deg")
    x, y = np.meshgrid(GRID * SCALE, GRID * SCALE)
    grid = sky.spherical_offsets_by(-x * u.arcsec, y * u.arcsec)
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    simple = simple_coverage(cache, sky)
    early, a = inspect_year(cache, sky, grid, "2011")
    late, b = inspect_year(cache, sky, grid, "2012")
    observations, releases = released_mosaics(cache, grid)
    for year, scandels in (("2011", a), ("2012", b)):
        for channel, image in scandels.items():
            observations[channel][year] = image
    values, images = difference(observations)
    if images:
        plot(values, images)
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=coords,
            position_source_sha256=text_digest(SOURCE),
            simple_epoch_coverage=simple,
            scandels=[early, late],
            released_mosaics=releases,
            difference=values,
            scope="ROI-level PBCD photometry with coverage masks, observed pixels, and local control scatter. "
            "Finite aperture and confusion allow no assignment of blended flux to the source alone.",
        ),
    )
    for ch, measurements in values.items():
        for year, row in measurements.items():
            print(ch, year, row["differential_njy"], row["scatter_njy"], len(row["controls_njy"]))


if __name__ == "__main__":
    main()

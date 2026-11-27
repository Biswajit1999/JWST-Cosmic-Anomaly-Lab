"""Measure the JADES 103607 position in Euclid DR1 H-band visit stacks."""

import json
from datetime import datetime
from urllib.parse import quote, urlencode
from xml.etree import ElementTree

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from check_jades_103607_followup import signed_margin_arcsec
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.dates import DateFormatter, DayLocator, HourLocator
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from measure_jades_103607_hst_history import SOURCE, difference
from read_jwst_i2d_cutout import MastImageSections, end_of_block
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

COLLECTION = "euclid_dr1_DpdNirStackedFrame"
SIA = "https://irsa.ipac.caltech.edu/SIA"
TAG = "module_b_001_euclid_dr1_h_history"
SCALE = 0.30
GRID = np.arange(-26, 27)
CENTER = int(len(GRID) // 2)
RADIUS = 0.75
ANNULUS = (1.2, 1.8)
POSITIONS = range(7, len(GRID) - 6, 6)


class EuclidSections(MastImageSections):
    """Read pinned SCI, RMS and data-quality stripes of an IRSA Euclid image."""

    range_survey = "NASA IRSA Euclid DR1"
    range_cache_prefix = "cross_detector/euclid_dr1_stripes"
    rejected_quality_bits = 0xFFFFFFFF

    def __init__(self, cache, cloud):
        self.cache = cache
        self.filename = cloud["key"].rsplit("/", 1)[-1]
        if cloud["bucket_name"] != "nasa-irsa-euclid-dr1" or cloud["region"] != "us-east-1":
            raise ValueError("Unexpected IRSA Euclid DR1 image bucket")
        self.url = f"https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/{quote(cloud['key'])}"
        with cache.session.head(self.url, timeout=(15, 40)) as response:
            response.raise_for_status()
            self.etag = response.headers["ETag"]
            self.size = int(response.headers["Content-Length"])
        self.ranges = []

    def planes(self, band="H"):
        primary, offset = self.header(0)
        if primary["NAXIS"] != 0 or primary["FILTER"] != band:
            raise ValueError(f"Expected {band}-band NISP image in {self.filename}")
        planes = {}
        for name in ("SCI", "RMS", "DQ"):
            header, start = self.header(offset)
            if header["EXTNAME"] != name or header["NAXIS"] != 2 or abs(header["BITPIX"]) != 32:
                raise ValueError(f"Unexpected Euclid {name} image in {self.filename}")
            planes[name] = header, start
            offset = start + end_of_block(4 * header["NAXIS1"] * header["NAXIS2"])
        return primary, planes

    def quality_stripe(self, entry, bounds):
        header, offset = entry
        x0, x1, y0, y1 = bounds
        zero = header.get("BZERO", 0)
        if header["BITPIX"] != 32 or zero not in (0, 2**31) or header.get("BSCALE", 1) != 1:
            raise ValueError("Unexpected Euclid DQ encoding")
        width = header["NAXIS1"]
        path = self.read(offset + 4 * width * y0, offset + 4 * width * y1)
        dq = np.memmap(path, dtype=">i4", mode="r", shape=(y1 - y0, width))
        flags = dq[:, x0:x1].astype(np.int64) + zero
        return (flags & self.rejected_quality_bits) == 0


def archive_rows(cache, sky, band="H"):
    query = dict(
        POS=f"CIRCLE {sky.ra.deg:.10f} {sky.dec.deg:.10f} 0.0001",
        COLLECTION=COLLECTION,
        MAXREC=500,
    )
    url = f"{SIA}?{urlencode(query)}"
    path = cache.fetch(
        url,
        f"metadata/{TAG}.vot",
        survey="NASA IRSA Euclid DR1",
        version="SIA2 NISP visit stacks",
        kind="observed_visit_inventory",
    )
    statuses = [
        item.attrib.get("value")
        for item in ElementTree.parse(path).findall(".//{*}INFO")
        if item.attrib.get("name") == "QUERY_STATUS"
    ]
    if statuses != ["OK"]:
        raise ValueError(f"Unexpected SIA query status {statuses}")
    table = Table.read(path, format="votable", use_names_over_ids=True)
    visits = []
    for row in table:
        if row["dataproduct_subtype"] != "science" or row["energy_bandpassname"] != band:
            continue
        polygon = str(row["s_region"]).replace("POLYGON ICRS ", "POLYGON ")
        inside, edge = signed_margin_arcsec(polygon, sky)
        start, end = Time([row["t_min"], row["t_max"]], format="mjd", scale="utc").isot
        visits.append(
            dict(
                obs_id=str(row["obs_id"]),
                archive_access_url=str(row["access_url"]),
                start_date_utc=start,
                end_date_utc=end,
                polygon_contains=inside,
                polygon_edge_arcsec=edge,
                nominal_exposure_s=float(row["t_exptime"]),
                s_region=str(row["s_region"]),
                cloud=json.loads(row["cloud_access"])["aws"],
            )
        )
    if not visits or len({row["obs_id"] for row in visits}) != len(visits):
        raise ValueError("Empty or nonunique H-band visit inventory")
    return sorted(visits, key=lambda row: row["start_date_utc"]), len(table), receipt(path, url)


def read_cutout(cache, row, grid, band="H"):
    if not row["polygon_contains"]:
        row["status"] = "outside footprint"
        return None
    remote = EuclidSections(cache, row["cloud"])
    primary, planes = remote.planes(band)
    header, _ = planes["SCI"]
    wcs = WCS(header, fix=False).celestial
    x, y = wcs.world_to_pixel(grid)
    for plane, _ in list(planes.values())[1:]:
        other_x, other_y = WCS(plane, fix=False).celestial.world_to_pixel(grid[CENTER, CENTER])
        if not np.allclose(
            [other_x, other_y], [x[CENTER, CENTER], y[CENTER, CENTER]], atol=1e-5, rtol=0
        ):
            raise ValueError("WCS center disagrees among Euclid image planes")
        if [plane["NAXIS1"], plane["NAXIS2"]] != [header["NAXIS1"], header["NAXIS2"]]:
            raise ValueError("Dimensions disagree among Euclid image planes")
    x0, x1 = int(np.floor(x.min())) - 1, int(np.ceil(x.max())) + 2
    y0, y1 = int(np.floor(y.min())) - 1, int(np.ceil(y.max())) + 2
    row.update(etag=remote.etag, object_bytes=remote.size, file_ranges=remote.ranges)
    if not (0 <= x0 < x1 <= header["NAXIS1"] and 0 <= y0 < y1 <= header["NAXIS2"]):
        row["status"] = "cutout crosses image edge"
        return None
    bounds = (x0, x1, y0, y1)
    # DR1 stacked SCI pixels already include PHOTFNU calibration (Polenta et al., C.3).
    native_area = float(abs(np.linalg.det(wcs.pixel_scale_matrix)) * 3600**2)
    calibration = 1000 * SCALE**2 / native_area
    cut = reproject_cutout(remote, planes, bounds, x, y) * calibration
    row.update(
        resampled_njy_per_input_microjy=calibration,
        native_pixel_area_arcsec2=native_area,
        header_photfnu=primary["PHOTFNU"],
        cutout_valid_fraction=float(np.isfinite(cut).mean()),
        header_exposure_start=primary["DATE_STA"],
        header_exposure_end=primary["DATE_END"],
    )
    if not np.isfinite(aperture(cut, CENTER, CENTER)):
        row["status"] = "masked target aperture or local annulus"
        return None
    row["status"] = "measured"
    return cut


def reproject_cutout(remote, planes, bounds, x, y):
    x0, _, y0, _ = bounds
    data = remote.stripe(planes["SCI"], bounds)
    rms = remote.stripe(planes["RMS"], bounds)
    quality = remote.quality_stripe(planes["DQ"], bounds)
    valid = np.isfinite(data) & np.isfinite(rms) & (rms > 0) & quality
    return map_coordinates(
        np.where(valid, data, np.nan), [y - y0, x - x0], order=1, cval=np.nan, prefilter=False
    )


def aperture(image, x, y, radius=RADIUS):
    yy, xx = np.ogrid[: len(GRID), : len(GRID)]
    distance = np.hypot(xx - x, yy - y) * SCALE
    core = image[distance <= radius]
    ring = image[(distance >= ANNULUS[0]) & (distance <= ANNULUS[1])]
    if not np.isfinite(core).all() or np.isfinite(ring).sum() < 0.8 * len(ring):
        return np.nan
    return float(core.sum() - len(core) * np.nanmedian(ring))


def photometry(image, radius=RADIUS):
    sky = aperture(image, CENTER, CENTER, radius)
    grid = [(x, y) for x in POSITIONS for y in POSITIONS if np.hypot(x - CENTER, y - CENTER) >= 9]
    controls = np.asarray([aperture(image, x, y, radius) for x, y in grid])
    valid = np.flatnonzero(np.isfinite(controls))
    controls = controls[valid]
    if len(controls) < 20 or not np.isfinite(sky):
        return None
    bias = float(np.median(controls))
    scatter = float(1.4826 * np.median(np.abs(controls - bias)))
    return dict(
        signed_aperture_excess_njy=sky - bias,
        empirical_scatter_njy=scatter,
        control_median_njy=bias,
        control_site_count=len(controls),
        residual_over_control_scatter=(sky - bias) / scatter,
        control_grid_indices=valid.tolist(),
        signed_control_residuals_njy=(controls - bias).tolist(),
    )


def visit_outlier(rows, differences, images, campaign_stacks):
    scored = [i for i, row in enumerate(rows) if row["delta"] is not None]
    peak = max(scored, key=lambda i: rows[i]["delta"]["residual_over_control_scatter"])
    top = rows[peak]
    target_score = top["delta"]["residual_over_control_scatter"]
    control_scores = [
        dict(
            zip(
                row["delta"]["control_grid_indices"],
                np.asarray(row["delta"]["signed_control_residuals_njy"])
                / row["delta"]["empirical_scatter_njy"],
                strict=True,
            )
        )
        for row in (rows[i] for i in scored)
    ]
    common = sorted(set.intersection(*(set(row) for row in control_scores)))
    peaks = [max(row[index] for row in control_scores) for index in common]
    peak_grid = control_scores[scored.index(peak)]
    partial = []
    for index in sorted(peak_grid):
        available = [row[index] for row in control_scores if index in row]
        if len(available) >= 0.8 * len(scored):
            partial.append(dict(grid_index=index, visits=len(available), max_score=max(available)))
    same_month = top["start_date_utc"][:7]
    indices = [
        i
        for i, row in enumerate(rows)
        if i != peak and row["start_date_utc"].startswith(same_month)
    ]
    same_month_reference = np.ma.median(np.ma.masked_invalid(images[indices]), axis=0).filled(
        np.nan
    )
    other_campaign = next(stack for month, stack in campaign_stacks.items() if month != same_month)
    comparisons = {
        "all_other_visits": differences[peak],
        "other_visits_in_same_campaign": images[peak] - same_month_reference,
        "other_campaign": images[peak] - other_campaign,
    }
    return dict(
        obs_id=top["obs_id"],
        start_date_utc=top["start_date_utc"],
        residual_over_control_scatter=target_score,
        visits_tested=len(scored),
        common_control_sites=len(common),
        common_control_grid_indices=common,
        control_max_scores=peaks,
        partial_controls_including_peak_visit=partial,
        partial_controls_meeting_or_exceeding_target=sum(
            bool(row["max_score"] >= target_score) for row in partial
        ),
        controls_meeting_or_exceeding_target=sum(bool(score >= target_score) for score in peaks),
        aperture_and_reference_checks={
            name: {str(radius): photometry(image, radius) for radius in (0.60, 0.75, 0.90)}
            for name, image in comparisons.items()
        },
        scope="Largest signed score at the target and at each fixed control site across "
        "the same visit scan. Partial controls have coverage in at least 80% of visits, "
        "including the peak visit, and may sample fewer opportunities for a noise peak. "
        "The ranks are descriptive because the references and control apertures share data. "
        "The aperture and reference variants reuse the peak visit.",
    ), differences[peak]


def plot(visits, campaigns, panels):
    fig = Figure(figsize=(10, 6.4), layout="constrained")
    outer = fig.add_gridspec(2, 1, height_ratios=[1.25, 1])
    windows = outer[0].subgridspec(1, len(campaigns))
    top_axes = [fig.add_subplot(windows[0, 0])]
    top_axes.extend(
        fig.add_subplot(windows[0, i], sharey=top_axes[0]) for i in range(1, len(campaigns))
    )
    for i, (month, axis) in enumerate(zip(campaigns, top_axes, strict=True)):
        rows = [r for r in visits if r["start_date_utc"].startswith(month)]
        valid = [r for r in rows if r["delta"] is not None]
        axis.errorbar(
            [datetime.fromisoformat(r["start_date_utc"]) for r in valid],
            [r["delta"]["signed_aperture_excess_njy"] for r in valid],
            yerr=[r["delta"]["empirical_scatter_njy"] for r in valid],
            fmt="o",
            ms=3.3,
            capsize=2,
            color="#265e75",
        )
        axis.axhline(0, color=".6", linewidth=0.8)
        if i == 0:
            axis.xaxis.set_major_locator(DayLocator())
            axis.xaxis.set_major_formatter(DateFormatter("%d %b"))
        else:
            axis.xaxis.set_major_locator(HourLocator(interval=6))
            axis.xaxis.set_major_formatter(DateFormatter("%d %b\n%H:%M"))
        axis.set(title=month, xlabel="Euclid H-band visit date (UTC)")
        if not i:
            axis.set_ylabel("Signed difference (nJy)")
    observed = outer[1].subgridspec(1, len(panels))
    extent = (-(CENTER + 0.5) * SCALE, (CENTER + 0.5) * SCALE) * 2
    for i, (label, panel) in enumerate(panels.items()):
        axis = fig.add_subplot(observed[0, i])
        limit = float(np.nanpercentile(np.abs(panel), 98))
        axis.imshow(panel, origin="lower", extent=extent, cmap="RdBu_r", vmin=-limit, vmax=limit)
        axis.add_patch(Circle((0, 0), RADIUS, fill=False, color="#272727"))
        axis.set(xlim=(-5.1, 5.1), ylim=(-5.1, 5.1), title=label, xlabel="East offset (arcsec)")
        if not i:
            axis.set_ylabel("North offset (arcsec)")
    fig.suptitle("JADES 103607 | observed Euclid DR1 H-band visits and campaign stacks")
    save_figure(fig, FIGURES / TAG)


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    offset = (GRID * SCALE / 3600) * sky.ra.unit
    east, north = np.meshgrid(offset, offset)
    grid = sky.spherical_offsets_by(east, north)
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    inventory, archive_count, archive = archive_rows(cache, sky)
    images, measured = [], []
    for row in inventory:
        image = read_cutout(cache, row, grid)
        if image is not None:
            images.append(image)
            measured.append(row)
        print(row["obs_id"], row["status"], row.get("cutout_valid_fraction"), flush=True)
    images = np.asarray(images)
    if len(images) < 3:
        raise ValueError("Insufficient usable Euclid visits")
    differences = [difference(images, i) for i in range(len(images))]
    for row, delta in zip(measured, differences, strict=True):
        row["delta"] = photometry(delta)
    months = sorted({row["start_date_utc"][:7] for row in measured})
    stacks = {
        month: np.ma.median(
            np.ma.masked_invalid(images[[r["start_date_utc"].startswith(month) for r in measured]]),
            axis=0,
        ).filled(np.nan)
        for month in months
    }
    if len(months) != 2:
        raise ValueError("Expected the two archived Euclid H-band campaigns")
    early, late = months
    delta = stacks[late] - stacks[early]
    outlier, outlier_image = visit_outlier(measured, differences, images, stacks)
    plot(
        measured,
        months,
        {
            **stacks,
            f"{late} minus {early}": delta,
            f"{outlier['start_date_utc'][:10]} excess": outlier_image,
        },
    )
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            position_source_sha256=text_digest(SOURCE),
            archive=archive,
            sia_rows_after_position_cut=archive_count,
            h_visits_screened=len(inventory),
            h_visits_with_usable_pixels=len(measured),
            h_visits_measured=sum(row["delta"] is not None for row in measured),
            aperture_radius_arcsec=RADIUS,
            annulus_arcsec=ANNULUS,
            visit_inventory=inventory,
            strongest_positive_visit=outlier,
            campaign_stack_difference=dict(
                comparison=f"{late} minus {early}", measurement=photometry(delta)
            ),
            scope="Signed photometry from observed Euclid DR1 H-band visit stacks, using only "
            "positive finite RMS and zero-flag data-quality pixels. Per-visit measurements "
            "use leave-one-out references, and monthly stacks use the median. Local control "
            "scatter includes noise and unresolved light. No PSF matching or aperture "
            "correction, and no constraint outside the returned visits.",
        ),
    )


if __name__ == "__main__":
    main()

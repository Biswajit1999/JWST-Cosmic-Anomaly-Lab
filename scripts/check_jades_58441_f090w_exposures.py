"""Trace the F090W streak near GOODS-N 58441 in calibrated detector exposures."""

import json

import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from audit_jades_58441_imaging import CONTROL_RADIUS, TARGET
from check_jwst_transient_precursor import contains
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from matplotlib.colors import SymLogNorm
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import pixel_scale
from read_jwst_i2d_cutout import S3, MastImageSections, end_of_block
from scipy.ndimage import map_coordinates
from vet_jwst_module_b import archive_query

from cosmic_anomaly_lab.imaging import aperture_measurement
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

HALF_ARCSEC = 2.1
# JWST pixel DQ DO_NOT_USE and SATURATED.
REJECT = 1 | 2
STRICT_REJECT = REJECT | 4 | 8 | 16
STEM = "jades_58441_f090w_exposures"


class CalibratedCutout(MastImageSections):
    """Read exact science, error and DQ stripes from a public detector _cal.fits."""

    def __init__(self, cache, filename):
        self.cache = cache
        self.filename = filename
        self.url = f"{S3}/{filename[:7]}/{filename[:13]}/{filename}"
        response = cache.session.head(self.url, timeout=(15, 30))
        response.raise_for_status()
        self.etag = response.headers["ETag"]
        self.size = int(response.headers["Content-Length"])
        self.ranges = []

    def extract(self, source):
        primary, start = self.header(0)
        if (
            primary["FILENAME"] != self.filename
            or primary["INSTRUME"] != "NIRCAM"
            or primary["FILTER"] != "F090W"
        ):
            raise ValueError(f"Detector version or identity mismatch: {self.filename}")
        extensions = {}
        for _ in range(10):
            h, offset = self.header(start)
            if h["XTENSION"] != "IMAGE" or h["NAXIS"] != 2:
                raise ValueError(f"Missing image planes in {self.filename}")
            if h["EXTNAME"] in ("SCI", "ERR", "DQ"):
                extensions[h["EXTNAME"]] = h, offset
            if len(extensions) == 3:
                break
            start = offset + end_of_block(h["NAXIS1"] * h["NAXIS2"] * abs(h["BITPIX"]) // 8)
        science, _ = extensions["SCI"]
        if science["BUNIT"] != "MJy/sr":
            raise ValueError(f"Unexpected science units in {self.filename}")
        wcs = WCS(science, fix=False).celestial
        x, y = map(float, wcs.world_to_pixel(source))
        scale = pixel_scale(wcs)
        half = int(np.ceil(HALF_ARCSEC / scale))
        ix, iy = round(x), round(y)
        nx, ny = science["NAXIS1"], science["NAXIS2"]
        if min(ix, iy, nx - ix - 1, ny - iy - 1) < half:
            raise ValueError("Insufficient detector coverage at target")
        x0, x1, y0, y1 = (ix - half, ix + half + 1, iy - half, iy + half + 1)
        data, error = (self.stripe(extensions[plane], (x0, x1, y0, y1)) for plane in ("SCI", "ERR"))
        dq_h, dq_offset = extensions["DQ"]
        if dq_h["BITPIX"] != 32 or dq_h.get("BSCALE", 1) != 1 or dq_h["NAXIS1"] != nx:
            raise ValueError("Unexpected detector DQ layout")
        a = self.read(dq_offset + y0 * 4 * nx, dq_offset + y1 * 4 * nx)
        raw = np.memmap(a, dtype=">i4", mode="r", shape=(y1 - y0, nx))[:, x0:x1]
        flags = (np.array(raw, dtype=np.int64) + dq_h.get("BZERO", 0)).astype(np.uint32)
        valid = (flags & REJECT == 0) & np.isfinite(data) & np.isfinite(error) & (error > 0)
        weights = np.divide(1.0, error.astype(float) ** 2, out=np.zeros(error.shape), where=valid)
        jy_per_unit = science["PIXAR_SR"] * 1e6
        meta = dict(
            file=self.filename,
            url=self.url,
            upstream_etag=self.etag,
            object_bytes=self.size,
            detector=primary["DETECTOR"],
            date_utc=primary["DATE-BEG"],
            start_mjd=float(primary["EXPSTART"]),
            exposure_seconds=float(primary["EFFEXPTM"]),
            pipeline=primary["CAL_VER"],
            crds_context=primary["CRDS_CTX"],
            ranges=self.ranges,
            native_pixel_scale_arcsec=scale,
            detector_xy=[x, y],
            cutout_bounds_x0_x1_y0_y1=[x0, x1, y0, y1],
        )
        shifted_wcs = wcs.slice((slice(y0, y1), slice(x0, x1)))
        return meta, data, weights, flags, shifted_wcs, jy_per_unit


def measure(data, weights, wcs, scale, centre, controls):
    points = {}
    for label, pos in (("target", centre), *controls.items()):
        try:
            fit = aperture_measurement(
                data,
                weights,
                wcs,
                pos.ra.deg,
                pos.dec.deg,
                jy_per_unit=scale,
                radius_arcsec=CONTROL_RADIUS,
            )
        except ValueError as exc:
            points[label] = dict(status="UNUSABLE", reason=str(exc))
        else:
            points[label] = dict(
                status="MEASURED",
                flux_njy=fit["aperture_flux_uJy"] * 1000,
                err_njy=fit["diagonal_error_uJy"] * 1000,
                background_njy_pixel=fit["background_uJy_pixel"] * 1000,
                aperture_pixels=fit["n_aperture_pixels"],
            )
    along = [v for k, v in points.items() if k.startswith("along")]
    across = [v for k, v in points.items() if k.startswith("across")]
    contrast = None
    if all(v["status"] == "MEASURED" for v in along + across):
        contrast = dict(
            flux_njy=float(
                np.mean([v["flux_njy"] for v in along]) - np.mean([v["flux_njy"] for v in across])
            ),
            diagonal_error_njy=float(
                np.linalg.norm([v["err_njy"] for v in along + across]) / len(along)
            ),
        )
    return points, contrast


def plot(results, images):
    fig = Figure(figsize=(10.8, 8.5), layout="constrained")
    grid = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.5])
    representatives = {}
    for row in results:
        contrast = row.get("directional_contrast")
        program = row["exposure"][:7]
        if contrast:
            previous = representatives.get(program)
            if previous is None or abs(contrast["flux_njy"]) > abs(
                previous["directional_contrast"]["flux_njy"]
            ):
                representatives[program] = row
    for i, (program, row) in enumerate(representatives.items()):
        ax = fig.add_subplot(grid[:2, i])
        ax.imshow(
            images[row["exposure"]],
            origin="lower",
            extent=(-1.1, 1.1, -1.1, 1.1),
            cmap="gray",
            norm=SymLogNorm(linthresh=0.3, vmin=-2, vmax=12),
        )
        ax.scatter(0, 0, facecolor="none", edgecolor="#48dbf6", s=120)
        ax.set(
            xlabel="East offset (arcsec)", title=f"PID {int(program[2:])} · {row['date_utc'][:10]}"
        )
        if i == 0:
            ax.set_ylabel("North offset (arcsec)")
    ax = fig.add_subplot(grid[2, :])
    labelled_streak = labelled_target = False
    for i, row in enumerate(results):
        contrasts, target = row.get("directional_contrast"), row.get("apertures", {}).get("target")
        if contrasts:
            ax.errorbar(
                i,
                contrasts["flux_njy"],
                yerr=contrasts["diagonal_error_njy"],
                fmt="o",
                capsize=2,
                color="#9e5632",
                label="Along minus across" if not labelled_streak else None,
            )
            labelled_streak = True
        if target and target["status"] == "MEASURED":
            ax.errorbar(
                i,
                target["flux_njy"],
                yerr=target["err_njy"],
                fmt="x",
                capsize=2,
                color="#246f81",
                label="Target aperture" if not labelled_target else None,
            )
            labelled_target = True
    ax.axhline(0, lw=0.7, color="gray")
    if labelled_target or labelled_streak:
        ax.legend(frameon=False)
    ax.set(
        xlabel="Detector exposure, in time order",
        ylabel="F090W flux (nJy)",
        xticks=range(len(results)),
        xticklabels=range(1, len(results) + 1),
    )
    fig.suptitle("GOODS-N 58441 · F090W across independent JWST exposures", fontsize=14)
    save_figure(fig, ROOT / f"figures/{STEM}")


def main():
    cache = DownloadCache(ROOT, 28_000 * 1024**2)
    target = next(r for r in json.loads(LEDGER.read_text())["records"] if r["dr4_id"] == TARGET)
    source = SkyCoord(target["ra"], target["dec"], unit="deg")
    prior_path = OUT / "58441_optical_ir_context.json"
    prior = json.loads(prior_path.read_text())
    image = ROOT / prior["input_images"]["path"]
    if digest(image) != prior["input_images"]["sha256"]:
        raise ValueError("Reference image checksum changed")
    with fits.open(image) as hdus:
        wcs = WCS(hdus[0].header).celestial
    controls = {
        f"{axis}{i}": wcs.pixel_to_world(p["x_pixel"], p["y_pixel"])
        for axis, points in prior["bands"][0]["streak_control"]["directional_apertures"].items()
        for i, p in enumerate(points)
    }
    table, inventory = archive_query(
        cache,
        "exposures",
        dict(ra=source.ra.deg, dec=source.dec.deg, arcmin=0.025, output="csv"),
        "goods_n_58441_exposure_inventory.csv",
    )
    matches = table[(table.extension == "rate") & table["filter"].eq("F090W-CLEAR")]
    matches = matches[matches.footprint.map(lambda footprint: contains(footprint, source))]
    if matches.empty:
        raise ValueError("No usable exposure metadata")
    results, images = [], {}
    radius = np.linspace(-1.1, 1.1, 75)
    east, north = np.meshgrid(radius, radius)
    grid = source.spherical_offsets_by(east * u.arcsec, north * u.arcsec)
    for row in matches.sort_values(["expstart", "file"]).itertuples():
        name = row.file + "_cal.fits"
        entry = dict(
            exposure=row.file,
            inventory_mjd=float(row.expstart),
            inventory_date_utc=Time(row.expstart, format="mjd").isot,
        )
        try:
            meta, data, weights, flags, wcs, scale = CalibratedCutout(cache, name).extract(source)
            if meta["detector"] != row.detector or abs(meta["start_mjd"] - row.expstart) > 0.005:
                raise ValueError(f"Detector or timestamp conflicts with exposure inventory: {name}")
            points, contrast = measure(data, weights, wcs, scale, source, controls)
            strict_weights = np.where(flags & STRICT_REJECT == 0, weights, 0)
            strict_points, strict_contrast = measure(
                data, strict_weights, wcs, scale, source, controls
            )
            x, y = wcs.world_to_pixel(grid)
            bg = points["target"].get("background_njy_pixel", np.nan)
            im = np.where(weights > 0, data * scale * 1e9 - bg, np.nan)
            images[row.file] = map_coordinates(im, [y, x], order=1, mode="constant", cval=np.nan)
            entry.update(
                meta,
                status="MEASURED",
                apertures=points,
                directional_contrast=contrast,
                strict_dq=dict(apertures=strict_points, directional_contrast=strict_contrast),
                valid_pixel_fraction=float((weights > 0).mean()),
                strict_valid_pixel_fraction=float((strict_weights > 0).mean()),
            )
        except (requests.RequestException, ValueError, KeyError) as exc:
            entry.update(status="UNAVAILABLE", reason=str(exc))
        results.append(entry)
        print(
            row.file,
            entry["status"],
            entry.get("directional_contrast", entry.get("reason")),
            flush=True,
        )
    summaries = []
    for visit in dict.fromkeys(row["exposure"][:13] for row in results):
        group = [row for row in results if row["exposure"].startswith(visit)]
        target_flux = [
            r["apertures"]["target"]["flux_njy"]
            for r in group
            if r.get("apertures", {}).get("target", {}).get("status") == "MEASURED"
        ]
        streak = [
            r["directional_contrast"]["flux_njy"] for r in group if r.get("directional_contrast")
        ]
        summaries.append(
            dict(
                visit=visit,
                exposures=len(group),
                targets=len(target_flux),
                median_target_flux_njy=float(np.median(target_flux)) if target_flux else None,
                target_flux_range_njy=[min(target_flux), max(target_flux)] if target_flux else None,
                contrasts=len(streak),
                median_directional_contrast_njy=float(np.median(streak)) if streak else None,
                contrast_range_njy=[min(streak), max(streak)] if streak else None,
            )
        )
    if images:
        plot(results, images)
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            inventory=inventory,
            reference_sha256=text_digest(prior_path),
            position_icrs_deg=[source.ra.deg, source.dec.deg],
            control_positions_icrs_deg={k: [p.ra.deg, p.dec.deg] for k, p in controls.items()},
            apertures_radius_arcsec=CONTROL_RADIUS,
            rejected_dq_bits=REJECT,
            strict_rejected_dq_bits=STRICT_REJECT,
            exposures=results,
            visit_summaries=summaries,
            scope="Exact public calibrated detector pixels. Same sky apertures across exposures. "
            "Formal ERR-plane uncertainties omit shared sky/flat-field systematics. "
            "North-up display uses interpolation; measurements use native pixels.",
        ),
    )


if __name__ == "__main__":
    main()

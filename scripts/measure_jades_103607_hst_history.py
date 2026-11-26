"""Measure historical F160W variability at JADES 103607 in CANDELS epochs."""

import argparse
import json
import re
from io import BytesIO

import numpy as np
from astropy.constants import c
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from read_jwst_i2d_cutout import MastImageSections
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

BASE = "https://archive.stsci.edu/hlsps/candels/goods-s"
SOURCE = OUT / "module_b_001_native_november_fade.json"
STEM = "module_b_001_hst_f160w_history"
EPOCHS = (
    "gsd01",
    "gsd02",
    "gsw03",
    "gsd04",
    "gsw05",
    "gsd06",
    "gsw07",
    "gsd08",
    "gsd10",
    "gsd11",
    "gsd12",
    "gsd14",
    "gsd15",
)
# The four final public products omit observation dates in their primary headers.
DOCUMENTED_DATES = {
    "gsd11": ("2011-09-10", "2011-09-22"),
    "gsd12": ("2011-11-02", "2011-11-07"),
    "gsd14": ("2011-12-22", "2011-12-30"),
    "gsd15": ("2012-02-14", "2012-02-20"),
}
SCALE = 0.06
RADIUS = 0.24
ANNULUS = (0.48, 0.72)
GRID = np.arange(-70, 71)


class CandelsSections(MastImageSections):
    """Read the primary image of a released CANDELS mosaic without fetching its full FITS."""

    def __init__(self, cache, epoch, product, band="F160W"):
        self.cache = cache
        self.filename = f"hlsp_candels_hst_wfc3_{epoch}_{band.lower()}_v0.5_{product}.fits"
        self.url = f"{BASE}/{epoch}/v0.5/{self.filename}"
        with cache.session.head(self.url, timeout=(15, 40)) as response:
            response.raise_for_status()
            self.etag = response.headers["ETag"]
            self.size = int(response.headers["Content-Length"])
        self.ranges = []

    def primary(self):
        with BytesIO(self.read(0, 230400).read_bytes()) as raw:
            header = fits.Header.fromfile(raw)
            return header, raw.tell()


def positions(header, grid):
    # The drizzled canvases have a TAN projection. Older headers retain SIP
    # keywords from detector processing even though their output CTYPE is TAN.
    header = header.copy()
    for key in set(header):
        if re.fullmatch(r"(A|B|AP|BP)_(ORDER|DMAX|[0-9]+_[0-9]+)", key):
            del header[key]
    wcs = WCS(header, fix=False)
    if tuple(wcs.wcs.ctype) != ("RA---TAN", "DEC--TAN"):
        raise ValueError("Unexpected CANDELS output projection")
    x, y = wcs.wcs_world2pix(grid.ra.deg, grid.dec.deg, 0)
    return np.asarray(x), np.asarray(y)


def read_epoch(cache, epoch, grid, band="F160W"):
    raw = {k: CandelsSections(cache, epoch, k, band) for k in ("wht", "drz")}
    wht_h, wht_start = raw["wht"].primary()
    x, y = positions(wht_h, grid)
    m = len(GRID) // 2
    row_width = 4 * wht_h["NAXIS1"]
    px, py = int(round(x[m, m])), int(round(y[m, m]))
    origin = wht_start + row_width * py + px * 4
    central = float(np.fromfile(raw["wht"].read(origin, origin + 4), dtype=">f4")[0])
    record = {
        "epoch": epoch,
        "header_start_date": wht_h.get("DATE-OBS"),
        "documented_window_utc": DOCUMENTED_DATES.get(epoch),
        "date_documentation": f"{BASE}/{epoch}/v0.5/hlsp_candels_hst_acs-wfc3_{epoch}_readme_v0.5.pdf"
        if epoch in DOCUMENTED_DATES
        else None,
        "center_weight": central,
        "ranges": raw["wht"].ranges,
    }
    if not central > 0:
        return record, None
    sci_h, sci_start = raw["drz"].primary()
    sx, sy = positions(sci_h, grid)
    if not np.allclose([x, y], [sx, sy], rtol=0, atol=0.001):
        raise ValueError(f"Science and weight WCS disagree in {epoch}")
    rectangle = (
        int(np.floor(x.min())) - 2,
        int(np.ceil(x.max())) + 3,
        int(np.floor(y.min())) - 2,
        int(np.ceil(y.max())) + 3,
    )
    xp, _, yp, _ = rectangle
    wht = raw["wht"].stripe((wht_h, wht_start), rectangle)
    counts = raw["drz"].stripe((sci_h, sci_start), rectangle)
    valid = np.isfinite(counts) & (wht > 0)
    factor_njy = sci_h["PHOTFLAM"] * sci_h["PHOTPLAM"] ** 2 / c.to_value("Angstrom/s") * 1e32
    cut = (
        map_coordinates(
            np.where(valid, counts, np.nan), [sy - yp, sx - xp], order=1, prefilter=False
        )
        * factor_njy
    )
    record.update(
        dates_from_sci_header={k: sci_h.get(k) for k in ("DATE-OBS", "DATE-END", "MJD-OBS")},
        calibration_njy_per_electron_second=factor_njy,
        ranges=raw["wht"].ranges + raw["drz"].ranges,
    )
    return record, cut


def aperture(data, x, y):
    yy, xx = np.ogrid[: len(GRID), : len(GRID)]
    r = np.hypot(xx - x, yy - y) * SCALE
    core, sky = data[r <= RADIUS], data[(r >= ANNULUS[0]) & (r < ANNULUS[1])]
    if not np.isfinite(core).all() or np.isfinite(sky).sum() < 0.8 * len(sky):
        return np.nan
    return float(core.sum() - len(core) * np.nanmedian(sky))


def control_centers():
    n = len(GRID) // 2
    return [
        (x, y)
        for x in range(14, len(GRID) - 14, 15)
        for y in range(14, len(GRID) - 14, 15)
        if np.hypot(x - n, y - n) > 16
    ]


def difference(image_stack, index):
    reference = np.ma.median(
        np.ma.masked_invalid(np.delete(image_stack, index, axis=0)), axis=0
    ).filled(np.nan)
    return image_stack[index] - reference


def measure(records, images):
    center = len(GRID) // 2
    image_stack = np.array(images)
    for i, (record, sample) in enumerate(zip(records, images, strict=True)):
        diff = difference(image_stack, i)
        controls = np.asarray([aperture(diff, x, y) for x, y in control_centers()])
        controls = controls[np.isfinite(controls)]
        target = aperture(diff, center, center)
        flux = aperture(sample, center, center)
        if len(controls) < 20 or not np.isfinite(target + flux):
            raise ValueError(f"Local photometry at {record['epoch']} has insufficient coverage")
        med = float(np.median(controls))
        shifted = [
            aperture(diff, center + dx, center + dy) - med for dx in (-1, 0, 1) for dy in (-1, 0, 1)
        ]
        finite_shifted = np.asarray(shifted)[np.isfinite(shifted)]
        record.update(
            centering_check_one_pixel_njy=[
                float(finite_shifted.min()),
                float(finite_shifted.max()),
            ],
            raw_aperture_njy=flux,
            target_minus_leave_one_out_njy=target,
            differential_njy=target - med,
            control_median_njy=med,
            controls_njy=controls.tolist(),
            empirical_scatter_njy=float(1.4826 * np.median(abs(controls - med))),
        )


def epoch_start(record):
    if record["documented_window_utc"]:
        return record["documented_window_utc"][0]
    return record["header_start_date"]


def figure(records, images, band, tag):
    fig = Figure(figsize=(10, 5.6), layout="constrained")
    layout = fig.add_gridspec(2, 3, height_ratios=[1.4, 1])
    timeline = fig.add_subplot(layout[0, :])
    stamps = [fig.add_subplot(layout[1, i]) for i in range(3)]
    dates = [Time(epoch_start(r)).decimalyear for r in records]
    timeline.errorbar(
        dates,
        [r["differential_njy"] for r in records],
        yerr=[r["empirical_scatter_njy"] for r in records],
        fmt="o",
        ms=4,
        capsize=3,
        color="#2e697c",
    )
    timeline.axhline(0, color=".53", lw=0.7)
    timeline.set(xlabel="HST epoch start year", ylabel=f"{band} aperture excess (nJy)")
    center = len(GRID) // 2
    selection = [0, len(records) // 2, len(records) - 1]
    lim = np.nanpercentile(
        np.asarray(images)[:, center - 31 : center + 32, center - 31 : center + 32], 99
    )
    for ax, i in zip(stamps, selection, strict=True):
        cut = images[i][center - 31 : center + 32, center - 31 : center + 32]
        ax.imshow(
            cut,
            origin="lower",
            cmap="magma",
            vmin=-0.3,
            vmax=max(lim, 0.5),
            extent=(-1.89, 1.89, -1.89, 1.89),
        )
        ax.scatter(0, 0, facecolors="none", edgecolors="#60d6de", s=34, lw=0.9)
        ax.set(title=epoch_start(records[i]), xlabel="East offset (arcsec)")
    stamps[0].set_ylabel("North offset (arcsec)")
    fig.suptitle("JADES 103607 · historical HST/CANDELS observations", fontsize=13)
    save_figure(fig, FIGURES / tag)


def load_series(cache, sky, band):
    east, north = np.meshgrid(*(GRID * SCALE / 3600,) * 2)
    grid = sky.spherical_offsets_by(east * sky.ra.unit, north * sky.dec.unit)
    coverage, images = [], []
    for epoch in EPOCHS:
        row, image = read_epoch(cache, epoch, grid, band)
        coverage.append(row)
        if image is not None:
            images.append(image)
        print(band, epoch, row["header_start_date"], row["center_weight"], flush=True)
    return coverage, np.asarray(images)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--band", choices=("F160W", "F125W"), default="F160W")
    band = parser.parse_args().band
    tag = STEM.replace("f160w", band.lower())
    position = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    coverage, images = load_series(cache, sky, band)
    records = [record for record in coverage if record["center_weight"] > 0]
    measure(records, images)
    write_json(
        OUT / f"{tag}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            position_source_sha256=text_digest(SOURCE),
            survey="HST CANDELS GOODS-S individual epoch mosaics",
            band=band,
            aperture_radius_arcsec=RADIUS,
            sky_annulus_arcsec=ANNULUS,
            coverage=coverage,
            scope="Same-band differential rates and local empirical scatter. Unmodelled PSF differences, shared "
            "references, filter differences relative to JWST and gaps in cadence limit recurrence conclusions.",
        ),
    )
    figure(records, images, band, tag)


if __name__ == "__main__":
    main()

"""Compare the first five GOODS-S ACS F850LP epochs at JADES 103607."""

import json
from warnings import catch_warnings, simplefilter

import numpy as np
from astropy import units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs import WCS
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

SOURCE = OUT / "module_b_001_native_november_fade.json"
BASE = "https://archive.stsci.edu/pub/hlsp/goods/v05"
STEM = "module_b_001_hst_acs_goods_epochs"
SCALE = 0.05
PIXEL = np.arange(-80, 81)
RADIUS = 0.2
SKY_RADII = (0.40, 0.65)
# The full set of 77 public F850LP tile headers was screened with the released
# TAN WCS. Prefer the tile furthest inside the detector at every epoch.
TILES = {1: (3, 2), 2: (25,), 3: (3, 2), 4: (25, 26), 5: (3, 2)}


def stamp(cache, ep, tile, grid):
    products, paths = [], {}
    for kind in ("wht", "sci"):
        name = f"h_s{ep}z{tile:02}{kind}v05_img.fits.gz"
        url = f"{BASE}/{name}"
        paths[kind] = cache.fetch(
            url,
            f"cross_detector/{name}",
            survey="HST GOODS-S",
            version="GOODS ACS v0.5",
            kind="observed_epoch_fits",
        )
        products.append(receipt(paths[kind], url))
    with (
        fits.open(paths["sci"], memmap=False) as science,
        fits.open(paths["wht"], memmap=False) as weights,
        catch_warnings(),
    ):
        simplefilter("ignore", VerifyWarning)
        counts, weight = science[0].data, weights[0].data
        sci = science[0].header
        wcs = WCS(sci, fix=False).celestial
        scale = abs(np.linalg.det(wcs.pixel_scale_matrix)) ** 0.5 * 3600
        if abs(scale - SCALE) > 0.001 or tuple(wcs.wcs.ctype) != ("RA---TAN", "DEC--TAN"):
            raise ValueError("Unexpected GOODS v0.5 science WCS")
        if counts.shape != weight.shape or int(sci["EXPTIME"]) != 1:
            raise ValueError("GOODS products do not match count-rate tile format")
        # Legacy weight WCS cards can precede the science astrometric correction.
        # Verify the native array association from the shared missing-data mask.
        zeros = counts == 0
        match = float(np.mean(weight[zeros] <= 0))
        if match < 0.97:
            raise ValueError("Native science gaps do not agree with the weight image")
        x, y = wcs.wcs_world2pix(grid.ra.deg, grid.dec.deg, 0)
        xmin, xmax = int(np.floor(x.min())) - 2, int(np.ceil(x.max())) + 3
        ymin, ymax = int(np.floor(y.min())) - 2, int(np.ceil(y.max())) + 3
        n = len(PIXEL) // 2
        if not (0 <= xmin < xmax <= counts.shape[1] and 0 <= ymin < ymax <= counts.shape[0]):
            return dict(tile=tile, status="local grid crosses tile edge", products=products), None
        if not weight[round(y[n, n]), round(x[n, n])] > 0:
            return dict(tile=tile, status="zero central weight", products=products), None
        region = np.s_[ymin:ymax, xmin:xmax]
        prepared = np.where(weight[region] > 0, counts[region], np.nan).astype(float)
        factor = sci["PHOTFLAM"] * sci["PHOTPLAM"] ** 2 / c.to_value("Angstrom / s") * 1e32
        output = map_coordinates(prepared, [y - ymin, x - xmin], order=1, prefilter=False) * factor
        return dict(
            tile=tile,
            status="covered",
            products=products,
            first_input_date=sci["DATE-OBS"],
            combined_exposure_s=sci.get("TEXPTIME"),
            calibration_njy_per_count_second=factor,
            fraction_of_zero_science_pixels_with_zero_weight=match,
        ), output


def aperture(image, x, y):
    yy, xx = np.ogrid[: len(PIXEL), : len(PIXEL)]
    r = np.hypot(xx - x, yy - y) * SCALE
    core = image[r <= RADIUS]
    sky = image[(r >= SKY_RADII[0]) & (r <= SKY_RADII[1])]
    if not np.isfinite(core).all() or np.isfinite(sky).sum() < 0.8 * len(sky):
        return np.nan
    return float(core.sum() - len(core) * np.nanmedian(sky))


def compare(records, images):
    n = len(PIXEL) // 2
    for i, record in enumerate(records):
        reference = np.ma.median(np.ma.masked_invalid(np.delete(images, i, axis=0)), axis=0).filled(
            np.nan
        )
        diff = images[i] - reference
        controls = [
            aperture(diff, x, y)
            for x in range(15, 146, 14)
            for y in range(15, 146, 14)
            if np.hypot(x - n, y - n) >= 20
        ]
        controls = np.asarray(controls)
        controls = controls[np.isfinite(controls)]
        target = aperture(diff, n, n)
        if len(controls) < 25 or not np.isfinite(target):
            raise ValueError(f"ACS epoch {record['epoch']} has insufficient valid photometry")
        center = float(np.median(controls))
        record.update(
            blended_aperture_njy=aperture(images[i], n, n),
            excess_njy=target - center,
            control_median_njy=center,
            empirical_scatter_njy=float(1.4826 * np.median(abs(controls - center))),
            controls_njy=controls.tolist(),
        )


def plot(records, images):
    fig = Figure(figsize=(10, 5.0), layout="constrained")
    layout = fig.add_gridspec(2, len(records), height_ratios=(1.55, 1))
    graph = fig.add_subplot(layout[0, :])
    graph.errorbar(
        range(len(records)),
        [row["excess_njy"] for row in records],
        yerr=[row["empirical_scatter_njy"] for row in records],
        color="#317386",
        fmt="o",
        ms=4,
        capsize=3,
    )
    graph.axhline(0, color=".52", lw=0.8)
    graph.set(
        xticks=range(len(records)),
        xticklabels=[r["selected"]["first_input_date"] for r in records],
        ylabel="F850LP aperture excess (nJy)",
    )
    n = len(PIXEL) // 2
    ex = 32
    v = np.nanpercentile(images[:, n - ex : n + ex + 1, n - ex : n + ex + 1], 99.3)
    for i, image in enumerate(images):
        ax = fig.add_subplot(layout[1, i])
        edge = (ex + 0.5) * SCALE
        ax.imshow(
            image[n - ex : n + ex + 1, n - ex : n + ex + 1],
            origin="lower",
            cmap="magma",
            extent=(-edge, edge, -edge, edge),
            vmin=-v / 3,
            vmax=v,
        )
        ax.add_patch(Circle((0, 0), RADIUS, fill=False, ec="#74d6dc", lw=1))
        ax.set(xticks=[], yticks=[], title=f"Epoch {i + 1}")
    fig.suptitle("JADES 103607 · 2002–2003 Hubble ACS imaging", fontsize=13)
    save_figure(fig, FIGURES / STEM)


def main():
    position = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    east, north = np.meshgrid(PIXEL * SCALE, PIXEL * SCALE)
    grid = sky.spherical_offsets_by(east * u.arcsec, north * u.arcsec)
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    records, images = [], []
    for ep, options in TILES.items():
        attempted = []
        for candidate in options:
            row, data = stamp(cache, ep, candidate, grid)
            attempted.append(row)
            n = len(PIXEL) // 2
            if data is not None and np.isfinite(aperture(data, n, n)):
                break
        else:
            raise ValueError(f"No GOODS ACS tile with usable target photometry for epoch {ep}")
        records.append(dict(epoch=ep, selected=row, attempted=attempted))
        images.append(data)
        print(ep, row["tile"], row["first_input_date"], flush=True)
    stack = np.asarray(images)
    compare(records, stack)
    plot(records, stack)
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            position_source_sha256=text_digest(SOURCE),
            aperture_radius_arcsec=RADIUS,
            background_annulus_arcsec=list(SKY_RADII),
            epochs=records,
            scope="Five early same-band GOODS epochs. F850LP and the JWST near-infrared bands sample different "
            "wavelengths. No point-source aperture correction or PSF matching was applied.",
        ),
    )


if __name__ == "__main__":
    main()

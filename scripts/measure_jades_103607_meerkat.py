"""Measure historical MeerKAT continuum at the JWST transient position."""

import json

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from read_jwst_i2d_cutout import MastImageSections
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

SOURCE = OUT / "module_b_001_native_november_fade.json"
TAG = "module_b_001_meerkat"
DATA = "https://archive-gw-1.kat.ac.za/public/repository/10.48479/7msw-r692/data"
PREFIX = "MIGHTEE_Continuum_DR1_CDFS-DEEP_5p5arcsec_"
OBSERVATION_DATES = ("2019-12-12", "2020-10-17")
CATALOGUE_RADIUS_ARCSEC = 60
CUTOUT_RADIUS_ARCSEC = 83


class SaraoSections(MastImageSections):
    """Read only the needed uncompressed FITS rows of a published MeerKAT mosaic."""

    range_survey = "SARAO MIGHTEE DR1"
    range_cache_prefix = "cross_detector/sarao_stripes"

    def __init__(self, cache, suffix):
        self.cache = cache
        self.filename = PREFIX + suffix
        self.url = f"{DATA}/{self.filename}"
        with cache.session.head(self.url, timeout=(15, 30)) as response:
            response.raise_for_status()
            self.etag = response.headers["ETag"]
            self.size = int(response.headers["Content-Length"])
        self.ranges = []

    def local_image(self, sky):
        hdr, start = self.header(0)
        # The released RMS header uses an upper-case spelling on a singleton axis.
        if hdr.get("CUNIT4") == "HZ":
            hdr["CUNIT4"] = "Hz"
        wcs = WCS(hdr, fix=False, naxis=2)
        x, y = map(float, wcs.world_to_pixel(sky))
        scale = float(np.sqrt(np.prod(proj_plane_pixel_scales(wcs))) * 3600)
        radius = int(np.ceil(CUTOUT_RADIUS_ARCSEC / scale))
        x0, x1, y0, y1 = (
            round(x) - radius,
            round(x) + radius + 1,
            round(y) - radius,
            round(y) + radius + 1,
        )
        if min(x0, y0) < 0 or x1 > hdr["NAXIS1"] or y1 > hdr["NAXIS2"]:
            raise ValueError(f"Cutout outside {self.filename}")
        image_header = hdr.copy()
        extra_axes = [hdr[f"NAXIS{i}"] for i in range(3, hdr["NAXIS"] + 1)]
        if any(size != 1 for size in extra_axes):
            raise ValueError(f"Expected a single image plane in {self.filename}")
        image_header["NAXIS"] = 2
        image = self.stripe((image_header, start), (x0, x1, y0, y1))
        return image, hdr, wcs, (x - x0, y - y0), (x0, y0), scale


def sample(image, x, y):
    return float(map_coordinates(image, [[y], [x]], order=1, prefilter=False)[0])


def nearest_catalogue(cache, sky):
    url = f"{DATA}/{PREFIX}I_v1.1_FinalCatalogue.srl.fits"
    path = cache.fetch(
        url,
        f"cross_detector/{PREFIX}I_v1.1_FinalCatalogue.srl.fits",
        survey="SARAO MIGHTEE DR1",
        version="CDFS-DEEP 5.5 arcsec v1.1",
        kind="observed_source_catalogue",
    )
    cat = Table.read(path, format="fits")
    positions = SkyCoord(cat["RA"], cat["DEC"], unit="deg")
    sep = sky.separation(positions).arcsec
    keep = np.flatnonzero(sep <= CATALOGUE_RADIUS_ARCSEC)
    nearest = int(np.argmin(sep))
    rows = []
    for i in sorted(set(keep) | {nearest}, key=lambda j: sep[j]):
        row = cat[i]
        rows.append(
            dict(
                catalogue_id=int(row["Source_id"]),
                separation_arcsec=float(sep[i]),
                icrs_deg=[float(row["RA"]), float(row["DEC"])],
                peak_ujy_beam=float(row["Peak_flux"]) * 1e6,
                peak_error_ujy_beam=float(row["E_Peak_flux"]) * 1e6,
                total_ujy=float(row["Total_flux"]) * 1e6,
            )
        )
    return rows, len(cat), receipt(path, url)


def plot(image, rms, center, scale, neighbours, sky, beam):
    fig = Figure(figsize=(8.4, 4.6), layout="constrained")
    axes = fig.subplots(1, 2, sharex=True, sharey=True)
    x, y = center
    ny, nx = image.shape
    extent = (
        (-x - 0.5) * scale,
        (nx - x - 0.5) * scale,
        (-y - 0.5) * scale,
        (ny - y - 0.5) * scale,
    )
    radii = (-38, 38)
    image_scale = float(np.nanmedian(rms)) * 1e6
    for axis, values, title, lim, label in zip(
        axes,
        (image * 1e6, image / rms),
        ("Stokes I", "Map divided by local RMS"),
        (image_scale, 1),
        ("µJy/beam", "RMS units"),
        strict=True,
    ):
        panel = axis.imshow(
            values, origin="lower", extent=extent, cmap="RdBu_r", vmin=-5 * lim, vmax=5 * lim
        )
        axis.scatter(0, 0, marker="+", color="#191919", s=75, label="JWST transient")
        for row in neighbours:
            coord = SkyCoord(*row["icrs_deg"], unit="deg")
            east, north = sky.spherical_offsets_to(coord)
            if radii[0] < -east.arcsec < radii[1] and radii[0] < north.arcsec < radii[1]:
                axis.scatter(
                    -east.arcsec, north.arcsec, s=50, facecolors="none", edgecolors="#5b284d"
                )
        axis.add_patch(
            Circle((-31, -31), radius=beam / 2, fill=False, color="#333333", linewidth=1)
        )
        axis.set(xlim=radii, ylim=radii, title=title)
        fig.colorbar(panel, ax=axis, shrink=0.83, label=label)
    axes[0].legend(loc="upper right", frameon=True, fontsize=8)
    fig.supxlabel("West offset (arcsec)")
    fig.supylabel("North offset (arcsec)")
    fig.suptitle("JADES 103607 | observed MeerKAT MIGHTEE CDFS-DEEP")
    save_figure(fig, FIGURES / TAG)


def main():
    position = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    products = {
        "intensity": SaraoSections(cache, "I_v1.1.fits"),
        "rms": SaraoSections(cache, "rms_v1.1.fits"),
        "effective_frequency": SaraoSections(cache, "effective-freq_v1.1.fits"),
    }
    image, header, wcs, center, origin, scale = products["intensity"].local_image(sky)
    rms, _, rms_wcs, rms_center, rms_origin, _ = products["rms"].local_image(sky)
    freq, freq_header, _, freq_center, _, _ = products["effective_frequency"].local_image(sky)
    if header["BUNIT"] != "Jy/beam" or freq_header["BUNIT"] != "Hz":
        raise ValueError("Unexpected MIGHTEE map units")
    check_x = [origin[0], origin[0] + image.shape[1] - 1]
    check_y = [origin[1], origin[1] + image.shape[0] - 1]
    discrepancy = (
        wcs.pixel_to_world(check_x, check_y)
        .separation(rms_wcs.pixel_to_world(check_x, check_y))
        .arcsec
    )
    max_discrepancy = float(np.max(discrepancy))
    if max_discrepancy > scale:
        raise ValueError("Radio image and noise-map coordinates differ by more than one pixel")
    iy, ix = np.indices(image.shape, dtype=float)
    locations = wcs.pixel_to_world(ix + origin[0], iy + origin[1])
    rms_x, rms_y = rms_wcs.world_to_pixel(locations)
    rms_display = map_coordinates(
        rms,
        [rms_y - rms_origin[1], rms_x - rms_origin[0]],
        order=1,
        mode="constant",
        cval=np.nan,
        prefilter=False,
    )
    peak = sample(image, *center) * 1e6
    noise = sample(rms, *rms_center) * 1e6
    effective_ghz = sample(freq, *freq_center) / 1e9
    if not np.isfinite([peak, noise, effective_ghz]).all() or noise <= 0:
        raise ValueError("Invalid local MIGHTEE measurement")
    beam = float(header["BMAJ"] * 3600)
    x, y = center
    offsets = np.arange(-70, 71, 10)
    values = np.asarray(
        [
            sample(image, x + dx / scale, y + dy / scale) * 1e6
            for dx in offsets
            for dy in offsets
            if 10 <= np.hypot(dx, dy) <= 70
        ]
    )
    values = values[np.isfinite(values)]
    if len(values) < 40:
        raise ValueError("Insufficient local radio control measurements")
    median = float(np.median(values))
    scatter = float(1.4826 * np.median(np.abs(values - median)))
    neighbours, ncat, catalogue_receipt = nearest_catalogue(cache, sky)
    plot(image, rms_display, center, scale, neighbours, sky, beam)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            position_sha256=text_digest(SOURCE),
            contributing_date_range_utc=OBSERVATION_DATES,
            observing_log="https://arxiv.org/pdf/2411.04958#page=22",
            noise_ujy_per_beam=noise,
            image_value_ujy_per_beam=peak,
            image_over_released_noise=peak / noise,
            max_native_grid_discrepancy_arcsec=max_discrepancy,
            beam_fwhm_arcsec=beam,
            effective_frequency_ghz=effective_ghz,
            control_median_ujy_per_beam=median,
            control_scatter_ujy_per_beam=scatter,
            control_samples=len(values),
            five_scatter_benchmark_ujy_per_beam=5 * max(scatter, noise),
            catalogue_rows=ncat,
            search_radius_arcsec=CATALOGUE_RADIUS_ARCSEC,
            catalogue_neighbours=neighbours,
            catalogue_archive=catalogue_receipt,
            image_products={
                role: dict(etag=obj.etag, object_bytes=obj.size, range_receipts=obj.ranges)
                for role, obj in products.items()
            },
            scope="Fixed-position Stokes I surface brightness from a single coadd. Spatial control "
            "measurements supplement the released RMS map. The beam can blend multiple sources, "
            "and the coadd gives no single-epoch sensitivity or spectral classification.",
        ),
    )


if __name__ == "__main__":
    main()

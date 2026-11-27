"""Measure fixed-position radio continuum across the public VLA image epochs."""

import json

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from download_jades_103607_vla import CATALOGUE, IMAGE
from download_jades_103607_vla import TAG as INVENTORY
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
from measure_jades_103607_alma import CONTEXT, beam_model, mad, sample
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_vla_radio_history"
MIN_CONTROLS = 50


def verified(item):
    path = ROOT / item["path"]
    if digest(path) != item["sha256"]:
        raise ValueError(f"VLA checksum mismatch: {path.name}")
    return path


def header_beam(header, is_nvss):
    if is_nvss:
        if header["SURVEY"] != "NVSS" or header["TELESCOP"] != "VLA":
            raise ValueError("Unexpected survey in SkyView radio cutout")
        return (45.0, 45.0, 0.0), 1.4
    if header["CTYPE3"] != "FREQ" or header["CTYPE4"] != "STOKES" or header["CRVAL4"] != 1:
        raise ValueError("Expected a Stokes I continuum map")
    beam = (header["BMAJ"] * 3600, header["BMIN"] * 3600, np.deg2rad(float(header["BPA"])))
    frequency = (header["CRVAL3"] * u.Unit(header["CUNIT3"])).to_value(u.GHz)
    return beam, frequency


def local_plane(item, sky, half_arcsec):
    with fits.open(verified(item), memmap=True) as hdus:
        hdu = hdus[0]
        header = hdu.header.copy()
        wcs = WCS(header, fix=False).celestial
        x, y = map(float, wcs.world_to_pixel(sky))
        radius = int(np.ceil(half_arcsec / np.min(np.linalg.svd(wcs.pixel_scale_matrix, compute_uv=False)) / 3600))
        x0, y0 = max(0, round(x)-radius), max(0, round(y)-radius)
        region = (slice(y0, round(y)+radius+1), slice(x0, round(x)+radius+1))
        if any(size != 1 for size in hdu.shape[:-2]):
            raise ValueError("Expected a single VLA sky plane")
        axes = (0,) * (len(hdu.shape)-2)
        image = hdu.data[axes + region].astype(float)
        if image.size == 0:
            raise ValueError("VLA cutout has no pixels at this position")
        if header["BUNIT"].lower() != "jy/beam":
            raise ValueError("Unexpected continuum map flux-density units")
    return image * 1e6, header, wcs.slice(region)


def catalogue_context(product, mosaic, sky):
    table = Table.read(verified(product), format="fits")
    coord = SkyCoord(table["RA"], table["DEC"])
    # The publisher's Table 1 specifies mJy/beam for these columns. The FITS
    # TUNIT says Jy/beam, so check the numeric scale against the released mosaic.
    with fits.open(verified(mosaic), memmap=True) as hdus:
        hdu = hdus[0]
        wcs = WCS(hdu.header, fix=False).celestial
        row, col = wcs.world_to_array_index(coord)
        peaks = np.asarray(table["Peak_flux"], dtype=float)
        in_bounds = (col >= 0) & (col < hdu.shape[-1]) & (row >= 0) & (row < hdu.shape[-2]) & (peaks >= 0.5)
        pixels = hdu.data[0, 0, row[in_bounds], col[in_bounds]]
    ratios = pixels / peaks[in_bounds]
    usable = np.isfinite(ratios)
    if np.count_nonzero(usable) < 50 or not 0.0008 < np.median(ratios[usable]) < 0.0012:
        raise ValueError("Published W-CDFS catalogue scale disagrees with its released Stokes I image")
    units = dict(
        catalogue_fits_peak_unit=str(table["Peak_flux"].unit), catalogue_fits_rms_unit=str(table["Isl_rms"].unit),
        publication_peak_and_rms_unit="mJy/beam", publication_table="https://www.mdpi.com/2218-1997/12/2/38",
        finite_matched_image_samples=int(usable.sum()),
        median_image_jy_per_beam_to_catalogue_numeric_peak=float(np.median(ratios[usable])),
    )
    distance = coord.separation(sky).arcsec
    indices = np.argsort(distance)
    closest = []
    for i in indices[:8]:
        closest.append(dict(
            id=str(table["Cmpt"][i]), position_icrs_deg=[float(coord[i].ra.deg), float(coord[i].dec.deg)],
            separation_arcsec=float(distance[i]),
            peak_ujy_beam=float(table["Peak_flux"][i]) * 1000,
            island_rms_ujy_beam=float(table["Isl_rms"][i]) * 1000,
        ))
    return closest, coord[distance <= 1100], len(table), units


def measure(item, sky, neighbor, cat, name, half_arcsec, plot_radius, rms_file=None, nvss=False):
    image, header, wcs = local_plane(item, sky, half_arcsec)
    beam, frequency = header_beam(header, nvss)
    matrix = wcs.pixel_scale_matrix * 3600
    cx, cy = map(float, wcs.world_to_pixel(sky))
    nx, ny = map(float, wcs.world_to_pixel(neighbor))
    target = beam_model(image, (cx, cy), matrix, beam)
    nearest = beam_model(image, (nx, ny), matrix, beam)
    if not np.isfinite(target):
        raise ValueError(f"Missing target beam support in {name}")
    grid = np.linspace(-0.81 * half_arcsec, 0.81 * half_arcsec, 19)
    offsets = np.asarray(
        [(x, y) for x in grid for y in grid if 2 * beam[0] < np.hypot(x, y) <= 0.81 * half_arcsec]
    )
    if len(offsets) < MIN_CONTROLS or grid[1] - grid[0] < beam[0]:
        raise ValueError(f"Control grid fails separation or size checks in {name}")
    places = sky.spherical_offsets_by(offsets[:, 0] * u.arcsec, offsets[:, 1] * u.arcsec)
    separated = np.ones(len(places), dtype=bool)
    if len(cat):
        _, sep, _ = places.match_to_catalog_sky(cat)
        separated &= sep.arcsec >= 2 * beam[0]
    x, y = wcs.world_to_pixel(places)
    fitted = np.array([beam_model(image, xy, matrix, beam) for xy in zip(x[separated], y[separated], strict=True)])
    finite = np.isfinite(fitted)
    values = fitted[finite]
    if len(values) < MIN_CONTROLS:
        raise ValueError(f"Insufficient local blank control apertures in {name}")
    baseline = float(np.median(values))
    scatter = mad(values)
    inner = np.hypot(*offsets[separated][finite].T) <= 0.46 * half_arcsec
    inner_scatter = mad(values[inner]) if np.count_nonzero(inner) >= 25 else None
    published_rms = None
    if rms_file:
        rms, _, rms_wcs = local_plane(rms_file, sky, half_arcsec)
        rx, ry = rms_wcs.world_to_pixel(SkyCoord([sky]))
        if np.isfinite(noise := sample(rms, (rx, ry))[0]) and noise > 0:
            published_rms = float(noise)
    error = max(scatter, inner_scatter or 0.0, published_rms or 0.0)
    if error <= 0:
        raise ValueError(f"Zero local radio noise in {name}")
    measured = float(target - baseline)
    score = measured / error
    summary = dict(
        name=name, target_signed_flux_ujy=measured, local_error_ujy=error, residual_score=score,
        neighbor_signed_flux_ujy=float(nearest - baseline) if np.isfinite(nearest) else None,
        neighbor_beam_fwhm_ratio=float(neighbor.separation(sky).arcsec / beam[0]),
        frequency_ghz=float(frequency), restoring_beam_arcsec=list(beam[:2]),
        restoring_beam_pa_deg=float(np.rad2deg(beam[2])), control_count=len(values),
        beam_fit_scatter_ujy=scatter, published_rms_ujy_beam=published_rms,
        inner_control_count=int(np.count_nonzero(inner)), inner_beam_fit_scatter_ujy=inner_scatter,
        five_scatter_depth_ujy=5 * error,
        recovery_90pct_5scatter_ujy=float(5 * error - np.quantile(values - baseline, 0.10, method="lower")),
        flux_increment_to_cross_5scatter_on_target_ujy=max(0, 5 * error - measured),
        fits_start_utc=header.get("DATE-OBS"), fits_end_utc=header.get("DATE-END"),
        science=item, rms=rms_file,
    )
    return summary, (image, wcs, plot_radius)


def plot(rows, sky, neighbor):
    row_count = (len(rows) + 2) // 3
    fig = Figure(figsize=(12, 3.5 * row_count + 0.35), layout="constrained")
    axes = list(fig.subplots(row_count, 3, squeeze=False).flat)
    for axis in axes[len(rows):]:
        axis.set_visible(False)
    n_east, n_north = sky.spherical_offsets_to(neighbor)
    for axis, (result, (image, wcs, radius)) in zip(axes[:len(rows)], rows, strict=True):
        line = np.linspace(-radius, radius, 181)
        east, north = np.meshgrid(line, line)
        grid = sky.spherical_offsets_by(east * u.arcsec, north * u.arcsec)
        density = sample(image, wcs.world_to_pixel(grid)) / result["local_error_ujy"]
        picture = axis.imshow(
            density, extent=(-radius, radius, -radius, radius), origin="lower",
            cmap="RdBu_r", vmin=-5, vmax=5,
        )
        axis.scatter(0, 0, marker="+", s=74, color="#262626")
        if radius < 100:
            axis.scatter(n_east.arcsec, n_north.arcsec, s=36, facecolors="none", edgecolor="#347354")
        major, minor = result["restoring_beam_arcsec"]
        axis.add_patch(Ellipse((radius * 0.73, -radius * 0.73), major, minor,
                               angle=90-result["restoring_beam_pa_deg"], fill=False, edgecolor="#444444"))
        datum = result["fits_start_utc"]
        day = datum[:10] if datum else "1990s mosaic"
        depth = result["five_scatter_depth_ujy"]
        axis.set(
            title=f"{result['name']} · {day}\nPoint-source sensitivity {depth:,.0f} µJy",
            xlim=(radius, -radius), xlabel="East offset (arcsec)",
        )
        fig.colorbar(picture, ax=axis, label="Local noise units", shrink=0.79)
    fig.supylabel("North offset (arcsec)")
    fig.suptitle("VLA images at the JADES 103607 transient position", fontsize=16)
    save_figure(fig, FIGURES / TAG)


def main():
    path = OUT / f"{INVENTORY}.json"
    inv = json.loads(path.read_text(encoding="utf-8"))
    sky = SkyCoord(*inv["position_icrs_deg"], unit="deg")
    rows = json.loads(CONTEXT.read_text(encoding="utf-8"))["simbad"]["closest"]
    other = next(row for row in rows if row["main_id"].endswith("103483"))
    neighbor = SkyCoord(other["ra"], other["dec"], unit="deg")
    wide_products = inv["w_cdfs"]["products"]
    closest, coords, count, units = catalogue_context(wide_products[CATALOGUE], wide_products[IMAGE], sky)
    cells = [
        measure(inv["nvss"]["science"], sky, neighbor, coords, "NVSS", 730, 270, nvss=True),
        measure(wide_products[IMAGE], sky, neighbor, coords, "W-CDFS", 220, 80),
    ]
    for product in inv["vlass"]["products"]:
        campaign = product["observation_id"].split(".")[:2]
        name = ".".join(campaign)
        item = measure(product["science"], sky, neighbor, coords, name, 106, 31, rms_file=product.get("rms"))
        timestamp = Time(item[0]["fits_start_utc"])
        reference = product["archive_mjd_range"][0]
        if reference is not None and abs((timestamp - Time(float(reference), format="mjd")).sec) > 2 * 3600:
            raise ValueError(f"VLASS archive and FITS times disagree: {product['observation_id']}")
        cells.append(item)
    cells.sort(key=lambda pair: pair[0]["fits_start_utc"] or "")
    plot(cells, sky, neighbor)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__), inventory_sha256=text_digest(path),
        context_sha256=text_digest(CONTEXT), position_icrs_deg=inv["position_icrs_deg"],
        neighbor_icrs_deg=[float(neighbor.ra.deg), float(neighbor.dec.deg)],
        neighbor_separation_arcsec=float(sky.separation(neighbor).arcsec),
        w_cdfs_catalogue_rows=count, nearest_w_cdfs_catalogue_components=closest, w_cdfs_catalogue_units=units,
        measurements=[result for result, _ in cells], event_archive=inv["event_vla_pointing_search"],
        methodology="Fixed-position unit-peak restoring-beam fits with a local planar background "
        "on science images in Jy per beam. Control positions exclude nearby W-CDFS catalogue "
        "components and are separated by at least one restoring-beam major FWHM. Use the larger "
        "of the local and inner scaled MAD of controls and the supplied VLASS RMS map value when available. "
        "The 90 percent recovery amplitude adds a matching compact-source contribution to the "
        "measured control residuals. NVSS uses the survey beam; the map is reprojected by SkyView.",
        limitation="The W-CDFS and NVSS mosaics combine exposures from their survey intervals. "
        "The VLASS Quick Look epoch 1.1 has documented pointing and flux-scale issues. "
        "The broad NVSS and W-CDFS beams blend the transient coordinate and nearby galaxy. "
        "An independent counterpart during the 2026 JWST visits requires contemporaneous target data.",
    ))


if __name__ == "__main__":
    main()

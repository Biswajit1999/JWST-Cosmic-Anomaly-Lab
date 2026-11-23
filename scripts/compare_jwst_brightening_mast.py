"""Cross-check the GOODS-S brightening with the official STScI stage-3 reduction."""

import json

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_area
from astroquery.mast import Observations
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from measure_jwst_brightening_offset import (
    BANDS,
    fit_core,
    local_alignment,
    match_pair,
    patch,
    sky_offsets,
)
from measure_jwst_epoch_candidate import RADII, photometry, pixel_scale, show_stamp
from measure_jwst_module_b_brightening import ANNULUS, FILTERS
from read_jwst_i2d_cutout import MastImageSections
from scipy.ndimage import binary_erosion, map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

PROGRAMS = ("08060", "12577")
STORAGE_MIB = 19_500
REGION_ARCSEC = 16.5


def official_associations():
    Observations.TIMEOUT = 45
    result = {}
    for program in PROGRAMS:
        selection = dict(obs_collection="JWST", proposal_id=str(int(program)))
        table = Observations.query_criteria(**selection)
        for filt in FILTERS:
            obs_id = f"jw{program}-o001_t001_nircam_clear-{filt.lower()}"
            (r,) = [row for row in table if str(row["obs_id"]) == obs_id and int(row["calib_level"]) == 3]
            product = str(r["dataURL"])
            filename = obs_id+"_i2d.fits"
            if r["dataRights"] != "PUBLIC" or product != "mast:JWST/product/"+filename:
                raise ValueError(f"Missing public stage-3 product for {obs_id}")
            result[program, filt] = dict(filename=filename, collection_query=selection, obs_id=obs_id,
                                         obsid=str(r["obsid"]), data_uri=product, archive_start_mjd=float(r["t_min"]))
    return result


def common_grid(template, newer):
    older_image, older_wcs, _, _, ox, oy = template
    new_image, new_wcs, _, _, nx, ny = newer
    yy, xx = np.indices(older_image.shape, dtype=float)
    ra, dec = older_wcs.pixel_to_world_values(xx+ox, yy+oy)
    new_x, new_y = new_wcs.world_to_pixel_values(ra, dec)
    coordinates = [new_y-ny, new_x-nx]
    mask = np.isfinite(new_image)
    safe = binary_erosion(mask, iterations=2)
    valid = map_coordinates(safe.astype(float), coordinates, order=1, mode="constant") >= .999
    mapped = map_coordinates(np.where(mask, new_image, 0.), coordinates, order=3, mode="constant")
    mapped[~valid] = np.nan
    scale = float(proj_plane_pixel_area(older_wcs)*np.deg2rad(1.)**2*1e15)
    return older_image*scale, mapped*scale


def measure_band(old, new, band):
    _, wcs, x, y, x0, y0 = old
    scale = pixel_scale(wcs)
    references_radius = BANDS.get(band, .22)
    old_image, new_image = common_grid(old, new)
    registration = local_alignment(old_image, new_image, x, y, scale, references_radius)
    template, registered, difference = match_pair(old_image, new_image, registration, scale, 0.)
    raw = new_image-old_image
    measurements = {label: photometry(data, x, y, scale, RADII, ANNULUS)
                    for label, data in (("older", old_image), ("newer_on_old_grid", new_image),
                                        ("unmatched_difference", raw), ("matched_difference", difference))}
    astrometry = None
    if band in BANDS:
        radius = BANDS[band]
        host, excess = (fit_core(*patch(im, x, y, scale, radius), radius) for im in (template, difference))
        if host is None or excess is None:
            raise ValueError(f"Unconstrained host or excess core in {band}")
        delta = sky_offsets(wcs, (x0, y0), x, y, scale, host, excess)
        astrometry = dict(host_core=host, excess_core=excess, excess_minus_host_east_north_arcsec=delta,
                          projected_offset_arcsec=float(np.linalg.norm(delta)))
    record = dict(registration=registration, aperture_photometry=measurements, astrometry=astrometry)
    return record, (template, registered, difference, x, y, scale, astrometry)


def plot(data, views, dawn):
    fig = Figure(figsize=(10.6, 9.1), layout="constrained")
    grid = fig.add_gridspec(4, 3, height_ratios=(1.4, 1, 1, 1))
    ax = fig.add_subplot(grid[0, :])
    for dx, colour, label, rows in (
        (-.075, "#245d83", "STScI stage 3", [data[b]["aperture_photometry"]["matched_difference"][-1] for b in FILTERS]),
        (.075, "#ac5437", "DAWN", [dawn[b]["unit_scale_apertures"][-1] for b in FILTERS]),
    ):
        ax.errorbar(np.arange(len(FILTERS))+dx, [v["flux_njy"] for v in rows],
                    yerr=[v["local_mad_sigma_njy"] for v in rows], color=colour, label=label,
                    marker="o", linestyle="none", capsize=3)
    ax.axhline(0, lw=.7, color=".5")
    ax.grid(axis="y", alpha=.2)
    ax.legend(loc="upper right", frameon=False)
    ax.set(ylabel="Change in 0.24″ aperture (nJy)", xticks=range(len(FILTERS)), xticklabels=FILTERS)
    for row, band in enumerate(BANDS, 1):
        old, new, residual, x, y, scale, pos = views[band]
        display = .51
        pair = np.concatenate([patch(d, x, y, scale, display)[2].ravel() for d in (old, new)])
        bottom, top = np.nanpercentile(pair, [4, 98.5])
        limit = np.nanpercentile(abs(patch(residual, x, y, scale, display)[2]), 98)
        for col, (im, name, cmap, bounds) in enumerate((
            (old, "Oct 2025", "gray", (bottom, top)), (new, "Oct 2026", "gray", (bottom, top)),
            (residual, "2026 minus 2025", "RdBu_r", (-limit, limit)),
        )):
            ax = fig.add_subplot(grid[row, col])
            show_stamp(ax, im, x, y, scale, bounds, cmap, .24, half_arcsec=display)
            ax.set_title(f"{band} · {name}", fontsize=9)
            if col == 2:
                for key, mark, colour in (("host_core", "+", "#158ba8"), ("excess_core", "x", "#b31972")):
                    ax.scatter(*pos[key]["offset_xy_arcsec"], marker=mark, s=44, color=colour, linewidths=1.2)
    fig.suptitle("JADES 103607 · Two reductions of the same JWST observations\n"
                 "STScI stage-3 images · cyan + quiescent core · magenta × excess", fontsize=12)
    save_figure(fig, FIGURES/"module_b_001_mast_comparison")


def main():
    input_path = OUT/"module_b_001_six_filter_change.json"
    evidence = json.loads(input_path.read_text())
    sky = SkyCoord(*evidence["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    associations = official_associations()
    data, products, views = {}, [], {}
    for band in FILTERS:
        images = []
        for program in PROGRAMS:
            record = associations[program, band]
            image, receipt = MastImageSections(cache, record["filename"]).cutout(sky, REGION_ARCSEC)
            if receipt["filter"] != band or abs(receipt["exposure_start_mjd"]-record["archive_start_mjd"]) > 1/1440:
                raise ValueError(f"Archive and FITS metadata disagree for {record['filename']}")
            products.append(dict(**receipt, mast_association=record))
            images.append(image)
            print(program, band, receipt["pipeline_version"], receipt["date_utc"], flush=True)
        data[band], views[band] = measure_band(*images, band)
        print(band, data[band]["aperture_photometry"]["matched_difference"][-1], flush=True)
    plot(data, views, evidence["difference_photometry"])
    write_json(OUT/"module_b_001_mast_comparison.json", dict(
        executed_utc=utcnow(), candidate=evidence["candidate"], icrs_deg=evidence["icrs_deg"],
        storage_cap_mib=STORAGE_MIB, half_cutout_arcsec=REGION_ARCSEC,
        official_mosaics=products, bands=data,
        dawn_photometry_sha256=text_digest(input_path),
        cutout_reader_sha256=text_digest(ROOT/"scripts/read_jwst_i2d_cutout.py"),
        relative_astrometry_script_sha256=text_digest(ROOT/"scripts/measure_jwst_brightening_offset.py"),
        script_sha256=text_digest(__file__),
        scope="STScI and DAWN perform separate resampling and image calibration of shared exposures."
              " The range receipts hash the downloaded sections, and ETags identify the upstream FITS versions.",
    ))


if __name__ == "__main__":
    main()

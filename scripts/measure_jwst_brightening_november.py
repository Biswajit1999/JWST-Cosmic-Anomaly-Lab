"""Measure the new November JWST epoch of the JADES 103607 brightening."""

import json

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_area
from astroquery.mast import Observations
from compare_jwst_brightening_mast import REGION_ARCSEC, common_grid
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from measure_jwst_brightening_offset import BANDS, fit_core, local_alignment, match_pair, patch
from measure_jwst_epoch_candidate import RADII, photometry, pixel_scale, show_stamp
from measure_jwst_module_b_brightening import ANNULUS
from read_jwst_i2d_cutout import MastImageSections
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

FILTERS = ("F150W", "F200W", "F277W", "F356W")
EPOCHS = ("Oct 2025", "Oct 2026", "Nov 2026")
ASSOCIATIONS = ("jw08060-o001_t001", "jw12577-o001_t001", "jw12577-o004_t004")
PREFIX = "module_b_001_november_photometry"
STORAGE_MIB = 25_500


def public_mosaics():
    Observations.TIMEOUT = 55
    query = dict(obs_collection="JWST", proposal_id="12577")
    rows = Observations.query_criteria(**query)
    catalogue = {}
    for row in rows:
        if not str(row["obs_id"]).startswith(ASSOCIATIONS[-1]) or int(row["calib_level"]) != 3:
            continue
        name = str(row["obs_id"])+"_i2d.fits"
        uri = str(row["dataURL"])
        if row["dataRights"] != "PUBLIC" or uri != "mast:JWST/product/"+name:
            raise ValueError(f"Expected public stage-3 data for {name}")
        catalogue[name] = dict(data_uri=uri, obsid=str(row["obsid"]), start_mjd=float(row["t_min"]), query=query)
    return catalogue


def aperture_measurement(stamp):
    image, wcs, x, y, _, _ = stamp
    n_jy_per_pixel = float(proj_plane_pixel_area(wcs)*np.deg2rad(1.)**2*1e15)
    return photometry(image*n_jy_per_pixel, x, y, pixel_scale(wcs), RADII, ANNULUS)


def difference(old, new, band, sign):
    _, wcs, x, y, x0, y0 = old
    scale = pixel_scale(wcs)
    a, b = common_grid(old, new)
    alignment = local_alignment(a, b, x, y, scale, BANDS.get(band, .22))
    template, registered, residual = match_pair(a, b, alignment, scale, 0.)
    measurements = photometry(residual, x, y, scale, RADII, ANNULUS)
    if measurements[-1]["flux_njy"] is None:
        raise ValueError(f"No usable matched target photometry in {band}")
    core = None
    if band in BANDS:
        radius = BANDS[band]
        core = fit_core(*patch(sign*residual, x, y, scale, radius), radius)
        if core is not None:
            dx, dy = np.array(core["offset_xy_arcsec"])/scale
            core["icrs_deg"] = list(map(float, wcs.pixel_to_world_values(x+x0+dx, y+y0+dy)))
    measured = dict(alignment=alignment, apertures=measurements, core_sign=sign, core=core)
    return measured, (template, registered, residual, x, y, scale)


def plot(data, views):
    fig = Figure(figsize=(10, 8.7), layout="constrained")
    grid = fig.add_gridspec(4, 3, height_ratios=(1.55, 1, 1, 1))
    ax = fig.add_subplot(grid[0, :])
    for j, (band, colour) in enumerate(zip(FILTERS, ("#006990", "#bc4c22", "#71519a", "#42744a"), strict=True)):
        points = [data[band]["matched_differences"][label]["apertures"][-1]
                  for label in ("october_2026_minus_2025", "november_2026_minus_2025")]
        ax.errorbar(np.arange(2)+(j-1.5)*.055, [p["flux_njy"] for p in points],
                    yerr=[p["local_mad_sigma_njy"] for p in points], fmt="o", capsize=3,
                    label=band, color=colour)
    ax.set(xticks=[0, 1], xticklabels=["Oct 2026", "Nov 2026"],
           ylabel="Excess over Oct 2025 (nJy)\n0.24″ radius")
    ax.axhline(0., color=".65", lw=.7)
    ax.legend(loc="best", ncols=2, frameon=False)
    for j, band in enumerate(("F150W", "F200W", "F277W"), start=1):
        oct_image, nov_image, residual, x, y, scale = views[band]
        half = .52
        pixels = np.concatenate([patch(im, x, y, scale, half)[2].ravel() for im in (oct_image, nov_image)])
        low, high = np.nanpercentile(pixels, (4, 98.5))
        limit = np.nanpercentile(np.abs(patch(residual, x, y, scale, half)[2]), 98)
        for i, (im, name, bounds, cmap) in enumerate((
            (oct_image, "Oct 2026", (low, high), "gray"),
            (nov_image, "Nov 2026", (low, high), "gray"),
            (residual, "Nov − Oct 2026", (-limit, limit), "RdBu_r"),
        )):
            ax = fig.add_subplot(grid[j, i])
            show_stamp(ax, im, x, y, scale, bounds, cmap, .24, half_arcsec=half)
            ax.set_title(f"{band} · {name}", fontsize=9)
    fig.suptitle("JADES 103607 · November 2026 JWST follow-up\n"
                 "STScI NIRCam mosaics · registered and PSF-matched", fontsize=12)
    save_figure(fig, FIGURES/PREFIX)


def main():
    prior = OUT/"module_b_001_mast_comparison.json"
    target = json.loads(prior.read_text())
    sky = SkyCoord(*target["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    new_catalogue = public_mosaics()
    verified_earlier = {row["filename"]: row["mast_association"] for row in target["official_mosaics"]}
    records, data, views = [], {}, {}
    for band in FILTERS:
        stamps = []
        for i, (epoch, assoc) in enumerate(zip(EPOCHS, ASSOCIATIONS, strict=True)):
            name = f"{assoc}_nircam_clear-{band.lower()}_i2d.fits"
            association = new_catalogue[name] if i == 2 else verified_earlier[name]
            stamp, receipt = MastImageSections(cache, name).cutout(sky, REGION_ARCSEC)
            start = association.get("archive_start_mjd", association.get("start_mjd"))
            if receipt["filter"] != band or abs(receipt["exposure_start_mjd"]-start) > 1/1440:
                raise ValueError(f"FITS and MAST metadata disagree for {name}")
            records.append(dict(**receipt, epoch=epoch, mast_association=association))
            stamps.append(stamp)
            print(epoch, band, receipt["date_utc"], flush=True)
        cases = (
            ("october_2026_minus_2025", 0, 1, 1),
            ("november_2026_minus_2025", 0, 2, 1),
            ("november_minus_october_2026", 1, 2, -1),
        )
        results = {}
        for label, a, b, sign in cases:
            results[label], scene = difference(stamps[a], stamps[b], band, sign)
            print(band, label, results[label]["apertures"][-1], flush=True)
        rising, fading = (results[cases[i][0]]["core"] for i in (0, 2))
        separation = None
        if rising is not None and fading is not None:
            first, second = (SkyCoord(*core["icrs_deg"], unit="deg") for core in (rising, fading))
            separation = float(first.separation(second).arcsec)
            print(band, "rise/fade core separation in arcsec", separation, flush=True)
        data[band] = dict(rise_fade_core_separation_arcsec=separation, direct_apertures={e: aperture_measurement(s) for e, s in zip(EPOCHS, stamps, strict=True)},
                          matched_differences=results)
        views[band] = scene
    plot(data, views)
    write_json(OUT/f"{PREFIX}.json", dict(
        executed_utc=utcnow(), candidate=target["candidate"], icrs_deg=target["icrs_deg"],
        official_mosaics=records, bands=data, previous_comparison_sha256=text_digest(prior),
        script_sha256=text_digest(__file__), storage_cap_mib=STORAGE_MIB,
        cutout_reader_sha256=text_digest(ROOT/"scripts/read_jwst_i2d_cutout.py"),
        pair_registration_sha256=text_digest(ROOT/"scripts/measure_jwst_brightening_offset.py"),
        scope="Four released filters of public GO 12577 observation 004."
              " Identical fixed sky position across epochs with local registration and PSF matching."
              " Error estimates are local aperture scatter, with shared data across comparisons.",
    ))


if __name__ == "__main__":
    main()

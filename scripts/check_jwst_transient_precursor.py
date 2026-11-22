"""Search local exposure footprints and measure two MIRI epochs before the flare."""

import ast
import json
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from audit_jwst_transient_wfss import MAST
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_transient_grism_geometry import STORAGE_MIB, fetch
from matplotlib.figure import Figure
from matplotlib.path import Path
from measure_jwst_epoch_candidate import aperture, pixel_scale, show_stamp, verified_path
from scipy.ndimage import map_coordinates
from search_jwst_cross_detector import FIGURES, check_same_grid, window

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

FILTERS = ("F1000W", "F2100W")
VISITS = ("jw01207001001", "jw05407040001")
FWHM = (.356, .685)
RADII = ((.4, .55, .75), (.65, .85, 1.05))
ANNULUS = ((.85, 1.25), (1.15, 1.65))
EARLIEST_BRIGHT = "2025-09-30"
LAST_CHECKED_FAINT = "2024-01-01"


def query(cache, endpoint, params, name):
    url = API+"/"+endpoint+"?"+urlencode(params)
    receipt = download(cache, url, "transient_precursor_"+name+".csv", "local_coverage_listing")
    return pd.read_csv(verified_path(receipt)), receipt


def contains(footprint, sky):
    """Test corners in a local tangent plane with east and north axes."""
    points = np.asarray(ast.literal_eval(footprint), dtype=float)
    if points.shape != (4, 2):
        raise ValueError(f"Expected rectangular footprint, got shape {points.shape}")
    e, n = sky.spherical_offsets_to(SkyCoord(points[:, 0], points[:, 1], unit="deg"))
    return bool(Path(np.column_stack([e.arcsec, n.arcsec])).contains_point((0, 0)))


def coverage(cache, sky):
    ra, dec = sky.ra.deg, sky.dec.deg
    exposures, exp_receipt = query(cache, "exposures", dict(ra=ra, dec=dec, arcmin=.008, output="csv"), "exposures")
    selected = exposures.loc[exposures.file.str.startswith("jw") & (exposures.extension == "rate")].copy()
    selected = selected[selected.footprint.map(lambda poly: contains(poly, sky))]
    selected["date_utc"] = Time(selected.expstart.to_numpy(), format="mjd").isot
    selected["visit"] = selected.file.str.slice(0, 13)
    rows = []
    for (visit, band, detector), group in selected.sort_values("expstart").groupby(
            ["visit", "filter", "detector"], sort=False):
        rows.append(dict(visit=visit, band=band, detector=detector,
                         first_exposure_utc=group.date_utc.iloc[0], last_exposure_utc=group.date_utc.iloc[-1],
                         native_exposure_count=len(group), total_exposure_seconds=float(group.exptime.sum())))

    shutters, msa_receipt = query(cache, "nirspec_slits", dict(coords=f"{ra},{dec}", size=15), "msa")
    shutter_centers = SkyCoord(shutters.ra, shutters.dec, unit="deg")
    shutters["distance"] = sky.separation(shutter_centers).arcsec
    shutters["contains_target"] = shutters.footprint.map(lambda poly: contains(poly, sky))
    science = shutters[shutters.valid]
    active = shutters[shutters.root.str.startswith("jw08060")]
    sources = active.loc[active.is_source & (active.source_id > 0)].drop_duplicates("source_id")
    source_sep = sky.separation(SkyCoord(sources.srcra, sources.srcdec, unit="deg")).arcsec
    m = int(np.argmin(source_sep))
    nearest = science.loc[science.distance.idxmin()]
    msa = dict(query_radius_arcsec=15, rows=len(shutters), valid_rows=len(science),
               exact_footprint_rows=int(shutters.contains_target.sum()),
               valid_exact_footprint_rows=int(science.contains_target.sum()),
               nearest_valid_shutter=dict(file=nearest.root, separation_arcsec=float(nearest.distance),
                                         date_utc=Time(nearest.expstart, format="mjd").isot),
               go8060_rows=len(active), go8060_nearest_shutter_arcsec=float(active.distance.min()),
               go8060_nearest_listed_target=dict(source_id=int(sources.iloc[m].source_id),
                                                 separation_arcsec=float(source_sep[m])))
    interval = [r for r in rows if LAST_CHECKED_FAINT < r["first_exposure_utc"][:10] < EARLIEST_BRIGHT]
    return dict(exposures=exp_receipt, nirspec_shutters=msa_receipt, matching_exposure_groups=rows,
                bracketing_dates_utc=[LAST_CHECKED_FAINT, EARLIEST_BRIGHT],
                intermediate_geometric_groups=interval, msa=msa)


def miri_pair(cache, position, band):
    entries, receipt = query(cache, "assoc_mosaic", dict(
        coords=f"{position[0]:.10f},{position[1]:.10f},.02", filters=band, output="csv"), band.lower())
    products = []
    for visit in VISITS:
        (row,) = list(entries.query("visit == @visit and tile == 672 and status == 2").itertuples())
        product = dict(visit=visit, band=band, date_utc=row.dateobs, assoc_name=row.assoc_name,
                       native_exposure_count=int(row.nexp))
        for kind in ("sci", "wht"):
            url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
            product[kind] = download(cache, url, url.rsplit("/", 1)[-1], "observed_miri_epoch_"+kind)
        products.append(product)
    return products, receipt


def data_windows(records, sky):
    arrays, grids = [], []
    for row in records:
        sci, wht = (window(row, kind, sky, 20) for kind in ("sci", "wht"))
        check_same_grid(sci, wht)
        field, wcs, origin, scale = sci
        arrays.append(np.where(np.isfinite(wht[0]) & (wht[0] > 0), field, np.nan))
        grids.append(sci)
    check_same_grid(*grids)
    xy = np.array(wcs.world_to_pixel_values(*sky))-origin
    return arrays, xy, scale, wcs.slice(np.s_[origin[1]:origin[1]+field.shape[0], origin[0]:origin[0]+field.shape[1]])


def coverage_stats(valid, xy, scale):
    y, x = np.nonzero(valid)
    distance = np.hypot(x-xy[0], y-xy[1])*scale
    ix, iy = (round(float(v)) for v in xy)
    central = bool(valid[iy, ix]) if 0 <= iy < valid.shape[0] and 0 <= ix < valid.shape[1] else False
    return dict(target_pixel_usable=central, usable_pixels_within_1_arcsec=int((distance <= 1).sum()),
                nearest_usable_pixel_arcsec=float(distance.min()) if len(distance) else None)


def official_coverage(cache, band, sky):
    filename = f"jw05407-o040_t085_miri_{band.lower()}_i2d.fits"
    url = MAST+urlencode(dict(uri="mast:JWST/product/"+filename))
    receipt = fetch(cache, url, filename, "observed_miri_coverage_audit",
                    version="GO 5407 public calibrated direct mosaic")
    with fits.open(verified_path(receipt)) as h:
        sci, weight = h["SCI"], h["WHT"]
        if h[0].header["FILTER"] != band or sci.header["BUNIT"] != "MJy/sr":
            raise ValueError(f"Unexpected filter or science calibration in {filename}")
        wcs = WCS(sci.header).celestial
        xy = wcs.world_to_pixel_values(*sky)
        valid = (weight.data > 0) & np.isfinite(sci.data)
        metrics = coverage_stats(valid, xy, pixel_scale(wcs))
    return dict(product=receipt, **metrics), valid, wcs


def aperture_series(images, xy, scale, radii, annulus):
    oldest, newest = images
    difference = newest-oldest
    step = max(1, round(3.4/scale))
    grid = np.arange(-int(15/scale), int(15/scale)+1, step)
    offsets = [(dx, dy) for dy in grid for dx in grid if 3 <= np.hypot(dx, dy)*scale <= 15]
    results = []
    for radius in radii:
        def measure(image, x, y):
            return aperture(image, x, y, scale, radius, annulus)

        old_flux, new_flux, diff_flux = (measure(image, *xy)/1000 for image in (*images, difference))
        row = dict(radius_arcsec=radius)
        for label, image, flux in zip(("2022", "2024", "difference"), (*images, difference),
                                      (old_flux, new_flux, diff_flux), strict=True):
            blank = np.array([measure(image, *(xy+[dx, dy]))/1000 for dx, dy in offsets])
            blank = blank[np.isfinite(blank)]
            noise = 1.4826*np.median(abs(blank-np.median(blank))) if len(blank) >= 20 else None
            row[label] = dict(flux_microjy=float(flux) if np.isfinite(flux) else None,
                              local_mad_microjy=float(noise) if noise is not None else None,
                              control_apertures=len(blank))
        results.append(row)
    return results


def plot(stamps):
    fig = Figure(figsize=(10, 7.4), layout="constrained")
    for row, (band, old_image, masks, xy, scale, measured, coverage) in zip(fig.subplots(2, 3), stamps, strict=True):
        x, y = map(round, xy)
        half = round(3/scale)
        core = np.s_[y-half:y+half+1, x-half:x+half+1]
        upper = float(np.nanpercentile(old_image[core], 98.5))
        radius = measured["radius_arcsec"]
        show_stamp(row[0], old_image, *xy, scale, (-.18*upper, upper), "gray", radius, half_arcsec=3)
        row[0].set_title(f"{band} · 2022 science image", fontsize=10)
        flux, noise = [measured["2022"][k] for k in ("flux_microjy", "local_mad_microjy")]
        row[0].set_xlabel(f"{radius:.2f}″ aperture   {flux:+.2f} ± {noise:.2f} µJy", fontsize=9)
        for ax, mask, title, info in zip(row[1:], masks, ("DJA", "MAST"), coverage, strict=True):
            show_stamp(ax, mask, *xy, scale, (0, 1), "gray", 1, half_arcsec=16)
            ax.set_title(title, fontsize=10)
            distance = info["nearest_usable_pixel_arcsec"]
            ax.set_xlabel(f"Nearest usable pixel {distance:.1f}″ away", fontsize=9)
    fig.suptitle("GOODS-S transient · mid-infrared archive coverage\n"
                 "2022 MIRI crop, 6″ across · 2024 weight maps, 32″ across", fontsize=13)
    save_figure(fig, FIGURES/"field_b_transient_precursor")


def main():
    source_path = OUT/"field_b_transient_vetting.json"
    transient = json.loads(source_path.read_text())
    position = transient["position_degrees"]
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    footprints = coverage(cache, SkyCoord(*position, unit="deg"))
    report, stamps = {}, []
    for band, fwhm, radii, annulus in zip(FILTERS, FWHM, RADII, ANNULUS, strict=True):
        records, listing = miri_pair(cache, position, band)
        images, xy, scale, wcs = data_windows(records, position)
        measurements = aperture_series(images, xy, scale, radii, annulus)
        mask = np.isfinite(images[1])
        dja = coverage_stats(mask, xy, scale)
        mast, official_mask, official_wcs = official_coverage(cache, band, position)
        y, x = np.indices(mask.shape, dtype=float)
        mx, my = official_wcs.world_to_pixel_values(*wcs.pixel_to_world_values(x, y))
        aligned = map_coordinates(official_mask, [my, mx], order=0, mode="constant", cval=0)
        report[band] = dict(exposures=records, listing=listing, approximate_psf_fwhm_arcsec=fwhm,
                            pixel_scale_arcsec=scale, sky_annulus_arcsec=list(annulus), apertures=measurements,
                            dja_2024_coverage=dja, mast_2024_coverage=mast)
        print(band, measurements, flush=True)
        stamps.append((band, images[0], (mask, aligned), xy, scale, measurements[1], (dja, mast)))
    plot(stamps)
    write_json(OUT/"field_b_transient_precursor.json", dict(
        executed_utc=utcnow(), candidate=transient["candidate"], position_degrees=position,
        coverage=footprints, miri=report, source_sha256=text_digest(source_path), script_sha256=text_digest(__file__),
        scope="Signed native-grid aperture differences and empirical local scatter; no PSF deblending or aperture correction",
    ))


if __name__ == "__main__":
    main()

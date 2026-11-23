"""Measure the three-filter Module B leads in native frames and other visits."""

import json
from datetime import datetime, timezone
from io import StringIO
from urllib.parse import urlencode
from zipfile import ZipFile

import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from audit_jwst_transient_wfss import mast_product
from crossmatch_jwst_epoch_candidate import published_positions
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_transient_grism_geometry import fetch
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import (
    audit_detector_target,
    centroids,
    load_dithers,
    pixel_scale,
    show_stamp,
    summarize_detector,
    verified_path,
)
from scipy.ndimage import maximum_filter
from search_jwst_cross_detector import window
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

CACHE_MIB = 15700
BANDS = ("F200W", "F356W", "F444W")
EXTRA = {
    "001": {
        "F200W": [("jw01180010001", "NRCB4"), ("jw01180012001", "NRCB4"), ("jw05324027001", "NRCA2")],
        "F356W": [("jw01180012001", "NRCBLONG"), ("jw05324027001", "NRCALONG")],
    },
    "002": {
        "F200W": [("jw01180018001", "NRCB4"), ("jw01180016001", "NRCB4"), ("jw01180222001", "NRCB3"),
                  ("jw07336002001", "NRCB3"), ("jw08060004001", "NRCB3"), ("jw05324026001", "NRCB4"),
                  ("jw05324027001", "NRCA3")],
        "F356W": [("jw04549004001", "NRCALONG"), ("jw08060004001", "NRCBLONG"),
                  ("jw05324027001", "NRCALONG")],
    },
}


def archive_query(cache, endpoint, params, filename):
    url = API+"/"+endpoint+"?"+urlencode(params)
    receipt = download(cache, url, filename, "position_archive_query")
    return pd.read_csv(verified_path(receipt)), receipt


def extra_products(cache, candidate, center):
    products, receipts = [], {}
    for band, selections in EXTRA[candidate.pointing].items():
        filt = band+"-CLEAR"
        rows, receipt = archive_query(cache, "assoc_mosaic", dict(
            coords=f"{center[0]:.8f},{center[1]:.8f},.02", filters=filt, output="csv"),
            f"module_b_{candidate.pointing}_{band}_epochs.csv")
        receipts[band] = receipt
        for visit, detector in selections:
            (row,) = list(rows.query("visit == @visit and detector == @detector and tile == 672 and status == 2").itertuples())
            item = dict(band=band, visit=visit, detector=detector, dateobs=row.dateobs, mjd=float(row.expstart),
                        assoc_name=row.assoc_name)
            for kind in ("sci", "wht"):
                url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
                item[kind] = download(cache, url, url.rsplit("/", 1)[-1], "observed_module_b_independent_epoch")
            products.append(item)
    return products, receipts


def native_products(cache, candidate, band, product, center):
    listing, receipt = archive_query(cache, "exposures", dict(
        ra=center[0], dec=center[1], arcmin=.025, filters=band+"-CLEAR", output="csv"),
        f"module_b_{candidate.pointing}_{band}_native.csv")
    records = []
    for row in listing.loc[(listing.assoc == product["assoc_name"]) & (listing.extension == "rate")].sort_values(
            "expstart").itertuples():
        name = row.file+"_cal.fits"
        relpath = "data/raw/cross_detector/"+name
        cached = next((r for r in reversed(cache.records) if r["path"] == relpath), None)
        if cached is not None and (ROOT/relpath).exists():
            result = {k: cached[k] for k in ("url", "path", "sha256")}
            verified_path(result)
        else:
            cloud = f"https://stpubdata.s3.us-east-1.amazonaws.com/jwst/public/{name[:7]}/{name[:13]}/{name}"
            try:
                result = fetch(cache, cloud, name, "observed_module_b_native", version="MAST public calibrated exposure")
            except requests.RequestException:
                _, result = mast_product(cache, name, "observed_module_b_native")
        records.append({**result, "file": name, "expstart_mjd": float(row.expstart)})
    if len(records) != product["exposure_count"]:
        raise ValueError(f"Expected {product['exposure_count']} exposures for {product['assoc_name']}, found {len(records)}")
    return dict(listing=receipt, products=records)


def native_alignment(product, receipt, center):
    image, wcs, origin, scale = window(product, "sci", center, 34)
    weight = window(product, "wht", center, 34)[0]
    x, y = np.array(wcs.world_to_pixel_values(*center))-origin
    yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
    distance = np.hypot(xx-x, yy-y)*scale
    clean = np.where((weight > 0) & np.isfinite(image), image, -np.inf)
    isolated = clean == maximum_filter(clean, size=2*round(.24/scale)+1)
    py, px = np.nonzero(isolated & (distance > 2) & (distance < 31) & (clean > 3))
    strongest = np.argsort(clean[py, px])[-100:]
    positions = np.column_stack((px[strongest], py[strongest]))
    centers = centroids(np.where(np.isfinite(clean), clean, np.nan), positions, core=4.5)
    finite = centers[np.isfinite(centers).all(axis=1)]+origin
    references = wcs.pixel_to_world_values(*finite.T)
    frames = load_dithers(receipt, references)
    rows = audit_detector_target(frames, center)
    for row, frame in zip(rows, frames, strict=True):
        x, y = row["forced_center_pixels_xy"]
        target = centroids(frame["clean"], [(x, y)], core=3.5)[0]
        if np.isfinite(target).all():
            coord = frame["wcs"].pixel_to_world_values(*(target-frame["translation"]))
            row["clean_centroid_ra_dec"] = list(map(float, coord))
            row["clean_centroid_offset_arcsec"] = float(np.hypot(*(target-[x, y]))*frame["scale"])
        else:
            row["clean_centroid_ra_dec"] = row["clean_centroid_offset_arcsec"] = None
    return frames, rows


def tns_catalogue(cache):
    day = datetime.now(timezone.utc).strftime("%Y%m%d")
    url = "https://www.wis-tns.org/sites/default/files/tns_public_objects/tns_public_objects.csv.zip"
    original = cache.session.headers["User-Agent"]
    cache.session.headers["User-Agent"] = 'tns_marker{"tns_id":0,"type":"guest","name":"JWST Cosmic Anomaly Lab"}'
    try:
        receipt = download(cache, url, f"tns_public_objects_snapshot_{day}.csv.zip", "public_tns_daily_catalogue")
    finally:
        cache.session.headers["User-Agent"] = original
    with ZipFile(verified_path(receipt)) as archive, archive.open("tns_public_objects.csv") as stream:
        snapshot = stream.readline().decode("ascii").strip()
        cols = ["objid", "name_prefix", "name", "ra", "declination", "type", "discoverydate", "time_received"]
        table = pd.read_csv(stream, usecols=cols, encoding="utf-8", encoding_errors="replace")
    if len(table) < 150_000 or not np.isfinite(table.ra+table.declination).all():
        raise ValueError("The TNS catalogue is incomplete or has invalid coordinates")
    return {**receipt, "snapshot_label_utc": snapshot, "catalogue_rows": len(table)}, table


def tns_match(tns, sky):
    receipt, table = tns
    separation = sky.separation(SkyCoord(table.ra, table.declination, unit="deg")).arcsec
    near = table.assign(separation_arcsec=separation).nsmallest(3, "separation_arcsec")
    return dict(**receipt, radius_arcsec=5, matches=int((separation <= 5).sum()),
                nearest=json.loads(near.to_json(orient="records")))


def catalogue_checks(cache, candidate, center, tns):
    sky = SkyCoord(*center, unit="deg")
    tag = "module_b_"+candidate.pointing

    old = json.loads((OUT/"field_b_transient_vetting.json").read_text())
    published = next(v["jades_2022_2023_transients"] for v in old.values()
                     if isinstance(v, dict) and "jades_2022_2023_transients" in v)
    ra, dec = zip(*published_positions(verified_path(published).read_text(encoding="utf-8")), strict=True)
    separation = sky.separation(SkyCoord(list(ra), list(dec), unit=("hourangle", "deg"))).arcsec
    survey = {k: published[k] for k in ("url", "path", "sha256")}
    survey.update(matches=int((separation <= 5).sum()), radius_arcsec=5, searched_positions=len(separation),
                  nearest_arcsec=float(separation.min()))

    url = API+"/nirspec_extractions?"+urlencode(dict(coords=f"{center[0]:.9f},{center[1]:.9f}",
                                                  size=2, output="csv"))
    nirspec = download(cache, url, tag+"_nirspec.csv", "public_position_spectrum_search")
    body = verified_path(nirspec).read_text(encoding="utf-8")
    valid_table = body.startswith("root,file,")
    parsed = pd.read_csv(verified_path(nirspec)) if valid_table else pd.DataFrame()
    nirspec.update(rows=len(parsed) if valid_table or body.lstrip().startswith("Nothing found for") else None)
    if len(parsed):
        dist = sky.separation(SkyCoord(parsed.ra, parsed.dec, unit="deg")).arcsec
        keys = ("root", "file", "srcid", "z", "grade", "comment")
        nearby = parsed.loc[:, keys].assign(separation_arcsec=dist).sort_values("separation_arcsec").head(6)
        nirspec["nearest"] = json.loads(nearby.to_json(orient="records"))
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"+urlencode({
        "-source": "V/159/gssample", "-c": f"{center[0]:.6f} {center[1]:.6f}", "-c.rs": 3, "-c.eq": "J2000",
        "-out.max": 30, "-out": "NIRCamID,RA_ICRS,DE_ICRS,KronRad,zphot,zphotl68,zphotu68",
    })
    jades = download(cache, url, tag+"_jades_dr2_position.tsv", "public_position_photometric_catalogue")
    body = verified_path(jades).read_text(encoding="utf-8")
    jades.update(radius_arcsec=3, rows=None)
    if "NIRCamID\tRA_ICRS\tDE_ICRS" in body:
        table = pd.read_csv(StringIO(body), comment="#", sep="\t")
        table = table.apply(pd.to_numeric, errors="coerce").dropna(subset=["NIRCamID", "RA_ICRS", "DE_ICRS"])
        dist = sky.separation(SkyCoord(table.RA_ICRS, table.DE_ICRS, unit="deg")).arcsec
        nearest = table.assign(separation_arcsec=dist).sort_values("separation_arcsec").head(5)
        jades.update(rows=len(table), nearest=json.loads(nearest.to_json(orient="records")))
    return dict(tns=tns_match(tns, sky), earlier_jades_transients=survey, nirspec=nirspec, jades_dr2=jades)


def plot_candidate(candidate, data, bright_dithers, samples):
    tag = f"module_b_{candidate.pointing}"
    fig = Figure(figsize=(9.5, 8), layout="constrained")
    axes = fig.subplots(3, 3)
    for row, band in zip(axes, BANDS, strict=True):
        old, new = (data[epoch, band][0] for epoch in ("older", "newer"))
        o, ow, ox, oy, ox0, oy0 = old
        n, nw, nx, ny, nx0, ny0 = new
        scale = pixel_scale(ow)
        if o.shape != n.shape or np.hypot(ox-nx, oy-ny) > .03 or abs(pixel_scale(nw)-scale) > 1e-7:
            raise ValueError("Expected common native mosaic grid for a direct difference")
        test_ra_dec = ow.pixel_to_world_values(ox0+ox+100, oy0+oy)
        px, py = nw.world_to_pixel_values(*test_ra_dec)
        if np.hypot(px-(nx0+nx+100), py-(ny0+ny)) > .03:
            raise ValueError("Mosaic displacement prohibits pixel-wise difference")
        difference = n-o
        ix, iy, radius = round(ox), round(oy), round(.72/scale)
        cut = np.s_[iy-radius:iy+radius+1, ix-radius:ix+radius+1]
        floor, top = np.nanpercentile(np.concatenate([o[cut].ravel(), n[cut].ravel()]), [4, 98.5])
        limit = np.nanpercentile(abs(difference[cut]), 97.5)
        for ax, im, name, cmap, bounds in zip(row, (o, n, difference),
                ("2025", "2026", "2026 minus 2025"), ("gray", "gray", "RdBu_r"),
                ((floor, top), (floor, top), (-limit, limit)), strict=True):
            show_stamp(ax, im, ox, oy, scale, bounds, cmap, .24, half_arcsec=.78)
            ax.set_title(f"{band} · {name}", fontsize=9)
    fig.suptitle(f"GOODS-S · {candidate.ra:.6f}°, {candidate.dec:+.6f}°\n"
                 "Observed JWST images · common intensity scale per filter", fontsize=12)
    save_figure(fig, FIGURES/(tag+"_mosaics"))

    fig = Figure(figsize=(9.5, 6.5), layout="constrained")
    axes = fig.subplots(2, 3)
    for row, (band, (frames, records)) in zip(axes, bright_dithers.items(), strict=True):
        for i, (ax, frame, record) in enumerate(zip(row, frames, records, strict=True), 1):
            x, y = record["forced_center_pixels_xy"]
            image, scale = frame["image"], frame["scale"]
            h = round(.50/scale)
            ix, iy = round(x), round(y)
            bounds = tuple(np.nanpercentile(image[iy-h:iy+h+1, ix-h:ix+h+1], [4, 98.5]))
            show_stamp(ax, image, x, y, scale, bounds, "gray", .16, half_arcsec=.50)
            if record["unmasked_peak_pixels_xy"]:
                u, v = record["unmasked_peak_pixels_xy"]
                ax.scatter((u-x)*scale, (v-y)*scale, marker="+", color="#ff6990", s=55)
            ax.set_title(f"{band} · dither {i} · DQ {record['unmasked_peak_dq']}", fontsize=9)
    fig.suptitle(f"GOODS-S · {candidate.ra:.6f}°, {candidate.dec:+.6f}°\n"
                 "Bright epoch native exposures · cyan at sky position · pink at brightest pixel", fontsize=12)
    save_figure(fig, FIGURES/(tag+"_native"))

    fig = Figure(figsize=(9, 4.5), layout="constrained")
    ax = fig.subplots()
    for band, color in zip(BANDS, ("#1b6180", "#765396", "#923e29"), strict=True):
        obs = sorted((x for x in samples if x["band"] == band and x["apertures"][-1]["flux_njy"] is not None),
                     key=lambda x: x["mjd"])
        ax.errorbar(pd.to_datetime([x["dateobs"] for x in obs]), [x["apertures"][-1]["flux_njy"] for x in obs],
                    np.array([x["apertures"][-1]["local_mad_sigma_njy"] for x in obs], dtype=float),
                    fmt="o", color=color, label=band, capsize=2)
    ax.set(ylabel="Signed 0.24″ aperture flux (nJy)", xlabel="UTC date")
    ax.set_title(f"GOODS-S · {candidate.ra:.6f}°, {candidate.dec:+.6f}°")
    ax.legend()
    ax.grid(alpha=.2)
    save_figure(fig, FIGURES/(tag+"_lightcurve"))


def vet(cache, candidate, products, tns):
    center = candidate.ra, candidate.dec
    selected = [r for r in products if r["pointing"] == candidate.pointing
                and r["detector"] == (candidate.sw_detector if r["band"] == "F200W" else "NRCBLONG")]
    data = {(r["epoch"], r["band"]): vet_mosaic(r, center) for r in selected}
    bright = "older" if candidate.sign < 0 else "newer"
    native, dither_images = {}, {}
    for band in ("F200W", "F444W"):
        (record,) = [r for r in selected if (r["band"], r["epoch"]) == (band, bright)]
        receipt = native_products(cache, candidate, band, record, center)
        frames, measures = native_alignment(record, receipt, center)
        native[band] = dict(receipt=receipt, **summarize_detector(measures))
        dither_images[band] = frames, measures
    other, listings = extra_products(cache, candidate, center)
    measurements = [val[1] for val in data.values()]
    measurements.extend(vet_mosaic(record, center)[1] for record in other)
    compare = {}
    for band in BANDS:
        (old, w, x, y, x0, y0), _ = data["older", band]
        (new, *unused), _ = data["newer", band]
        residual = candidate.sign*(new-old)
        scale = pixel_scale(w)
        px, py = centroids(residual, [(x, y)], round(.14/scale), search=round(.14/scale))[0]
        if np.isfinite(px+py):
            ra, dec = w.pixel_to_world_values(px+x0, py+y0)
            compare[band] = dict(positive_residual_center_icrs_deg=[float(ra), float(dec)],
                                 offset_arcsec=float(np.hypot(x-px, y-py)*scale))
        else:
            compare[band] = None
    plot_candidate(candidate, data, dither_images, measurements)
    return dict(candidate=candidate.candidate, icrs_deg=list(center), sign=int(candidate.sign), native=native,
                extra_products=other, association_receipts=listings, measurements=measurements,
                direct_difference_centroids=compare, public_matches=catalogue_checks(cache, candidate, center, tns))


def main():
    matches = OUT/"module_b_matches.csv"
    products = OUT/"module_b_products.json"
    shortlist = pd.read_csv(matches, dtype={"pointing": str}).query("passes_three_filter_gate")
    if set(shortlist.pointing) != set(EXTRA) or len(shortlist) != len(EXTRA):
        raise ValueError("The fixed follow-up expects one lead in each pointing")
    records = json.loads(products.read_text())["images"]
    cache = DownloadCache(ROOT, CACHE_MIB*1024**2)
    tns = tns_catalogue(cache)
    results = []
    for candidate in shortlist.itertuples():
        results.append(vet(cache, candidate, records, tns))
        write_json(OUT/"module_b_vetting.json", {
            "executed_utc": utcnow(), "match_list_sha256": text_digest(matches), "products_sha256": text_digest(products),
            "storage_cap_mib": CACHE_MIB, "candidates": results, "script_sha256": text_digest(__file__),
            "original_goods_s_transient_tns": tns_match(tns, SkyCoord(53.1620359953, -27.7463839937, unit="deg")),
        })
        print(candidate.candidate, "done", flush=True)


if __name__ == "__main__":
    main()

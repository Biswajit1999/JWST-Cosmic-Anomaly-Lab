"""Measure the GOODS-S brightening in six filters, including three held-out bands."""

import json

import numpy as np
from download_jwst_cross_detector import OUT, ROOT, download
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import RADII, aperture, photometry, pixel_scale, show_stamp
from scipy.ndimage import maximum_filter
from search_jwst_cross_detector import check_same_grid
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic
from vet_jwst_module_b import archive_query

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

STORAGE_MIB = 18_000
FILTERS = ("F115W", "F150W", "F200W", "F277W", "F356W", "F444W")
HELD_OUT = ("F115W", "F150W", "F277W")
OLD, NEW = "jw08060001001", "jw12577001001"
ANNULUS = (.30, .48)


def held_out_products(cache, sky):
    products, listings = [], {}
    for band in HELD_OUT:
        short = band != "F277W"
        detector = "NRCB4" if short else "NRCBLONG"
        alternate = "NRCA2" if short else "NRCALONG"
        filt = band + "-CLEAR"
        rows, receipt = archive_query(cache, "assoc_mosaic", dict(
            coords=f"{sky[0]:.8f},{sky[1]:.8f},.02", filters=filt, output="csv"),
            f"module_b_001_{band}_heldout_epochs.csv")
        listings[band] = receipt
        choices = [(visit, detector) for visit in ("jw01180010001", "jw01180012001", OLD, NEW)]
        if band != "F115W":
            choices.append(("jw05324027001", alternate))
        for visit, det in choices:
            (row,) = rows.query("visit == @visit and detector == @det and tile == 672 and status == 2").itertuples()
            product = dict(band=band, visit=visit, detector=det, dateobs=row.dateobs,
                           mjd=float(row.expstart), assoc_name=row.assoc_name, exposure_count=int(row.nexp))
            for kind in ("sci", "wht"):
                url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
                product[kind] = download(cache, url, url.rsplit("/", 1)[1], "observed_brightening_heldout_epoch")
            products.append(product)
            print(band, visit, det, flush=True)
    return products, listings


def local_flux_reference(old, new, valid, x, y, scale):
    yy, xx = np.ogrid[:old.shape[0], :old.shape[1]]
    distance = np.hypot(xx-x, yy-y)*scale
    old_valid = np.where(valid, old, -np.inf)
    peaks = maximum_filter(old_valid, size=2*round(.42/scale)+1) == old_valid
    py, px = np.nonzero(peaks & valid & (distance > 2) & (distance < 14))
    refs = []
    for i, j in zip(px, py, strict=True):
        first, second = (aperture(im, i, j, scale, .24, ANNULUS) for im in (old, new))
        if np.isfinite(first + second) and first >= 75 and second > 0:
            refs.append(dict(offset_arcsec=[float((i-x)*scale), float((j-y)*scale)],
                             old_njy=first, new_njy=second, new_to_old=second/first))
    if len(refs) < 8:
        raise ValueError(f"Too few bright local reference sources: {len(refs)}")
    ratios = np.array([ref["new_to_old"] for ref in refs])
    ratio = float(np.median(ratios))
    scatter = float(1.4826*np.median(abs(ratios-ratio)))
    return dict(median_gain=ratio, robust_scatter=scatter, references=refs)


def aligned_excess(old, new):
    o, ow, ox, oy, ox0, oy0 = old
    n, nw, nx, ny, nx0, ny0 = new
    scale = pixel_scale(ow)
    check_same_grid((o, ow, (ox0, oy0), scale), (n, nw, (nx0, ny0), pixel_scale(nw)))
    valid = np.isfinite(o) & np.isfinite(n)
    reference = local_flux_reference(o, n, valid, ox, oy, scale)
    difference = n - o
    scaled = n - reference["median_gain"]*o
    return (difference, nx, ny, scale), dict(
        unit_scale_apertures=photometry(difference, nx, ny, scale, RADII, ANNULUS),
        local_flux_scale=reference,
        local_scale_apertures=photometry(scaled, nx, ny, scale, RADII, ANNULUS),
    )


def plot(data, comparisons, stamps):
    fig = Figure(figsize=(11, 6.7), layout="constrained")
    grid = fig.add_gridspec(3, 3, height_ratios=(1.45, 1, 1))
    ax = fig.add_subplot(grid[0, :])
    points = [comparisons[band]["unit_scale_apertures"][-1] for band in FILTERS]
    for i, (band, val) in enumerate(zip(FILTERS, points, strict=True)):
        held = band in HELD_OUT
        label = ("Held-out filters" if held else "Search filters") if band in {"F115W", "F200W"} else None
        ax.errorbar(i, val["flux_njy"], yerr=val["local_mad_sigma_njy"], capsize=3, label=label,
                    marker="o" if held else "s", color="#176287" if held else "#a34d30")
    ax.axhline(0, lw=.7, color="#7c858a")
    ax.grid(axis="y", alpha=.2)
    ax.set(ylabel="Change in a 0.24″ aperture (nJy)", xticks=range(len(FILTERS)), xticklabels=FILTERS)
    ax.legend(loc="upper right", frameon=False)
    for i, band in enumerate(FILTERS):
        ax = fig.add_subplot(grid[1+i//3, i%3])
        residual, x, y, scale = stamps[band]
        half = round(.72/scale)
        ix, iy = round(x), round(y)
        region = residual[iy-half:iy+half+1, ix-half:ix+half+1]
        limit = np.nanpercentile(abs(region), 98)
        show_stamp(ax, residual, x, y, scale, (-limit, limit), "RdBu_r", .24, half_arcsec=.72)
        first, second = (data[visit, band][1]["dateobs"][:10] for visit in (OLD, NEW))
        ax.set_title(f"{band} · {first} to {second}", fontsize=9)
    fig.suptitle("JADES 103607 · JWST difference photometry and images\n"
                 "2026 minus 2025 · brighter in red · per-filter display stretch", fontsize=12)
    save_figure(fig, FIGURES/"module_b_001_six_filter_change")


def main():
    vetting_path = OUT/"module_b_vetting.json"
    search_products_path = OUT/"module_b_products.json"
    (target,) = (c for c in json.loads(vetting_path.read_text())["candidates"] if c["sign"] > 0)
    sky = target["icrs_deg"]
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    extra, listings = held_out_products(cache, sky)
    search = [p for p in json.loads(search_products_path.read_text())["images"] if p["pointing"] == "001"
              and p["detector"] == ("NRCB4" if p["band"] == "F200W" else "NRCBLONG")]
    data = {(p["visit"], p["band"]): vet_mosaic(p, sky) for p in [*search, *extra]}
    comparisons, stamps = {}, {}
    for band in FILTERS:
        stamps[band], comparisons[band] = aligned_excess(data[OLD, band][0], data[NEW, band][0])
    plot(data, comparisons, stamps)
    write_json(OUT/"module_b_001_six_filter_change.json", dict(
        executed_utc=utcnow(), candidate=target["candidate"], icrs_deg=sky, baseline_visit=OLD,
        bright_visit=NEW, held_out_filters=list(HELD_OUT), observed_apertures=[v[1] for v in data.values()],
        difference_photometry=comparisons, association_queries=listings, additional_products=extra,
        storage_cap_mib=STORAGE_MIB, search_products_sha256=text_digest(search_products_path),
        vetting_sha256=text_digest(vetting_path), script_sha256=text_digest(__file__),
        scope="Signed aperture fluxes without total-flux correction; errors are empirical local sky scatter."
              " Held-out filters were excluded from ranking and share the bright observing visit.",
    ))


if __name__ == "__main__":
    main()

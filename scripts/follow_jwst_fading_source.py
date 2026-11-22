"""Extend the GOODS-S fading lead with public F444W epochs and local controls."""

import re
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import (
    FIGURES,
    OUT,
    aperture,
    centroids,
    load_mosaic,
    photometry,
    pixel_scale,
    read_receipt,
    reference_sky,
    show_stamp,
)

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = "GS-JWST-53.195414-27.779328"
ARCHIVE_VISITS = ("jw01180007001", "jw06541003001", "jw06541005001")
SKY_ANNULUS = (.28, .44)
STORAGE_MIB = 6500


def download_archival(sky):
    ra, dec = sky
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    params = urlencode(dict(coords=f"{ra:.8f},{dec:.8f},0.025", filters="F444W-CLEAR", output="csv"))
    url = "https://grizli-cutout.herokuapp.com/assoc_mosaic?" + params

    def fetch(url, name, kind):
        path = cache.fetch(url, "epoch_search/"+name, survey="DAWN JWST Archive", version="association mosaic v7.0",
                           kind=kind)
        return dict(url=url, path=path.relative_to(ROOT).as_posix(), sha256=digest(path))

    listing = fetch(url, "fading_lead_f444w_associations.csv", "archive_association_metadata")
    records = []
    for visit in ARCHIVE_VISITS:
        rows = pd.read_csv(ROOT / listing["path"]).query("visit == @visit and tile == 672 and status == 2")
        if len(rows) != 1:
            raise ValueError(f"Expected one F444W association matching {visit}")
        row = rows.iloc[0]
        rec = dict(band=row["filter"], visit=visit, dateobs=row.dateobs, assoc_name=row.assoc_name,
                   exposure_count=int(row.nexp), mjd=float(row.expstart))
        for kind in ("sci", "wht"):
            url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
            rec[kind] = fetch(url, url.rsplit("/", 1)[1], "observed_epoch_"+kind)
        records.append(rec)
        print(visit, row.dateobs, flush=True)
    return records, listing


def catalogue_context(sky):
    ra, dec = sky
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    queries = {
        "simbad": "https://simbad.cfa.harvard.edu/simbad/sim-coo?" + urlencode({
            "Coord": f"{ra:.8f} {dec:.8f}", "CooFrame": "ICRS", "Radius": 5,
            "Radius.unit": "arcsec", "output.format": "ASCII"}),
        "tns": "https://www.wis-tns.org/search?" + urlencode({
            "ra": ra, "decl": dec, "radius": 5, "coords_unit": "arcsec", "format": "html"}),
    }
    result = {}
    for name, url in queries.items():
        try:
            path = cache.fetch(url, f"epoch_search/fading_lead_{name}.txt", survey=name, version="5 arcsec cone",
                               kind="public_position_query")
        except requests.RequestException as exc:
            result[name] = dict(url=url, status="UNAVAILABLE", error=str(exc))
            continue
        content = path.read_text(encoding="utf-8")
        metadata = dict(url=url, path=path.relative_to(ROOT).as_posix(), sha256=digest(path), radius_arcsec=5)
        if name == "simbad":
            rows = re.findall(r"^\s*\d+\s*\|\s*([\d.]+)\|([^|]+)\|([^|]+)", content, re.MULTILINE)
            count = re.search(r"Number of objects\s*:\s*(\d+)", content)
            if count and int(count[1]) == len(rows) and rows:
                nearest = min(rows, key=lambda row: float(row[0]))
                metadata.update(matches=len(rows), nearest_separation_arcsec=float(nearest[0]),
                                nearest_id=nearest[1].strip(), nearest_type=nearest[2].strip())
        elif 'class="no-results">No results found' in content:
            metadata["matches"] = 0
        result[name] = dict(status="PARSED" if "matches" in metadata else "UNVERIFIED", **metadata)
    return result


def measure(image_info, refs):
    image, wcs, x, y, x0, y0 = image_info
    scale = pixel_scale(wcs)
    rx, ry = wcs.world_to_pixel_values(*refs)
    control = np.array([aperture(image, cx-x0, cy-y0, scale, .24, SKY_ANNULUS)
                        for cx, cy in zip(rx, ry, strict=True)])
    flux = photometry(image, x, y, scale, (.24,), SKY_ANNULUS)[0]
    cx, cy = centroids(image, [(x, y)], core=4.5, search=5)[0]
    recentered = aperture(image, cx, cy, scale, .24, SKY_ANNULUS) if np.isfinite(cx+cy) else np.nan
    return flux, recentered, control


def figure_for(rows, mosaics):
    figure = Figure(figsize=(10.3, 7.4), layout="constrained")
    grid = figure.add_gridspec(3, 3, height_ratios=[1.65, 1, 1])
    chart = figure.add_subplot(grid[0, :])
    dates = [datetime.fromisoformat(r["utc_start"]) for r in rows]
    chart.errorbar(dates, np.array([r["scaled_fixed_flux_njy"] for r in rows], dtype=float),
                   yerr=np.array([r["scaled_blank_scatter_njy"] for r in rows], dtype=float),
                   fmt="o", color="#216680", label="Fixed aperture, scaled using local sources")
    chart.plot(dates, np.array([r["scaled_recentered_flux_njy"] for r in rows], dtype=float),
               "x", color="#b94f34", label="Recentered aperture, same scaling")
    unscaled = [i for i, r in enumerate(rows) if r["control_epoch_to_2025_median_ratio"] is None]
    if unscaled:
        chart.errorbar([dates[i] for i in unscaled],
                       np.array([rows[i]["raw_fixed"]["flux_njy"] for i in unscaled], dtype=float),
                       yerr=np.array([rows[i]["raw_fixed"]["local_mad_sigma_njy"] for i in unscaled], dtype=float),
                       fmt="D", mfc="white", color="#626a6d", label="Fixed aperture, insufficient controls")
    chart.set_ylabel("F444W aperture flux (nJy)")
    chart.grid(axis="y", alpha=.25)
    chart.legend(loc="best", fontsize=9)
    peaks = []
    for image, _, x, y, _, _ in mosaics:
        peaks.append(image[round(y)-5:round(y)+6, round(x)-5:round(x)+6].ravel())
    bounds = tuple(np.nanpercentile(np.concatenate(peaks), [3, 99]))
    for k, (record, (image, wcs, x, y, _, _)) in enumerate(zip(rows, mosaics, strict=True)):
        ax = figure.add_subplot(grid[1+k//3, k % 3])
        show_stamp(ax, image, x, y, pixel_scale(wcs), bounds, "gray", .24, half_arcsec=.55)
        ax.set_title(record["utc_start"][:10], fontsize=9)
    figure.suptitle("GOODS-S variability lead · RA 53.195414°, Dec −27.779328°\n"
                    "Observed F444W cutouts · 0.24″ aperture · common image stretch", fontsize=13)
    return figure


def main():
    fading = read_receipt("candidate_photometry.json")["fading_shortlist_2025"]
    lead = next(r for r in fading if r["candidate"] == CANDIDATE)
    sky = lead["ra"], lead["dec"]
    additional, association = download_archival(sky)
    initial = read_receipt("products.json")["images"] + read_receipt("target_followup_products.json")["images"]
    records = sorted(additional + [r for r in initial if r["band"] == "F444W-CLEAR"], key=lambda r: r["mjd"])
    mosaics = [load_mosaic(r, sky, half=425) for r in records]
    anchor = next(i for i, r in enumerate(records) if r["dateobs"].startswith("2025"))
    refs = reference_sky(mosaics[anchor])
    measurements = [measure(mosaic, refs) for mosaic in mosaics]
    base = measurements[anchor][2]
    rows = []
    for record, (raw, recentered, controls) in zip(records, measurements, strict=True):
        good = np.isfinite(controls+base) & (base >= 50) & (controls > 0)
        factors = controls[good]/base[good]
        if len(factors) < 6:
            factor = scatter = None
        else:
            factor = float(np.median(factors))
            scatter = float(1.4826*np.median(abs(factors-factor)))
        rows.append({
            "visit": record["visit"], "utc_start": record["dateobs"], "mjd": record["mjd"], "raw_fixed": raw,
            "raw_recentered_flux_njy": float(recentered) if np.isfinite(recentered) else None,
            "usable_control_sources": len(factors), "control_epoch_to_2025_median_ratio": factor,
            "control_ratio_mad_sigma": scatter,
            "scaled_fixed_flux_njy": raw["flux_njy"]/factor if factor and raw["flux_njy"] is not None else None,
            "scaled_recentered_flux_njy": float(recentered/factor) if factor and np.isfinite(recentered) else None,
            "scaled_blank_scatter_njy": raw["local_mad_sigma_njy"]/factor if factor and raw["local_mad_sigma_njy"] is not None else None,
        })
    write_json(OUT / "fading_lead_epochs.json", {
        "executed_utc": utcnow(), "candidate": CANDIDATE, "ra": sky[0], "dec": sky[1],
        "selection": "Existing source with F444W aperture brightening in 2025 relative to 2023 and 2026",
        "aperture_radius_arcsec": .24, "sky_annulus_arcsec": SKY_ANNULUS,
        "control_scaling": "Median same-aperture epoch/2025 ratio of nearby 2025 peaks with baseline flux >= 50 nJy",
        "uncertainties": "Scaled local blank-aperture scatter excludes control ratio scatter and PSF systematics",
        "archive_association_listing": association, "products": records, "epochs": rows,
        "catalogue_context": catalogue_context(sky),
        "candidate_receipt_sha256": text_digest(OUT / "candidate_photometry.json"),
        "script_sha256": text_digest(__file__),
    })
    save_figure(figure_for(rows, mosaics), FIGURES / "fading_lead_epochs")


if __name__ == "__main__":
    main()

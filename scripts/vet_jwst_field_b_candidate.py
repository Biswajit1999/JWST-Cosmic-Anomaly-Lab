"""Test the field B fading lead in independent epochs and native exposures."""

import json
import re
from datetime import datetime
from io import StringIO
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from crossmatch_jwst_epoch_candidate import published_positions
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_field_b_candidate import STORAGE_MIB
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import (
    audit_detector_target,
    centroids,
    load_dithers,
    photometry,
    pixel_scale,
    reference_sky,
    show_stamp,
)
from search_jwst_cross_detector import check_same_grid, window

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

FIGURES = ROOT / "figures/jwst_cross_detector"
STYLES = {"F200W": "#126785", "F356W": "#7d5598", "F444W": "#9a542a"}


def check_catalogues(sky):
    ra, dec = sky
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    queries = {
        "tns": ("https://www.wis-tns.org/search", dict(ra=f"{ra:.6f}", decl=f"{dec:.6f}", radius=5,
                                                           coords_unit="arcsec", format="html")),
        "simbad": ("https://simbad.cfa.harvard.edu/simbad/sim-coo",
                   {"Coord": f"{ra:.8f} {dec:.8f}", "CooFrame": "ICRS", "Radius": 5,
                    "Radius.unit": "arcsec", "output.format": "ASCII"}),
        "nirspec": (API+"/nirspec_extractions", dict(coords=f"{ra:.8f},{dec:.8f}", size=2, output="csv")),
        "jades_dr2": ("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv", {
            "-source": "V/159/gssample", "-c": f"{ra:.8f} {dec:.8f}", "-c.rs": 3, "-out.max": 50,
            "-out": "NIRCamID,RA_ICRS,DE_ICRS,KronRad,zphot,zphotl68,zphotu68",
        }),
    }
    results = {}
    for label, (base, params) in queries.items():
        url = base + "?" + urlencode(params)
        try:
            receipt = download(cache, url, f"field_b_candidate_{label}.txt", "public_coordinate_check")
        except requests.RequestException as exc:
            results[label] = dict(status="UNAVAILABLE", url=url, error=str(exc))
            continue
        response = (ROOT / receipt["path"]).read_text(encoding="utf-8")
        metadata = dict(status="RETRIEVED", **receipt)
        if label == "tns" and 'class="no-results">No results found' in response:
            metadata.update(matches=0, radius_arcsec=5)
        if label == "nirspec" and response.lstrip().startswith("Nothing found for"):
            metadata.update(matches=0, side_arcsec=4)
        if label == "simbad":
            rows = re.findall(r"^\s*\d+\s*\|\s*([\d.]+)\|([^|]+)\|([^|]+)", response, re.MULTILINE)
            count = re.search(r"Number of objects\s*:\s*(\d+)", response)
            if count and len(rows) == int(count[1]):
                metadata.update(matches=len(rows), radius_arcsec=5, objects=[
                    dict(separation_arcsec=float(distance), name=name.strip(), type=kind.strip())
                    for distance, name, kind in rows
                ])
        if label == "jades_dr2" and "NIRCamID\tRA_ICRS\tDE_ICRS" in response:
            table = pd.read_csv(StringIO(response), comment="#", sep="\t")
            table["NIRCamID"] = pd.to_numeric(table.NIRCamID, errors="coerce")
            table = table[table.NIRCamID.notna()].astype(float)
            separation = SkyCoord(*sky, unit="deg").separation(
                SkyCoord(table.RA_ICRS, table.DE_ICRS, unit="deg")).arcsec
            sources = table.assign(separation_arcsec=separation).sort_values("separation_arcsec")
            metadata.update(matches=len(sources), radius_arcsec=3, sources=[
                dict(nircam_id=int(row.NIRCamID), ra=row.RA_ICRS, dec=row.DE_ICRS,
                     separation_arcsec=row.separation_arcsec, kron_radius_arcsec=row.KronRad,
                     photo_z=row.zphot, photo_z_68=[row.zphotl68, row.zphotu68])
                for row in sources.replace({np.nan: None}).itertuples()
            ])
        results[label] = metadata
    paper = "https://arxiv.org/html/2406.05060v4"
    try:
        primary = cache.fetch(paper, "epoch_search/decoursey2024_arxiv_v4.html",
                              survey="DeCoursey et al.", version="arXiv v4", kind="public_novelty_check")
    except requests.RequestException as exc:
        results["jades_2022_2023_transients"] = dict(status="UNAVAILABLE", url=paper, error=str(exc))
    else:
        ras, decs = zip(*published_positions(primary.read_text(encoding="utf-8")), strict=True)
        separation = SkyCoord(*sky, unit="deg").separation(SkyCoord(list(ras), list(decs), unit=("hourangle", "deg"))).arcsec
        results["jades_2022_2023_transients"] = dict(
            url=paper, path=primary.relative_to(ROOT).as_posix(), sha256=digest(primary), searched_positions=len(separation),
            nearest_separation_arcsec=float(separation.min()), matches=int((separation <= 5).sum()), radius_arcsec=5)
    return results


def vet_mosaic(record, sky):
    science, weight = (window(record, kind, sky, 16.5) for kind in ("sci", "wht"))
    check_same_grid(science, weight)
    image, wcs, (x0, y0), scale = science
    image[~np.isfinite(weight[0]) | (weight[0] <= 0)] = np.nan
    x, y = np.array(wcs.world_to_pixel_values(*sky))-[x0, y0]
    mosaic = image, wcs, x, y, x0, y0
    band = record["band"].split("-")[0]
    measurements = photometry(image, x, y, scale, (.12, .18, .24), (.30, .48))
    center = centroids(image, [(x, y)], round(.13/scale), search=round(.14/scale))[0]
    measured_sky = [float(v) for v in wcs.pixel_to_world_values(*(center+[x0, y0]))] if np.isfinite(center).all() else None
    row = dict(band=band, detector=record["detector"], visit=record["visit"], dateobs=record["dateobs"],
               mjd=record["mjd"], apertures=measurements, local_centroid_ra_dec=measured_sky)
    return mosaic, row


def native_check(products, mosaic, sky):
    frames = load_dithers(products, reference_sky(mosaic))
    rows = audit_detector_target(frames, sky)
    for row, frame in zip(rows, frames, strict=True):
        x, y = row["forced_center_pixels_xy"]
        centroid = centroids(frame["clean"], [(x, y)], core=3.5)[0]
        if np.isfinite(centroid).all():
            ra, dec = frame["wcs"].pixel_to_world_values(*(centroid-frame["translation"]))
            row["clean_centroid_ra_dec"] = [float(ra), float(dec)]
            row["clean_centroid_offset_arcsec"] = float(np.hypot(*(centroid-[x, y]))*frame["scale"])
        else:
            row["clean_centroid_ra_dec"] = row["clean_centroid_offset_arcsec"] = None
    return frames, rows


def plot_epochs(measurements, stamps):
    figure = Figure(figsize=(11.3, 8.3), layout="constrained")
    grid = figure.add_gridspec(3, 4, height_ratios=(1.44, 1, 1))
    ax = figure.add_subplot(grid[0, :])
    for band, colour in STYLES.items():
        rows = [m for m in measurements if m["band"] == band and m["apertures"][-1]["flux_njy"] is not None
                and m["apertures"][-1]["local_mad_sigma_njy"] is not None]
        ax.errorbar([datetime.fromisoformat(r["dateobs"]) for r in rows], [r["apertures"][-1]["flux_njy"] for r in rows],
                    yerr=[r["apertures"][-1]["local_mad_sigma_njy"] for r in rows], color=colour,
                    marker="o", linestyle="none", capsize=2, label=band)
    ax.grid(axis="y", alpha=.27)
    ax.legend(ncols=3, loc="upper left")
    ax.set_ylabel("Aperture flux (nJy)")
    by_band = [m for m in measurements if m["band"] == "F200W"]
    if len(by_band) != 8:
        raise ValueError("Expected eight F200W observations for the epoch figure")
    bright = next(stamps["F200W", r["visit"]] for r in by_band if r["visit"] == "jw08060002001")
    image, _, x, y, _, _ = bright
    bright_core = image[round(y)-7:round(y)+8, round(x)-7:round(x)+8]
    bounds = (-.04, float(np.nanpercentile(bright_core, 98)))
    for i, row in enumerate(by_band):
        target_ax = figure.add_subplot(grid[1+i//4, i%4])
        image, wcs, x, y, _, _ = stamps["F200W", row["visit"]]
        show_stamp(target_ax, image, x, y, pixel_scale(wcs), bounds, "gray", .24, half_arcsec=.76)
        target_ax.set_title(row["dateobs"][:10] + " · " + row["detector"], fontsize=9)
    figure.suptitle("GOODS-S fading transient candidate · RA 53.162036°, Dec −27.746384°\n"
                    "Public JWST measurements · 0.24″ aperture · observed F200W cutouts", fontsize=13)
    save_figure(figure, FIGURES / "field_b_transient_epochs")


def plot_native(groups):
    figure = Figure(figsize=(9.0, 6.0), layout="constrained")
    axes = figure.subplots(2, 3)
    for row, (band, (frames, records)) in zip(axes, groups.items(), strict=True):
        for i, (ax, frame, record) in enumerate(zip(row, frames, records, strict=True), 1):
            x, y = record["forced_center_pixels_xy"]
            scale = frame["scale"]
            image = frame["image"]
            radius = round(.58/scale)
            ix, iy = round(x), round(y)
            local = image[iy-radius:iy+radius+1, ix-radius:ix+radius+1]
            bounds = tuple(np.nanpercentile(local, [5, 98.5]))
            show_stamp(ax, image, x, y, scale, bounds, "gray", .16, half_arcsec=.58)
            det_x, det_y = record["unmasked_peak_pixels_xy"]
            ax.scatter((det_x-x)*scale, (det_y-y)*scale, marker="+", color="#ff819c", s=50, linewidths=1)
            ax.set_title(f"{band} · dither {i} · peak ({det_x}, {det_y})", fontsize=9)
    figure.suptitle("Bright epoch · six native NIRCam detector exposures\n"
                    "Cyan aperture at the same sky position · pink cross at the native peak", fontsize=13)
    save_figure(figure, FIGURES / "field_b_transient_native")


def main():
    inpath = OUT / "field_b_candidate_products.json"
    data = json.loads(inpath.read_text())
    sky = data["ra"], data["dec"]
    results, stamps = [], {}
    rows = sorted(data["survey_visits"]+data["older_epochs"], key=lambda r: r["mjd"])
    for record in rows:
        mosaic, row = vet_mosaic(record, sky)
        results.append(row)
        stamps[row["band"], row["visit"]] = mosaic
        print(row["dateobs"], row["band"], row["apertures"][-1]["flux_njy"], flush=True)
    frames = {band: native_check(data["fitted_exposures"][band], stamps[band, "jw08060002001"], sky)
              for band in ("F200W", "F444W")}
    plot_epochs(results, stamps)
    plot_native(frames)
    write_json(OUT / "field_b_transient_vetting.json", {
        "executed_utc": utcnow(), "candidate": data["candidate"], "position_degrees": sky,
        "observed_aperture_flux": results, "photometry_correction": "Signed local background removed, no aperture correction",
        "native_bright_visit": {b: records for b, (_, records) in frames.items()},
        "public_prior_checks": check_catalogues(sky), "input_sha256": text_digest(inpath),
        "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

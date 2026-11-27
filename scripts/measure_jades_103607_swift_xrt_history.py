"""Measure Swift XRT archive coverage and forced count rates at JADES 103607."""

import json
import warnings
from collections import Counter, defaultdict
from datetime import datetime
from urllib.parse import urlencode
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import FITSFixedWarning
from check_jades_103607_independent_2026 import HEASARC
from check_jades_103607_swift import (
    XRT_BACKGROUND_ARCSEC,
    XRT_BANDS,
    XRT_SOURCE_ARCSEC,
    matching_xrt_timing,
    measure_xrt,
    position_on_image,
)
from download_jades_103607_swift_history import PLAN
from download_jades_103607_swift_history import TAG as DOWNLOAD_TAG
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.colors import LogNorm
from matplotlib.dates import DateFormatter, YearLocator
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from measure_jwst_epoch_candidate import verified_path
from scipy.stats import binom, chi2
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_swift_xrt_history"
PRIOR = OUT / "goods_s_prior_records.json"
SXPS_API = "https://www.swift.ac.uk/API/main.php"
SXPS_BANDS = ("Total", "Soft", "Medium", "Hard")
PAIR_ARCSEC = 6.0
NATIVE_BANDS = (*XRT_BANDS, (0.3, 1), (1, 2), (2, 10))
SXPS_TO_NATIVE_BAND = dict(Total=(0.3, 10), Soft=(0.3, 1), Medium=(1, 2), Hard=(2, 10))


def historical_neighbors(point):
    prior = json.loads(PRIOR.read_text(encoding="utf-8"))
    chandra = prior["literature_catalogues"]["luo_2017"]
    matches = [
        record["nearest"] for key, record in chandra["positions"].items()
        if key.startswith("GS-001-B")
    ]
    if len(matches) != 1:
        raise ValueError("Expected one Chandra neighbour list for GS-001-B")
    neighbors = []
    for row in matches[0]:
        coord = SkyCoord(row["RAJ2000"], row["DEJ2000"], frame="fk5", unit="deg").icrs
        offset = point.spherical_offsets_to(SkyCoord(coord.ra, coord.dec, frame="icrs"))
        neighbors.append(dict(
            **row, offset_east_north_arcsec=[angle.arcsec for angle in offset],
            position_icrs_deg=[coord.ra.deg, coord.dec.deg],
        ))
    if not (2*PAIR_ARCSEC < neighbors[0]["separation_arcsec"] < XRT_SOURCE_ARCSEC-PAIR_ARCSEC):
        raise ValueError("Paired 6 arcsec apertures must be disjoint and inside the main aperture")
    catalog = pd.read_csv(verified_path(chandra["product"]), comment="#", sep="\t")
    nearest = neighbors[0]
    same = (
        pd.to_numeric(catalog["RAJ2000"], errors="coerce").eq(nearest["RAJ2000"])
        & pd.to_numeric(catalog["DEJ2000"], errors="coerce").eq(nearest["DEJ2000"])
    )
    if same.sum() != 1:
        raise ValueError("Chandra catalog does not reproduce the stored nearest position")
    entry = catalog.loc[same].iloc[0]
    nearest["chandra_full_band_aperture_corrected_counts"] = float(entry["FB"])
    nearest["chandra_full_band_flux_erg_cm2_s"] = float(entry["FFB"])
    return neighbors, chandra["product"]


def exposure_coverage(path, point, neighbor):
    with fits.open(path) as hdus:
        hdu = hdus[0]
        try:
            x, y, scale = position_on_image(hdu, point)
        except ValueError:
            return dict(status="position outside map array")
        px = np.arange(max(0, int(x)-30), min(hdu.data.shape[1], int(x)+31))
        py = np.arange(max(0, int(y)-30), min(hdu.data.shape[0], int(y)+31))
        pixels = np.asarray(hdu.data[np.ix_(py, px)], float)
        radius = np.hypot(px[None, :]-x, py[:, None]-y) * scale
        core = radius <= XRT_SOURCE_ARCSEC
        ring = (radius >= XRT_BACKGROUND_ARCSEC[0]) & (radius <= XRT_BACKGROUND_ARCSEC[1])
        valid = (pixels > 0) & np.isfinite(pixels)
        point_exposure = float(hdu.data[round(y), round(x)])
        metadata = dict(
            point_exposure_s=point_exposure if np.isfinite(point_exposure) else None,
            source_positive_exposure_fraction=float(valid[core].mean()),
            background_positive_exposure_fraction=float(valid[ring].mean()),
        )
        status = "full exposure"
        if metadata["source_positive_exposure_fraction"] == 0:
            status = "no exposure in source aperture"
        elif not valid[core].all():
            status = "partial exposure in source aperture"
        elif not valid[ring].any():
            status = "no exposure in background annulus"
        if status != "full exposure":
            return dict(status=status, **metadata)
        nx, ny, _ = position_on_image(hdu, neighbor)
        nr = np.hypot(px[None, :]-nx, py[:, None]-ny) * scale
        metadata["paired_aperture_mean_exposure_s"] = [
            float(pixels[radius <= PAIR_ARCSEC].mean()),
            float(pixels[nr <= PAIR_ARCSEC].mean()),
        ]
        return dict(status=status, **metadata)


def select_events(products):
    candidates = [label for label in PLAN if label.startswith("events_")
                  and products[label]["status"] == "available"]
    if not candidates:
        return None, "PC cleaned event file unavailable"
    map_header = fits.getheader(ROOT / products["exposure_map"]["path"])
    matched = [
        label for label in candidates if matching_xrt_timing(
            fits.getheader(ROOT / products[label]["path"], "EVENTS"), map_header
        )
    ]
    if len(matched) != 1:
        return None, "no unique PC event file matches exposure map timing"
    return PLAN[matched[0]][1], "matching file"


def measure_history(cache, ledger, point, neighbors):
    neighbor = SkyCoord(*neighbors[0]["position_icrs_deg"], unit="deg", frame="icrs")
    near_x, near_y = neighbors[0]["offset_east_north_arcsec"]
    rows, event_positions = [], {}
    for observation in ledger["observation_rows"]:
        prod = observation["products"]
        row = dict(
            obsid=observation["obsid"], nominal_start_utc=observation["start_utc"],
            pointing_separation_arcmin=observation["pointing_separation_arcmin"],
        )
        rows.append(row)
        if not prod:
            row["status"] = "no nominal XRT exposure in master log"
            continue
        if prod["exposure_map"]["status"] != "available":
            row["status"] = "PC exposure map unavailable"
            continue
        try:
            coverage = exposure_coverage(ROOT / prod["exposure_map"]["path"], point, neighbor)
            row["map_coverage"] = coverage
            if coverage["status"] != "full exposure":
                row["status"] = coverage["status"]
                continue
            suffix, status = select_events(prod)
            if suffix is None:
                row["status"] = status
                continue
            xrt, offsets = measure_xrt(cache, observation, point, event_suffix=suffix, bands=NATIVE_BANDS)
            row.update(status="measured", xrt=xrt, event_file_suffix=suffix)
            tests = []
            for band in xrt["bands"]:
                n, b = band["source_counts"], band["annulus_counts"]
                alpha = xrt["background_exposure_area_scale"]
                tests.append(float(binom.sf(n-1, n+b, alpha/(1+alpha))))
            row["conditional_on_off_tails"] = tests
            event_positions[observation["obsid"]] = offsets
            east, north = offsets
            row["paired_aperture_total_counts"] = [
                int(np.count_nonzero(np.hypot(east, north) <= PAIR_ARCSEC)),
                int(np.count_nonzero(np.hypot(east-near_x, north-near_y) <= PAIR_ARCSEC)),
            ]
        except (KeyError, OSError, ValueError) as exc:
            row.update(status="PC product validation failed", detail=str(exc))
    return rows, event_positions


def summarize_bands(rows):
    exposure = sum(r["xrt"]["mean_source_aperture_exposure_s"] for r in rows)
    summary = []
    for index, first in enumerate(rows[0]["xrt"]["bands"]):
        distribution = np.ones(1)
        n = b = 0
        expectation = 0.0
        for row in rows:
            band = row["xrt"]["bands"][index]
            source, ring = band["source_counts"], band["annulus_counts"]
            alpha = row["xrt"]["background_exposure_area_scale"]
            n += source
            b += ring
            expectation += band["estimated_aperture_background_counts"]
            if source + ring:
                distribution = np.convolve(
                    distribution, binom.pmf(np.arange(source+ring+1), source+ring, alpha/(1+alpha))
                )
        summary.append(dict(
            energy_kev=first["energy_kev_lower_inclusive_upper_exclusive"],
            source_counts=n, annulus_counts=b, estimated_background_counts=expectation,
            signed_excess_counts=n-expectation,
            conditional_background_tail=min(1.0, float(distribution[n:].sum())),
            gross_95_percent_rate_per_s_upper=float(chi2.ppf(.95, 2*(n+1)) / (2*exposure)),
        ))
    return dict(visits=len(rows), summed_aperture_exposure_s=exposure, bands=summary)


def sxps_search(cache, coords, rows):
    payload = dict(
        APIFunc="getSXPSUL", APIVersion="1.0.5", catalogue="LSXPS", RA=coords[0], Dec=coords[1],
        bands=[band.lower() for band in SXPS_BANDS], whichData="all", timeFormat="MJD",
        sigma=3, detectionsAsRates=True, skipDetections=False, useCatValues=False,
    )
    path = cache.fetch(
        SXPS_API, f"metadata/{TAG}_lsxps.json", survey="UK Swift Science Data Centre",
        version="LSXPS getSXPSUL forced photometry", kind="observed_rate_and_sensitivity",
        json_payload=payload,
    )
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("OK") != 1 or not isinstance(data.get("ULData"), list):
        raise ValueError(f"Unexpected LSXPS upper-limit response: {sorted(data)}")
    by_obs = {row["obsid"]: row for row in rows}
    datasets, flags = [], []
    seen = set()
    for value in data["ULData"]:
        obs = value["ObsID"]
        if len(obs) != 11 or not obs.isdigit() or obs in seen:
            raise ValueError(f"Unexpected LSXPS observation identifier {obs}")
        seen.add(obs)
        src = dict(
            obsid=obs, stack=obs.startswith("1"),
            start_mjd=float(value["StartTime"]), stop_mjd=float(value["StopTime"]),
            source_exposure_s=float(value["SourceExposure"]),
            native_status=by_obs.get(obs, {}).get("status"), bands={},
        )
        if not (src["start_mjd"] <= src["stop_mjd"] and src["source_exposure_s"] > 0):
            raise ValueError(f"Invalid date or exposure in LSXPS result {obs}")
        for band in SXPS_BANDS:
            def opt(key):
                number = value.get(f"{band}_{key}")
                return None if number is None or not np.isfinite(float(number)) else float(number)

            detected = opt("IsDetected")
            if detected not in (0, 1, None):
                raise ValueError(f"Unexpected LSXPS detection flag {detected}")
            measure = dict(
                detected=None if detected is None else bool(detected),
                upper_limit_3sigma_per_s=opt("UpperLimit"),
                rate_per_s=opt("Rate"), rate_error_minus=opt("RateNeg"), rate_error_plus=opt("RatePos"),
                counts=opt("Counts"), background_counts=opt("BGCounts"), correction_factor=opt("CorrectionFactor"),
            )
            src["bands"][band] = measure
            if measure["detected"]:
                native_row = by_obs.get(obs, {})
                native = None
                if native_row.get("status") == "measured":
                    index = NATIVE_BANDS.index(SXPS_TO_NATIVE_BAND[band])
                    native = dict(
                        **native_row["xrt"]["bands"][index],
                        conditional_background_tail=native_row["conditional_on_off_tails"][index],
                    )
                flags.append(dict(
                    obsid=obs, band=band, stack=src["stack"],
                    start_utc=Time(src["start_mjd"], format="mjd", scale="utc").isot,
                    native_same_band=native, **measure,
                ))
        datasets.append(src)
    return dict(
        query=payload, receipt=receipt(path, SXPS_API), datasets=datasets, detection_flags=flags,
        nearby_catalogued_sources=data.get("DetData", []),
    )


def sxps_neighbors(cache, point):
    query = (
        "SELECT TOP 100 name,ra,dec,error_radius,count_rate_fb,rate_fb_flag FROM swiftlsxps "
        f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{point.ra.deg:.12f},"
        f"{point.dec.deg:.12f},0.03))=1"
    )
    url = f"{HEASARC}?{urlencode(dict(REQUEST='doQuery', FORMAT='votable', LANG='ADQL', QUERY=query))}"
    path = cache.fetch(
        url, f"metadata/{TAG}_lsxps_cone.vot", survey="NASA HEASARC",
        version="Swift LSXPS catalogue positional query", kind="observed_source_catalogue",
    )
    statuses = [node.get("value") for node in ElementTree.parse(path).findall(".//{*}INFO")
                if node.get("name") == "QUERY_STATUS"]
    table = Table.read(path, format="votable", use_names_over_ids=True)
    if statuses != ["OK"] or len(table) >= 100:
        raise ValueError("Unexpected or truncated LSXPS catalogue cone")
    table["separation_arcsec"] = point.separation(
        SkyCoord(table["ra"], table["dec"], unit="deg"),
    ).arcsec
    table.sort("separation_arcsec")
    return dict(
        radius_arcsec=108, receipt=receipt(path, url),
        rows=json.loads(table.to_pandas().to_json(orient="records", double_precision=12)),
    )


def paired_counts(measured):
    counts = np.array([r["paired_aperture_total_counts"] for r in measured]).sum(axis=0)
    times = np.array([r["map_coverage"]["paired_aperture_mean_exposure_s"] for r in measured]).sum(axis=0)
    return dict(
        aperture_radius_arcsec=PAIR_ARCSEC, order=["JWST position", "nearest 7 Ms Chandra position"],
        counts=counts.tolist(), summed_mean_exposure_s=times.tolist(),
    )


def plot(measured, positions, neighbors, peak):
    fig = Figure(figsize=(11.3, 8.3), layout="constrained")
    spec = fig.add_gridspec(2, 2, height_ratios=(1, 1.24))
    top = fig.add_subplot(spec[0, :])
    dates = [datetime.fromisoformat(row["xrt"]["date_obs_utc"]) for row in measured]
    limits = [row["xrt"]["bands"][0]["gross_95_percent_aperture_rate_per_s_upper"] for row in measured]
    top.scatter(dates, limits, marker="v", s=16, color="#17687c", alpha=.68)
    top.set(yscale="log", ylabel="95% gross aperture rate bound (counts/s)",
            title="Swift photon-counting observations at JADES 103607")
    top.xaxis.set_major_locator(YearLocator(3))
    top.xaxis.set_major_formatter(DateFormatter("%Y"))
    top.grid(alpha=.2)
    all_events = np.concatenate([np.asarray(positions[r["obsid"]]) for r in measured], axis=1)
    stacked = fig.add_subplot(spec[1, 0])
    single = fig.add_subplot(spec[1, 1])
    for axis, title, offsets in (
        (stacked, "All usable visits, 0.3–10 keV", all_events),
        (single, f"Strongest nominal visit, {peak['xrt']['date_obs_utc'][:10]}", positions[peak["obsid"]]),
    ):
        hist, xedge, yedge = np.histogram2d(*offsets, bins=np.linspace(-70, 70, 41))
        vmax = max(2, float(hist.max()))
        img = axis.imshow(
            np.ma.masked_less(hist.T, 1), cmap="Blues", origin="lower",
            extent=(xedge[0], xedge[-1], yedge[0], yedge[-1]),
            interpolation="nearest", norm=LogNorm(vmin=1, vmax=vmax),
        )
        axis.add_patch(Circle((0, 0), XRT_SOURCE_ARCSEC, fill=False, color="#98472a", lw=1.35))
        for radius in XRT_BACKGROUND_ARCSEC:
            axis.add_patch(Circle((0, 0), radius, fill=False, color="#808080", linestyle=":", lw=.85))
        axis.scatter(0, 0, marker="+", s=60, color="#98472a", label="JWST position")
        axis.scatter(*np.asarray([row["offset_east_north_arcsec"] for row in neighbors]).T,
                     marker="x", s=35, color="#476e43", label="Nearby Chandra sources")
        axis.set(title=title, xlabel="East (arcsec)", ylabel="North (arcsec)")
        fig.colorbar(img, ax=axis, label="Observed events per bin", shrink=.72)
    stacked.legend(loc="upper right", fontsize=8)
    FIGURES.mkdir(parents=True, exist_ok=True)
    save_figure(fig, FIGURES / TAG)


def main():
    path = OUT / f"{DOWNLOAD_TAG}.json"
    ledger = json.loads(path.read_text(encoding="utf-8"))
    coords = ledger["position_icrs_deg"]
    point = SkyCoord(*coords, unit="deg", frame="icrs")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    neighbors, chandra_catalog = historical_neighbors(point)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*RADECSYS.*", category=FITSFixedWarning)
        rows, positions = measure_history(cache, ledger, point, neighbors)
    measured = [row for row in rows if row["status"] == "measured"]
    if not measured:
        raise ValueError("No usable Swift PC coverage found")
    counts = Counter(row["status"] for row in rows)
    print(dict(counts), flush=True)
    bands = summarize_bands(measured)
    print(bands, flush=True)
    groups = defaultdict(list)
    for row in measured:
        groups[row["xrt"]["date_obs_utc"][:4]].append(row)
    tests = [(tail, i, row) for row in measured for i, tail in enumerate(row["conditional_on_off_tails"])]
    tail, index, peak = min(tests, key=lambda item: item[0])
    sxps = sxps_search(cache, coords, rows)
    sources = sxps_neighbors(cache, point)
    print("LSXPS dataset count:", len(sxps["datasets"]), "flags:", sxps["detection_flags"], flush=True)
    plot(measured, positions, neighbors, peak)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        product_ledger_sha256=text_digest(path), chandra_prior_results_sha256=text_digest(PRIOR),
        position_icrs_deg=coords, observation_status_counts=dict(counts), native_observations=rows,
        native_summed_measurement=bands,
        native_by_year={key: summarize_bands(group) for key, group in sorted(groups.items())},
        native_best_individual=dict(
            obsid=peak["obsid"], utc=peak["xrt"]["date_obs_utc"], band=peak["xrt"]["bands"][index],
            conditional_tail=tail, tests=len(tests), bonferroni_p=min(1.0, len(tests)*tail),
        ),
        nearest_chandra_sources=neighbors, target_and_nearest_neighbor=paired_counts(measured),
        chandra_catalogue_receipt=chandra_catalog,
        lsxps=sxps,
        lsxps_catalogue_cone=sources,
        scope="Native screened XRT observations and LSXPS rates can share the same photons. "
        "The native conditional on/off test accounts for each observation's source-to-background exposure ratio. "
        "The aggregate and paired-aperture counts can include nearby sources. Native rate bounds omit "
        "PSF and spectral response corrections. The 3-sigma LSXPS stack rates have overlapping time coverage.",
    ))


if __name__ == "__main__":
    main()

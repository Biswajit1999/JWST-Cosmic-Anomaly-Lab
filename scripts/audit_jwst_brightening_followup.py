"""Audit follow-up coverage and measure two earlier visits of JADES 103607."""

import json
from datetime import datetime, timezone

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.mast import Observations
from check_jwst_transient_precursor import contains
from download_jwst_cross_detector import OUT, ROOT, download
from measure_jwst_module_b_brightening import STORAGE_MIB
from vet_jwst_field_b_candidate import vet_mosaic
from vet_jwst_module_b import archive_query

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

BANDS = ("F150W", "F200W", "F277W", "F356W")
EARLIER = ("jw08060004001", "jw05324026001")
BRIGHT_VISIT = "jw12577001001"


def mast_followup(sky):
    Observations.TIMEOUT = 40
    selection = dict(obs_collection="JWST", proposal_id="12577")
    table = Observations.query_criteria(**selection)
    if not len(table):
        raise ValueError("MAST returned no programme data")
    visits = {}
    for visit, obs in (("jw12577003001", "jw12577-o003"), ("jw12577004001", "jw12577-o004")):
        rows = []
        for row in table:
            if not row["obs_id"].startswith((visit, obs)):
                continue
            footprint = str(row["s_region"])
            if footprint.startswith("POLYGON "):
                points = np.asarray(footprint.split()[1:], dtype=float).reshape(-1, 2)
                intersects = contains(str(points.tolist()), sky) if len(points) == 4 else None
            else:
                intersects = None
            t_min = row["t_min"]
            start = float(t_min) if not np.ma.is_masked(t_min) and np.isfinite(t_min) else None
            rows.append(dict(obs_id=str(row["obs_id"]), filters=str(row["filters"]),
                             start_mjd=start, calibration_level=int(row["calib_level"]),
                             rights=str(row["dataRights"]),
                             data_url=None if np.ma.is_masked(row["dataURL"]) else str(row["dataURL"]),
                             s_region=footprint, predicted_polygon_contains_position=intersects))
        if not rows:
            raise ValueError(f"Programme visit missing in MAST result: {visit}")
        visits[visit] = rows
    return dict(queried_utc=utcnow(), query=selection, programme_rows=len(table), visits=visits,
                scope="Untimed entries at calibration level -1 describe planned footprints.")


def exposure_coverage(cache, sky, day):
    data, receipt = archive_query(cache, "exposures", dict(
        ra=sky.ra.deg, dec=sky.dec.deg, arcmin=.013, output="csv"), f"module_b_001_coverage_{day}.csv")
    native = data[(data.extension == "rate") & (data.expstart >= Time("2025-10-09").mjd)].copy()
    native = native[native.footprint.map(lambda p: contains(p, sky))]
    native["visit"] = native.file.str.slice(0, 13)
    matching = []
    for (visit, detector, filt), rows in native.sort_values("expstart").groupby(["visit", "detector", "filter"]):
        matching.append(dict(visit=visit, detector=detector, filter=filt,
                             earliest_start_utc=Time(rows.expstart.min(), format="mjd").isot,
                             exposure_files=rows.file.tolist()))
    if not any(row["visit"] == BRIGHT_VISIT for row in matching):
        raise ValueError("The archive query failed to recover the known bright visit")
    return dict(query=receipt, matching_groups=matching)


def earlier_photometry(cache, sky, coverage):
    centre = [sky.ra.deg, sky.dec.deg]
    measurements, products, queries = [], [], []
    for band in BANDS:
        filt = band + "-CLEAR"
        covered = [r for r in coverage["matching_groups"] if r["filter"] == filt and r["visit"] in EARLIER]
        if len(covered) != 2:
            raise ValueError(f"Expected two earlier geometrically covered visits in {band}")
        table, query = archive_query(cache, "assoc_mosaic", dict(
            coords=f"{centre[0]:.8f},{centre[1]:.8f},.02", filters=filt, output="csv"),
            f"module_b_001_{band}_intermediate_mosaics.csv")
        queries.append(query)
        for ep in covered:
            visit, detector = ep["visit"], ep["detector"]
            (row,) = table.query("visit == @visit and detector == @detector and tile == 672 and status == 2").itertuples()
            product = dict(band=band, visit=visit, detector=detector, dateobs=row.dateobs,
                           mjd=float(row.expstart), assoc_name=row.assoc_name, exposure_count=int(row.nexp))
            for kind in ("sci", "wht"):
                url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
                product[kind] = download(cache, url, url.rsplit("/", 1)[1], "observed_brightening_intermediate_epoch")
            products.append(product)
            _, obs = vet_mosaic(product, centre)
            measurements.append(obs)
            print(band, visit, detector, obs["apertures"][-1]["flux_njy"], flush=True)
    return dict(epoch_products=products, association_queries=queries, measurements=measurements)


def main():
    path = OUT/"module_b_001_six_filter_change.json"
    lead = json.loads(path.read_text())
    sky = SkyCoord(*lead["icrs_deg"], unit="deg")
    day = datetime.now(timezone.utc).strftime("%Y%m%d")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    mast = mast_followup(sky)
    coverage = exposure_coverage(cache, sky, day)
    prior = earlier_photometry(cache, sky, coverage)
    write_json(OUT/"module_b_001_followup_coverage.json", dict(
        executed_utc=utcnow(), candidate=lead["candidate"], icrs_deg=lead["icrs_deg"],
        mast=mast, dawn=coverage, earlier_visits=prior,
        six_filter_sha256=text_digest(path), script_sha256=text_digest(__file__),
    ))


if __name__ == "__main__":
    main()

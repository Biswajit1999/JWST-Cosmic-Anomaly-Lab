"""Audit the public Einstein Probe archive at the JADES 103607 transient position."""

import json
from urllib.parse import urlencode
from xml.etree import ElementTree

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_einstein_probe_coverage"
POSITION = OUT / "module_b_001_native_november_fade.json"
TAP = "https://emds.esac.esa.int/service/tap/sync"
START, STOP = "2026-10-01", "2026-11-27"
CENTER_RADIUS = 3.0
STORAGE_MIB = 34_200
OBSCORE = "einsteinprobe.obscore"


def query(cache, label, adql):
    url = TAP + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable", QUERY=adql))
    path = cache.fetch(url, f"einstein_probe/{label}.vot", survey="Einstein Probe",
                       version="ESA EPSA public", kind="observation_metadata")
    statuses = [
        element.attrib["value"] for element in ElementTree.parse(path).iter()
        if element.attrib.get("name") == "QUERY_STATUS"
    ]
    if statuses != ["OK"]:
        raise ValueError(f"EPSA {label} TAP status: {statuses}")
    return Table.read(path, format="votable"), receipt(path, url)


def to_utc(value):
    return Time(value, format="mjd", scale="utc").isot


def main():
    position = json.loads(POSITION.read_text(encoding="utf-8"))["icrs_deg"]
    ra, dec = position
    start, stop = Time([START, STOP], scale="utc").mjd
    window = f"t_max >= {start} AND t_min < {stop}"
    footprint = f"CONTAINS(POINT('ICRS',{ra},{dec}),s_region)=1"
    center = f"CONTAINS(POINT('ICRS',s_ra,s_dec),CIRCLE('ICRS',{ra},{dec},{CENTER_RADIUS}))=1"
    cache = DownloadCache(ROOT, STORAGE_MIB * 1024**2)
    results, products = {}, {}
    for label, adql in [
        ("obscore_epoch_range", "SELECT instrument_name, COUNT(*) AS n, MIN(t_min) AS lo, "
         f"MAX(t_max) AS hi FROM {OBSCORE} GROUP BY instrument_name"),
        ("wxt_product_size", "SELECT COUNT(*) AS n FROM einsteinprobe.wxt_product"),
        ("geometry", f"SELECT obs_id,instrument_name,t_min,t_max,s_region FROM {OBSCORE} "
         f"WHERE {footprint}"),
        ("centers", f"SELECT obs_id,instrument_name,t_min,t_max,s_ra,s_dec,s_fov,s_region,"
         f"t_exptime FROM {OBSCORE} WHERE {center}"),
        ("event_interval", f"SELECT obs_id,instrument_name,t_min,t_max,s_ra,s_dec,s_region "
         f"FROM {OBSCORE} WHERE {window} AND ({footprint} OR {center})"),
    ]:
        results[label], products[label] = query(cache, label, adql)
    centers = results["centers"]
    distances = SkyCoord(*position, unit="deg").separation(
        SkyCoord(centers["s_ra"], centers["s_dec"], unit="deg")
    ).arcsec
    nearby = []
    for obs_id in np.unique(centers["obs_id"]):
        indices = np.flatnonzero(centers["obs_id"] == obs_id)
        index = indices[distances[indices].argmin()]
        row = centers[index]
        nearby.append(dict(
            observation_id=str(obs_id), instrument=str(row["instrument_name"]),
            center_icrs_deg=[float(row["s_ra"]), float(row["s_dec"])],
            center_distance_arcsec=float(distances[index]),
            start_utc=to_utc(float(row["t_min"])), stop_utc=to_utc(float(row["t_max"])),
            footprint_stcs=None if np.ma.is_masked(row["s_region"]) else str(row["s_region"]),
            metadata_product_rows=len(indices),
        ))
    nearby.sort(key=lambda row: row["center_distance_arcsec"])
    archive_date_ranges = [
        dict(instrument=str(row["instrument_name"]), product_rows=int(row["n"]),
             first_utc=to_utc(float(row["lo"])), last_utc=to_utc(float(row["hi"])))
        for row in results["obscore_epoch_range"]
    ]
    hits = {
        label: sorted(set(map(str, results[label]["obs_id"])))
        for label in ("geometry", "event_interval")
    }
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        position_sha256=text_digest(POSITION), icrs_deg=position,
        event_interval_utc_inclusive_start_exclusive_end=[START, STOP],
        raw_storage_cap_mib=STORAGE_MIB, archive_queries=products,
        public_obscore_date_ranges=archive_date_ranges,
        public_wxt_product_table_rows=int(results["wxt_product_size"]["n"][0]),
        footprint_match_observation_ids=hits["geometry"],
        event_interval_footprint_or_center_observation_ids=hits["event_interval"],
        center_search_radius_deg=CENTER_RADIUS, nearest_public_observation_centers=nearby,
        limitations=[
            "The public ESA archive inventory does not test unpublished WXT survey exposure.",
            "An observation center can lie in the cone while the published footprint misses.",
            "The metadata search does not measure an X-ray flux or sensitivity.",
        ],
    ))


if __name__ == "__main__":
    main()

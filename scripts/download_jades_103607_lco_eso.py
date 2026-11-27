"""Inventory ESO and Las Cumbres observations and fetch public LCO FITS images."""

import csv
import json
from urllib.parse import urlencode, urlsplit

from astropy.coordinates import SkyCoord
from astropy.time import Time
from check_jades_103607_independent_2026 import END, SOURCE, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from requests import RequestException

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_lco_eso_inventory"
LCO = "https://archive-api.lco.global/frames/"
ESO = "https://archive.eso.org/tap_obs/sync"
BUDGET = 32_750 * 1024**2
ESO_MAX_ROWS = 2500
LCO_MAX_ROWS = 100


def lco_listing(cache, coord):
    query = dict(
        covers=f"POINT({coord.ra.deg:.12f} {coord.dec.deg:.12f})",
        limit=LCO_MAX_ROWS, reduction_level=91, configuration_type="EXPOSE", end=END,
    )
    url = LCO + "?" + urlencode(query)
    path = cache.fetch(
        url, f"metadata/{TAG}_lco.json", survey="Las Cumbres Observatory",
        version="public Archive API footprint query", kind="archive_metadata",
    )
    page = json.loads(path.read_text(encoding="utf-8"))
    rows = page["results"]
    # The archive's count field can differ from the number of returned rows.
    if page.get("next") or len(rows) >= LCO_MAX_ROWS:
        raise ValueError("LCO footprint query requires additional pages")
    if not rows or len({row["id"] for row in rows}) != len(rows):
        raise ValueError("LCO query returned no frames or duplicate frame IDs")
    return rows, dict(
        receipt=receipt(path, url), returned_rows=len(rows),
        archive_reported_count=page.get("count"), next_page_present=bool(page.get("next")),
    )


def fetch_asset(cache, row):
    keys = (
        "id", "basename", "observation_date", "public_date", "target_name", "proposal_id",
        "instrument_id", "site_id", "telescope_id", "exposure_time", "primary_optical_element",
        "area", "configuration_type", "observation_id", "reduction_level",
    )
    record = {key: row[key] for key in keys}
    stable = f"{LCO}{row['id']}/"
    version = max(row["version_set"], key=lambda item: item["id"])
    checksum = "md5:" + version["md5"]
    name = row["basename"] + version["extension"]
    destination = f"cross_detector/lco/{name}"
    local = cache.raw / destination
    signed_url = version["url"]
    try:
        if not local.is_file():
            # Resolve a fresh public URL on first acquisition, since list snapshots are cached.
            with cache.session.get(stable, timeout=(20, 60), headers={"Cache-Control": "no-cache"}) as result:
                result.raise_for_status()
                detail = result.json()
            if detail["id"] != row["id"]:
                raise ValueError("Public frame API returned another frame")
            resolved = next(item for item in detail["version_set"] if item["id"] == version["id"])
            if resolved["md5"] != version["md5"]:
                raise ValueError("Published checksum changed for the same LCO version")
            signed_url = resolved["url"]
        host = urlsplit(signed_url)
        if host.scheme != "https" or host.hostname != "archive-lco-global.s3.us-west-2.amazonaws.com":
            raise ValueError("Unexpected LCO public object store")
        path = cache.fetch(
            stable, destination, survey="Las Cumbres Observatory", version=f"BANZAI E91 {version['id']}",
            kind="observed_optical_image", checksum=checksum, signed_download_url=signed_url,
        )
    except RequestException as error:
        response = getattr(error, "response", None)
        status = response.status_code if response is not None else None
        return dict(**record, status="archive request failed", http_status=status, frame_url=stable)
    print(f"LCO frame {row['id']} {row['observation_date']} {path.stat().st_size / 1024**2:.1f} MiB", flush=True)
    return dict(
        **record, status="downloaded", archive=receipt(path, stable),
        version_id=version["id"], published_md5=version["md5"],
    )


def query_eso(cache, label, table, fields, conditions, order):
    query = f"SELECT {fields} FROM {table} WHERE {conditions} ORDER BY {order}"
    url = ESO + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=ESO_MAX_ROWS, QUERY=query))
    path = cache.fetch(
        url, f"metadata/{TAG}_eso_{label}.csv", survey="ESO Science Archive Facility",
        version="ESO TAP positional and date query", kind="archive_metadata",
    )
    with path.open(encoding="utf-8-sig", newline="") as stream:
        table = csv.DictReader(stream)
        if not set(fields.split(",")) <= set(table.fieldnames or []):
            raise ValueError(f"Unexpected ESO TAP return fields: {label}")
        rows = [{key: value or None for key, value in row.items()} for row in table]
    if len(rows) >= ESO_MAX_ROWS:
        raise ValueError("ESO query reached the record limit")
    return dict(receipt=receipt(path, url), returned_rows=len(rows), rows=rows)


def eso_inventory(cache, coord):
    start, end = Time([START, END]).mjd
    ra, dec = coord.ra.deg, coord.dec.deg
    point = f"POINT('ICRS',{ra:.12f},{dec:.12f})"
    date_raw = f"mjd_obs >= {start} AND mjd_obs < {end} AND dp_cat='SCIENCE'"
    match = f"1=INTERSECTS(s_region,{point})"
    raw_fields = (
        "dp_id,instrument,dp_tech,dp_type,ra,dec,mjd_obs,date_obs,release_date,"
        "s_region,exposure,filter_path,telescope,prog_id"
    )
    # Roughly one degree to each side, enlarged in RA to account for declination.
    rectangle = (
        f"ra BETWEEN {ra - 1.2:.9f} AND {ra + 1.2:.9f} "
        f"AND dec BETWEEN {dec - 1.05:.9f} AND {dec + 1.05:.9f}"
    )
    return dict(
        raw_within_center_rectangle=query_eso(
            cache, "raw_centers", "dbo.raw", raw_fields, f"{date_raw} AND {rectangle}", "mjd_obs,dp_id",
        ),
        raw_regions_intersecting_target=query_eso(
            cache, "raw_intersects", "dbo.raw", raw_fields, f"{date_raw} AND {match}", "mjd_obs,dp_id",
        ),
        reduced_regions_intersecting_target=query_eso(
            cache, "phase3_intersects", "ivoa.ObsCore",
            "dp_id,s_ra,s_dec,instrument_name,dataproduct_type,obs_collection,t_min,t_max,s_region",
            f"t_min < {end} AND t_max >= {start} AND {match}", "t_min,dp_id",
        ),
        center_rectangle_icrs_deg=[ra - 1.2, ra + 1.2, dec - 1.05, dec + 1.05],
        scope="The rectangle uses ESO raw pointing centers and does not prove target coverage. "
        "Published raw regions can consist only of a pointing center; the separate "
        "target-intersection results are for the geometry actually stored in the archive.",
    )


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    coord = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, BUDGET)
    eso = eso_inventory(cache, coord)
    rows, listing = lco_listing(cache, coord)
    frames = [fetch_asset(cache, row) for row in sorted(rows, key=lambda row: row["observation_date"])]
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        position_source_sha256=text_digest(SOURCE), position_icrs_deg=position,
        jwst_event_window_utc=[START, END], lco_listing=listing, lco_frames=frames,
        eso=eso, raw_cache_budget_bytes=BUDGET,
        scope="LCO query selects publicly returned final reduced EXPOSE frames whose archived "
        "geometry covers the target. The exposure list alone supplies no photometric limit. "
        "Only stable frame URLs and published checksums are exposed in the results.",
    ))


if __name__ == "__main__":
    main()

"""Inventory ALMA footprints and download selected public continuum FITS cutouts."""

import json
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode, urlsplit

import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from check_jades_103607_independent_2026 import END, SOURCE, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_alma_inventory"
ESO = "https://almascience.eso.org"
NRAO = "https://almascience.nrao.edu"
NAOJ = "https://almascience.nao.ac.jp"
BUDGET = 33_250 * 1024**2
RADIUS_DEG = 0.005555556
MAX_PB_BYTES = 12 * 1024**2
MAX_ROWS = 2000


def tap_query(cache, label, where):
    fields = (
        "obs_id,member_ous_uid,proposal_id,target_name,t_min,t_max,s_ra,s_dec,"
        "band_list,is_mosaic,s_resolution,em_min,em_max,s_region,data_rights,obs_release_date,qa2_passed"
    )
    query = f"SELECT {fields} FROM ivoa.obscore WHERE {where} ORDER BY t_min,member_ous_uid,obs_id"
    url = ESO + "/tap/sync?" + urlencode(dict(
        REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=MAX_ROWS, QUERY=query,
    ))
    path = cache.fetch(
        url, f"metadata/{TAG}_{label}.csv", survey="ALMA Science Archive",
        version="ObsCore local region search", kind="archive_coverage_metadata",
    )
    frame = pd.read_csv(path)
    if set(fields.split(",")) != set(frame.columns) or len(frame) >= MAX_ROWS:
        raise ValueError(f"Unexpected or truncated ALMA TAP result: {label}")
    return frame, receipt(path, url)


def fits_urls(row, links):
    kinds = {"science": ".mfs.I.pbcor.fits", "primary_beam": ".mfs.I.pb.fits.gz"}
    assets = {}
    for kind, suffix in kinds.items():
        matches = [entry for entry in links if str(entry["access_url"]).endswith(suffix)]
        if len(matches) > 1:
            raise ValueError(f"Multiple {kind} continuum image matches for {row.obs_id}")
        if matches:
            entry = matches[0]
            assets[kind] = dict(
                filename=str(entry["access_url"]).rsplit("/", 1)[-1],
                listed_url=str(entry["access_url"]), listed_bytes=int(entry["content_length"]),
            )
    return assets


def requests_for_links(frame):
    tasks = []
    for row in frame.itertuples(index=False):
        url = ESO + "/datalink/sync?" + urlencode(dict(ID=row.obs_id))
        filename = f"metadata/alma/{row.obs_id.removeprefix('uid://').replace('/', '_')}.vot"
        tasks.append((row.obs_id, url, filename, "ALMA DataLink", "continuum_product_links"))
    return tasks


def prepared(cache, task):
    _, url, filename, _, _ = task
    if (cache.raw / filename).is_file():
        return task, None, None
    try:
        remote = requests.get(url, stream=True, timeout=(12, 34))
        if not remote.ok:
            status = remote.status_code
            remote.close()
            return task, None, str(status)
    except requests.RequestException as exc:
        return task, None, type(exc).__name__
    return task, remote, None


def bounded_downloads(cache, tasks):
    # Only the main thread commits to DownloadCache. Each wave opens up to four public GETs.
    done = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for start in range(0, len(tasks), 4):
            batch = list(pool.map(lambda spec: prepared(cache, spec), tasks[start:start + 4]))
            try:
                for (key, url, filename, version, kind), response, error in batch:
                    if error is not None:
                        done[key] = dict(status="archive request failed", url=url, error=error)
                        continue
                    path = cache.fetch(
                        url, filename, survey="ALMA Science Archive",
                        version=version, kind=kind, response=response,
                    )
                    done[key] = dict(status="downloaded", **receipt(path, url))
            finally:
                for _, response, _ in batch:
                    if response is not None:
                        response.close()
    return done


def science_request(record, coord, server):
    asset = record["assets"]["science"]
    parameters = urlencode(dict(
        ID=asset["filename"], POS=f"CIRCLE {coord.ra.deg:.12f} {coord.dec.deg:.12f} {RADIUS_DEG}",
    ))
    url = server + "/soda/sync?" + parameters
    rel = f"cross_detector/alma/stamps/{asset['filename']}.{urlsplit(server).hostname}.fits"
    return (record["obs_id"], url, rel, "ALMA SODA continuum FITS cutout", "observed_continuum_image")


def products(cache, data, coord):
    pages = bounded_downloads(cache, requests_for_links(data))
    rows = []
    missing = pd.isna(data)
    dictionaries = json.loads(data.drop(columns="s_region").to_json(orient="records"))
    for row, values in zip(data.itertuples(index=False), dictionaries, strict=True):
        page = pages[row.obs_id]
        assets = {}
        if page["status"] == "downloaded":
            table = Table.read(ROOT / page["path"], format="votable")
            if row.obs_id not in set(table["ID"]):
                raise ValueError(f"DataLink identifier mismatch for {row.obs_id}")
            assets = fits_urls(row, table)
        rows.append(dict(**values, datalink=page, assets=assets))
    if not rows or missing["s_region"].any():
        raise ValueError("Incomplete public ALMA footprint inventory")
    ready = [row for row in rows if "science" in row["assets"]]
    science = bounded_downloads(cache, [science_request(row, coord, ESO) for row in ready])
    fallback = [row for row in ready if science[row["obs_id"]]["status"] != "downloaded"]
    if fallback:
        alternate = bounded_downloads(cache, [science_request(row, coord, NAOJ) for row in fallback])
        for row in fallback:
            alternate[row["obs_id"]]["eso_request"] = science[row["obs_id"]]
        science.update(alternate)
    pb_tasks = []
    for row in ready:
        asset = row["assets"].get("primary_beam")
        if asset is None or asset["listed_bytes"] > MAX_PB_BYTES:
            continue
        name = asset["filename"]
        url = NRAO + "/dataPortal/" + name
        pb_tasks.append((row["obs_id"], url, f"cross_detector/alma/pb/{name}",
                         "ALMA continuum primary beam", "primary_beam_map"))
    beams = bounded_downloads(cache, pb_tasks)
    for row in rows:
        obs_id = row["obs_id"]
        row["science_cutout"] = science.get(obs_id)
        row["primary_beam_image"] = beams.get(obs_id)
    return rows


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, BUDGET)
    point = f"POINT('ICRS',{sky.ra.deg:.12f},{sky.dec.deg:.12f})"
    data, coverage = tap_query(cache, "target_position", f"1=INTERSECTS(s_region,{point})")
    start, end = Time([START, END]).mjd
    date = f"t_min < {end} AND t_max >= {start}"
    nearby, local_search = tap_query(
        cache, "event_pointing_centers", f"{date} AND "
        f"s_ra BETWEEN {sky.ra.deg - 1.2:.9f} AND {sky.ra.deg + 1.2:.9f} AND "
        f"s_dec BETWEEN {sky.dec.deg - 1:.9f} AND {sky.dec.deg + 1:.9f}",
    )
    centers = SkyCoord(nearby["s_ra"], nearby["s_dec"], unit="deg")
    nearby["center_separation_arcsec"] = centers.separation(sky).arcsec
    event = data[(data["t_min"] < end) & (data["t_max"] >= start)]
    usable = data[(data["data_rights"] == "Public") & (data["qa2_passed"] == "T")]
    print(f"ALMA {len(data)} target-footprint rows, {len(event)} during the event, "
          f"{usable['member_ous_uid'].nunique()} public QA2-passed historical datasets.", flush=True)
    examined = products(cache, usable, sky)
    print(f"ALMA downloaded {sum((r['science_cutout'] or {}).get('status') == 'downloaded' for r in examined)} "
          "continuum cutouts.", flush=True)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__), position_source_sha256=text_digest(SOURCE),
        position_icrs_deg=position, archive_position_query=coverage, archive_event_center_query=local_search,
        event_window_utc=[START, END], covering_spectral_windows=len(data),
        event_footprint_rows=json.loads(event.drop(columns="s_region").to_json(orient="records")),
        event_center_rows=json.loads(nearby.drop(columns="s_region").to_json(orient="records")),
        public_images=examined, cutout_query_radius_arcsec=RADIUS_DEG * 3600,
        max_primary_beam_product_bytes=MAX_PB_BYTES, raw_cache_budget_bytes=BUDGET,
        scope="Archive footprints and pointing centers describe observation metadata. The completeness "
        "of the pixel assessment is recorded per spectral-window continuum cutout. Products are "
        "not concatenated with image cubes. Primary-beam files above the per-product allocation "
        "are omitted and qualified limits in those windows are measured directly from the delivered "
        "primary-beam-corrected images.",
    ))


if __name__ == "__main__":
    main()

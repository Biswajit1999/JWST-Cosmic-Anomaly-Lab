"""Cache archival Swift XRT photon-counting products for the JADES 103607 sightline."""

import json
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from functools import partial

import requests
from check_jades_103607_swift import TAG as COVERAGE_TAG
from check_jades_103607_swift import product, product_location
from download_jwst_cross_detector import OUT, ROOT

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_swift_xrt_products"
PLAN = {
    "events_w2": ("xrt/event", "xpcw2po_cl.evt.gz"),
    "events_w3": ("xrt/event", "xpcw3po_cl.evt.gz"),
    "exposure_map": ("xrt/products", "xpc_ex.img.gz"),
}
BATCH = 12
WORKERS = 6


def open_get(url, headers):
    for attempt in range(2):
        try:
            # Receive the small FITS product before the serialized cache write.
            # Idle streaming connections closed early on some archive transfers.
            return requests.get(url, headers=headers, timeout=(15, 45)), None
        except requests.RequestException as exc:
            if attempt == 1:
                return None, str(exc)


def acquire_batch(cache, executor, rows, existing, recent_missing):
    jobs = []
    for row in rows:
        for label, parts in PLAN.items():
            url, relative = product_location(row, *parts)
            cached = (url, f"data/raw/{relative}") in existing and (cache.raw / relative).is_file()
            if not cached and url in recent_missing:
                row["products"][label] = recent_missing[url]
                continue
            jobs.append((row, label, parts, url, cached))
    pending_urls = [item[3] for item in jobs if not item[4]]
    get = partial(open_get, headers={"User-Agent": cache.session.headers["User-Agent"]})
    with ExitStack() as responses:
        opened = {}
        for url, result in zip(pending_urls, executor.map(get, pending_urls), strict=True):
            response, error = result
            if response is not None:
                responses.enter_context(response)
            opened[url] = (response, error)
        for row, label, parts, url, cached in jobs:
            response, error = (None, None) if cached else opened[url]
            if error is not None:
                row["products"][label] = dict(status="request failure", url=url, detail=error)
                continue
            try:
                _, source = product(cache, row, *parts, response=response)
            except requests.RequestException as exc:
                status = exc.response.status_code if exc.response is not None else None
                row["products"][label] = dict(
                    status="archive GET failed", http_status=status, url=url, detail=str(exc),
                    checked_utc=utcnow(),
                )
            else:
                row["products"][label] = dict(status="available", **source)


def main():
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    path = OUT / f"{COVERAGE_TAG}.json"
    coverage = json.loads(path.read_text(encoding="utf-8"))
    rows = [dict(**row, products={}) for row in coverage["history_rows"]]
    if len({row["obsid"] for row in rows}) != len(rows):
        raise ValueError("Unexpected duplicate Swift observation IDs")
    candidates = [row for row in rows if (row["xrt_exposure"] or 0) > 0]
    existing = {(record["url"], record["path"]) for record in cache.records}
    destination = OUT / f"{TAG}.json"
    recent_missing = {}
    if destination.exists():
        previous = json.loads(destination.read_text(encoding="utf-8"))
        today = utcnow()[:10]
        for row in previous["observation_rows"]:
            for value in row["products"].values():
                date = value.get("checked_utc", previous["executed_utc"])
                if value.get("http_status") == 404 and date[:10] == today:
                    recent_missing[value["url"]] = value | {"checked_utc": date}
    # Network GETs overlap, while the manifest and storage budget still have one writer.
    with ThreadPoolExecutor(max_workers=WORKERS) as executor:
        for offset in range(0, len(candidates), BATCH):
            batch = candidates[offset:offset+BATCH]
            acquire_batch(cache, executor, batch, existing, recent_missing)
            complete = sum(
                row["products"]["exposure_map"]["status"] == "available"
                and any(row["products"][label]["status"] == "available" for label in ("events_w2", "events_w3"))
                for row in batch
            )
            print(f"{offset+len(batch)}/{len(candidates)} Swift candidates checked, "
                  f"{complete}/{len(batch)} complete PC pairs in this batch", flush=True)
    write_json(destination, dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        coverage_result_sha256=text_digest(path), position_icrs_deg=coverage["position_icrs_deg"],
        candidates_with_nonzero_xrt_master_exposure=len(candidates), observation_rows=rows,
        scope="Public Swift photon-counting cleaned event files and matching pipeline exposure maps. "
        "Both standard w2 and w3 event products were checked for every positive nominal XRT exposure "
        "in the 0.5-degree master-log query. Same-day HTTP 404 results can be reused. "
        "Availability alone does not establish coverage at the target. No Windowed Timing or BAT products.",
    ))


if __name__ == "__main__":
    main()

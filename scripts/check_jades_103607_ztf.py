"""Check public ZTF image coverage and archive depth metadata at JADES 103607."""

import argparse
import csv
import json
from collections import Counter
from datetime import datetime
from urllib.parse import urlencode

from astropy.time import Time
from check_jades_103607_independent_2026 import END, SOURCE, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_ztf_coverage"
IBE = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci"
FINK = "https://api.ztf.fink-portal.org/api/v1/conesearch"
COLUMNS = ("pid", "expid", "obsjd", "obsdate", "ipac_pub_date", "filtercode", "maglimit", "infobits")


def query(cache, position, tag, *, window=None, surrounding=False, refresh=False):
    params = dict(
        POS=",".join(map(str, position)), COLUMNS=",".join(COLUMNS),
        INTERSECT="OVERLAPS" if surrounding else "CENTER", CT="csv",
    )
    if surrounding:
        params["SIZE"] = "0.25"
    if window:
        lo, hi = Time(window, scale="utc").jd
        params["WHERE"] = f"obsjd >= {lo:.1f} AND obsjd < {hi:.1f}"
    url = f"{IBE}?{urlencode(params)}"
    path = cache.fetch(
        url, f"metadata/{TAG}_{tag}.csv",
        survey="NASA IRSA Zwicky Transient Facility", version="IBE calibrated science quadrants",
        kind="observed_detector_footprint_metadata", refresh=refresh,
    )
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not set(COLUMNS).issubset(reader.fieldnames):
            raise ValueError("Unexpected IRSA ZTF IBE metadata columns")
        rows = sorted(reader, key=lambda row: float(row["obsjd"]))
    if len({row["pid"] for row in rows}) != len(rows):
        raise ValueError(f"IRSA returned duplicate processed quadrants in {tag}")
    counts = {row["expid"] for row in rows}
    source = receipt(path, url)
    accessed = next(
        record["accessed_utc"] for record in reversed(cache.records)
        if record["path"] == source["path"] and record["url"] == source["url"]
    )
    return rows, dict(
        **source, accessed_utc=accessed, matching_quadrants=len(rows),
        unique_exposures=len(counts), size_deg=0.25 if surrounding else None, utc_window=window,
    )


def summary(rows):
    good = [row for row in rows if row["infobits"] == "0" and 0 < float(row["maglimit"]) < 30]
    dates = sorted({row["obsdate"][:10] for row in rows})
    by_band = {}
    for band in sorted({row["filtercode"] for row in good}):
        limits = [float(row["maglimit"]) for row in good if row["filtercode"] == band]
        by_band[band] = dict(
            deepest_image_maglimit=max(limits), shallowest_image_maglimit=min(limits),
            zero_infobits_quadrants=len(limits),
        )
    return dict(
        utc_days=dates, quadrant_filters=dict(Counter(row["filtercode"] for row in rows)),
        zero_infobits_exposure_limit_metadata=by_band,
    )


def fink(cache, position, radius, refresh):
    params = {"ra": position[0], "dec": position[1], "radius": radius, "output-format": "json"}
    if radius > 5:
        params.update(
            startdate=f"{START} 00:00:00",
            window=(datetime.fromisoformat(END) - datetime.fromisoformat(START)).days,
        )
    url = f"{FINK}?{urlencode(params)}"
    path = cache.fetch(
        url, f"metadata/{TAG}_fink_{radius}.json",
        survey="Fink ZTF broker", version="Cone search",
        kind="observed_alert_metadata", refresh=refresh,
    )
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list) or any("i:objectId" not in row for row in data):
        raise ValueError("Unexpected Fink ZTF cone-search response")
    return dict(
        **receipt(path, url), radius_arcsec=radius, query=params,
        object_ids=sorted({row["i:objectId"] for row in data}),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="Fetch current public archive results.")
    args = parser.parse_args()
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    all_rows, history = query(cache, position, "all", refresh=args.refresh)
    recent, event = query(cache, position, "2026_center", window=(START, END), refresh=args.refresh)
    _, nearby = query(cache, position, "2026_box", window=(START, END), surrounding=True, refresh=args.refresh)
    _, euclid = query(cache, position, "2024_euclid", window=("2024-08-05", "2024-08-08"), refresh=args.refresh)
    alerts = {label: fink(cache, position, radius, args.refresh)
              for radius, label in ((5, "target_full_history"), (120, "surrounding_2026_first_detection"))}
    lo, hi = Time([START, END], scale="utc").jd
    if {row["pid"] for row in recent} != {row["pid"] for row in all_rows if lo <= float(row["obsjd"]) < hi}:
        raise ValueError("ZTF targeted and complete history date cuts disagree")
    published = [
        row for row in all_rows if row["ipac_pub_date"] not in ("", "null")
        and datetime.fromisoformat(row["ipac_pub_date"]) <= datetime.fromisoformat(history["accessed_utc"])
    ]
    report = dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        position_icrs_deg=position, position_source_sha256=text_digest(SOURCE),
        complete_center_query=history, recent_center_query=event, recent_overlapping_box_query=nearby,
        around_euclid_h_excess_query=euclid, released_matching_quadrants=len(published),
        public_alert_searches=alerts,
        most_recent_returned_quadrant=all_rows[-1] if all_rows else None,
        available_history=summary(published),
        scope="Anonymous IBE calibrated-image metadata for geometrically covering quadrants. "
        "Release dates indicate catalog metadata availability as of the query. Exposure-level "
        "limiting magnitudes do not measure flux or sensitivity at the target position. "
        "Recent public metadata gaps cannot exclude observations not yet visible in the archive. "
        "The 120-arcsec broker query restricts first detections and the 5-arcsec query spans "
        "the broker history. The alerts do not give target exposure depth.",
    )
    print(f"{len(all_rows)} position matches, {len(published)} listed as released. "
          f"Event window {len(recent)}, expanded footprint window {nearby['matching_quadrants']}."
          f" Euclid bracket {euclid['matching_quadrants']}.", flush=True)
    write_json(OUT / f"{TAG}.json", report)


if __name__ == "__main__":
    main()

"""Check additional 2025 follow-up visits against the November 2026 release."""

from itertools import product

import pandas as pd
from download_jwst_cross_detector import OUT, ROOT
from search_jwst_november_catalogues import (
    BANDS,
    LATER,
    STORAGE_MIB,
    compare,
    cross_band,
    read_catalogue,
)

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

EARLIER = tuple(f"08060-{i:03d}" for i in range(3, 7))
STEM = "november_goods_s_prior_visits"


def main():
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    data, receipts = {}, []
    for visit, band in product((*EARLIER, *LATER), BANDS):
        data[visit, band], description = read_catalogue(cache, visit, band)
        receipts.append(description)
        print(visit, band, description["usable"], flush=True)
    diagnostics, leads = [], []
    for old, new in product(EARLIER, LATER):
        bands = {}
        for band in BANDS:
            bands[band], summary = compare(data[old, band], data[new, band])
            diagnostics.append(dict(older=old, newer=new, band=band, **summary))
        matches = cross_band(bands)
        leads.extend(dict(older=old, newer=new, **match) for match in matches)
        print(old, new, "cross-channel leads", len(matches), flush=True)
    cols = ("older", "newer", "ra", "dec", "sign", "bands", "rank")
    table = pd.DataFrame(leads) if leads else pd.DataFrame(columns=cols)
    path = OUT/f"{STEM}_matches.csv"
    table.sort_values(["bands", "rank"], ascending=False).to_csv(path, index=False)
    write_json(OUT/f"{STEM}.json", dict(
        executed_utc=utcnow(), archive_products=receipts, comparisons=diagnostics,
        pair_level_candidates=len(table), distinct_followup_visits=len(EARLIER),
        ranking_method_sha256=text_digest(ROOT/"scripts/search_jwst_november_catalogues.py"),
        script_sha256=text_digest(__file__), candidates_sha256=text_digest(path),
        scope="Matched source-catalogue comparisons between November/December 2025 and November 2026.",
    ))


if __name__ == "__main__":
    main()

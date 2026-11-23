"""Search both public GOODS-S pointings recorded by NIRCam module B."""

import json
from collections import Counter

import pandas as pd
from download_jwst_cross_detector import OUT, ROOT, download
from measure_jwst_epoch_candidate import verified_path
from search_jwst_cross_detector import (
    DIRECT_RADIUS,
    MATCH_RADIUS,
    MIN_FRAC_CHANGE,
    MIN_SCORE,
    pairs,
    search_pair,
)

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

STORAGE_MIB = 14100
POINTINGS = {
    "001": {"older": "jw08060001001", "newer": "jw12577001001"},
    "002": {"older": "jw08060002001", "newer": "jw12577002001"},
}
DETECTORS = {
    "F444W": ("NRCBLONG",),
    "F356W": ("NRCBLONG",),
    "F200W": ("NRCB1", "NRCB2", "NRCB3", "NRCB4"),
}
LEDGER = OUT / "module_b_products.json"
SEARCH = OUT / "module_b_search.json"


def retrieve():
    listings_path = OUT / "mosaics.json"
    listings = json.loads(listings_path.read_text())["associations"]
    cache = DownloadCache(ROOT, STORAGE_MIB * 1024**2)
    products = []
    for band, detectors in DETECTORS.items():
        (receipt,) = [r for r in listings if r["path"].endswith(f"/b_{band}_associations.csv")]
        entries = pd.read_csv(verified_path(receipt))
        for pointing, visits in POINTINGS.items():
            for detector in detectors:
                for epoch, visit in visits.items():
                    (row,) = list(entries.query(
                        "visit == @visit and tile == 672 and status == 2 and detector == @detector"
                    ).itertuples())
                    record = dict(pointing=pointing, band=band, detector=detector, epoch=epoch, visit=visit,
                                  detector_center=[float(row.ra), float(row.dec)], dateobs=row.dateobs,
                                  mjd=float(row.expstart), assoc_name=row.assoc_name, exposure_count=int(row.nexp))
                    for kind in ("sci", "wht"):
                        url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
                        record[kind] = download(cache, url, url.rsplit("/", 1)[-1], "observed_module_b_"+kind)
                    products.append(record)
                    print("Retrieved", pointing, visit, band, detector, flush=True)
            write_json(LEDGER, {
                "executed_utc": utcnow(), "images": products, "association_ledger_sha256": text_digest(listings_path),
                "storage_cap_mib": STORAGE_MIB, "script_sha256": text_digest(__file__),
            })
    return products


def run_pointing(pointing, products):
    reference = None
    peaks, diagnostics = [], []
    for band, detectors in DETECTORS.items():
        for detector in detectors:
            images = {r["epoch"]: r for r in products
                      if (r["pointing"], r["band"], r["detector"]) == (pointing, band, detector)}
            results, data, metrics = search_pair(images, band, reference=reference if band == "F200W" else None)
            if band == "F444W":
                reference = {k: data[k] for k in ("covered", "wcs", "origin", "scale")}
            del data
            peaks.append(results)
            diagnostics.append(metrics)
    singles = pd.concat(peaks, ignore_index=True)
    matches = pairs(singles)
    if not matches.empty:
        matches["candidate"] = [f"GS-{pointing}-B-{x.ra:.6f}{x.dec:+.6f}" for x in matches.itertuples()]
        confirms = matches.get("f356w_aperture_gate", pd.Series(False, index=matches.index))
        matches["passes_three_filter_gate"] = matches.aperture_gate & confirms.eq(True)
    summary = {
        "ranked_peaks": len(singles), "cross_channel_pairs": len(matches),
        "pair_signs": dict(Counter(matches.sign)) if not matches.empty else {},
        "passes_pair_aperture_gate": int(matches.aperture_gate.sum()) if not matches.empty else 0,
        "passes_three_filter_gate": int(matches.passes_three_filter_gate.sum()) if not matches.empty else 0,
        "summed_chip_intersection_f200w_f444w_arcmin2": sum(
            r["overlap_with_f444w_arcmin2"] for r in diagnostics if r["band"] == "F200W"),
        "diagnostics": diagnostics,
    }
    singles.insert(0, "pointing", pointing)
    if not matches.empty:
        matches.insert(0, "pointing", pointing)
    print("Completed", pointing, {k: v for k, v in summary.items() if k != "diagnostics"}, flush=True)
    return singles, matches, summary


def main():
    products = retrieve()
    all_peaks, all_matches, summaries = [], [], {}
    for pointing in POINTINGS:
        peaks, matches, summary = run_pointing(pointing, products)
        all_peaks.append(peaks)
        all_matches.append(matches)
        summaries[pointing] = summary
        pd.concat(all_peaks, ignore_index=True).to_csv(OUT / "module_b_peaks.csv", index=False)
        pd.concat(all_matches, ignore_index=True).to_csv(OUT / "module_b_matches.csv", index=False)
        write_json(SEARCH, {
            "executed_utc": utcnow(), "pointings": summaries, "aperture_radius_arcsec": DIRECT_RADIUS,
            "minimum_score": MIN_SCORE, "match_radius_arcsec": MATCH_RADIUS,
            "minimum_fraction_change": MIN_FRAC_CHANGE, "products_sha256": text_digest(LEDGER),
            "search_module_sha256": text_digest(ROOT / "scripts/search_jwst_cross_detector.py"),
            "diagnostic_module_sha256": text_digest(ROOT / "src/cosmic_anomaly_lab/epoch_search.py"),
            "script_sha256": text_digest(__file__),
        })


if __name__ == "__main__":
    main()

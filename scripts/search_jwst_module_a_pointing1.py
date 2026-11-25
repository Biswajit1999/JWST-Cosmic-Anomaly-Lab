"""Search the remaining Module A footprint of the two October GOODS-S epochs."""

import json

import pandas as pd
from check_jwst_module_b_coincidences import counts, shifted_trials
from download_jwst_cross_detector import CONTROL_VISITS, DETECTORS, OUT, ROOT, download
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

STORAGE_MIB = 24_000
PREFIX = "module_a_001"
PRODUCTS = OUT/f"{PREFIX}_products.json"
PEAKS = OUT/f"{PREFIX}_peaks.csv"
MATCHES = OUT/f"{PREFIX}_matches.csv"
RESULTS = OUT/f"{PREFIX}_search.json"


def retrieve(cache):
    source = OUT/"mosaics.json"
    archive = json.loads(source.read_text())["associations"]
    images, metadata = [], []
    for band, detectors in DETECTORS.items():
        (listing,) = [r for r in archive if r["path"].endswith(f"/b_{band}_associations.csv")]
        rows = pd.read_csv(verified_path(listing))
        metadata.append(listing)
        for detector in detectors:
            for era, visit in CONTROL_VISITS.items():
                (row,) = rows.query("visit == @visit and tile == 672 and status == 2 and detector == @detector").itertuples()
                record = dict(pointing="001", band=band, detector=detector, epoch=era, visit=visit,
                              detector_center=[float(row.ra), float(row.dec)], dateobs=row.dateobs,
                              mjd=float(row.expstart), assoc_name=row.assoc_name, exposure_count=int(row.nexp))
                for kind in ("sci", "wht"):
                    url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
                    record[kind] = download(cache, url, url.rsplit("/", 1)[1], f"observed_module_a_001_{kind}")
                images.append(record)
                print("Verified", visit, band, detector, flush=True)
        write_json(PRODUCTS, dict(executed_utc=utcnow(), images=images, association_queries=metadata,
                                  discovery_ledger_sha256=text_digest(source), storage_cap_mib=STORAGE_MIB,
                                  script_sha256=text_digest(__file__)))
    return images


def scan(products):
    reference = None
    peaks, diagnostics = [], []
    for band, detectors in DETECTORS.items():
        for detector in detectors:
            records = {r["epoch"]: r for r in products if (r["band"], r["detector"]) == (band, detector)}
            items, data, fit = search_pair(records, band, reference=reference if band == "F200W" else None)
            if band == "F444W":
                reference = {name: data[name] for name in ("covered", "wcs", "origin", "scale")}
            peaks.append(items)
            diagnostics.append(fit)
            del data
    return pd.concat(peaks, ignore_index=True), diagnostics


def main():
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    products = retrieve(cache)
    peaks, diagnostics = scan(products)
    match = pairs(peaks)
    if not match.empty:
        match["candidate"] = [f"GS-001-A-{r.ra:.6f}{r.dec:+.6f}" for r in match.itertuples()]
        crosscheck = match.get("f356w_aperture_gate", pd.Series(False, index=match.index))
        match["passes_three_filter_gate"] = match.aperture_gate & crosscheck.eq(True)
    peaks.to_csv(PEAKS, index=False)
    match.to_csv(MATCHES, index=False)
    observed = counts(peaks)
    shifted = shifted_trials(peaks)
    searched = sum(fit["overlap_with_f444w_arcmin2"] for fit in diagnostics if fit["band"]=="F200W")
    write_json(RESULTS, dict(
        executed_utc=utcnow(), products_sha256=text_digest(PRODUCTS),
        peaks_sha256=text_digest(PEAKS), matches_sha256=text_digest(MATCHES),
        ranked_peaks=len(peaks), common_area_arcmin2=searched, observed=observed, shifted_trials=shifted,
        pair_diagnostics=diagnostics, fixed_gates=dict(minimum_ranking_score=MIN_SCORE,
            match_radius_arcsec=MATCH_RADIUS, aperture_radius_arcsec=DIRECT_RADIUS,
            minimum_fractional_change=MIN_FRAC_CHANGE),
        search_method_sha256=text_digest(ROOT/"scripts/search_jwst_cross_detector.py"),
        diagnostic_method_sha256=text_digest(ROOT/"src/cosmic_anomaly_lab/epoch_search.py"),
        shifted_control_sha256=text_digest(ROOT/"scripts/check_jwst_module_b_coincidences.py"),
        script_sha256=text_digest(__file__),
        scope="Visits 001 of GO 8060 and 12577. Areas are sums across short-wave detector intersections."
              " Shifted matches check incidental overlap without calibrating a global false-alarm rate.",
    ))
    print("Module A pointing 1", len(peaks), searched, observed, flush=True)


if __name__ == "__main__":
    main()

"""Verify the two Module A pointing 001 leads in observed images and prior records."""

import json
from zipfile import ZipFile

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from audit_jwst_goods_s_prior_records import CATALOGUES, catalogue, proposal
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from measure_jwst_brightening_offset import BANDS, local_alignment, match_pair, patch
from measure_jwst_epoch_candidate import RADII, photometry, pixel_scale, show_stamp, verified_path
from measure_jwst_module_b_brightening import ANNULUS
from search_jwst_cross_detector import check_same_grid
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic
from vet_jwst_module_b import tns_match

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

FILTERS = ("F200W", "F356W", "F444W")
PREFIX = "module_a_001"
PLAN_DATE = "20261123"
STORAGE_MIB = 25_500


def inspect_images(candidate, products):
    center = (candidate.ra, candidate.dec)
    results, views = {}, {}
    for band in FILTERS:
        detector = candidate.sw_detector if band == "F200W" else "NRCALONG"
        records = {p["epoch"]: p for p in products if p["band"] == band and p["detector"] == detector}
        old, before = vet_mosaic(records["older"], center)
        new, after = vet_mosaic(records["newer"], center)
        o, wcs, x, y, x0, y0 = old
        n, nwcs, _, _, nx0, ny0 = new
        scale = pixel_scale(wcs)
        check_same_grid((o, wcs, (x0, y0), scale), (n, nwcs, (nx0, ny0), pixel_scale(nwcs)))
        try:
            registration = local_alignment(o, n, x, y, scale, BANDS.get(band, .24))
        except ValueError as exc:
            if not str(exc).startswith(("Too few bright local reference sources", "Insufficient local core matches")):
                raise
            registration = dict(method="unmatched_common_grid", reason=str(exc))
            o_match, n_match, delta = o, n, n-o
        else:
            registration["method"] = "local_translation_gain_and_core_width"
            o_match, n_match, delta = match_pair(o, n, registration, scale, 0.)
        signed = photometry(delta, x, y, scale, RADII, ANNULUS)
        results[band] = dict(detector=detector, older=before, newer=after,
                             registration=registration, difference_apertures=signed)
        views[band] = (o_match, n_match, delta, x, y, scale)
        print(candidate.candidate, band, signed[-1], flush=True)
    return results, views


def plot(candidate, views, results):
    fig = Figure(figsize=(8.1, 7.9), layout="constrained")
    axes = fig.subplots(3, 3)
    for row, band in zip(axes, FILTERS, strict=True):
        old, new, difference, x, y, scale = views[band]
        matched = results[band]["registration"]["method"] != "unmatched_common_grid"
        half = .67
        flux = np.concatenate([patch(d, x, y, scale, half)[2].ravel() for d in (old, new)])
        low, high = np.nanpercentile(flux, (4, 98.5))
        span = np.nanpercentile(abs(patch(difference, x, y, scale, half)[2]), 98)
        for ax, (data, epoch, cmap, bounds) in zip(row, (
            (old, "Oct 2025", "gray", (low, high)),
            (new, "Oct 2026", "gray", (low, high)),
            (difference, "2026 − 2025", "RdBu_r", (-span, span)),
        ), strict=True):
            show_stamp(ax, data, x, y, scale, bounds, cmap, .24, half_arcsec=half)
            suffix = " · unmatched" if not matched and epoch == "2026 − 2025" else ""
            ax.set_title(f"{band} · {epoch}{suffix}", fontsize=10)
    fig.suptitle(f"GOODS-S · {candidate.ra:.6f}°, {candidate.dec:.6f}°\n"
                 "Observed JWST mosaics · shared astrometric grid", fontsize=12)
    save_figure(fig, FIGURES/f"{PREFIX}_{candidate.sw_detector.lower()}_vetting")


def public_registry():
    prior = json.loads((OUT/"module_b_001_release_20261124.json").read_text())
    receipt = {k: v for k, v in prior["tns"].items()
               if k in ("url", "path", "sha256", "snapshot_label_utc", "catalogue_rows")}
    with ZipFile(verified_path(receipt)) as archive, archive.open("tns_public_objects.csv") as stream:
        snapshot = stream.readline().decode("ascii").strip()
        columns = ["objid", "name_prefix", "name", "ra", "declination", "type", "discoverydate", "time_received"]
        table = pd.read_csv(stream, usecols=columns, encoding="utf-8", encoding_errors="replace")
    if snapshot != receipt["snapshot_label_utc"] or len(table) != receipt["catalogue_rows"]:
        raise ValueError("TNS snapshot differs from the verified archive receipt")
    return receipt, table


def main():
    ledger, matchfile = (OUT/f"{PREFIX}_{name}" for name in ("products.json", "matches.csv"))
    products = json.loads(ledger.read_text())["images"]
    candidates = pd.read_csv(matchfile).query("passes_three_filter_gate")
    if len(candidates) != 2:
        raise ValueError("Unexpected three-filter shortlist size")
    centres = {r.candidate: SkyCoord(r.ra, r.dec, unit="deg") for r in candidates.itertuples()}
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    plans = [proposal(cache, number, centres, PLAN_DATE) for number in (8060, 12577)]
    published = {name: catalogue(cache, name, spec, centres) for name, spec in CATALOGUES.items()}
    tns = public_registry()
    inspected = []
    for row in candidates.itertuples():
        print(row.candidate, [(p["program"], [(c["name"], c["positions"][row.candidate]["matches"])
                                               for c in p["catalogues"] if c["positions"][row.candidate]["matches"]])
                             for p in plans], flush=True)
        results, views = inspect_images(row, products)
        plot(row, views, results)
        inspected.append(dict(candidate=row.candidate, icrs_deg=[row.ra, row.dec], sign=int(row.sign), bands=results,
                              tns=tns_match(tns, centres[row.candidate])))
    write_json(OUT/f"{PREFIX}_vetting.json", dict(
        executed_utc=utcnow(), candidates=inspected, public_plans=plans, published_catalogues=published,
        products_sha256=text_digest(ledger), matches_sha256=text_digest(matchfile), script_sha256=text_digest(__file__),
        scope="Two preselected three-filter leads. Photometry uses mosaics from the same two epochs as the search;"
              " local alignment accounts approximately for nearby zero-point, centroid and core-width differences.",
    ))


if __name__ == "__main__":
    main()

"""Search the November GOODS-S release with archival NIRCam source catalogues."""

import argparse
from itertools import product

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import fetch
from measure_jwst_epoch_candidate import verified_path
from read_jwst_i2d_cutout import S3

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

EARLIER = ("08060-001", "08060-002", "12577-001", "12577-002")
LATER = ("12577-003", "12577-004")
BANDS = ("F150W", "F200W", "F277W", "F356W")
SHORT = ("F150W", "F200W")
LONG = ("F277W", "F356W")
PREFIX = "november_goods_s_catalogue"
OCTOBER_PREFIX = "october_goods_s_catalogue"
OCTOBER_PAIRS = tuple(zip(EARLIER[:2], EARLIER[2:], strict=True))
MATCH_ARCSEC = .13
ISOLATION_ARCSEC = .28
MIN_FRACTION = .18
MIN_SCORE = 3.
STRONG_SCORE = 6.
STORAGE_MIB = 28_000


def product_name(visit, band, suffix):
    program, obs = visit.split("-")
    return f"jw{program}-o{obs}_t{obs}_nircam_clear-{band.lower()}_{suffix}"


def coordinates(frame):
    return SkyCoord(frame.ra.to_numpy(), frame.dec.to_numpy(), unit="deg")


def read_catalogue(cache, visit, band):
    program, obs = visit.split("-")
    name = product_name(visit, band, "cat.ecsv")
    url = f"{S3}/jw{program}/L3/t/o{obs}/{name}"
    receipt = fetch(cache, url, name, "observed_november_catalogue",
                    survey="MAST public S3", version="JWST stage-3 source catalogue")
    q = QTable.read(verified_path(receipt), format="ascii.ecsv")
    sky = q["sky_centroid"].icrs
    f = pd.DataFrame(dict(ra=sky.ra.deg, dec=sky.dec.deg, source_id=np.asarray(q["label"]),
                          flux=q["aper70_flux"].to_value(u.nJy), err=q["aper70_flux_err"].to_value(u.nJy)))
    ok = np.isfinite(f[["ra", "dec", "flux", "err"]]).all(axis=1) & (f.err > 0)
    f = f[ok].reset_index(drop=True)
    if len(f) < 100:
        raise ValueError(f"Too few usable rows in {name}")
    f["isolation_arcsec"] = coordinates(f).match_to_catalog_sky(coordinates(f), nthneighbor=2)[1].arcsec
    return f, dict(product=receipt, visit=visit, band=band, rows=len(q), usable=len(f),
                   calibration_pipeline=q.meta["version"]["jwst"],
                   aperture_pixels=float(q.meta["aperture_params"]["aperture_radii"][-1]))


def pairs(old, new):
    a, b = coordinates(old), coordinates(new)
    nearest, distance, _ = a.match_to_catalog_sky(b)
    reverse = b.match_to_catalog_sky(a)[0]
    use = (distance.arcsec <= MATCH_ARCSEC) & (reverse[nearest] == np.arange(len(old)))
    first, second = old[use].reset_index(drop=True), new.iloc[nearest[use]].reset_index(drop=True)
    table = second[["ra", "dec"]].copy()
    for tag, values in (("old", first), ("new", second)):
        for key in ("flux", "err", "isolation_arcsec", "source_id"):
            table[f"{tag}_{key}"] = values[key]
    table["position_difference_arcsec"] = distance.arcsec[use]
    return table


def robust_mad(x):
    return float(1.4826*np.median(abs(x-np.median(x))))


def compare(before, after):
    all_matches = pairs(before, after)
    table = all_matches[(all_matches.old_isolation_arcsec >= ISOLATION_ARCSEC)
                        & (all_matches.new_isolation_arcsec >= ISOLATION_ARCSEC)].copy()
    stable = table[(table.old_flux >= 200) & (table.old_flux <= 8000) & (table.new_flux > 0)
                   & (table.old_flux/table.old_err >= 20) & (table.new_flux/table.new_err >= 20)]
    stats = dict(mutual_matches=len(all_matches), isolated_pairs=len(table), reference_count=len(stable))
    if len(stable) < 24:
        return table.iloc[:0], dict(**stats, insufficient_overlap=True)
    ratio = np.log(stable.new_flux.to_numpy()/stable.old_flux.to_numpy())
    gain = float(np.exp(np.median(ratio)))
    scatter = robust_mad(ratio)
    delta = table.new_flux-gain*table.old_flux
    brighter = np.maximum(table.new_flux, gain*table.old_flux)
    sigma = np.sqrt(table.new_err**2+(gain*table.old_err)**2+(scatter*np.maximum(brighter, 0))**2)
    # Catalog uncertainties can underestimate drizzled pixel correlation.
    # Calibrate the standardized residual with the central population.
    normalized = delta/sigma
    valid = (table.old_flux > 0) & (table.new_flux > 0)
    width = max(robust_mad(normalized[valid]), 1.)
    table["delta_njy"] = delta
    table["score"] = normalized/width
    table["fraction"] = delta/np.where(brighter > 0, brighter, np.nan)
    table["sign"] = np.sign(delta).astype(int)
    keep = valid & (abs(table.score) >= MIN_SCORE) & (abs(table.fraction) >= MIN_FRACTION)
    return table[keep].reset_index(drop=True), dict(**stats, gain=gain, reference_log_ratio_mad=scatter,
                                                     normalization=width, shortlist=int(keep.sum()))


def cross_band(tables):
    picks = []
    for sign in (-1, 1):
        by_band = {b: t[t.sign == sign].reset_index(drop=True) for b, t in tables.items() if len(t)}
        for primary in SHORT:
            if primary not in by_band:
                continue
            anchors = by_band[primary]
            match = {}
            for band, frame in by_band.items():
                if band != primary and len(frame):
                    match[band] = coordinates(anchors).match_to_catalog_sky(coordinates(frame))[:2]
            for i, anchor in enumerate(anchors.itertuples(index=False)):
                measured = {primary: anchors.iloc[i]}
                for band, (nearest, distance) in match.items():
                    if distance[i].arcsec <= MATCH_ARCSEC:
                        measured[band] = by_band[band].iloc[nearest[i]]
                # Anchor twice only if no changed F200W detection is available.
                if primary == "F150W" and "F200W" in measured:
                    continue
                if not any(b in measured for b in LONG) or max(abs(m.score) for m in measured.values()) < STRONG_SCORE:
                    continue
                row = dict(ra=anchor.ra, dec=anchor.dec, sign=sign, bands=len(measured), anchor_band=primary,
                           rank=min(abs(m.score) for m in measured.values()))
                for b, m in measured.items():
                    for field in ("score", "fraction", "old_flux", "new_flux", "delta_njy",
                                  "old_source_id", "new_source_id", "position_difference_arcsec"):
                        row[f"{b.lower()}_{field}"] = m[field]
                picks.append(row)
    return picks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--october", action="store_true", help="Compare the homologous October visits.")
    args = parser.parse_args()
    chosen = OCTOBER_PAIRS if args.october else tuple(product(EARLIER, LATER))
    stem = OCTOBER_PREFIX if args.october else PREFIX
    visits = tuple(dict.fromkeys(v for pair in chosen for v in pair))
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    archive, receipts = {}, []
    for visit, band in product(visits, BANDS):
        archive[visit, band], receipt = read_catalogue(cache, visit, band)
        receipts.append(receipt)
        print(visit, band, receipt["usable"], flush=True)
    leads, comparisons, excesses = [], [], []
    for old, new in chosen:
        tables = {}
        for band in BANDS:
            tables[band], stats = compare(archive[old, band], archive[new, band])
            comparisons.append(dict(older=old, newer=new, band=band, **stats))
            excesses.append(tables[band].assign(older=old, newer=new, band=band))
        accepted = cross_band(tables)
        print(old, new, len(accepted), sum(r["bands"] >= 3 for r in accepted), flush=True)
        leads.extend(dict(older=old, newer=new, **lead) for lead in accepted)
    matches = OUT/f"{stem}_matches.csv"
    columns = ["older", "newer", "ra", "dec", "sign", "bands", "anchor_band", "rank"]
    collected = pd.DataFrame(leads) if leads else pd.DataFrame(columns=columns)
    collected.sort_values(["bands", "rank"], ascending=False).to_csv(matches, index=False)
    measurements = OUT/f"{stem}_residuals.csv"
    pd.concat(excesses, ignore_index=True).to_csv(measurements, index=False)
    write_json(OUT/f"{stem}.json", dict(
        executed_utc=utcnow(), products=receipts, comparisons=comparisons, pair_level_candidates=len(leads),
        matches_sha256=text_digest(matches), residuals_sha256=text_digest(measurements),
        script_sha256=text_digest(__file__),
        gates=dict(match_arcsec=MATCH_ARCSEC, isolation_arcsec=ISOLATION_ARCSEC, fractional_change=MIN_FRACTION,
                   ranking_score=MIN_SCORE, at_least_one_band_score=STRONG_SCORE),
        scope="Mutual catalogue matches in the listed two-visit comparisons. Both epochs need positive"
              " measured flux. Ranking gates use empirical scatter and require short- and long-wave agreement."
              " Pair-level entries can duplicate sky positions. Image validation is required.",
    ))


if __name__ == "__main__":
    main()

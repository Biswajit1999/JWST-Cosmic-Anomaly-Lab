"""Search a second GOODS-S footprint for changes on distinct NIRCam detectors."""

import json
from collections import Counter

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from download_jwst_cross_detector import DETECTORS, OUT, ROOT
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import aperture, pixel_scale, show_stamp, verified_path

from cosmic_anomaly_lab.epoch_search import aperture_score, fit_registration_residual, local_peaks
from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

FIGURES = ROOT / "figures/jwst_cross_detector"
SEARCH_RADII = {"F200W": .10, "F356W": .14, "F444W": .16}
DIRECT_RADIUS = .24
MIN_SCORE = 4.5
MATCH_RADIUS = .14
MIN_FRAC_CHANGE = .18


def window(record, kind, center, half):
    with fits.open(verified_path(record[kind]), memmap=False) as h:
        header = h[0].header
        wcs = WCS(header).celestial
        scale = pixel_scale(wcs)
        x, y = wcs.world_to_pixel_values(*center)
        radius = int(np.ceil(half/scale))
        x0, y0 = round(float(x))-radius, round(float(y))-radius
        x1, y1 = x0+2*radius+1, y0+2*radius+1
        xa, xb = max(0, x0), min(header["NAXIS1"], x1)
        ya, yb = max(0, y0), min(header["NAXIS2"], y1)
        if xa >= xb or ya >= yb:
            raise ValueError(f"Requested region misses {record['assoc_name']}")
        out = np.zeros((2*radius+1, 2*radius+1), dtype=np.float32)
        sub = out[ya-y0:yb-y0, xa-x0:xb-x0]
        sub[:] = h[0].data[ya:yb, xa:xb]
        unit = float(header["PHOTFNU"])*1e9
        sub *= unit if kind == "sci" else 1/unit**2
    return out, wcs, np.array([x0, y0]), scale


def check_same_grid(a, b):
    arr, wa, (ax, ay), _ = a
    other, wb, (bx, by), _ = b
    if arr.shape != other.shape:
        raise ValueError("Epoch windows have different dimensions")
    coords = np.array([[0, 0], [arr.shape[1]-1, arr.shape[0]-1], [0, arr.shape[0]-1]])
    world = wa.pixel_to_world_values(coords[:, 0]+ax, coords[:, 1]+ay)
    tx, ty = wb.world_to_pixel_values(*world)
    shift = np.max(np.hypot(tx-bx-coords[:, 0], ty-by-coords[:, 1]))
    if shift >= .03:
        raise ValueError(f"Epoch pixel grid displacement {shift:.4f}")


def sky_coverage(reference, shortwave):
    """Map long-wave validity onto the short-wave pixels using the common TAN grid."""
    yy, xx = np.nonzero(shortwave["covered"])
    aw, bw = shortwave["wcs"].wcs, reference["wcs"].wcs
    if list(aw.ctype) != list(bw.ctype) or not np.allclose(aw.crval, bw.crval, rtol=0, atol=1e-8):
        raise ValueError("Different projection centres in coverage calculation")
    transform = np.linalg.solve(reference["wcs"].pixel_scale_matrix, shortwave["wcs"].pixel_scale_matrix)
    xy = np.column_stack((xx, yy))
    fitted = ((xy+shortwave["origin"]+1-aw.crpix) @ transform.T + bw.crpix-1-reference["origin"])
    rx, ry = np.rint(fitted).astype(int).T
    coverage = reference["covered"]
    good = (rx >= 0) & (ry >= 0) & (rx < coverage.shape[1]) & (ry < coverage.shape[0])
    return int(coverage[ry[good], rx[good]].sum())*shortwave["scale"]**2/3600


def search_pair(records, band, reference=None):
    older, newer = (records[epoch] for epoch in ("older", "newer"))
    center = newer["detector_center"]
    half = 45 if band == "F200W" else 87
    entries = [window(r, kind, center, half) for kind in ("sci", "wht") for r in (older, newer)]
    for other in entries[1:]:
        check_same_grid(entries[0], other)
    old, new, wo, wn = (e[0] for e in entries)
    _, wcs, origin, scale = entries[0]
    valid = np.isfinite(old+new+wo+wn) & (wo > 0) & (wn > 0)
    old[~valid] = new[~valid] = 0
    difference, safe, fit = fit_registration_residual(old, new, valid)
    _, score, _ = aperture_score(difference, wo, wn, safe, SEARCH_RADII[band]/scale,
                                  block=round(7.2/scale))
    old[~valid] = new[~valid] = np.nan
    difference[~safe] = np.nan
    del entries, wo, wn, valid, safe
    covered = np.isfinite(score)
    results = []
    for sign in (1, -1):
        maxima = local_peaks(score, sign, MIN_SCORE, round(.22/scale))
        for y, x in maxima:
            if min(y, x, score.shape[0]-y-1, score.shape[1]-x-1) < round(.55/scale):
                continue
            before, after = [aperture(im, x, y, scale, DIRECT_RADIUS, (.30, .48)) for im in (old, new)]
            if not np.isfinite(before+after):
                continue
            ra, dec = wcs.pixel_to_world_values(x+origin[0], y+origin[1])
            brightest = max(before, after)
            results.append(dict(band=band, detector=newer["detector"], sign=sign, ra=float(ra), dec=float(dec),
                                score=float(score[y, x]), old_njy=before, new_njy=after,
                                consistent_aperture_change=(brightest > 0 and sign*(after-before) >=
                                                            MIN_FRAC_CHANGE*brightest)))
    data = dict(difference=difference, covered=covered, wcs=wcs, origin=origin, scale=scale,
                older=old, newer=new)
    fit.update(band=band, detector=newer["detector"], retained_peaks=len(results),
               field_area_arcmin2=int(covered.sum())*scale**2/3600,
               overlap_with_f444w_arcmin2=sky_coverage(reference, data) if reference is not None else None)
    print(band, newer["detector"], len(results), fit["overlap_with_f444w_arcmin2"], flush=True)
    return pd.DataFrame(results), data, fit


def pairs(peaks):
    def coords(rows):
        return SkyCoord(rows.ra.to_numpy(), rows.dec.to_numpy(), unit="deg")

    matches = []
    for sign in (1, -1):
        by_band = {b: peaks[(peaks.band == b) & (peaks.sign == sign)].reset_index(drop=True)
                   for b in SEARCH_RADII}
        blue, red, middle = (by_band[b] for b in ("F200W", "F444W", "F356W"))
        if blue.empty or red.empty:
            continue
        red_ids, dist, _ = coords(blue).match_to_catalog_sky(coords(red))
        middle_sky = coords(middle) if not middle.empty else None
        for i in np.nonzero(dist.arcsec <= MATCH_RADIUS)[0]:
            sw, lw = blue.iloc[i], red.iloc[red_ids[i]]
            candidate_id = f"GS-B-{sw.ra:.6f}{sw.dec:+.6f}"
            row = dict(candidate=candidate_id, ra=sw.ra, dec=sw.dec, sign=sign, sw_detector=sw.detector,
                       offset_arcsec=dist.arcsec[i], sw_score=sw.score, lw_score=lw.score,
                       minimum_score=min(abs(sw.score), abs(lw.score)), sw_old_njy=sw.old_njy,
                       sw_new_njy=sw.new_njy, lw_old_njy=lw.old_njy, lw_new_njy=lw.new_njy,
                       aperture_gate=bool(sw.consistent_aperture_change and lw.consistent_aperture_change
                                          and max(abs(sw.score), abs(lw.score)) >= 6))
            if middle_sky is not None:
                mid_dist = middle_sky.separation(SkyCoord(sw.ra, sw.dec, unit="deg")).arcsec
                m = int(np.argmin(mid_dist))
                if mid_dist[m] <= MATCH_RADIUS:
                    mid = middle.iloc[m]
                    row.update(f356w_score=mid.score, f356w_old_njy=mid.old_njy, f356w_new_njy=mid.new_njy,
                               f356w_aperture_gate=bool(mid.consistent_aperture_change))
            matches.append(row)
    return pd.DataFrame(matches).sort_values("minimum_score", ascending=False) if matches else pd.DataFrame()


def atlas(table, arrays):
    selected = table[table.aperture_gate].head(8)
    if selected.empty:
        return
    figure = Figure(figsize=(10.7, 2.3*len(selected)+.7), layout="constrained")
    axes = figure.subplots(len(selected), 6, squeeze=False)
    for row, panels in zip(selected.itertuples(), axes, strict=True):
        for band, detector, subaxes in [("F200W", row.sw_detector, panels[:3]), ("F444W", "NRCALONG", panels[3:])]:
            group = arrays[band, detector]
            x, y = np.array(group["wcs"].world_to_pixel_values(row.ra, row.dec))-group["origin"]
            scale = group["scale"]
            half = round(.5/scale)
            cut = np.s_[round(y)-half:round(y)+half+1, round(x)-half:round(x)+half+1]
            old, new, diff = (group[key] for key in ("older", "newer", "difference"))
            lim = np.nanpercentile(abs(diff[cut]), 97.5)
            bounds = tuple(np.nanpercentile(np.r_[old[cut].ravel(), new[cut].ravel()], [6, 98.5]))
            for ax, name, image in zip(subaxes, ("2025", "2026", "residual"), (old, new, diff), strict=True):
                show_stamp(ax, image, x, y, scale, (-lim, lim) if name == "residual" else bounds,
                           "RdBu_r" if name == "residual" else "gray", .24, half_arcsec=.5)
                ax.set_title(f"{band} · {name}", fontsize=8)
        panels[0].set_ylabel(f"{row.ra:.6f}°\n{row.dec:+.6f}°", fontsize=8)
    figure.suptitle("GOODS-S field B · cross-detector candidates\nObserved JWST cutouts · fixed 0.24″ apertures", fontsize=13)
    save_figure(figure, FIGURES / "field_b_candidates")


def main():
    receipt = OUT / "mosaics.json"
    products = json.loads(receipt.read_text())["images"]
    results, arrays, metrics = [], {}, []
    for band, detectors in DETECTORS.items():
        for detector in detectors:
            records = {r["epoch"]: r for r in products
                       if r["group"] == "b" and r["detector"] == detector and r["band"] == band+"-CLEAR"}
            peaks, data, fit = search_pair(records, band, reference=arrays.get(("F444W", "NRCALONG")))
            results.append(peaks)
            arrays[band, detector] = data
            metrics.append(fit)
    peaks = pd.concat(results, ignore_index=True)
    results = pairs(peaks)
    peaks.to_csv(OUT / "field_b_peaks.csv", index=False)
    results.to_csv(OUT / "field_b_matches.csv", index=False)
    if results.empty:
        signed_counts, shortlisted = {}, 0
    else:
        signed_counts, shortlisted = dict(Counter(results.sign)), int(results.aperture_gate.sum())
        atlas(results, arrays)
    write_json(OUT / "field_b_search.json", {
        "executed_utc": utcnow(), "bands": list(DETECTORS), "sw_lw_matched_peaks": len(results),
        "sign_counts": signed_counts, "aperture_shortlist": shortlisted,
        "sum_chip_common_sw_lw_area_arcmin2": sum(r["overlap_with_f444w_arcmin2"] for r in metrics if r["band"] == "F200W"),
        "aperture_radius_arcsec": DIRECT_RADIUS, "match_radius_arcsec": MATCH_RADIUS,
        "gate_minimum_fractional_change": MIN_FRAC_CHANGE, "local_ranking_score_minimum": MIN_SCORE,
        "pair_diagnostics": metrics, "mosaics_sha256": text_digest(receipt),
        "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

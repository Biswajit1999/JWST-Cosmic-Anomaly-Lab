"""Measure blind November catalogue leads in real JWST science images."""

import argparse
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from audit_jwst_goods_s_prior_records import proposal
from compare_jwst_brightening_mast import REGION_ARCSEC, common_grid
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from measure_jwst_brightening_offset import BANDS as CORE_RADII
from measure_jwst_brightening_offset import fit_core, local_alignment, match_pair, patch
from measure_jwst_epoch_candidate import RADII, photometry, pixel_scale, show_stamp
from measure_jwst_module_b_brightening import ANNULUS
from read_jwst_i2d_cutout import MastImageSections
from search_jwst_november_catalogues import BANDS, OCTOBER_PREFIX, PREFIX, STORAGE_MIB, product_name
from search_jwst_november_prior_visits import STEM as PRIOR_VISITS_PREFIX
from search_jwst_november_unmatched import STEM as UNMATCHED_PREFIX
from vet_jwst_field_b_candidate import FIGURES
from vet_jwst_module_b import tns_catalogue, tns_match

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

DEDUP_ARCSEC = .16
PREVIOUS = {
    "JADES 103607 brightening": (53.18866983298658, -27.812539402767552),
    "GO 8060 transient 2510038": (53.162036, -27.746384),
    "GO 8060 transient 2510132": (53.147252, -27.790612),
    "GO 8060 transient 2510131": (53.187419, -27.756589),
    "GO 8060 transient 2510042": (53.19541375266448, -27.77932804913642),
}


def shortlist(matches, maximum):
    if matches.empty:
        return []
    kept = []
    for row in matches.itertuples(index=False):
        sky = SkyCoord(row.ra, row.dec, unit="deg")
        if any(sky.separation(other).arcsec <= DEDUP_ARCSEC for other, _ in kept):
            continue
        name = f"GSN-{row.ra:.6f}{row.dec:+.6f}"
        kept.append((sky, dict(candidate=name, icrs_deg=[row.ra, row.dec], sign=int(row.sign),
                               catalogue_bands=int(row.bands), catalogue_rank=float(row.rank),
                               older=row.older, newer=row.newer)))
    return kept[:maximum]


def fit_difference(template, new, band, sign):
    a, b = common_grid(template, new)
    _, wcs, x, y, x0, y0 = template
    scale = pixel_scale(wcs)
    radius = CORE_RADII.get(band, .24)
    try:
        reg = local_alignment(a, b, x, y, scale, radius)
    except ValueError as exc:
        if not str(exc).startswith(("Too few bright local reference sources", "Insufficient local core matches")):
            raise
        before, after, difference = a, b, b-a
        reg = dict(method="wcs_only", reason=str(exc))
    else:
        before, after, difference = match_pair(a, b, reg, scale, 0.)
        reg["method"] = "local_shift_gain_and_core_width"
    profile = photometry(difference, x, y, scale, RADII, ANNULUS)
    core = fit_core(*patch(sign*difference, x, y, scale, radius), radius)
    if core is not None:
        dx, dy = np.array(core["offset_xy_arcsec"])/scale
        core["icrs_deg"] = list(map(float, wcs.pixel_to_world_values(x+x0+dx, y+y0+dy)))
    pre = photometry(before, x, y, scale, RADII[-1:], ANNULUS)[0]
    post = photometry(after, x, y, scale, RADII[-1:], ANNULUS)[0]
    return dict(alignment=reg, change=profile, core=core, before=pre, after=post), (before, after, difference, x, y, scale)


def figure(lead, views):
    fig = Figure(figsize=(8.8, 8.3), layout="constrained")
    axes = fig.subplots(len(views), 3, squeeze=False)
    for row, (band, view) in zip(axes, views.items(), strict=True):
        before, after, difference, x, y, scale = view
        half = 1.35 if lead["catalogue_bands"] <= 3 else .55
        pixels = np.r_[patch(before, x, y, scale, half)[2].ravel(), patch(after, x, y, scale, half)[2].ravel()]
        if not np.isfinite(pixels).any():
            continue
        lo, hi = np.nanpercentile(pixels, (4, 98.5))
        values = patch(difference, x, y, scale, half)[2]
        limit = np.nanpercentile(abs(values), 98.5) if np.isfinite(values).any() else 1.
        for axis, (im, label, bounds, cmap) in zip(row, (
            (before, lead["older"], (lo, hi), "gray"),
            (after, lead["newer"], (lo, hi), "gray"),
            (difference, "New minus old", (-limit, limit), "RdBu_r"),
        ), strict=True):
            show_stamp(axis, im, x, y, scale, bounds, cmap, .24, half_arcsec=half)
            axis.set_title(f"{band} · {label}", fontsize=9)
    fig.suptitle(f"{lead['candidate']} · Observed STScI NIRCam mosaics", fontsize=12)
    stem = lead["image_basename"]
    save_figure(fig, FIGURES/stem)


def inspect(cache, sky, lead):
    views, bands = {}, {}
    for band in BANDS:
        images, products = [], []
        for visit in (lead["older"], lead["newer"]):
            name = product_name(visit, band, "i2d.fits")
            try:
                stamp, commit = MastImageSections(cache, name).cutout(sky, REGION_ARCSEC)
            except ValueError as exc:
                if not str(exc).startswith("Target lies outside the mosaic"):
                    raise
                break
            images.append(stamp)
            products.append(commit)
        if len(images) != 2:
            bands[band] = dict(products=products, status="outside_image")
            continue
        result, view = fit_difference(*images, band, lead["sign"])
        result["products"] = products
        bands[band], views[band] = result, view
        measurement = result["change"][-1]
        print(lead["candidate"], band, measurement["flux_njy"], measurement["local_mad_sigma_njy"], flush=True)
    if views:
        figure(lead, views)
    return bands


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=6)
    parser.add_argument("--source", choices=(PREFIX, OCTOBER_PREFIX, UNMATCHED_PREFIX, PRIOR_VISITS_PREFIX), default=PREFIX)
    args = parser.parse_args()
    source = OUT/f"{args.source}_matches.csv"
    candidates = shortlist(pd.read_csv(source), args.limit)
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    centre_lookup = {lead["candidate"]: sky for sky, lead in candidates}
    day = datetime.now(timezone.utc).strftime("%Y%m%d")
    plans = [proposal(cache, program, centre_lookup, day) for program in (8060, 12577)]
    tns = tns_catalogue(cache)
    results = []
    summary = OUT/f"{args.source}_vetting.json"
    for sky, lead in candidates:
        prior = {name: sky.separation(SkyCoord(*coords, unit="deg")).arcsec for name, coords in PREVIOUS.items()}
        lead["prior_project_sources_within_1arcsec"] = {name: float(value) for name, value in prior.items()
                                                      if value <= 1.}
        name = lead["candidate"]
        lead["plan_matches_within_1arcsec"] = [
            dict(program=plan["program"], catalogue=cat["name"], label_field=cat["label_field"], matches=matches)
            for plan in plans for cat in plan["catalogues"] if (matches := cat["positions"][name]["matches"])
        ]
        lead["tns"] = tns_match(tns, sky)
        print(lead["candidate"], lead["catalogue_bands"], lead["older"], lead["newer"],
              lead["prior_project_sources_within_1arcsec"], flush=True)
        lead["image_basename"] = f"{args.source}_{name.replace('.', 'p')}"
        lead["bands"] = inspect(cache, sky, lead)
        results.append(lead)
        write_json(summary, dict(
            executed_utc=utcnow(), planned=len(candidates), candidates=results,
            observation_plans=[dict(program=p["program"], product=p["product"]) for p in plans],
            shortlist_sha256=text_digest(source), script_sha256=text_digest(__file__),
            scope="Image-domain checks of the highest ranked unique sky positions. This scope does not"
                  " establish completeness, publication priority, or a physical class.",
        ))


if __name__ == "__main__":
    main()

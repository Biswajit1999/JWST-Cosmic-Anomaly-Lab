"""Find strong multi-band November detections absent from earlier source lists."""

import json
from itertools import product

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.table import QTable
from astroquery.mast import Observations
from download_jwst_cross_detector import OUT
from matplotlib.path import Path
from measure_jwst_epoch_candidate import verified_path
from search_jwst_november_catalogues import BANDS, EARLIER, LATER, LONG, PREFIX, product_name

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json

STEM = "november_goods_s_unmatched"
MIN_FLUX_NJY = 25.
MIN_APERTURE_SNR = 12.
SAME_EPOCH_RADIUS_ARCSEC = .13
OPPOSITE_EPOCH_RADIUS_ARCSEC = .20
EDGE_MARGIN_ARCSEC = 2.


def footprint(catalogue, visit, band):
    obs = product_name(visit, band, "i2d.fits").removesuffix("_i2d.fits")
    row = next(r for r in catalogue if r["obs_id"] == obs and int(r["calib_level"]) == 3)
    region = str(row["s_region"])
    if row["dataRights"] != "PUBLIC" or not region.startswith("POLYGON "):
        raise ValueError(f"Expected public stage-3 polygons for {obs}")
    return dict(obs_id=obs, obsid=str(row["obsid"]), s_region=region,
                data_uri=str(row["dataURL"]), start_mjd=float(row["t_min"]))


def inside(region, sky):
    ra0, de0 = sky.ra.deg[0], sky.dec.deg[0]
    scale = np.array([3600*np.cos(np.deg2rad(de0)), 3600.])
    zero = np.array([ra0, de0])
    points = (np.column_stack((sky.ra.deg, sky.dec.deg))-zero)*scale
    covered = np.zeros(len(sky), dtype=bool)
    for piece in region.strip().split("POLYGON")[1:]:
        vertices = (np.asarray(piece.split(), dtype=float).reshape(-1, 2)-zero)*scale
        x, y = vertices.T
        area2 = x@np.roll(y, -1)-y@np.roll(x, -1)
        counterclockwise = vertices if area2 > 0 else vertices[::-1]
        covered |= Path(counterclockwise).contains_points(points, radius=-2*EDGE_MARGIN_ARCSEC)
    return covered


def tables(ledger):
    catalogues = {}
    for entry in ledger["products"]:
        q = QTable.read(verified_path(entry["product"]), format="ascii.ecsv")
        coords = q["sky_centroid"].icrs
        ok = np.isfinite(coords.ra.deg+coords.dec.deg)
        flux = q["aper70_flux"].to_value(u.nJy)
        error = q["aper70_flux_err"].to_value(u.nJy)
        snr = np.divide(flux, error, out=np.full(len(flux), np.nan), where=error > 0)
        strong = ok & (flux >= MIN_FLUX_NJY) & (snr >= MIN_APERTURE_SNR)
        catalogues[entry["visit"], entry["band"]] = dict(
            all_coordinates=coords[ok], positive_coordinates=coords[strong], positive_snr=snr[strong],
            labels=np.asarray(q["label"])[strong],
        )
    return catalogues


def scan(catalogues, coverage, visit, comparisons):
    anchor = catalogues[visit, "F200W"]
    sky = anchor["positive_coordinates"]
    corroborates, no_counterpart, nearest_snr = {}, {}, {}
    for band in BANDS:
        lookup = catalogues[visit, band]
        point_id, distance, _ = sky.match_to_catalog_sky(lookup["positive_coordinates"])
        corroborates[band] = distance.arcsec <= SAME_EPOCH_RADIUS_ARCSEC
        nearest_snr[band] = lookup["positive_snr"][point_id]
        other_distances = [sky.match_to_catalog_sky(catalogues[other, band]["all_coordinates"])[1].arcsec
                           for other in comparisons]
        no_counterpart[band] = np.min(other_distances, axis=0) > OPPOSITE_EPOCH_RADIUS_ARCSEC
    cross_filter = {band: corroborates[band] & no_counterpart[band] for band in BANDS}
    initial = cross_filter["F200W"] & (np.sum(list(cross_filter.values()), axis=0) >= 3)
    rows = []
    for other in comparisons:
        if other.startswith("08060"):
            continue
        test = {b: cross_filter[b] & inside(coverage[other, b]["s_region"], sky) for b in BANDS}
        possible = test["F200W"] & (test[LONG[0]] | test[LONG[1]]) & (np.sum(list(test.values()), axis=0) >= 3)
        for i in np.flatnonzero(possible):
            bands = [b for b in BANDS if test[b][i]]
            new_source = visit in LATER
            old, new = (other, visit) if new_source else (visit, other)
            rows.append(dict(ra=float(sky.ra.deg[i]), dec=float(sky.dec.deg[i]), source_id=int(anchor["labels"][i]),
                             older=old, newer=new, sign=1 if new_source else -1, bands=len(bands),
                             rank=float(min(nearest_snr[b][i] for b in bands)), filters=",".join(bands)))
    return rows, dict(active_visit=visit, positive_f200w_anchors=len(sky),
                      before_geometry=int(initial.sum()), retained_pair_level=len(rows))


def main():
    prior = OUT/f"{PREFIX}.json"
    cats = tables(json.loads(prior.read_text()))
    Observations.TIMEOUT = 40
    archive = Observations.query_criteria(obs_collection="JWST", proposal_id="12577")
    coverage = {(v, b): footprint(archive, v, b) for v, b in product((*EARLIER[2:], *LATER), BANDS)}
    rows, diag = [], []
    for visit in (*EARLIER[2:], *LATER):
        other = EARLIER if visit in LATER else LATER
        matches, diagnostics = scan(cats, coverage, visit, other)
        print(diagnostics, flush=True)
        rows.extend(matches)
        diag.append(diagnostics)
    fields = ["ra", "dec", "source_id", "older", "newer", "sign", "bands", "rank", "filters"]
    path = OUT/f"{STEM}_matches.csv"
    pd.DataFrame(rows, columns=fields).sort_values(["bands", "rank"], ascending=False).to_csv(path, index=False)
    write_json(OUT/f"{STEM}.json", dict(
        executed_utc=utcnow(), inputs_sha256=text_digest(prior), products=list(coverage.values()),
        script_sha256=text_digest(__file__), output_sha256=text_digest(path), survivor_counts=diag,
        scope="Strong F200W anchors, with coeval detections in at least two other filters. Crossmatch"
              " absence uses the complete finite-position catalogue, including rows without valid flux."
              " MAST polygon overlap requires approximately two arcseconds of interior margin."
              " Geometric coverage alone does not establish useful pixel weight.",
        minimum_catalogue_flux_njy=MIN_FLUX_NJY, minimum_formal_aperture_snr=MIN_APERTURE_SNR,
        same_epoch_radius_arcsec=SAME_EPOCH_RADIUS_ARCSEC, opposite_epoch_radius_arcsec=OPPOSITE_EPOCH_RADIUS_ARCSEC,
        interior_margin_arcsec=EDGE_MARGIN_ARCSEC,
    ))


if __name__ == "__main__":
    main()

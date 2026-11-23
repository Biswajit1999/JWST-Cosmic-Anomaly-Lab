"""Verify the provenance of nearby catalogue redshifts for JADES 103607."""

import json
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from compare_jwst_brightening_mast import STORAGE_MIB
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import fetch
from measure_jwst_epoch_candidate import verified_path

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
ZFOURGE_ID = 14409


def tap_csv(cache, url, query, name, survey):
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query, MAXREC=4000)
    receipt = fetch(cache, url+"?"+urlencode(params), name, "public_redshift_catalogue",
                    survey=survey, version="TAP catalogue results")
    return pd.read_csv(verified_path(receipt)), receipt


def position_matches(table, sky, ra, dec, radius=1.0):
    if not np.isfinite(table[ra] + table[dec]).all():
        raise ValueError("Non-finite catalogue coordinates")
    distance = sky.separation(SkyCoord(table[ra], table[dec], unit="deg")).arcsec
    ranked = table.assign(separation_arcsec=distance).sort_values("separation_arcsec")
    return dict(rows=len(table), radius_arcsec=radius, matches=int((distance <= radius).sum()),
                closest=json.loads(ranked.head(8).to_json(orient="records", double_precision=10)))


def simbad_nearby(cache, sky):
    query = ("SELECT main_id,otype,ra,dec,rvz_redshift,rvz_bibcode FROM basic "
             f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{sky.ra.deg:.10f},{sky.dec.deg:.10f},"
             "0.002777777778))=1")
    data, receipt = tap_csv(cache, "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync", query,
                           "module_b_001_simbad_10arcsec.csv", "SIMBAD")
    return dict(product=receipt, **position_matches(data, sky, "ra", "dec"))


def vizier_entry(cache, catalog, fields, id_column, identifier, label, survey, sky):
    query = {"-source": catalog, "-out.max": 2, "-out": fields, id_column: identifier}
    receipt = fetch(cache, VIZIER+"?"+urlencode(query), f"module_b_001_{label}.tsv",
                    "published_redshift_provenance", survey=survey, version=catalog)
    data = pd.read_csv(verified_path(receipt), sep="\t", comment="#")
    selected = data.loc[pd.to_numeric(data[id_column], errors="coerce") == identifier]
    numeric = selected.apply(pd.to_numeric, errors="coerce")
    (row,) = json.loads(numeric.to_json(orient="records", double_precision=10))
    if row.get("zspec") is not None and row["zspec"] < 0:
        row["zspec"] = None
    if "RAJ2000" in row:
        coordinate = SkyCoord(row["RAJ2000"], row["DEJ2000"], unit="deg")
        row["separation_from_search_arcsec"] = float(sky.separation(coordinate).arcsec)
    return dict(product=receipt, row=row)


def zfourge_entry(cache, sky):
    source = "J/ApJ/830/51/"
    requests = {
        "zf_cdfs": "Seq,RAJ2000,DEJ2000,SNR,Nghb,zspec",
        "cdfsz": "Seq,zspec,zpk,zsp,zp,l68,u68",
        "cdfsxmatch": "Seq,[SWM2014],RAJ2000,DEJ2000,Sep,irAGN,radAGN,xAGN",
    }
    return {catalog: vizier_entry(cache, source+catalog, fields, "Seq", ZFOURGE_ID, f"zfourge_{catalog}",
                                  "ZFOURGE / CDS VizieR", sky) for catalog, fields in requests.items()}


def muse_catalogue(cache, sky):
    query = "SELECT ID,RAJ2000,DEJ2000,z,e_z,zconf,DataSet FROM safcat.AMUSED_MAIN_SOURCE_CAT_V1"
    data, receipt = tap_csv(cache, "https://archive.eso.org/tap_cat/sync", query, "module_b_001_muse_udf_dr2.csv",
                           "ESO Archive / MUSE HUDF DR2")
    if len(data) != 2221:
        raise ValueError("The MUSE UDF release query did not return the complete catalogue")
    return dict(product=receipt, **position_matches(data, sky, "RAJ2000", "DEJ2000"))


def main():
    prior = OUT/"module_b_001_six_filter_change.json"
    record = json.loads(prior.read_text())
    sky = SkyCoord(*record["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    simbad = simbad_nearby(cache, sky)
    zfourge = zfourge_entry(cache, sky)
    linked_id = zfourge["cdfsxmatch"]["row"]["[SWM2014]"]
    hduv = vizier_entry(cache, "J/ApJS/237/12/goodss", "ID,RAJ2000,DEJ2000,zbest,zpeak,zspec",
                        "ID", linked_id, "hduv_goodss", "HDUV / CDS VizieR", sky)
    write_json(OUT/"module_b_001_redshift_context.json", dict(
        executed_utc=utcnow(), candidate=record["candidate"], icrs_deg=record["icrs_deg"],
        simbad=simbad, zfourge=zfourge, hduv=hduv, muse_hudf_dr2=muse_catalogue(cache, sky),
        six_filter_sha256=text_digest(prior), script_sha256=text_digest(__file__),
        scope="Positional neighbours remain separate. ZFOURGE zpk and zp are photometric estimates.",
    ))


if __name__ == "__main__":
    main()

"""Check new spectroscopic catalogues and planned JWST visits at JADES 103607."""

import json
from io import StringIO
from urllib.parse import urlencode
from warnings import catch_warnings, simplefilter
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astroquery.mast import Observations
from audit_jwst_brightening_followup import mast_footprint_contains
from audit_jwst_brightening_redshifts import position_matches
from check_jwst_transient_precursor import contains
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from vet_jwst_field_b_candidate import API, FIGURES
from vet_jwst_module_b import archive_query

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_followup_20261126"
STSCI = "https://www.stsci.edu/jwst-program-info/visits/?download=&program=12577"
VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
SMILES_REV = "ca43c5e7ad769121b3dd9dc672fe5417e2aa7cd5"
SMILES_FILE = "hlsp_smiles_jwst_nirspec_goodss_multi_v1.0_catalog.fits"
VUDS = "https://vo.lam.fr/aspic_vuds/q/cone_ecdfs/scs.xml"
CONE_ARCSEC = 60
MATCH_ARCSEC = 1
PLANNED = (5, 6, 9, 10)


def cone(cache, sky, label, catalog, columns, numeric):
    url = (
        VIZIER
        + "?"
        + urlencode(
            {
                "-source": catalog,
                "-c": f"{sky.ra.deg:.9f} {sky.dec.deg:.9f}",
                "-c.rs": CONE_ARCSEC,
                "-out.max": 500,
                "-out": ",".join(columns),
            }
        )
    )
    path = cache.fetch(
        url,
        f"metadata/{TAG}_{label}.tsv",
        survey="CDS VizieR",
        version=catalog,
        kind="published_spectroscopic_catalogue",
    )
    body = path.read_text(encoding="utf-8")
    if f"#Name: {catalog}" not in body or "#Truncated" in body:
        raise ValueError(f"Unexpected VizieR response for {catalog}")
    data = pd.read_csv(StringIO(body), sep="\t", comment="#")
    for k in ("RAJ2000", "DEJ2000", *numeric):
        data[k] = pd.to_numeric(data[k], errors="coerce")
    data = data.dropna(subset=["RAJ2000", "DEJ2000"]).copy()
    return dict(
        product=receipt(path, url),
        query_radius_arcsec=CONE_ARCSEC,
        **position_matches(data, sky, "RAJ2000", "DEJ2000", MATCH_ARCSEC),
    )


def vuds_cone(cache, sky):
    url = (
        VUDS + "?" + urlencode(dict(RA=sky.ra.deg, DEC=sky.dec.deg, SR=CONE_ARCSEC / 3600, VERB=1))
    )
    path = cache.fetch(
        url,
        f"metadata/{TAG}_vuds.xml",
        survey="VUDS / CeSAM",
        version="VUDS DR1 ECDFS cone service",
        kind="published_spectroscopic_catalogue",
    )
    with catch_warnings():
        simplefilter("ignore", u.UnitsWarning)
        data = Table.read(path, format="votable")[
            ["ident", "alpha", "delta", "z", "zflags"]
        ].to_pandas()
    return dict(
        product=receipt(path, url),
        query_radius_arcsec=CONE_ARCSEC,
        **position_matches(data, sky, "alpha", "delta", MATCH_ARCSEC),
    )


def smiles_catalogue(cache, sky):
    url = f"https://raw.githubusercontent.com/ydzhuastro/smiles_dr2/{SMILES_REV}/{SMILES_FILE}"
    path = cache.fetch(
        url,
        f"metadata/{SMILES_FILE}",
        survey="SMILES",
        version=f"DR2 release {SMILES_REV}",
        kind="published_spectroscopic_catalogue",
    )
    data = Table.read(path, hdu="Z_SPEC")[["NIRCAM_ID", "RA", "DEC", "Z_SPEC", "Q"]].to_pandas()
    data["Q"] = data["Q"].str.decode("ascii")
    if len(data) != 166:
        raise ValueError("Unexpected size of the SMILES DR2 catalogue")
    return dict(
        product=receipt(path, url), **position_matches(data, sky, "RA", "DEC", MATCH_ARCSEC)
    )


def dawn_nirspec(cache, sky):
    url = (
        API
        + "/nirspec_extractions?"
        + urlencode(dict(coords=f"{sky.ra.deg:.10f},{sky.dec.deg:.10f}", size=3, output="csv"))
    )
    path = cache.fetch(
        url,
        f"metadata/{TAG}_dja_nirspec.csv",
        survey="DAWN JWST Archive",
        version="position query size=3",
        kind="public_position_spectrum_search",
    )
    data = pd.read_csv(path)[["root", "file", "srcid", "ra", "dec", "z", "grade", "grating"]]
    return dict(
        product=receipt(path, url), **position_matches(data, sky, "ra", "dec", MATCH_ARCSEC)
    )


def signed_margin_arcsec(footprint, sky):
    inside = mast_footprint_contains(footprint, sky)
    if inside is None:
        return None, None
    distance = []
    for piece in footprint.strip().split("POLYGON")[1:]:
        vertices = np.asarray(piece.split(), float).reshape(-1, 2)
        east, north = sky.spherical_offsets_to(SkyCoord(vertices[:, 0], vertices[:, 1], unit="deg"))
        start = np.column_stack((east.arcsec, north.arcsec))
        direction = np.roll(start, -1, axis=0) - start
        length = (direction**2).sum(axis=1)
        good = length > 0
        shift = np.clip(-np.einsum("ij,ij->i", start[good], direction[good]) / length[good], 0, 1)
        closest = start[good] + shift[:, None] * direction[good]
        distance.append(np.linalg.norm(closest, axis=1).min())
    signed = float(min(distance)) * (1 if inside else -1)
    return inside, signed


def visit_check(cache, sky):
    path = cache.fetch(
        STSCI,
        f"metadata/{TAG}_stsci.xml",
        survey="STScI",
        version="JWST GO 12577 visit status",
        kind="public_observing_schedule",
    )
    xml = ElementTree.parse(path).getroot()
    if xml.attrib.get("id") != "12577":
        raise ValueError("Unexpected STScI proposal number")
    schedule = {
        int(item.attrib["observation"]): {
            "status": item.findtext("status"),
            "configuration": item.findtext("configuration"),
            "target": item.findtext("target"),
            "plan_window": item.findtext("planWindow"),
            "actual_start_utc": item.findtext("startTime"),
            "actual_end_utc": item.findtext("endTime"),
        }
        for item in xml.findall("visit")
        if int(item.attrib["visit"]) == 1
    }
    Observations.TIMEOUT = 50
    selection = dict(obs_collection="JWST", proposal_id="12577")
    table = Observations.query_criteria(**selection)
    rows = []
    for r in table:
        obs_id = str(r["obs_id"])
        for group in sorted(schedule):
            if not obs_id.startswith((f"jw12577{group:03d}001", f"jw12577-o{group:03d}_")):
                continue
            footprint = str(r["s_region"])
            hit, margin = signed_margin_arcsec(footprint, sky)
            url = None if np.ma.is_masked(r["dataURL"]) else str(r["dataURL"])
            t = r["t_min"]
            mjd = float(t) if not np.ma.is_masked(t) and np.isfinite(t) else None
            rows.append(
                dict(
                    observation=group,
                    id=obs_id,
                    level=int(r["calib_level"]),
                    filter=str(r["filters"]),
                    where_polygon=footprint,
                    contains=hit,
                    edge_margin_arcsec=margin,
                    data_url=url,
                    archive_start_mjd=mjd,
                )
            )
    if not all(any(row["observation"] == i for row in rows) for i in (1, 4, *PLANNED)):
        raise ValueError("Programme query omitted an expected archive or planning visit")
    return dict(
        schedule_source=receipt(path, STSCI),
        schedule_generated=xml.findtext("reportTime"),
        schedule=schedule,
        mast_selection=selection,
        mast_total_rows=len(table),
        mast_matching_programme_rows=rows,
    )


def figure(sky, visits, catalogues):
    fig = Figure(figsize=(10, 4.7), layout="constrained")
    sky_ax, cat_ax = fig.subplots(1, 2)
    for obs, colour in ((4, "#367c8b"), (5, "#ad6440"), (6, "#76679b")):
        rows = [r for r in visits["mast_matching_programme_rows"] if r["observation"] == obs]
        r = next(t for t in rows if t["level"] == 3) if obs == 4 else rows[0]
        for i, poly in enumerate(r["where_polygon"].split("POLYGON")[1:]):
            points = np.asarray(poly.split(), float).reshape(-1, 2)
            east, north = sky.spherical_offsets_to(SkyCoord(points[:, 0], points[:, 1], unit="deg"))
            sky_ax.plot(
                np.r_[east.arcsec, east.arcsec[0]],
                np.r_[north.arcsec, north.arcsec[0]],
                color=colour,
                label=f"Visit {obs:03d}" if i == 0 else None,
                ls="-" if obs == 4 else "--",
                lw=1.3,
            )
    sky_ax.plot(0, 0, "*", ms=10, color="#141e23", label="Transient")
    sky_ax.set(
        xlabel="East offset (arcsec)",
        ylabel="North offset (arcsec)",
        title="Observed and planned JWST footprints",
    )
    sky_ax.set_aspect("equal")
    sky_ax.legend(fontsize=8, frameon=False)
    for i, (name, record) in enumerate(catalogues.items()):
        d = [
            r["separation_arcsec"]
            for r in record["closest"]
            if r["separation_arcsec"] <= CONE_ARCSEC
        ]
        if d:
            cat_ax.scatter(d, [i] * len(d), color="#367c8b", s=21)
    cat_ax.axvline(1, color="#a74636", lw=1, ls="--")
    cat_ax.set(
        xlim=(0, CONE_ARCSEC),
        xlabel="Separation from transient (arcsec)",
        yticks=range(len(catalogues)),
        yticklabels=list(catalogues),
        title="Nearest published spectroscopic entries",
    )
    fig.suptitle("JADES 103607 · archive and host checks", fontsize=13)
    save_figure(fig, FIGURES / TAG)


def main():
    old = OUT / "module_b_001_native_november_fade.json"
    source = json.loads(old.read_text())
    sky = SkyCoord(*source["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    cats = {
        "ACES": cone(
            cache,
            sky,
            "aces",
            "J/MNRAS/425/2116/table2",
            ("Seq", "Mask", "RAJ2000", "DEJ2000", "z", "zhel", "Q"),
            ("z", "zhel", "Q"),
        ),
        "VANDELS": cone(
            cache,
            sky,
            "vandels",
            "J/A+A/647/A150/cdfs",
            ("Name", "RAJ2000", "DEJ2000", "zsp", "q_zsp", "zph"),
            ("zsp", "zph"),
        ),
        "VUDS DR1": vuds_cone(cache, sky),
        "SMILES DR2": smiles_catalogue(cache, sky),
        "DAWN NIRSpec": dawn_nirspec(cache, sky),
    }
    dja, dja_receipt = archive_query(
        cache,
        "exposures",
        dict(ra=sky.ra.deg, dec=sky.dec.deg, arcmin=0.025, output="csv"),
        f"{TAG}_dja.csv",
    )
    rate = dja.loc[dja.extension.eq("rate")]
    rate = rate[rate.footprint.map(lambda p: contains(p, sky))]
    later = rate[rate.expstart > Time("2026-11-24T20:50:11").mjd]
    visits = visit_check(cache, sky)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            native_analysis_sha256=text_digest(old),
            icrs_deg=source["icrs_deg"],
            spectroscopic_catalogues=cats,
            jwst_go12577=visits,
            dawn=dict(
                receipt=dja_receipt,
                covering_rate_exposures=len(rate),
                latest_rate_start_mjd=float(rate.expstart.max()),
                later_than_november_visit=json.loads(
                    later[["file", "filter", "expstart"]].to_json(orient="records")
                ),
            ),
            scope="Search results are limited to the named catalogues and DAWN inventory. MAST future polygons "
            "are planning metadata. An MSA footprint does not document a slit assignment.",
        ),
    )
    figure(sky, visits, cats)
    for name, r in cats.items():
        print(name, r["rows"], r["matches"], r["closest"][:1])
    for obs in (4, *PLANNED):
        hits = [
            r["edge_margin_arcsec"]
            for r in visits["mast_matching_programme_rows"]
            if r["observation"] == obs and r["edge_margin_arcsec"] is not None
        ]
        print(
            obs,
            len(hits),
            min(hits, default=None),
            max(hits, default=None),
            visits["schedule"][obs],
        )


if __name__ == "__main__":
    main()

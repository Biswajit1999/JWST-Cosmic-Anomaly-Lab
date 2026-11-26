"""Record independent archives checked during the 2026 JADES 103607 event."""

import json
from urllib.parse import urlencode
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astroquery.mast import Observations
from check_jades_103607_followup import signed_margin_arcsec
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

SOURCE = OUT / "module_b_001_native_november_fade.json"
STEM = "module_b_001_independent_2026_coverage"
START = "2026-10-01"
END = "2026-11-27"
HEASARC = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
XSA = "https://nxsa.esac.esa.int/tap-server/tap/sync"


def footprint(region, sky):
    pieces = region.split()
    if pieces and pieces[0].upper() == "CIRCLE":
        values = pieces[2:] if pieces[1].upper() in {"ICRS", "J2000"} else pieces[1:]
        if len(values) == 3:
            ra, dec, radius = map(float, values)
            separation = sky.separation(SkyCoord(ra, dec, unit="deg")).deg
            margin = float((radius - separation) * 3600)
            return margin >= 0, margin
    return signed_margin_arcsec(region, sky)


def optional(row, key):
    value = row[key]
    if np.ma.is_masked(value):
        return None
    value = value.item() if isinstance(value, np.generic) else value
    return None if isinstance(value, float) and not np.isfinite(value) else value


def hst(sky, bounds):
    Observations.TIMEOUT = 50
    table = Observations.query_criteria(
        obs_collection="HST", coordinates=sky, radius=4 * u.arcmin, t_min=bounds
    )
    result = []
    for row in table:
        region = str(row["s_region"])
        contains, margin = footprint(region, sky)
        result.append(
            dict(
                id=optional(row, "obs_id"),
                proposal=optional(row, "proposal_id"),
                instrument=optional(row, "instrument_name"),
                filter=optional(row, "filters"),
                start_mjd=optional(row, "t_min"),
                end_mjd=optional(row, "t_max"),
                dataproduct_type=optional(row, "dataproduct_type"),
                rights=optional(row, "dataRights"),
                url=optional(row, "dataURL"),
                region=region,
                contains=contains,
                signed_edge_margin_arcsec=margin,
            )
        )
    return sorted(result, key=lambda r: (r["start_mjd"], r["id"]))


def high_energy(cache, mission, sky, bounds):
    q = (
        f"SELECT TOP 100 obsid,time,ra,dec,name FROM {mission} "
        f"WHERE CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{sky.ra.deg:.12f},{sky.dec.deg:.12f},0.5))=1 "
        f"AND time BETWEEN {bounds[0]} AND {bounds[1]} ORDER BY time"
    )
    url = HEASARC + "?" + urlencode(dict(REQUEST="doQuery", FORMAT="votable", LANG="ADQL", QUERY=q))
    path = cache.fetch(
        url,
        f"metadata/{STEM}_{mission}.vot",
        survey="HEASARC",
        version="TAP public mission log position and date query",
        kind="archive_metadata",
    )
    root = ElementTree.parse(path).getroot()
    statuses = [
        info.get("value")
        for info in root.findall(".//{*}INFO")
        if info.get("name") == "QUERY_STATUS"
    ]
    if statuses != ["OK"]:
        raise ValueError(f"{mission} did not return an untruncated archive log: {statuses}")
    table = Table.read(path, format="votable", use_names_over_ids=True)
    return dict(
        product=receipt(path, url), rows=json.loads(table.to_pandas().to_json(orient="records"))
    )


def xmm_observation_log(cache, sky):
    # Search a center rectangle wider than the half-degree cone at this declination.
    q = (
        "SELECT TOP 100 observation_id,ra,dec,start_utc FROM xsa.v_all_observations "
        f"WHERE ra BETWEEN {sky.ra.deg - 0.58:.9f} AND {sky.ra.deg + 0.58:.9f} "
        f"AND dec BETWEEN {sky.dec.deg - 0.51:.9f} AND {sky.dec.deg + 0.51:.9f} "
        f"AND start_utc BETWEEN '{START}' AND '{END}' ORDER BY start_utc"
    )
    url = XSA + "?" + urlencode(dict(REQUEST="doQuery", FORMAT="csv", LANG="ADQL", QUERY=q))
    path = cache.fetch(
        url,
        f"metadata/{STEM}_xsa.csv",
        survey="XMM Science Archive",
        version="ESA TAP position and date query",
        kind="archive_metadata",
    )
    table = pd.read_csv(path, dtype={"observation_id": str})
    if len(table) >= 100 or not {"observation_id", "ra", "dec", "start_utc"} <= set(table.columns):
        raise ValueError("XMM log response has unexpected schema or is truncated")
    if len(table):
        sep = sky.separation(SkyCoord(table.ra, table.dec, unit="deg")).deg
        table = table.assign(separation_deg=sep).loc[sep <= 0.5]
    return dict(product=receipt(path, url), rows=json.loads(table.to_json(orient="records")))


def main():
    coords = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*coords, unit="deg")
    bounds = [float(Time(d).mjd) for d in (START, END)]
    observations = hst(sky, bounds)
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    cxc = high_energy(cache, "chanmaster", sky, bounds)
    xmm = xmm_observation_log(cache, sky)
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_source_sha256=text_digest(SOURCE),
            position_icrs_deg=coords,
            start_filter_window_mjd=bounds,
            mast_hst_center_query_arcmin=4,
            mast_hst=observations,
            xray_pointing_cone_deg=0.5,
            chanmaster=cxc,
            xsa=xmm,
            scope="HST published coverage polygons and slit circles. Chandra and XMM pointings are checked by "
            "center in a 0.5-degree cone. Unsearched facilities and survey slews are outside this query.",
        ),
    )
    print(
        "HST",
        len(observations),
        "covering",
        [r["id"] for r in observations if r["contains"]],
        "unknown",
        [r["id"] for r in observations if r["contains"] is None],
    )
    print(
        "Examples at closest polygon or circle edge",
        max(
            (
                r["signed_edge_margin_arcsec"]
                for r in observations
                if r["signed_edge_margin_arcsec"] is not None
            ),
            default=None,
        ),
    )
    print("Chandra", len(cxc["rows"]), "XMM-Newton", len(xmm["rows"]))


if __name__ == "__main__":
    main()

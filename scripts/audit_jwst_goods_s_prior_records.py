"""Check GOODS-S variability positions against public observing plans and catalogues."""

import hashlib
import json
from datetime import datetime, timezone
from io import StringIO
from urllib.parse import urlencode
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import fetch
from measure_jwst_epoch_candidate import verified_path

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

RADIUS = 1.0
CATALOGUES = {
    "pouliasis_2019": ("vizier.cfa.harvard.edu", "J/MNRAS/487/4285/table3", 113, ("ID", "Grade", "z")),
    "villforth_2010": ("vizier.idia.ac.za", "J/ApJ/723/737/table5", 155, ("GOODS", "Var")),
    "sarajedini_2011": ("vizier.idia.ac.za", "J/ApJ/731/97/table1", 85, ("Name", "Signi", "z")),
    "lyu_2024": ("vizier.cfa.harvard.edu", "J/ApJ/966/229/table2", 217, ("ID", "JADES", "Sample", "z")),
    "luo_2017": ("vizier.cfa.harvard.edu", "J/ApJS/228/2/sources", 1055, ("OType", "errPos", "zspec", "zF")),
}


def distances(data, centres, ra, dec):
    sources = SkyCoord(data[ra], data[dec], unit="deg")
    return {name: sky.separation(sources).arcsec for name, sky in centres.items()}


def positional_rows(data, separation, columns):
    nearest = data.iloc[np.argsort(separation)[:3]].loc[:, columns].copy()
    nearest["separation_arcsec"] = separation[np.argsort(separation)[:3]]
    return json.loads(nearest.to_json(orient="records", double_precision=10))


def catalogue(cache, name, specification, centres):
    host, table, expected_count, columns = specification
    params = {"-source": table, "-out.max": "unlimited"}
    if name in {"pouliasis_2019", "lyu_2024"}:
        params["-out.all"] = 1
    url = f"https://{host}/viz-bin/asu-tsv?{urlencode(params)}"
    receipt = fetch(cache, url, f"goods_s_{name}_full.tsv", "published_counterpart_catalogue",
                    survey="CDS VizieR", version=table)
    data = pd.read_csv(verified_path(receipt), comment="#", sep="\t")
    ra, dec = "RAJ2000", "DEJ2000"
    data[[ra, dec]] = data[[ra, dec]].apply(pd.to_numeric, errors="coerce")
    data = data.dropna(subset=[ra, dec])
    if len(data) != expected_count or not np.isfinite(data[[ra, dec]].to_numpy()).all():
        raise ValueError(f"{table}: expected {expected_count} valid positions, found {len(data)}")
    coordinates = {}
    for label, sep in distances(data, centres, ra, dec).items():
        coordinates[label] = dict(matches=int((sep <= RADIUS).sum()),
                                  nearest=positional_rows(data, sep, [ra, dec, *columns]))
    return dict(product=receipt, rows=len(data), radius_arcsec=RADIUS, positions=coordinates)


def plan_catalogue(node, centres, plans):
    name = node.findtext("./{*}Name")
    csv = node.findtext("./{*}CatalogAsCsv")
    header, *lines = csv.strip().splitlines()
    data = pd.read_csv(StringIO("\n".join([header.lstrip("#"), *[line for line in lines if not line.startswith("#")]])))
    if data[["ID", "RA", "DEC"]].isna().any().any():
        raise ValueError(f"Incomplete APT target coordinates in {name}")
    label = next((field for field in data if field.startswith("[LABEL] - ")), None)
    columns = ["ID", "RA", "DEC", "Redshift", "Weight"] + ([label] if label else [])
    assigned = {plan: configs for plan, configs in plans.items() if configs["catalog"]["name"] == name}
    positions = {}
    for key, sep in distances(data, centres, "RA", "DEC").items():
        nearest = positional_rows(data, sep, columns)
        in_radius = data.iloc[np.flatnonzero(sep <= RADIUS)].loc[:, columns].copy()
        in_radius["separation_arcsec"] = sep[sep <= RADIUS]
        matches = json.loads(in_radius.sort_values("separation_arcsec").to_json(
            orient="records", double_precision=10))
        for row in matches:
            row["configurations"] = [
                dict(plan=plan, configuration=c["name"], role=role)
                for plan, p in assigned.items() for c in p["configs"]
                for role in ("primary", "filler") if row["ID"] in c[f"{role}Ids"]
            ]
        positions[key] = dict(nearest=nearest, matches=matches)
    return dict(name=name, rows=len(data), label_field=label,
                csv_sha256=hashlib.sha256(csv.encode()).hexdigest(), positions=positions)


def proposal(cache, number, centres, day):
    url = f"https://www.stsci.edu/jwst-program-info/download/jwst/apt/{number}/"
    receipt = fetch(cache, url, f"goods_s_{number}_{day}.aptx", "public_observing_plan",
                    survey="STScI JWST observing programmes", version=f"Public APT {number}")
    with ZipFile(verified_path(receipt)) as archive:
        files = [name for name in archive.namelist() if name.endswith(".json") and not name.startswith("edu.")]
        plans = {name: json.loads(archive.read(name)) for name in files}
        root = ElementTree.fromstring(archive.read(f"{number}.xml"))
        content = [plan_catalogue(node, centres, plans) for node in root.iter()
                   if node.find("./{*}CatalogAsCsv") is not None]
    if len(content) < 4 or not plans:
        raise ValueError(f"Incomplete APT inventory for {number}")
    return dict(program=number, product=receipt, radius_arcsec=RADIUS, catalogues=content)


def main():
    source = OUT/"module_b_vetting.json"
    original = OUT/"field_b_transient_ancillary.json"
    leads = json.loads(source.read_text())["candidates"]
    centres = {lead["candidate"]: SkyCoord(*lead["icrs_deg"], unit="deg") for lead in leads}
    row = json.loads(original.read_text())
    centres[row["candidate"]] = SkyCoord(*row["position_degrees"], unit="deg")
    cache = DownloadCache(ROOT, 18_000*1024**2)
    day = datetime.now(timezone.utc).strftime("%Y%m%d")
    plans = [proposal(cache, n, centres, day) for n in (8060, 12577)]
    catalogues = {name: catalogue(cache, name, spec, centres) for name, spec in CATALOGUES.items()}
    write_json(OUT/"goods_s_prior_records.json", dict(
        executed_utc=utcnow(), positions_icrs_deg={key: [v.ra.deg, v.dec.deg] for key, v in centres.items()},
        public_apt=plans, literature_catalogues=catalogues,
        vetting_sha256=text_digest(source), original_position_sha256=text_digest(original),
        script_sha256=text_digest(__file__),
        scope="Public planning snapshots and selected historical catalogues. A planning row is neither an observed"
              " spectrum nor evidence of recognition of variability unless its label or assignment supports that inference.",
    ))


if __name__ == "__main__":
    main()

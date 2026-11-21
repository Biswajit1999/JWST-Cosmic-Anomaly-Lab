"""Restore the fixed long-period search sample from public Kepler listings."""

import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlencode

import requests
import yaml

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/kepler_long_period"
CONFIG_PATH = ROOT / "configs/kepler_long_period.yaml"
config = yaml.safe_load(CONFIG_PATH.read_text())
cache = DownloadCache(ROOT, config["raw_storage_budget_mib"] * 1024**2)
base = "https://archive.stsci.edu/missions/kepler/lightcurves/"


def fetch(url, name, kind, survey="Kepler"):
    path = cache.fetch(url, f"exoplanets/{name}", survey=survey, version="public archive", kind=kind)
    return dict(url=url, path=path.relative_to(ROOT).as_posix(), sha256=digest(path))


original = json.loads((ROOT / "results/exoplanet_pilot/selection.json").read_text())
listing = fetch(base + config["search_group"] + "/", "kepler_0114.html", "directory")
eligible = {int(v) for v in re.findall(r'href="(\d{9})/"', (ROOT / listing["path"]).read_text())}
previous = original["search_targets"]
ranked = sorted(eligible.difference(previous, [original["known_control"]]),
                key=lambda kic: hashlib.sha256(f"{config['selection_seed']}:{kic}".encode()).hexdigest())
roster = ([dict(kic=config["control_kic"], role="confirmed_control")]
          + [dict(kic=kic, role="previous_short_period_sample") for kic in previous]
          + [dict(kic=kic, role="new_search_sample") for kic in ranked[:config["new_targets"]]])
write_json(OUT / "selection.json", dict(executed_utc=utcnow(), roster=roster,
    listing=listing, config_sha256=text_digest(CONFIG_PATH), before_lightcurve_search=True))
results = []
for star in roster:
    kic = star["kic"]
    padded = f"{kic:09d}"
    url = f"{base}{padded[:4]}/{padded}/"
    products = []
    directory = None
    errors = []
    try:
        directory = fetch(url, f"kic{kic}_directory.html", "directory")
    except (OSError, requests.RequestException) as exc:
        errors.append(dict(url=url, error=str(exc)))
    if directory:
        files = sorted(set(re.findall(r'href="(kplr\d+-\d+_llc.fits)"',
                                      (ROOT / directory["path"]).read_text())))
        for name in files:
            try:
                products.append(fetch(url + name, name, "observed_lightcurve"))
            except (OSError, requests.RequestException) as exc:
                errors.append(dict(url=url + name, error=str(exc)))
    results.append(dict(**star, directory=directory, products=products, errors=errors))
    write_json(OUT / "products.json", dict(executed_utc=utcnow(), targets=results))
    print(kic, star["role"], len(products), len(errors), flush=True)

# This census is inspected after the candidate ephemerides are fitted.
ids = ",".join(str(star["kic"]) for star in roster)
queries = {
    "koi": f"select kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_time0bk,koi_duration from cumulative where kepid in ({ids})",
    "tce": f"select kepid,tce_plnt_num,tce_period,tce_time0bk from q1_q17_dr25_tce where kepid in ({ids})",
}
catalogues = {}
for table, adql in queries.items():
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urlencode(dict(query=adql, format="json"))
    record = fetch(url, f"kepler_long_period_{table}.json", "catalogue_snapshot", "NASA Exoplanet Archive")
    rows = json.loads((ROOT / record["path"]).read_text())
    if not isinstance(rows, list) or not any(row["kepid"] == config["control_kic"] for row in rows):
        raise ValueError(f"Unexpected {table} query response")
    catalogues[table] = dict(**record, rows=len(rows))
write_json(OUT / "catalogues.json", dict(executed_utc=utcnow(), **catalogues))

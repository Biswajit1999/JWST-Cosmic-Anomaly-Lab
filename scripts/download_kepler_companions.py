"""Freeze a host sample from the pinned KOI table and restore public light curves."""

import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlencode

import requests
import yaml

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_companions"
config_path = root / "configs/kepler_companions.yaml"
config = yaml.safe_load(config_path.read_text())
cache = DownloadCache(root, config["raw_storage_budget_mib"] * 1024**2)


def fetch(url, name, kind, survey="Kepler"):
    path = cache.fetch(url, f"exoplanets/{name}", survey=survey, version="public archive", kind=kind)
    return dict(url=url, path=path.relative_to(root).as_posix(), sha256=digest(path))


expected = json.loads((root / "results/exoplanet_pilot/summary.json").read_text())["koi_sha256"]
catalogue = next(r for r in cache.records if r["path"] == "data/raw/exoplanets/koi_cumulative.json"
                 and r["sha256"] == expected)
cached = fetch(catalogue["url"], "koi_cumulative.json", "catalogue_snapshot", "NASA Exoplanet Archive")
assert cached["sha256"] == expected, "Pinned KOI snapshot has changed"
by_host = defaultdict(list)
for record in json.loads((root / cached["path"]).read_text()):
    by_host[record["kepid"]].append(record)
controls = config["positive_controls"]
excluded = set(controls) | {config["alias_control_kic"]}
eligible = []
for kic, rows in by_host.items():
    confirmed = [r for r in rows if r["koi_disposition"] == "CONFIRMED" and r["kepler_name"]]
    if len(confirmed) != 1 or kic in excluded:
        continue
    r = confirmed[0]
    if (r["koi_srad"] is not None and config["stellar_radius_solar"][0] <= r["koi_srad"] <= config["stellar_radius_solar"][1]
            and r["koi_period"] is not None and config["known_period_days"][0] <= r["koi_period"] <= config["known_period_days"][1]
            and r["koi_time0bk"] is not None and r["koi_duration"] and r["koi_duration"] > 0):
        eligible.append(kic)
eligible.sort(key=lambda kic: hashlib.sha256(f"{config['selection_seed']}:{kic}".encode()).hexdigest())
selection = [(kic, "withheld_planet_control", name) for kic, name in controls.items()]
selection += [(config["alias_control_kic"], "known_alias_control", None)]
selection += [(kic, "search_host", None) for kic in eligible[:config["single_transit_hosts"]]]
roster = []
for kic, role, withheld in selection:
    confirmed = [r for r in by_host[kic] if r["koi_disposition"] == "CONFIRMED"]
    masks = [r for r in confirmed if r["kepler_name"] != withheld]
    assert all(r["koi_period"] and r["koi_time0bk"] is not None and r["koi_duration"] for r in masks)
    roster.append(dict(kic=kic, role=role, withheld_planet=withheld, masks=masks,
                       cumulative_koi=by_host[kic]))
write_json(out / "selection.json", dict(executed_utc=utcnow(), config_sha256=text_digest(config_path),
    source_catalogue=cached, eligible_single_confirmed_transit_hosts=len(eligible), targets=roster))

products = []
for star in roster:
    kic = star["kic"]
    key = f"{kic:09d}"
    url = f"https://archive.stsci.edu/missions/kepler/lightcurves/{key[:4]}/{key}/"
    files, failures = [], []
    try:
        listing = fetch(url, f"kic{kic}_directory.html", "directory")
        names = sorted(set(re.findall(r'href="(kplr\d+-\d+_llc.fits)"',
                                      (root / listing["path"]).read_text())))
        for name in names:
            try:
                files.append(fetch(url+name, name, "observed_lightcurve"))
            except (requests.RequestException, OSError) as exc:
                failures.append(dict(url=url+name, error=str(exc)))
    except (requests.RequestException, OSError) as exc:
        listing = None
        failures.append(dict(url=url, error=str(exc)))
    products.append(dict(kic=kic, listing=listing, products=files, failures=failures))
    write_json(out / "products.json", dict(executed_utc=utcnow(), targets=products))
    print(kic, star["role"], len(files), len(failures), flush=True)

ids = ",".join(str(r["kic"]) for r in roster)
for name, query in {
    "tce": f"select kepid,tce_plnt_num,tce_period,tce_time0bk from q1_q17_dr25_tce where kepid in ({ids})",
    "recent_koi": f"select kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_time0bk,koi_duration from cumulative where kepid in ({ids})",
}.items():
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urlencode(dict(query=query, format="json"))
    path = fetch(url, f"kepler_companion_{name}.json", "catalogue_snapshot", "NASA Exoplanet Archive")
    rows = json.loads((root / path["path"]).read_text())
    assert rows and all(isinstance(r["kepid"], int) for r in rows)
    write_json(out / (name+"_receipt.json"), dict(executed_utc=utcnow(), **path, rows=len(rows)))

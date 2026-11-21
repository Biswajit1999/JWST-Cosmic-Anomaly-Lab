"""Retrieve a bounded TESS replication sample and further Kepler quarters."""

import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlencode

from cosmic_anomaly_lab.provenance import DownloadCache, digest, utcnow, write_json

parser = argparse.ArgumentParser()
parser.add_argument("mission", choices=["kepler", "tess"])
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
out = root / "results" / (args.mission + "_followup")
cache = DownloadCache(root, 4 * 1024**3)


def fetch(url, name, kind, survey):
    path = cache.fetch(url, f"exoplanets/{name}", survey=survey, version="public archive",
                       kind=kind)
    return dict(url=url, path=path.relative_to(root).as_posix(), sha256=digest(path))


products = []
if args.mission == "kepler":
    kid = 11446424
    base = f"https://archive.stsci.edu/missions/kepler/lightcurves/0114/{kid:09d}/"
    listing = fetch(base, f"kic{kid}_directory.html", "directory", "Kepler")
    names = sorted(set(re.findall(r'href="(kplr\d+-\d+_llc.fits)"',
                                  (root / listing["path"]).read_text())))
    assert len(names) >= 2
    for name in names:
        products.append(fetch(base + name, name, "observed_lightcurve", "Kepler"))
    tpf_base = base.replace("lightcurves", "target_pixel_files")
    pixel_list = fetch(tpf_base, f"kic{kid}_tpf_directory.html", "directory", "Kepler")
    available = set(re.findall(r'href="(kplr\d+-\d+_lpd-targ.fits.gz)"',
                               (root / pixel_list["path"]).read_text()))
    old = json.loads((root / "results/exoplanet_pilot/summary.json").read_text())
    catalogue_input = next(r for r in cache.records if r["path"] == "data/raw/exoplanets/koi_cumulative.json"
                           and r["sha256"] == old["koi_sha256"])
    catalogue = fetch(catalogue_input["url"], "koi_cumulative.json", "catalogue_snapshot",
                      "NASA Exoplanet Archive")
    if catalogue["sha256"] != old["koi_sha256"]:
        raise ValueError("The live KOI table no longer matches the pinned catalogue snapshot")
    ancillary = dict(directory=listing, pixel_directory=pixel_list,
                     listed_target_pixel_files=sorted(available), catalogue=catalogue)
else:
    control = 100100827
    listings, urls = [], {}
    for sector in [2, 3]:
        url = f"https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{sector}_lc.sh"
        rec = fetch(url, f"tess_sector_{sector}_lc.sh", "download_manifest", "TESS")
        listings.append(rec)
        found = re.findall(r'(-s\d{4}-(\d{16})-\d+-s_lc\.fits)\s+(https://\S+)',
                           (root / rec["path"]).read_text())
        urls[sector] = {int(tic): link for _, tic, link in found}
    common = urls[2].keys() & urls[3].keys()
    assert control in common
    ranked = sorted(common - {control},
                    key=lambda tic: hashlib.sha256(f"tess-s02-s03-v1:{tic}".encode()).hexdigest())
    ids = [control, *ranked[:16]]
    write_json(out / "selection.json", dict(control=control, search_targets=ranked[:16],
               eligible_intersection=len(common), selection="SHA256 ascending with fixed seed tess-s02-s03-v1",
               chosen_before_light_curve_inspection=True, listings=listings))
    for tic in ids:
        for sector in [2, 3]:
            url = urls[sector][tic]
            name = url.rsplit("/", 1)[-1]
            products.append(dict(tic=tic, sector=sector, role="control" if tic == control else "search",
                            **fetch(url, name, "observed_lightcurve", "TESS")))
            print(tic, sector, flush=True)
    def query(statement, name):
        url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urlencode(
            {"query": statement, "format": "json"})
        return fetch(url, name, "catalogue_snapshot", "NASA Exoplanet Archive")

    id_list = ",".join(map(str, ids))
    names = ",".join(f"'TIC {tic}'" for tic in ids)
    ancillary = dict(
        toi=query(f"select tid,toi,tfopwg_disp,pl_orbper from toi where tid in ({id_list})", "tess_followup_toi.json"),
        planets=query(f"select pl_name,tic_id,pl_orbper from pscomppars where tic_id in ({names})",
                      "tess_followup_planets.json"),
    )
write_json(out / "products.json", dict(executed_utc=utcnow(), products=products, ancillary=ancillary))

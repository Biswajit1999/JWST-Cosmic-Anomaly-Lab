"""Finite public Kepler search sample and published JWST atmosphere benchmark."""

import json
import re
from pathlib import Path
from urllib.parse import urlencode

from cosmic_anomaly_lab.provenance import DownloadCache, digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
out = root / "results/exoplanet_pilot"
cache = DownloadCache(root, budget_bytes=4096*1024**2)
base = "https://archive.stsci.edu/missions/kepler/lightcurves/0114/"
listing = cache.fetch(base, "exoplanets/kepler_0114.html", survey="Kepler", version="live-2026-11-16", kind="directory")
ids = sorted(set(map(int, re.findall(r'href="(\d{9})/"', listing.read_text()))))
control = 11446443  # TrES-2: known positive control, excluded from discovery count.
selected = sorted([i for i in ids if i != control], key=lambda i: (abs(i-control), i))[:12]
write_json(out / "selection.json", dict(executed_utc=utcnow(), known_control=control,
    search_targets=selected, method="Twelve nearest numerical KIC IDs to control in fixed 0114 directory; not angular neighbours",
    quarters_requested=[2, 3], selection_before_lightcurve_or_KOI_inspection=True,
    directory_sha256=digest(listing)))
products = []
for kid in [control] + selected:
    url = base + f"{kid:09d}/"
    directory = cache.fetch(url, f"exoplanets/kic{kid}_directory.html", survey="Kepler", version="live-2026-11-16", kind="directory")
    links = re.findall(r'href="([^"]+_llc.fits)"', directory.read_text())
    for quarter, stamp in [(2, "2009259160929"), (3, "2009350155506")]:
        record = dict(kic=kid, role="known_control" if kid == control else "search", quarter=quarter)
        candidates = [name for name in links if stamp in name]
        if len(candidates) != 1:
            record.update(status="UNAVAILABLE", reason="No unique listed long-cadence product for fixed quarter")
        else:
            path = cache.fetch(url+candidates[0], "exoplanets/"+candidates[0], survey="Kepler", version="released-long-cadence", kind="observed_lightcurve")
            record.update(status="DOWNLOADED", url=url+candidates[0], path=path.relative_to(root).as_posix(), sha256=digest(path))
        products.append(record)
        write_json(out / "products.json", dict(executed_utc=utcnow(), products=products))
        print(kid, quarter, record["status"], flush=True)
queries = {
    "confirmed_count": "select count(*) as n from pscomppars",
    "koi_cumulative": "select kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_time0bk,koi_duration,koi_depth,koi_prad,koi_srad,koi_steff,koi_smass,ra,dec from cumulative",
    "wasp39_parameters": "select pl_name,hostname,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_radj,pl_radjerr1,pl_radjerr2,pl_bmassj,pl_bmassjerr1,pl_bmassjerr2,pl_eqt,st_teff,st_rad,sy_dist from pscomppars where pl_name='WASP-39 b'",
}
for name, query in queries.items():
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urlencode(dict(query=query, format="json"))
    path = cache.fetch(url, f"exoplanets/{name}.json", survey="NASA Exoplanet Archive", version="live-2026-11-16", kind="catalogue_snapshot")
    json.loads(path.read_text())
    print(name, path.stat().st_size, flush=True)
metadata = cache.fetch("https://zenodo.org/api/records/7185300", "exoplanets/wasp39_zenodo.json", survey="JWST ERS", version="7185300", kind="release_metadata")
entry = next(f for f in json.loads(metadata.read_text())["files"] if f["key"] == "G395H_paper_data.zip")
path = cache.fetch(entry["links"]["self"], "exoplanets/G395H_paper_data.zip", survey="JWST ERS", version="7185300", kind="published_atmospheric_spectra_and_models", checksum=entry["checksum"])
write_json(out / "atmosphere_product.json", dict(path=path.relative_to(root).as_posix(), sha256=digest(path),
    url=entry["links"]["self"], published_checksum=entry["checksum"], doi="10.5281/zenodo.7185300",
    role="known_planet_benchmark_not_discovery"))
# The known-control pixel cube was selected from the archived directory listing.
pixel_record = out / "pixel_product.json"
if pixel_record.exists():
    p = json.loads(pixel_record.read_text())
    path = cache.fetch(p["url"], p["path"].removeprefix("data/raw/"), survey="Kepler",
                       version="released-Q2", kind="observed_target_pixel_cube")
    assert digest(path) == p["sha256"]

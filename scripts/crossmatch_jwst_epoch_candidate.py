"""Cross-check the brightening position against specific public lists."""

import json
import re
from pathlib import Path
from urllib.parse import urlencode

import requests
from astropy.coordinates import SkyCoord

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jwst_epoch_search"
RADIUS = 5


def published_positions(html):
    coordinates = []
    for table_name in ("S3.T4.2", "S5.T9.2"):
        pattern = rf'<table\b[^>]*\bid="{re.escape(table_name)}"[^>]*>(.*?)</table>'
        table = re.search(pattern, html, re.DOTALL)
        if table is None:
            raise ValueError(f"Missing published table {table_name}")
        positions = re.findall(r"\b(\d{2}:\d{2}:\d{2}\.\d+)\s*<[^>]+>\s*<[^>]+>\s*([+-]\d{2}:\d{2}:\d{2}\.\d+)",
                               table.group(1))
        if not positions:
            raise ValueError(f"No coordinates decoded from {table_name}")
        coordinates.extend(positions)
    return coordinates


def main():
    path = OUT / "novelty_queries.json"
    receipt = json.loads(path.read_text())
    ra, dec = receipt["ra"], receipt["dec"]
    sky = SkyCoord(ra, dec, unit="deg")
    cache = DownloadCache(ROOT, 6200 * 1024**2)
    queries = receipt["queries"]
    tns = queries["tns"]
    tns_path = ROOT / tns["path"]
    if digest(tns_path) != tns["sha256"]:
        raise ValueError("TNS response checksum changed")
    if 'class="no-results">No results found' in tns_path.read_text(encoding="utf-8"):
        tns.update(matches_within_radius=0, radius_arcsec=RADIUS)

    simbad = "https://simbad.cfa.harvard.edu/simbad/sim-coo?" + urlencode(
        {"Coord": f"{ra:.8f} {dec:.8f}", "CooFrame": "ICRS", "Radius": RADIUS,
         "Radius.unit": "arcsec", "output.format": "ASCII"})
    paper = "https://arxiv.org/html/2406.05060v4"
    for name, url, filename in [
        ("simbad_harvard", simbad, "simbad_harvard_target.txt"),
        ("jades_2022_2023_transients", paper, "decoursey2024_arxiv_v4.html"),
    ]:
        try:
            downloaded = cache.fetch(url, "epoch_search/"+filename, survey=name, version="fixed-coordinate comparison",
                                     kind="public_novelty_check")
        except requests.RequestException as exc:
            queries[name] = dict(status="UNAVAILABLE", url=url, error=str(exc))
            continue
        result = dict(url=url, path=downloaded.relative_to(ROOT).as_posix(),
                      sha256=digest(downloaded), radius_arcsec=RADIUS)
        content = downloaded.read_text(encoding="utf-8")
        if name == "simbad_harvard":
            count = re.search(r"Number of objects\s*:\s*(\d+)", content)
            objects = re.findall(r"^\s*\d+\|\s*([\d.]+)\|([^|]+)\|([^|]+)", content, re.MULTILINE)
            if count is None or int(count[1]) != len(objects):
                raise ValueError("SIMBAD result row count is unverified")
            result.update(matches_within_radius=len(objects),
                          objects=[dict(separation_arcsec=float(d), name=n.strip(), type=t.strip())
                                   for d, n, t in objects])
        else:
            positions = published_positions(content)
            ras, decs = zip(*positions, strict=True)
            separation = sky.separation(SkyCoord(list(ras), list(decs), unit=("hourangle", "deg"))).arcsec
            result.update(searched_positions=len(positions), nearest_separation_arcsec=float(separation.min()),
                          matches_within_radius=int((separation <= RADIUS).sum()), published_tables=[4, 9])
        queries[name] = dict(status="PARSED", **result)
    write_json(path, {**receipt, "expanded_utc": utcnow(), "crossmatch_script_sha256": text_digest(__file__)})


if __name__ == "__main__":
    main()

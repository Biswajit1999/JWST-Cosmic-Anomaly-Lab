"""Retrieve line-control spectra and the source-linked instrument reference."""

import re
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd

from cosmic_anomaly_lab.provenance import DownloadCache, digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
records = []
catalogue = pd.read_csv(root / "data/catalogues/observations.csv.gz")
row = catalogue[
    (catalogue.source_id == "CAL-82610407f3945338") & (catalogue.grating == "G235M")
].iloc[0]
path = cache.fetch(
    row.url, "spectra/" + row.file, survey="JADES", version="DJA-v4.4", kind="observed_spectrum"
)
records.append(dict(path=path.relative_to(root).as_posix(), sha256=digest(path)))
directory = "https://jades.herts.ac.uk/DR4/goods-s/spectra/f170lp-g235m/goods-s-mediumjwst/"
listing = cache.fetch(
    directory,
    "metadata/dr4_goods-s-mediumjwst_g235m.html",
    survey="JADES",
    version="DR4-listing-snapshot",
    kind="archive_directory",
)
matches = [
    href
    for href in re.findall(r'href="([^"]+)"', listing.read_text())
    if "-00159717_" in href and href.endswith("_x1d.fits")
]
assert len(matches) == 1
path = cache.fetch(
    urljoin(directory, matches[0]),
    "independent/" + matches[0],
    survey="JADES",
    version="DR4-v1.0",
    kind="independently_reduced_observed_spectrum",
)
records.append(dict(path=path.relative_to(root).as_posix(), sha256=digest(path)))
url = "https://jwst-docs.stsci.edu/files/216455634/216455636/1/1762452060535/jwst_nirspec_g395h_disp.fits"
path = cache.fetch(
    url,
    "instrument/g395h_disp.fits",
    survey="STScI",
    version="JDox-linked-reference-2026-11-15",
    kind="instrument_resolution_reference",
)
records.append(dict(path=path.relative_to(root).as_posix(), sha256=digest(path)))
write_json(
    root / "results/line_anchor/downloads.json", dict(executed_utc=utcnow(), products=records)
)
print(records)

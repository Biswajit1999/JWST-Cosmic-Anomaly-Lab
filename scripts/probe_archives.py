"""Retrieve small metadata first, then the explicitly budgeted DJA snapshot."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cosmic_anomaly_lab.provenance import DownloadCache, utcnow, write_json

root = Path(__file__).resolve().parents[1]
cache = DownloadCache(root)
products = [
    (
        "https://zenodo.org/records/15472354/files/dja_msaexp_emission_lines_v4.4.columns.csv",
        "dja_columns.csv",
        "DJA",
        "4.4",
        None,
    ),
    (
        "https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits",
        "jades_dr4.fits",
        "JADES",
        "DR4-1.2.1",
        None,
    ),
    (
        "https://drive.google.com/uc?export=download&id=10BQOEvwzcazLkJFmH4EzNWlwY70xejYi",
        "uncover_dr4_1.dat",
        "UNCOVER",
        "DR4.1",
        None,
    ),
    (
        "https://zenodo.org/records/15472354/files/dja_msaexp_emission_lines_v4.4.csv.gz",
        "dja_v4.4.csv.gz",
        "DJA",
        "4.4",
        "md5:721c58c446f4a4eb1570aef02c4ab47b",
    ),
]
outcomes = []
for url, name, survey, version, checksum in products:
    try:
        path = cache.fetch(
            url, name, survey=survey, version=version, kind="catalogue", checksum=checksum
        )
        result = dict(survey=survey, url=url, bytes=path.stat().st_size, status="downloaded")
    except Exception as exc:
        result = dict(survey=survey, url=url, status="failed", error=str(exc))
    print(result, flush=True)
    outcomes.append(result | {"accessed_utc": utcnow()})
    write_json(root / "data/manifests/archive_probes.json", outcomes)

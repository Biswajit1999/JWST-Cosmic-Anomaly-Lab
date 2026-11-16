"""Budgeted public HST/Spitzer photometry and Chandra catalogue retrieval."""

from pathlib import Path

from cosmic_anomaly_lab.provenance import DownloadCache, digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
products = []
for field, directory in [("goodsn", "GOODS-N"), ("goodss", "GOODS-S")]:
    base = "https://archive.stsci.edu/missions/hlsp/3d-hst/RELEASE_V4.0/Photometry/" + directory + "/"
    products.extend([
        (base + field + "_3dhst.v4.1.cats.tar.gz", f"ancillary/{field}_3dhst.v4.1.cats.tar.gz", "3D-HST", "v4.1"),
        (base + field + f"_3dhst.v4.1.cats/{field}_readme.v4.1.txt", f"ancillary/{field}_readme.v4.1.txt", "3D-HST", "v4.1"),
    ])
for field, journal, table in [("cdfn", "224/15", "table3.dat"), ("cdfs", "228/2", "table4.dat")]:
    for name in ["ReadMe", table]:
        products.append((f"https://cdsarc.cds.unistra.fr/ftp/J/ApJS/{journal}/{name}",
                         f"ancillary/{field}/{name}", "Chandra", journal))
records = []
for url, filename, survey, version in products:
    record = dict(url=url, survey=survey, version=version)
    try:
        path = cache.fetch(url, filename, survey=survey, version=version, kind="published_catalogue_or_documentation")
        record.update(status="DOWNLOADED", path=path.relative_to(root).as_posix(), sha256=digest(path))
    except Exception as exc:
        record.update(status="UNAVAILABLE", reason=str(exc))
    records.append(record)
    write_json(root / "results/multiwavelength/downloads.json", dict(executed_utc=utcnow(), products=records))
    print(filename, record["status"], record.get("reason", ""), flush=True)

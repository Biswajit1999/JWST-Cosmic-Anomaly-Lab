"""Retrieve the published Chandra main catalogues through NASA's TAP mirror."""

from pathlib import Path
from urllib.parse import urlencode

from astropy.table import Table

from cosmic_anomaly_lab.provenance import DownloadCache, digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
records = []
for table, sample, expected in [("cdfn2msnew", "Main", 683), ("chandfs7ms", "M", 1008)]:
    query = f"SELECT * FROM {table} WHERE source_sample='{sample}'"
    url = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable", QUERY=query, MAXREC=5000))
    record = dict(table=table, query=query, expected_rows=expected, url=url)
    try:
        path = cache.fetch(url, f"ancillary/{table}.xml", survey="Chandra",
                           version="HEASARC-main-2026-11-16", kind="published_catalogue")
        data = Table.read(path, format="votable").to_pandas()
        if len(data) != expected or not {"ra", "dec", "error_radius"}.issubset(data.columns):
            raise ValueError(f"Unexpected table shape: {data.shape}")
        record.update(status="VERIFIED", rows=len(data), path=path.relative_to(root).as_posix(),
                      sha256=digest(path))
    except Exception as exc:
        record.update(status="UNAVAILABLE", reason=str(exc))
    records.append(record)
    write_json(root / "results/multiwavelength/chandra_downloads.json",
               dict(executed_utc=utcnow(), products=records))
    print(table, record["status"], record.get("reason", ""), flush=True)

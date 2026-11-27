"""Retrieve TESS full-frame cutouts and planned footprints at JADES 103607."""

import json
import re
from datetime import datetime
from urllib.parse import urlencode

from bs4 import BeautifulSoup
from check_jades_103607_independent_2026 import END, SOURCE, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_tess_inventory"
CUTTER = "https://mast.stsci.edu/tesscut/api/v0.1"
POINTING = (
    "https://heasarc.gsfc.nasa.gov/wsgi-scripts/TESS/TESS-point_Web_Tool/"
    "TESS-point_Web_Tool/wtv_v2.0.py/RADec_result"
)
SCHEDULE = "https://tess.mit.edu/tess-year-9-observations/"


def rows(path):
    doc = BeautifulSoup(path.read_bytes(), "html.parser")
    return [
        [cell.get_text(" ", strip=True) for cell in tr.find_all(["th", "td"])]
        for tr in doc.find_all("tr")
    ]


def coverage_models(cache, position):
    form = {key: str(value) for key, value in position.items()}
    model = cache.fetch(
        POINTING, "metadata/module_b_001_tess_point.html", survey="NASA TESS Science Support Center",
        version="TESS-point v2.6 public pointing model", kind="planned_detector_footprints",
        form_payload=form,
    )
    headers = ["RA", "Dec", "EclipLong", "EclipLat", "Sector", "Camera", "CCD", "ColPix", "RowPix", "Sector Midpoint"]
    modelling = rows(model)
    if headers not in modelling:
        raise ValueError("Missing TESS-point coordinate result headings")
    predictions = [
        dict(sector=int(r[4]), camera=int(r[5]), ccd=int(r[6]),
             column=float(r[7]), row=float(r[8]), sector_midpoint=r[9])
        for r in modelling if len(r) == len(headers) and r[4].isdigit()
    ]
    if not predictions:
        raise ValueError("Unexpected empty TESS-point response")
    schedule = cache.fetch(
        SCHEDULE, "metadata/module_b_001_tess_year9.html", survey="TESS Science Office",
        version="Public Year 9 observing calendar", kind="planned_observation_dates",
    )
    date_ranges = {}
    for row in rows(schedule):
        if len(row) < 2 or not row[0].isdigit():
            continue
        dates = re.findall(r"\d{1,2}/\d{1,2}/\d{2,4}", row[1])
        if len(dates) != 2:
            continue
        fmt = "%m/%d/%y" if len(dates[0].rsplit("/", 1)[1]) == 2 else "%m/%d/%Y"
        date_ranges[int(row[0])] = [datetime.strptime(value, fmt).date().isoformat() for value in dates]
    if not {108, 109, 110, 111} <= date_ranges.keys():
        raise ValueError("Missing Year 9 observing dates around the JWST window")
    recent = [
        dict(sector=number, planned_start_utc=start, planned_end_utc=end,
             predicted_detectors=[p for p in predictions if p["sector"] == number])
        for number, (start, end) in sorted(date_ranges.items()) if start < END and end >= START
    ]
    return dict(predicted_covering_detectors=predictions, sectors_during_2026_window=recent,
                pointing_receipt=receipt(model, POINTING), schedule_receipt=receipt(schedule, SCHEDULE))


def main():
    coords = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    position = dict(ra=coords[0], dec=coords[1])
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    url = f"{CUTTER}/sector?{urlencode(position)}"
    listing = cache.fetch(
        url, "metadata/module_b_001_tess_sectors.json", survey="NASA TESS / MAST",
        version="TESSCut sectors 0.1", kind="observed_sector_coverage",
    )
    archive = json.loads(listing.read_text(encoding="utf-8"))
    fields = sorted(archive["results"], key=lambda row: int(row["sector"]))
    if not fields or len({int(row["sector"]) for row in fields}) != len(fields):
        raise ValueError("Expected one public TESSCut detector per covering sector")
    coverage = coverage_models(cache, position)
    products = []
    for row in fields:
        sector = int(row["sector"])
        size = 15 if sector == 108 else 13
        query = urlencode(dict(**position, x=size, y=size, sector=sector))
        remote = f"{CUTTER}/astrocut?{query}"
        path = cache.fetch(
            remote, f"cross_detector/tess_ffi/module_b_001_tess_{sector}_{size}x{size}.zip",
            survey="NASA TESS / MAST", version=f"TESSCut SPOC FFI {size}x{size}",
            kind="observed_ffi_cutout",
        )
        products.append(dict(
            sector=sector, camera=int(row["camera"]), ccd=int(row["ccd"]),
            size_pixels=size, cutout=receipt(path, remote),
        ))
        print(f"TESS sector {sector}: {path.stat().st_size / 1024**2:.1f} MiB", flush=True)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        position_source_sha256=text_digest(SOURCE), position_icrs_deg=coords,
        archive_sector_receipt=receipt(listing, url), tesscut_products=products,
        event_window_utc=[START, END], **coverage,
        scope="TESSCut identifies public SPOC FFI products at the position. "
        "NASA TESS-point and the TESS Science Office calendar give predicted footprints "
        "and planned dates. Neither alone establishes usable pixel data.",
    ))


if __name__ == "__main__":
    main()

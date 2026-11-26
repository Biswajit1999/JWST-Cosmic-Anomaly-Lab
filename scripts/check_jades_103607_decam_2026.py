"""Check DECam detector footprints and archive release dates around the JWST event."""

import argparse
import json
from collections import Counter
from math import atan2, isfinite

from astropy.coordinates import SkyCoord
from check_jades_103607_followup import signed_margin_arcsec
from check_jades_103607_independent_2026 import END, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

SOURCE = OUT / "module_b_001_native_november_fade.json"
TAG = "module_b_001_decam_2026_coverage"
URL = "https://astroarchive.noirlab.edu/api/adv_search/find/?rectype=hdu&limit=2000"
FILE_COLUMNS = (
    "archive_filename",
    "dateobs_center",
    "md5sum",
    "exposure",
    "ifilter",
    "proposal",
    "release_date",
    "depth",
    "seeing",
)
CHIP_COLUMNS = ("hdu_idx", "ra_min", "ra_max", "dec_min", "dec_max", "ra_center", "dec_center")
CORNERS = tuple(f"COR{i}{coordinate}1" for i in range(1, 5) for coordinate in ("RA", "DEC"))
CENTRE_BOX_HALF_WIDTH_DEG = 0.25


def query_body(ra, dec):
    return dict(
        outfields=(
            [f"file:{column}" for column in FILE_COLUMNS]
            + [f"hdu:{column}" for column in (*CHIP_COLUMNS, *CORNERS, "DETPOS")]
        ),
        search=[
            ["file:instrument", "decam"],
            ["file:dateobs_center", START, END],
            ["file:proc_type", "instcal"],
            ["file:prod_type", "image"],
            ["hdu:ra_center", ra - CENTRE_BOX_HALF_WIDTH_DEG, ra + CENTRE_BOX_HALF_WIDTH_DEG],
            ["hdu:dec_center", dec - CENTRE_BOX_HALF_WIDTH_DEG, dec + CENTRE_BOX_HALF_WIDTH_DEG],
        ],
    )


def decode_row(row, sky):
    metadata = {key: row[f"file:{key}"] for key in FILE_COLUMNS}
    chip = {key: row[f"hdu:{key}"] for key in CHIP_COLUMNS}
    vertices = [row.get(f"hdu:{key}") for key in CORNERS]
    chip["archived_corners_deg"] = vertices
    chip["detector_position"] = row.get("hdu:DETPOS")
    if all(value is not None and isfinite(float(value)) for value in vertices):
        positions = [(float(vertices[i]), float(vertices[i + 1])) for i in range(0, 8, 2)]
        mean_ra = sum(ra for ra, _ in positions) / len(positions)
        mean_dec = sum(dec for _, dec in positions) / len(positions)
        # FITS corner numbers follow detector pixel indices. Order them around the perimeter.
        positions.sort(key=lambda pos: atan2(pos[1] - mean_dec, pos[0] - mean_ra))
        corners = " ".join(f"{ra} {dec}" for ra, dec in positions)
        contains, edge = signed_margin_arcsec(f"POLYGON {corners}", sky)
    else:
        contains, edge = None, None
    chip.update(contains_target=contains, approximate_signed_edge_arcsec=edge)
    return dict(**metadata, **chip)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="Query the live archive again.")
    args = parser.parse_args()
    ra, dec = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(ra, dec, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    body = query_body(ra, dec)
    filename = f"metadata/{TAG}_hdu.json"
    path = cache.fetch(
        URL,
        filename,
        survey="NSF NOIRLab Astro Data Archive",
        version="Public advanced search HDU metadata",
        kind="observed_detector_footprint_metadata",
        json_payload=body,
        refresh=args.refresh,
    )
    header, *rows = json.loads(path.read_text(encoding="utf-8"))
    if (
        header.get("META", {}).get("endpoint") != "adv_search/find"
        or len(rows) >= 2000
        or header["PARAMETERS"]["json_payload"] != body
    ):
        raise ValueError("Unexpected or truncated NOIRLab chip response")
    matching = [decode_row(row, sky) for row in rows if row["file:dateobs_center"][:10] < END]
    matching.sort(key=lambda r: (r["dateobs_center"], r["md5sum"], r["hdu_idx"]))
    inside = [row for row in matching if row["contains_target"]]
    if len({(row["md5sum"], row["hdu_idx"]) for row in matching}) != len(matching):
        raise ValueError("Repeated image and detector product returned by NOIRLab")
    archive_key = f"data/raw/{filename}"
    accessed = next(
        record["accessed_utc"]
        for record in reversed(cache.records)
        if record["path"] == archive_key and record["url"] == URL
    )
    observation_dates = sorted({row["dateobs_center"][:10] for row in inside})
    pending = sorted(
        {
            row["release_date"]
            for row in inside
            if row["release_date"] and row["release_date"] > accessed[:10]
        }
    )
    depths = [row["depth"] for row in inside if row["depth"] is not None]
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=[ra, dec],
            position_sha256=text_digest(SOURCE),
            search_window_utc_end_exclusive=[START, END],
            archive=receipt(path, URL),
            archive_retrieved_utc=accessed,
            archive_request=body,
            chip_rows_screened=len(matching),
            unclassifiable_chip_rows=sum(row["contains_target"] is None for row in matching),
            containing_exposure_count=len({row["md5sum"] for row in inside}),
            containing_chips_per_night=dict(Counter(row["dateobs_center"][:10] for row in inside)),
            containing_observation_dates_utc=observation_dates,
            reported_available_exposures=len(
                {
                    row["md5sum"]
                    for row in inside
                    if row["release_date"] and row["release_date"] <= accessed[:10]
                }
            ),
            pending_release_date_range=[pending[0], pending[-1]] if pending else None,
            pipeline_depth_magnitude_range=[min(depths), max(depths)] if depths else None,
            containing_chip_rows=inside,
            screened_chip_rows=matching,
            scope="Polygons join four calibrated-chip corner coordinates from public metadata. "
            "Chip geometry can include masked pixels or imperfect corners. Exposure-level "
            "pipeline magnitude depths do not measure sensitivity at the target. "
            "No proprietary images or pixels were accessed.",
        ),
    )


if __name__ == "__main__":
    main()

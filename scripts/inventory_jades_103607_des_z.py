"""Inventory Dark Energy Survey z-band exposures at the JADES 103607 position."""

import csv
import json
import re
from collections import Counter, defaultdict
from datetime import date
from urllib.parse import parse_qs, urlencode, urlsplit
from xml.etree import ElementTree

from astropy.coordinates import SkyCoord
from astropy.table import Table
from check_jades_103607_decam_2026 import SOURCE, URL, decode_row, query_body
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_des_z_inventory"
SIA = "https://datalab.noirlab.edu/sia/des_dr2_se"
TAP = "https://datalab.noirlab.edu/tap/sync"
PROPOSAL = "2012B-0001"
BAND = "z DECam SDSS c0004 9260.0 1520.0"
WINDOW = (date(2013, 7, 1), date(2019, 7, 1))


def footprint_query(cache, sky, start, end):
    body = query_body(sky.ra.deg, sky.dec.deg)
    body["search"] = [
        criterion for criterion in body["search"] if criterion[0] != "file:dateobs_center"
    ] + [
        ["file:dateobs_center", str(start), str(end)],
        ["file:proposal", PROPOSAL],
        ["file:ifilter", BAND],
    ]
    body["outfields"] += ["file:EXPNUM", "file:OBJECT", "hdu:CCDNUM", "hdu:FWHM"]
    path = cache.fetch(
        URL,
        f"metadata/{TAG}/noirlab_{start}_{end}.json",
        survey="NSF NOIRLab Astro Data Archive",
        version="Historical DES DECam z-band CCD footprints",
        kind="observed_detector_footprint_metadata",
        json_payload=body,
    )
    header, *rows = json.loads(path.read_text(encoding="utf-8"))
    if header["META"]["endpoint"] != "adv_search/find" or header["PARAMETERS"]["json_payload"] != body:
        raise ValueError("Unexpected NOIRLab detector query response")
    if len(rows) >= 2000:
        mid = start + (end - start) // 2
        if mid == start:
            raise ValueError(f"One day of DECam z-band metadata reaches the archive row limit: {start}")
        for span in ((start, mid), (mid, end)):
            yield from footprint_query(cache, sky, *span)
        return
    matching = []
    unresolved = 0
    for row in rows:
        if not str(start) <= row["file:dateobs_center"][:10] < str(end):
            continue
        decoded = decode_row(row, sky)
        if decoded["contains_target"] is None:
            unresolved += 1
        if decoded["contains_target"]:
            matching.append(dict(
                **decoded, expnum=int(row["file:EXPNUM"]), ccdnum=int(row["hdu:CCDNUM"]),
                object=row["file:OBJECT"], archive_chip_fwhm_pixels=row["hdu:FWHM"],
            ))
    print(f"DES z {start} to {end}: {len(rows)} CCD rows, {len(matching)} containing target", flush=True)
    yield matching, dict(
        **receipt(path, URL), start=str(start), end_exclusive=str(end),
        archive_rows=len(rows), unresolved_footprints=unresolved,
    )


def lookup_cutouts(cache, sky):
    query_url = f"{SIA}?{urlencode(dict(POS=f'{sky.ra.deg},{sky.dec.deg}', SIZE='0.008,0.0075', VERB=2))}"
    path = cache.fetch(
        query_url, f"metadata/{TAG}/des_dr2_se.vot",
        survey="NSF NOIRLab Astro Data Lab", version="DES DR2 single-epoch image service",
        kind="observed_cutout_inventory",
    )
    status = [
        item.get("value") for item in ElementTree.parse(path).findall(".//{*}INFO")
        if item.get("name") == "QUERY_STATUS"
    ]
    if status != ["OK"]:
        raise ValueError(f"Unexpected DES SIA status {status}")
    table = Table.read(path, format="votable", use_names_over_ids=True)
    products = defaultdict(dict)
    for row in table:
        url = str(row["access_url"])
        parts = urlsplit(url)
        params = parse_qs(parts.query)
        if (parts.scheme, parts.netloc, parts.path) != ("https", "datalab.noirlab.edu", "/svc/cutout"):
            raise ValueError(f"Unexpected DES cutout service URL {url}")
        if params["col"] != ["des_dr2_se"]:
            raise ValueError(f"Unexpected DES cutout collection {url}")
        filename = params["siaRef"][0]
        match = re.fullmatch(r"D(\d{8})_([grizY])_c(\d{2})_r\d+p\d+_immasked\.fits\.fz", filename)
        ext = params["extn"][0]
        if match is None or ext not in ("1", "2"):
            raise ValueError(f"Unexpected DES image name or extension {url}")
        if match[2] != "z":
            continue
        channel = "science" if ext == "1" else "mask"
        if channel == "science" and str(row["obs_bandpass"]).strip() != BAND:
            raise ValueError(f"DES z-band science metadata disagrees with filename {filename}")
        candidates = products[(int(match[1]), int(match[3]), filename)]
        if channel in candidates and candidates[channel] != url:
            raise ValueError(f"Ambiguous cutout URL for {filename} {channel}")
        candidates[channel] = url
    return products, len(table), receipt(path, query_url)


def calibration(cache, keys):
    all_ids = sorted({expnum for expnum, _, _ in keys})
    ccds = ",".join(str(ccdnum) for ccdnum in sorted({c for _, c, _ in keys}))
    calibrations, inputs = {}, []
    for offset in range(0, len(all_ids), 80):
        subset = ",".join(map(str, all_ids[offset:offset+80]))
        query = (
            "SELECT expnum,ccdnum,mag_zero,sigma_mag_zero,flag,source,version,imagename "
            "FROM des_dr2.y6_gold_zeropoint "
            f"WHERE band='z' AND ccdnum IN ({ccds}) AND expnum IN ({subset})"
        )
        url = f"{TAP}?{urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', MAXREC=10000, QUERY=query))}"
        path = cache.fetch(
            url, f"metadata/{TAG}/fgcm_{offset:04d}.csv",
            survey="NSF NOIRLab Astro Data Lab", version="DES Y6 Gold FGCM y6a1_v2.1",
            kind="observed_photometric_calibration",
        )
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != [
                "expnum", "ccdnum", "mag_zero", "sigma_mag_zero", "flag", "source", "version", "imagename"
            ]:
                raise ValueError(f"Unexpected DES calibration response {path}")
            rows = list(reader)
        if len(rows) >= 10000:
            raise ValueError("DES calibration TAP query reached MAXREC")
        for row in rows:
            key = (int(row["expnum"]), int(row["ccdnum"]), row["imagename"] + ".fz")
            if key not in keys:
                continue
            value = dict(
                ab_zp=float(row["mag_zero"]), zeropoint_error_mag=float(row["sigma_mag_zero"]),
                calibration_flag=int(row["flag"]), calibration_source=row["source"],
                calibration_version=row["version"],
            )
            if key in calibrations and value != calibrations[key]:
                raise ValueError(f"Ambiguous FGCM calibration for {key}")
            calibrations[key] = value
        inputs.append(receipt(path, url))
    return calibrations, inputs


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    footprints, metadata = [], []
    boundaries = [date(year, month, 1) for year in range(2013, 2020) for month in (1, 4, 7, 10)]
    for start, end in zip(boundaries, boundaries[1:], strict=False):
        if start < WINDOW[0] or end > WINDOW[1]:
            continue
        for rows, source in footprint_query(cache, sky, start, end):
            footprints.extend(rows)
            metadata.append(source)
    cutouts, sia_rows, sia_receipt = lookup_cutouts(cache, sky)
    coverage = {(row["expnum"], row["ccdnum"]) for row in footprints}
    matching_keys = {key for key in cutouts if key[:2] in coverage}
    calibrations, tap_receipts = calibration(cache, matching_keys)
    products = [
        dict(expnum=exposure, ccdnum=ccd, filename=name,
             cutouts=cutouts[exposure, ccd, name], fgcm=calibrations.get((exposure, ccd, name)))
        for exposure, ccd, name in sorted(matching_keys)
    ]
    print(f"{len(coverage)} covering exposure/CCD pairs, {len(products)} matched cutout products, "
          f"{len(calibrations)} with FGCM", flush=True)
    print(f"Calibration flags: {Counter(row['calibration_flag'] for row in calibrations.values())}", flush=True)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        position_source_sha256=text_digest(SOURCE), position_icrs_deg=position,
        window_utc=[str(x) for x in WINDOW], proposal=PROPOSAL, filter=BAND,
        footprint_queries=metadata, cutout_inventory=sia_receipt, calibration_queries=tap_receipts,
        sia_row_count_all_filters=sia_rows, covering_exposure_ccd_pairs=len(coverage),
        covering_exposure_count=len({e for e, _ in coverage}),
        matching_cutout_products=len(products), matching_zeropoints=len(calibrations),
        footprints=sorted(footprints, key=lambda row: (row["dateobs_center"], row["md5sum"], row["hdu_idx"])),
        products=products,
        scope="Four archived CCD corner coordinates establish geometry. DES DR2 product identification "
        "uses archive exposure and CCD IDs. FGCM additionally requires the exact DESDM image basename. "
        "Nominal target coverage does not qualify pixel masks, seeing or point-source sensitivity.",
    ))


if __name__ == "__main__":
    main()

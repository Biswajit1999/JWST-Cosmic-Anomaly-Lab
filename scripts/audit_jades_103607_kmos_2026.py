"""Check the individual KMOS integral-field units between the two JWST visits."""

import json
import math
from collections import Counter
from datetime import datetime, timezone

from astropy.coordinates import FK5, SkyCoord
from astropy.io import fits
from astropy.time import Time
from bs4 import BeautifulSoup
from download_jades_103607_lco_eso import BUDGET
from download_jades_103607_lco_eso import TAG as INVENTORY
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_kmos_2026"
HEADER = "https://archive.eso.org/hdr?DpId="
PATROL_RADIUS_ARCSEC = 216.0
IFU_SIDE_ARCSEC = 2.8


def sexagesimal_number(value):
    # ESO stores HHMMSS.s or signed DDMMSS.s as a floating-point number.
    absolute = abs(float(value))
    main = int(absolute // 10000)
    minutes = int(absolute // 100) % 100
    seconds = absolute % 100
    if minutes >= 60 or seconds >= 60:
        raise ValueError(f"Invalid ESO sexagesimal number: {value}")
    return math.copysign(main + minutes / 60 + seconds / 3600, value)


def header_data(path, dataset):
    content = BeautifulSoup(path.read_bytes(), "html.parser").find("pre")
    if content is None:
        raise ValueError(f"ESO header text missing for {dataset}")
    cards = [line.strip() for line in content.get_text().splitlines() if line.strip()]
    header = fits.Header.fromstring("\n".join(cards), sep="\n")
    if header.get("ARCFILE") != dataset + ".fits" or header.get("INSTRUME") != "KMOS":
        raise ValueError(f"ESO header identity mismatch for {dataset}")
    return header


def arm_coordinates(header, target):
    arms = []
    for index in range(1, 25):
        ocs = f"ESO OCS ARM{index} "
        instrument_key = f"ESO INS ARM{index} LOCK"
        item = dict(arm=index, target_name=header.get(ocs + "NAME"), type=header.get(ocs + "TYPE"),
                    vignetted=header.get(ocs + "VIGNET"), locked=header.get(instrument_key))
        if ocs + "ALPHA" in header and ocs + "DELTA" in header:
            ra = 15 * sexagesimal_number(header[ocs + "ALPHA"])
            dec = sexagesimal_number(header[ocs + "DELTA"])
            sky = SkyCoord(ra, dec, unit="deg", frame=FK5(equinox=Time("J2000")))
            if item["type"] in ("O", "S"):
                matched = sky.icrs
                separation = float(matched.separation(target).arcsec)
                item.update(
                    icrs_deg=[float(matched.ra.deg), float(matched.dec.deg)],
                    separation_arcsec=separation,
                    inside_nominal_ifu_circumscribed_circle=(
                        separation <= (IFU_SIDE_ARCSEC / math.sqrt(2))
                    ),
                )
        else:
            if item["type"] in ("O", "S"):
                raise ValueError(f"Missing coordinates for assigned arm {index}")
        arms.append(item)
    if not any("icrs_deg" in arm for arm in arms):
        raise ValueError("No assigned KMOS arm centers in the ESO header")
    return arms


def measure(cache, item, target, center_sep):
    dataset = item["dp_id"]
    url = HEADER + dataset
    path = cache.fetch(
        url, f"metadata/{TAG}_{dataset.replace(':', '-')}.html", survey="ESO VLT KMOS",
        version="published raw FITS header", kind="detector_field_geometry",
    )
    header = header_data(path, dataset)
    arms = arm_coordinates(header, target)
    nearest = min((arm for arm in arms if "separation_arcsec" in arm), key=lambda arm: arm["separation_arcsec"])
    pointing = SkyCoord(header["RA"], header["DEC"], unit="deg", frame=FK5(equinox=Time("J2000")))
    return dict(
        dp_id=dataset, programme=item["prog_id"], center_separation_arcsec=center_sep,
        date_obs_utc=header["DATE-OBS"], exposure_seconds=header["EXPTIME"],
        target=header["OBJECT"], observation_block_id=header["ESO OBS ID"],
        exposure_sequence_number=header["ESO TPL EXPNO"],
        telescope_pointing_icrs_deg=[float(pointing.icrs.ra.deg), float(pointing.icrs.dec.deg)],
        nearest_arm=nearest, arms=arms, header_receipt=receipt(path, url),
    )


def main():
    path = OUT / f"{INVENTORY}.json"
    incoming = json.loads(path.read_text(encoding="utf-8"))
    position = incoming["position_icrs_deg"]
    target = SkyCoord(*position, unit="deg")
    rows = incoming["eso"]["raw_within_center_rectangle"]["rows"]
    cache = DownloadCache(ROOT, BUDGET)
    centers = []
    for row in (row for row in rows if row["instrument"] == "KMOS"):
        pointing = SkyCoord(float(row["ra"]), float(row["dec"]), unit="deg", frame=FK5(equinox=Time("J2000")))
        key = (row["prog_id"], row["date_obs"][:10])
        centers.append((key, row, float(pointing.icrs.separation(target).arcsec)))
    # Keep an additional 2 arcsec at the patrol boundary for the IFU diagonal and rounding.
    nights = {key for key, _, separation in centers if separation <= PATROL_RADIUS_ARCSEC + 2}
    # An observed programme-night includes the entire object and sky nod sequence.
    measured = [measure(cache, row, target, separation) for key, row, separation in centers if key in nights]
    other = [
        dict(dp_id=row["dp_id"], separation_arcsec=separation)
        for key, row, separation in centers if key not in nights
    ]
    if not measured:
        raise ValueError("No KMOS headers tested inside the patrol region")
    earliest = min(datetime.fromisoformat(row["date_obs_utc"]).replace(tzinfo=timezone.utc) for row in measured)
    last = max(datetime.fromisoformat(row["date_obs_utc"]).replace(tzinfo=timezone.utc) for row in measured)
    mode_counts = Counter(
        "object" if any(arm["type"] == "O" for arm in row["arms"]) else "sky only" for row in measured
    )
    nearest = min(measured, key=lambda row: row["nearest_arm"]["separation_arcsec"])
    best_distance = nearest["nearest_arm"]["separation_arcsec"]
    print(f"KMOS checked {len(measured)} headers, closest IFU center {best_distance:.2f} arcsec.", flush=True)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__), inventory_sha256=text_digest(path),
        position_icrs_deg=position, patrol_radius_arcsec=PATROL_RADIUS_ARCSEC,
        ifu_side_arcsec=IFU_SIDE_ARCSEC, first_pointing_utc=earliest.isoformat(), last_pointing_utc=last.isoformat(),
        mode_counts=dict(mode_counts), closest_frame=nearest["dp_id"], closest_arm=nearest["nearest_arm"],
        checked_programme_nights=[list(key) for key in sorted(nights)],
        checked_headers=measured, other_pointing_rows=other,
        scope="Tests assigned object and sky arm centers in public FITS headers, converted from "
        "FK5 J2000 to ICRS. The circumscribed circle bounds a 2.8 by 2.8 arcsec square "
        "at any rotation about its nominal center. These are geometric checks; arm coordinates "
        "are nominal and do not certify calibrated IFU spaxel coverage or a spectral sensitivity.",
    ))


if __name__ == "__main__":
    main()

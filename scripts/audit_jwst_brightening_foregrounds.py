"""Check roll-angle and optical foreground alternatives for the GOODS-S event."""

import json
from io import StringIO
from urllib.parse import urlencode

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import fetch
from measure_jwst_epoch_candidate import verified_path

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

STORAGE_MIB = 25_500
ARCSEC = 60


def position_angles(records):
    angles = {}
    for entry in records:
        start, science = entry["ranges"][:2]
        first = fits.Header.fromstring(verified_path(start).read_bytes().decode("ascii", errors="replace"))
        second = fits.Header.fromstring(verified_path(science).read_bytes().decode("ascii", errors="replace"))
        if first["FILENAME"] != entry["filename"] or second["EXTNAME"] != "SCI":
            raise ValueError("Unexpected stage-3 header order")
        angles.setdefault(entry["filter"], {})[entry["epoch"]] = dict(
            pa_v3_deg=float(second["PA_V3"]), science_header=science, midpoint_mjd=float(first["EXPMID"]),
        )
    return angles


def gaia(cache, sky, epoch):
    params = {
        "-source": "I/355/gaiadr3",
        "-c": f"{sky.ra.deg:.9f} {sky.dec.deg:.9f}", "-c.rs": ARCSEC, "-out.max": 200,
        "-out": "Source,RA_ICRS,DE_ICRS,Gmag,Plx,pmRA,pmDE",
    }
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"+urlencode(params)
    product = fetch(cache, url, "module_b_001_gaia_dr3.tsv", "public_foreground_catalogue",
                    survey="CDS VizieR", version="Gaia DR3 I/355/gaiadr3")
    text = verified_path(product).read_text(encoding="utf-8")
    if "#Name: I/355/gaiadr3" not in text or "#Truncated" in text:
        raise ValueError("Unexpected Gaia DR3 query response")
    rows = pd.read_csv(StringIO(text), comment="#", sep="\t", dtype=str)
    for field in ("RA_ICRS", "DE_ICRS", "Gmag", "Plx", "pmRA", "pmDE"):
        rows[field] = pd.to_numeric(rows[field], errors="coerce")
    rows = rows.dropna(subset=["RA_ICRS", "DE_ICRS"]).copy()
    coords = SkyCoord(rows.RA_ICRS, rows.DE_ICRS, unit="deg")
    rows["separation_arcsec_at_catalogue_epoch"] = sky.separation(coords).arcsec
    has_motion = rows[["pmRA", "pmDE"]].notna().all(axis=1).to_numpy()
    elapsed = (epoch - Time(2016., format="jyear")).to_value(u.yr)
    projected = coords[has_motion].spherical_offsets_by(
        rows.loc[has_motion, "pmRA"].to_numpy()*elapsed*u.mas,
        rows.loc[has_motion, "pmDE"].to_numpy()*elapsed*u.mas,
    )
    rows["separation_arcsec_propagated"] = np.nan
    rows.loc[has_motion, "separation_arcsec_propagated"] = sky.separation(projected).arcsec
    return dict(
        catalogue=product, radius_arcsec=ARCSEC, catalogue_epoch_jyear=2016.,
        propagated_epoch_utc=epoch.isot, rows=json.loads(rows.to_json(orient="records", double_precision=10)),
        limiting_scope="This cone includes only Gaia DR3 detections selected at catalogue coordinates."
                       " The propagation is linear and omits parallax and radial velocity.",
    )


def nearby_small_bodies(cache, sky, times):
    parameters = dict(COMMAND="'-170'", CENTER="'500@399'", EPHEM_TYPE="VECTORS",
                      TLIST="'"+",".join(f"{jd:.11f}" for jd in times.tdb.jd)+"'",
                      OUT_UNITS="KM-S", REF_PLANE="FRAME", REF_SYSTEM="ICRF",
                      VEC_TABLE=2, CSV_FORMAT="YES", OBJ_DATA="NO", format="json")
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(parameters)
    state_receipt = fetch(cache, url, "module_b_001_jwst_geocentric_states.json", "public_jwst_ephemeris",
                          survey="NASA JPL Horizons", version="ICRF Earth-centered vectors")
    payload = json.loads(verified_path(state_receipt).read_text())
    result = payload["result"]
    if "(-170)" not in result or "Earth (399)" not in result or "Reference frame : ICRF" not in result:
        raise ValueError("Unexpected observer state metadata")
    lines = result.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    ra = sky.ra.to_string(unit=u.hourangle, sep="-", precision=2, pad=True)
    dec = "M"+(-sky.dec).to_string(unit=u.deg, sep="-", precision=2, pad=True)
    searches = []
    for i, (when, line) in enumerate(zip(times, lines, strict=True), 1):
        date, _calendar, *vector = line.split(",")[:8]
        if abs(float(date)-when.tdb.jd) > 1/86400:
            raise ValueError("Horizons epoch differs from the requested time")
        state = [float(x) for x in vector]
        args = {"xobs": ",".join(f"{x:.12g}" for x in state), "obs-time": f"{when.jd:.11f}",
                "fov-ra-center": ra, "fov-dec-center": dec, "fov-ra-hwidth": .025, "fov-dec-hwidth": .025,
                "two-pass": "true", "mag-required": "true", "vmag-lim": 28.5}
        url = "https://ssd-api.jpl.nasa.gov/sb_ident.api?"+urlencode(args)
        receipt = fetch(cache, url, f"module_b_001_small_bodies_{i}.json", "public_small_body_search",
                        survey="NASA JPL Small-Body Identification", version="Webb state at October mosaic midpoint")
        data = json.loads(verified_path(receipt).read_text())
        if data["signature"]["version"] != "1.1" or data["n_second_pass"] != len(data.get("data_second_pass", [])):
            raise ValueError("Unexpected JPL identification result")
        searches.append(dict(utc=when.isot, utc_jd=float(when.jd), geocentric_state_km_km_s=state,
                             catalogue=receipt, first_pass_rows=data["n_first_pass"],
                             integrated_matches=data.get("data_second_pass", [])))
        print(when.isot, "integrated JPL matches", data["n_second_pass"], flush=True)
    return dict(horizons=state_receipt, field_half_width_degrees=.025, visual_magnitude_limit=28.5,
                searches=searches, scope="Known small bodies with magnitude parameters at two midpoints."
                " A visual magnitude cut does not set a completeness limit for near-infrared objects.")


def main():
    source = OUT/"module_b_001_november_photometry.json"
    report = json.loads(source.read_text())
    sky = SkyCoord(*report["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    angles = position_angles(report["official_mosaics"])
    gaia_rows = gaia(cache, sky, Time(angles["F200W"]["Oct 2026"]["midpoint_mjd"], format="mjd", scale="utc"))
    for filt, epochs in angles.items():
        first, second = (epochs[epoch]["pa_v3_deg"] for epoch in ("Oct 2026", "Nov 2026"))
        print(filt, "orientation change (deg)", (second-first+180)%360-180, flush=True)
    bright = Time([angles[band]["Oct 2026"]["midpoint_mjd"] for band in ("F200W", "F150W")],
                  format="mjd", scale="utc")
    solar_system = nearby_small_bodies(cache, sky, bright)
    write_json(OUT/"module_b_001_orientations_and_gaia.json", dict(
        executed_utc=utcnow(), icrs_deg=report["icrs_deg"], visit_orientations=angles, gaia=gaia_rows, small_body_identification=solar_system,
        november_sha256=text_digest(source), script_sha256=text_digest(__file__),
    ))


if __name__ == "__main__":
    main()

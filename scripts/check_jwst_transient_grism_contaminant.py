"""Test whether the shared detector trace includes a known galaxy's emission."""

import csv
import json
from urllib.parse import urlencode

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_transient_grism_geometry import STORAGE_MIB
from extract_jwst_transient_grism import extract_column, flat_strip, trace_for_column
from measure_jwst_epoch_candidate import verified_path

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

CATALOGUE_ID = 1178
PA_ALPHA_VACUUM_UM = 1.87561
LINE_HALF_WIDTH_UM = .008
SIDEBAND_UM = (.013, .028)


def feature(wave, flux, center):
    """Use the local median sideband and report the strongest finite line sample."""
    offset = abs(wave-center)
    good = np.isfinite(flux)
    line = np.flatnonzero(good & (offset <= LINE_HALF_WIDTH_UM))
    side = flux[good & (offset >= SIDEBAND_UM[0]) & (offset <= SIDEBAND_UM[1])]
    measurement = dict(line_samples=len(line), sideband_samples=len(side), peak_um=None, peak_above_sideband=None)
    if len(line) >= 3 and len(side) >= 5:
        peak = line[np.argmax(flux[line])]
        measurement.update(peak_um=float(wave[peak]), peak_above_sideband=float(flux[peak]-np.median(side)))
    return measurement


def nirspec_reference(coords):
    query = API+"/nirspec_extractions?"+urlencode({
        "coords": f"{coords[0]:.9f},{coords[1]:.9f}", "size": .6, "output": "csv",
    })
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    catalogue = download(cache, query, "transient_grism_contaminant_nirspec.csv", "nirspec_spatial_query")
    with verified_path(catalogue).open(newline="", encoding="utf-8") as stream:
        (entry,) = [row for row in csv.DictReader(stream)
                    if (row["srcid"], row["grating"], row["grade"]) == ("217321", "G235M_F170LP", "3")]
    distance = SkyCoord(*coords, unit="deg").separation(
        SkyCoord(float(entry["ra"]), float(entry["dec"]), unit="deg")).arcsec
    if distance > .1:
        raise ValueError("NIRSpec spectrum falls outside the direct catalogue match radius")
    spec = download(cache, f'https://s3.amazonaws.com/msaexp-nirspec/extractions/{entry["root"]}/{entry["file"]}',
                    entry["file"], "observed_nirspec_contaminant_reference")
    expected = PA_ALPHA_VACUUM_UM*(1+float(entry["z"]))
    with fits.open(verified_path(spec)) as h:
        d = h["SPEC1D"].data
        flux = np.where(np.isfinite(d["err"]) & (d["err"] > 0), d["flux"], np.nan)
        result = feature(d["wave"], flux, expected)
    return expected, dict(query=catalogue, product=spec, catalogue_entry=entry,
                          separation_arcsec=float(distance), feature_flux_unit="microJy", **result)


def inspect_frame(frame, header, models, flats, target, neighbour, expected):
    path = verified_path({"path": header["x1d_path"], "sha256": header["x1d_sha256"]})
    with fits.open(path) as h:
        ext = h["EXTRACT1D"]
        if ext.header["SPORDER"] != 1 or ext.columns["FLUX"].unit != "Jy":
            raise ValueError("Expected first-order GRISMR flux in Jy")
        (row,) = ext.data[ext.data["SOURCE_ID"] == CATALOGUE_ID]
        count = int(row["N_ALONGDISP"])
        wave = row["WAVELENGTH"][:count]
        flux = np.where(row["DQ"][:count] == 0, row["FLUX"][:count]*1e6, np.nan)
        detector_x = np.arange(count)+int(row["EXTRACT2D_XSTART"])-1
        measurement = feature(wave, flux, expected)
    cols, target_wave, target_y, other_y = trace_for_column(header, models, target, neighbour)
    projected = float(np.interp(np.interp(expected, wave, detector_x), cols, target_wave, left=np.nan, right=np.nan))
    if not np.isfinite(projected):
        raise ValueError("Expected the contaminant line to project within the transient wavelength range")
    select = abs(target_wave-projected) <= SIDEBAND_UM[1]
    pixels, rows = flat_strip(frame, flats, cols[select], target_y[select])
    raw = np.array([extract_column(values, rows, y, neighbour_y, 1) for values, y, neighbour_y in
                    zip(pixels.T, target_y[select], other_y[select], strict=True)])
    blend = feature(target_wave[select], raw*1000, projected)
    return dict(filename=frame["filename"], archived_x1d_feature_microjy=measurement,
                projected_line_on_transient_wavelength_um=projected, forced_blend_feature_mdn_per_s=blend)


def main():
    filenames = {
        "rates": "field_b_transient_grism_rates.json", "headers": "field_b_transient_grism_products.json",
        "positions": "field_b_transient_wfss.json", "overlap": "field_b_transient_grism_overlap.json",
    }
    paths = {key: OUT / value for key, value in filenames.items()}
    data = {key: json.loads(path.read_text()) for key, path in paths.items()}
    rates, headers, positions, overlap = (data[key] for key in filenames)
    (match,) = [row for row in overlap["overlap_results"][0]["overlapping_sources"]
                if row["catalogue_label"] == CATALOGUE_ID and row["order"] == 1]
    expected, nirspec = nirspec_reference(match["icrs"])
    records = {r["rate_filename"]: r for r in headers["rate_headers"]}
    nearest = {r["filename"].replace("_x1d.fits", "_rate.fits"): r["nearest"]["icrs_deg"]
               for r in positions["x1d_exposures"]}
    results = []
    for frame in rates["native_rate_frames"]:
        name = frame["filename"]
        results.append(inspect_frame(frame, records[name], headers["calibrations"], rates["calibrations"],
                                     overlap["position_degrees"], nearest[name], expected))
    write_json(OUT / "field_b_transient_grism_contaminant.json", {
        "executed_utc": utcnow(), "candidate": overlap["candidate"], "catalogued_source": match,
        "pa_alpha_rest_vacuum_um": PA_ALPHA_VACUUM_UM, "predicted_contaminant_pa_alpha_um": expected,
        "line_half_width_um": LINE_HALF_WIDTH_UM, "sideband_offsets_um": list(SIDEBAND_UM),
        "nirspec": nirspec, "grism": results,
        "scope": "Peak samples relative to fixed sidebands, without line fitting or transient deblending",
        "input_sha256": {key: text_digest(path) for key, path in paths.items()}, "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

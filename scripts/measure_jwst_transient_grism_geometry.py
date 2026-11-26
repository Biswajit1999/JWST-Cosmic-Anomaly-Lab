"""Map the transient onto the bright-epoch grism detector using CRDS models."""

import json
from contextlib import ExitStack
from io import BytesIO

import asdf
import numpy as np
from astropy.constants import c
from astropy.io import fits
from astropy.modeling import models
from astropy.wcs import WCS
from download_jwst_cross_detector import OUT, ROOT
from gwcs.geometry import CartesianToSpherical, SphericalToCartesian
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import (
    format_flux,
    photometry,
    pixel_scale,
    show_stamp,
    verified_path,
)
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

NOMINAL = (3.14, 3.98)
REFERENCE_PIXEL_EDGE = 4
C_M_S = c.to_value("m/s")


def pointing_header(path):
    raw = BytesIO(path.read_bytes())
    primary = fits.Header.fromfile(raw)
    if primary["NAXIS"] != 0:
        raise ValueError(f"Unexpected primary data array in {path}")
    science = fits.Header.fromfile(raw)
    if science["EXTNAME"] != "SCI":
        raise ValueError(f"Expected SCI as the first extension in {path}")
    return primary, science


def imaging_transform(science, distortion, offset):
    v2, v3 = science["V2_REF"], science["V3_REF"]
    scale = science["VA_SCALE"]
    aberration = ((models.Scale(scale) | models.Shift((1-scale)*v2))
                  & (models.Scale(scale) | models.Shift((1-scale)*v3)))
    angles = [v2/3600, -v3/3600, science["ROLL_REF"], science["DEC_REF"], -science["RA_REF"]]
    attitude = (models.Scale(1/3600) & models.Scale(1/3600)
                | SphericalToCartesian(wrap_lon_at=180)
                | models.RotationSequence3D(angles, axes_order="zyxyz")
                | CartesianToSpherical(wrap_lon_at=360))
    filter_shift = models.Shift(offset["column_offset"]) & models.Shift(offset["row_offset"])
    return filter_shift | distortion | aberration | attitude


def dispersed_trace(reference, order, x0, y0, velocity):
    """Use field-dependent trace coefficients and a dense monotone parameter grid."""
    t = np.linspace(0., 1., 8193)
    j = reference["orders"].index(order)
    wave = np.polynomial.polynomial.polyval(t, [float(p(x0, y0)) for p in reference["displ"][j]])
    x_models, y_models = reference["dispx"][j], reference["dispy"][j]
    if len(x_models) != 1 or x_models[0].n_inputs != 1:
        raise ValueError("Expected a GRISMR single along-dispersion polynomial")
    x = x0 + x_models[0](t)
    y = y0 + np.polynomial.polynomial.polyval(t, [float(p(x0, y0)) for p in y_models])
    dx = np.diff(x)
    if not ((np.all(dx > 0) or np.all(dx < 0)) and np.all(np.diff(wave) > 0)):
        raise ValueError("Expected monotone dispersion and wavelength across the reference trace")
    return x, y, wave/(1+velocity/C_M_S)


def validate_pointing(record, transform, dispersion, velocity):
    path = ROOT / record["x1d_path"]
    if digest(path) != record["x1d_sha256"]:
        raise ValueError(f"Archived x1d checksum mismatch for {path.name}")
    position_errors, wavelength_errors = [], []
    with fits.open(path) as h:
        for r in h["EXTRACT1D"].data:
            ra, dec, cat_x, cat_y = (float(r[name]) for name in (
                "SOURCE_RA", "SOURCE_DEC", "SOURCE_XPOS", "SOURCE_YPOS"))
            if not np.isfinite(ra+dec+cat_x+cat_y):
                continue
            x0, y0 = transform.inverse(ra, dec)
            position_errors.append(float(np.hypot(x0-cat_x, y0-cat_y)))
            wave = r["WAVELENGTH"][:int(r["N_ALONGDISP"])]
            k = np.arange(len(wave), dtype=float)+int(r["EXTRACT2D_XSTART"])-1
            x, _, model_wave = dispersed_trace(dispersion, 1, x0, y0, velocity)
            if x[0] > x[-1]:
                x, model_wave = x[::-1], model_wave[::-1]
            expected = np.interp(k, x, model_wave, left=np.nan, right=np.nan)
            valid = np.isfinite(wave+expected) & (wave > 0)
            wavelength_errors.extend(abs(wave[valid]-expected[valid]))
    if not position_errors or not wavelength_errors:
        raise ValueError(f"No usable dispersion validation rows in {path.name}")
    return dict(catalogue_positions=len(position_errors), position_error_median_pix=float(np.median(position_errors)),
                position_error_max_pix=float(np.max(position_errors)),
                wavelength_sample_count=len(wavelength_errors),
                absolute_wavelength_error_median_um=float(np.median(wavelength_errors)),
                absolute_wavelength_error_p95_um=float(np.percentile(wavelength_errors, 95)))


def coverage(wave, x, y, nx, ny, low, high):
    mask = ((wave >= low) & (wave <= high) & (x >= REFERENCE_PIXEL_EDGE) & (x <= nx-1-REFERENCE_PIXEL_EDGE)
            & (y >= REFERENCE_PIXEL_EDGE+3) & (y <= ny-1-REFERENCE_PIXEL_EDGE-3))
    if not mask.any():
        return dict(first_um=None, last_um=None, fraction_of_band=0.)
    first, last = wave[mask][[0, -1]]
    return dict(first_um=float(first), last_um=float(last), fraction_of_band=float((last-first)/(high-low)))


def evaluate_exposure(record, calibrations, sky):
    primary, science = pointing_header(verified_path(record))
    if (primary["FILTER"], primary["PUPIL"], primary["DETECTOR"]) != ("F356W", "GRISMR", "NRCALONG"):
        raise ValueError("Unexpected grism pointing configuration")
    with ExitStack() as stack:
        references = {key: stack.enter_context(asdf.open(verified_path(calibrations[name])))
                      for key, name in record["references"].items()}
        (offset,) = [item for item in references["R_FILOFF"]["filters"]
                     if (item["filter"], item["pupil"]) == ("F356W", "GRISMR")]
        transform = imaging_transform(science, references["R_DISTOR"]["model"], offset)
        x0, y0 = transform.inverse(*sky)
        velocity = science["VELOSYS"]
        x, y, wave = dispersed_trace(references["R_SPCWCS"], 1, x0, y0, velocity)
        validation = validate_pointing(record, transform, references["R_SPCWCS"], velocity)
        ranges = references["R_WAVRAN"]["wavelengthrange"]
        (_, _, low, high), = [item for item in ranges if item[:2] == [1, "F356W"]]
    nx, ny = primary["SUBSIZE1"], primary["SUBSIZE2"]
    return dict(
        filename=record["rate_filename"], date_utc=f'{primary["DATE-OBS"]}T{primary["TIME-OBS"]}',
        date_mjd=primary["EXPSTART"], estimated_direct_center_pix=[float(x0), float(y0)],
        nominal=coverage(wave, x, y, nx, ny, *NOMINAL), crds=coverage(wave, x, y, nx, ny, low, high),
        maximum_barycentric_window_um=[low, high], astronomical_pointing_validation=validation,
    )


def direct_cutout(record, sky):
    with fits.open(verified_path(record), memmap=True, lazy_load_hdus=True) as h:
        science = h["SCI"]
        if h[0].header["PUPIL"] != "CLEAR" or science.header["BUNIT"] != "MJy/sr":
            raise ValueError("Expected calibrated F356W direct imaging in MJy/sr")
        wcs = WCS(science.header).celestial
        x, y = map(float, wcs.world_to_pixel_values(*sky))
        scale = pixel_scale(wcs)
        half = int(np.ceil(8/scale))
        x0, y0 = round(x)-half, round(y)-half
        x1, y1 = x0+2*half+1, y0+2*half+1
        xa, ya = max(0, x0), max(0, y0)
        xb, yb = min(science.header["NAXIS1"], x1), min(science.header["NAXIS2"], y1)
        image = np.full((2*half+1, 2*half+1), np.nan, dtype=np.float64)
        if xb > xa and yb > ya:
            pixels = np.s_[ya:yb, xa:xb]
            values = science.data[pixels]*(float(science.header["PIXAR_SR"])*1e15)
            valid = np.isfinite(h["WHT"].data[pixels]) & (h["WHT"].data[pixels] > 0)
            image[ya-y0:yb-y0, xa-x0:xb-x0] = np.where(valid, values, np.nan)
        date = f'{h[0].header["DATE-OBS"]}T{h[0].header["TIME-OBS"]}'
    xx, yy = x-x0, y-y0
    return (image, xx, yy, scale), {
        "filename": record["path"].rsplit("/", 1)[-1], "date_utc": date, "scale_arcsec": scale,
        "apertures": photometry(image, xx, yy, scale, (.12, .18, .24), (.30, .48)),
    }


def plot_result(exposures, cutouts, measurements):
    fig = Figure(figsize=(10, 6), layout="constrained")
    grid = fig.add_gridspec(2, len(cutouts), height_ratios=(1.05, 1))
    geo = fig.add_subplot(grid[0, :])
    geo.axvspan(*NOMINAL, color="#e7e9f5", label="Nominal F356W passband")
    for i, row in enumerate(exposures, 1):
        limits = row["crds"]
        if limits["first_um"] is not None:
            geo.hlines(i, limits["first_um"], limits["last_um"], linewidth=3, color="#84858f")
        nominal = row["nominal"]
        if nominal["first_um"] is not None:
            geo.hlines(i, nominal["first_um"], nominal["last_um"], linewidth=3.5, color="#227a9b")
    geo.set(ylim=(.5, len(exposures)+.5), xlim=(3.10, 4.105), xlabel="Observed wavelength (µm)",
            ylabel="GRISMR exposures in time order", yticks=range(1, len(exposures)+1, 2))
    geo.legend(loc="upper left", fontsize=8)
    nearby = np.concatenate([im[round(y)-15:round(y)+16, round(x)-15:round(x)+16].ravel()
                             for im, x, y, _ in cutouts])
    near = nearby[np.isfinite(nearby)]
    upper = float(np.percentile(near, 99)) if len(near) else 1.
    for i, ((im, x, y, scale), row) in enumerate(zip(cutouts, measurements, strict=True)):
        ax = fig.add_subplot(grid[1, i])
        show_stamp(ax, im, x, y, scale, (-.10, upper), "gray", .24, half_arcsec=.64)
        ax.set_title(f'{row["date_utc"][:16]} · {format_flux(row["apertures"][-1])}', fontsize=8)
    fig.suptitle("GOODS-S transient · grism trace coverage and observed F356W images\n"
                 "Predicted detector coverage above · fixed-aperture imaging below", fontsize=12)
    save_figure(fig, FIGURES / "field_b_transient_grism_geometry")


def main():
    path = OUT / "field_b_transient_grism_products.json"
    persisted = json.loads(path.read_text())
    sky = persisted["position_degrees"]
    exposures = []
    for record in persisted["rate_headers"]:
        row = evaluate_exposure(record, persisted["calibrations"], sky)
        exposures.append(row)
        print(row["filename"], row["nominal"], row["astronomical_pointing_validation"]["position_error_max_pix"],
              flush=True)
    direct = [direct_cutout(record, sky) for record in persisted["direct_images"]]
    candidate = json.loads((OUT / "field_b_candidate_products.json").read_text())
    for product in (p for p in candidate["survey_visits"] if p["band"] == "F356W-CLEAR"):
        (im, wcs, x, y, _, _), row = vet_mosaic(product, sky)
        scale = pixel_scale(wcs)
        direct.append(((im, x, y, scale), {**row, "date_utc": row["dateobs"], "scale_arcsec": scale}))
    stamps, measurements = zip(*direct, strict=True)
    plot_result(exposures, stamps, measurements)
    write_json(OUT / "field_b_transient_grism_geometry.json", {
        "executed_utc": utcnow(), "candidate": persisted["candidate"], "position_degrees": sky,
        "nominal_passband_um": NOMINAL, "reference_pixel_edge": REFERENCE_PIXEL_EDGE,
        "exposure_coverage": exposures, "direct_f356w_photometry": measurements,
        "scope": "Predicted geometric trace inside detector science area; throughput and DQ not evaluated",
        "products_sha256": text_digest(path),
        "comparison_mosaics_sha256": text_digest(OUT / "field_b_candidate_products.json"),
        "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

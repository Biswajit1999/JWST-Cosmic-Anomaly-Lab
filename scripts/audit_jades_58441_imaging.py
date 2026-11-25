"""Measure the short-wavelength excess near GOODS-N 58441 in released JWST images."""

import json
from urllib.parse import urlencode

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.colors import AsinhNorm
from matplotlib.figure import Figure
from matplotlib.patches import Circle

from cosmic_anomaly_lab.imaging import aperture_measurement
from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import check_registration, save_figure

TARGET = "goods-n-mediumjwst_58441"
BANDS = ("F090W-CLEAR", "F115W-CLEAR", "F150W-CLEAR")
C_A_S = 2.99792458e18
RADII = (0.12, 0.15, 0.20, 0.30)
CONTROL_RADIUS = 0.15
CONTROL_OFFSETS = (-0.78, -0.45, 0.45, 0.78)
AXIS_ANNULUS = (0.30, 0.85)
THROUGHPUT = (
    "https://jwst-docs.stsci.edu/files/216457506/216457511/1/1762454197925/"
    "F090W_mean_system_throughput.txt"
)


def imaging(path, source):
    rows, panels = [], []
    with fits.open(path) as hdus:
        hdus.verify("exception")
        for band in BANDS:
            science, weight = (
                next(
                    h
                    for h in hdus
                    if h.header.get("FILTER") == band and h.header.get("EXTVER") == ext
                )
                for ext in ("SCI", "WHT")
            )
            wcs = WCS(science.header).celestial
            data = np.asarray(science.data, dtype=float)
            wht = np.asarray(weight.data, dtype=float)
            valid = np.isfinite(data) & np.isfinite(wht) & (wht > 0)
            conversion = science.header["PHOTFNU"]
            photometry = [
                aperture_measurement(
                    data,
                    wht,
                    wcs,
                    source.ra.deg,
                    source.dec.deg,
                    jy_per_unit=conversion,
                    radius_arcsec=radius,
                )
                for radius in RADII
            ]
            x0, y0 = wcs.world_to_pixel(source)
            yy, xx = np.indices(data.shape)
            distance = np.hypot(xx - x0, yy - y0) * photometry[0]["pixel_scale_arcsec"]
            background = photometry[0]["background_uJy_pixel"]
            calibrated = np.where(valid, data * conversion * 1e6 - background, np.nan)
            centres = []
            for radius in RADII:
                mask = valid & (distance <= radius)
                weight_pos = np.maximum(calibrated[mask], 0)
                xc = float(np.average(xx[mask], weights=weight_pos))
                yc = float(np.average(yy[mask], weights=weight_pos))
                pos = wcs.pixel_to_world(xc, yc)
                east, north = source.spherical_offsets_to(pos)
                centres.append(
                    dict(
                        radius_arcsec=radius,
                        x_pixel=xc,
                        y_pixel=yc,
                        east_mas=float(east.mas),
                        north_mas=float(north.mas),
                        ra=float(pos.ra.deg),
                        dec=float(pos.dec.deg),
                    )
                )
            registration = check_registration(panels[0]["wcs"], wcs, data.shape) if panels else 0.0
            panels.append(
                dict(
                    wcs=wcs,
                    data=calibrated,
                    distance=distance,
                    valid=valid,
                    science=data,
                    weight=wht,
                    centre=(x0, y0),
                )
            )
            rows.append(
                dict(
                    band=band,
                    apertures=photometry,
                    centres=centres,
                    photfnu=float(conversion),
                    wcs_grid_mismatch_pixels=registration,
                )
            )
    return rows, panels


def streak_controls(panels, rows):
    first = panels[0]
    x0, y0 = first["centre"]
    yy, xx = np.indices(first["data"].shape)
    lo, hi = AXIS_ANNULUS
    annulus = first["valid"] & (first["distance"] >= lo) & (first["distance"] <= hi)
    vectors = np.column_stack((xx[annulus] - x0, yy[annulus] - y0))
    weights = np.maximum(first["data"][annulus], 0)
    moment = (vectors.T * weights) @ vectors / weights.sum()
    eigenvalues, eigenvectors = np.linalg.eigh(moment)
    parallel = eigenvectors[:, -1]
    directions = {"along": parallel, "across": np.array([-parallel[1], parallel[0]])}
    for panel, row in zip(panels, rows, strict=True):
        scale = row["apertures"][1]["pixel_scale_arcsec"]
        controls = {}
        for label, direction in directions.items():
            points = []
            for offset in CONTROL_OFFSETS:
                x, y = panel["centre"] + direction * offset / scale
                position = panel["wcs"].pixel_to_world(x, y)
                result = aperture_measurement(
                    panel["science"],
                    panel["weight"],
                    panel["wcs"],
                    position.ra.deg,
                    position.dec.deg,
                    jy_per_unit=row["photfnu"],
                    radius_arcsec=CONTROL_RADIUS,
                )
                points.append(
                    dict(
                        offset_arcsec=offset,
                        x_pixel=float(x),
                        y_pixel=float(y),
                        flux_njy=result["aperture_flux_uJy"] * 1000,
                        diagonal_error_njy=result["diagonal_error_uJy"] * 1000,
                    )
                )
            controls[label] = points
        along = np.array([p["flux_njy"] for p in controls["along"]])
        paired_means = [float((along[0] + along[3]) / 2), float((along[1] + along[2]) / 2)]
        target_flux = row["apertures"][1]["aperture_flux_uJy"] * 1000
        row["streak_control"] = dict(
            directional_apertures=controls,
            along_pair_means_njy=paired_means,
            target_minus_along_pair_njy=[target_flux - flux for flux in paired_means],
        )
    return dict(
        major_axis_xy=parallel.tolist(),
        moment_eigenvalue_ratio=float(eigenvalues[-1] / eigenvalues[0]),
        axis_fit_annulus_arcsec=list(AXIS_ANNULUS),
        aperture_radius_arcsec=CONTROL_RADIUS,
    )


def line_budget(path, scaling):
    wave_micron, response = np.loadtxt(path, skiprows=1, unpack=True)
    if not np.all(np.diff(wave_micron) > 0) or np.any(response < 0):
        raise ValueError("Invalid NIRCam throughput curve")
    reference = json.loads(scaling.read_text())
    line_a = reference["oxygen_doublet_alternative"]["main_wavelength_a"]
    line_response = np.interp(line_a / 1e4, wave_micron, response)
    # Photon-weighted equivalent flat Fnu for a monochromatic, integrated line flux.
    bandpass = np.trapezoid(response / wave_micron, wave_micron)
    n_jy_per_flux_unit = 1e-20 * line_response * line_a / C_A_S / bandpass * 1e32
    records = []
    for row in reference["comparisons"]:
        for reduction, measurement in zip(
            ("JADES", "DAWN"), row["released_line_fluxes"], strict=True
        ):
            records.append(
                dict(
                    mode=row["mode"],
                    reduction=reduction,
                    line_equivalent_f090w_njy=float(measurement["flux"] * n_jy_per_flux_unit),
                    diagonal_error_njy=float(measurement["error"] * n_jy_per_flux_unit),
                )
            )
    return dict(
        wavelength_a=line_a,
        throughput=float(line_response),
        logarithmic_bandpass=float(bandpass),
        records=records,
    )


def draw(panels, rows, source):
    fig = Figure(figsize=(9.5, 8.5))
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.885, wspace=0.38, hspace=0.25)
    axes = fig.subplots(2, 2).ravel()
    for img, row, ax in zip(panels, rows, axes[:3], strict=True):
        high = float(np.nanpercentile(img["data"][img["distance"] < 1.35], 99.3))
        ax.imshow(
            img["data"],
            origin="lower",
            cmap="gray",
            norm=AsinhNorm(linear_width=high / 12, vmin=-high / 20, vmax=high),
        )
        xc, yc = img["wcs"].world_to_pixel(source)
        ax.set(
            xlim=(xc - 27, xc + 27),
            ylim=(yc - 27, yc + 27),
            title=row["band"].split("-")[0],
            xticks=[],
            yticks=[],
        )
        ax.add_patch(Circle((xc, yc), 3.0, fill=False, lw=1.2, edgecolor="#68e0ed"))
        ax.plot(xc - 22 + np.array([0, 20]), [yc - 22, yc - 22], color="white", lw=1.6)
        ax.text(xc - 22, yc - 20, "1″", color="white", fontsize=8)
        for label, color in (("along", "#f2b45e"), ("across", "#c3a9ef")):
            for point in row["streak_control"]["directional_apertures"][label]:
                ax.add_patch(
                    Circle(
                        (point["x_pixel"], point["y_pixel"]), 3.0, fill=False, lw=1, edgecolor=color
                    )
                )
    phot = axes[3]
    for row, color in zip(rows, ("#865a1b", "#267184", "#775da2"), strict=True):
        measurements = row["apertures"]
        phot.errorbar(
            RADII,
            [m["aperture_flux_uJy"] * 1000 for m in measurements],
            yerr=[m["diagonal_error_uJy"] * 1000 for m in measurements],
            marker="o",
            color=color,
            label=row["band"].split("-")[0],
            capsize=3,
        )
    residuals = rows[0]["streak_control"]["target_minus_along_pair_njy"]
    phot.scatter(
        [CONTROL_RADIUS] * len(residuals),
        residuals,
        marker="x",
        c="#222222",
        label="F090W after streak subtraction",
    )
    phot.set(xlabel="Aperture radius (arcsec)", ylabel="Aperture flux (nJy)", xticks=RADII)
    phot.legend(frameon=False, loc="upper left")
    fig.suptitle("GOODS-N 58441 · F090W streak and control apertures", fontsize=14, y=0.97)
    fig.text(
        0.5,
        0.924,
        "Cyan target  ·  Amber along streak  ·  Purple across streak",
        ha="center",
        fontsize=9,
    )
    save_figure(fig, ROOT / "figures/jades_58441_optical_ir_context")


def main():
    target = next(r for r in json.loads(LEDGER.read_text())["records"] if r["dr4_id"] == TARGET)
    source = SkyCoord(target["ra"], target["dec"], unit="deg")
    query = urlencode(
        dict(
            ra=source.ra.deg,
            dec=source.dec.deg,
            size=5,
            filters=("f775w-clear,f814w-clear,f090w-clear,f115w-clear,f150w-clear"),
            output="fits_weight",
        )
    )
    url = "https://grizli-cutout.herokuapp.com/thumb?" + query
    cache = DownloadCache(ROOT, 28_000 * 1024**2)
    cutout = cache.fetch(
        url,
        "images/jades_58441_optical_ir_context.fits",
        survey="DAWN JWST Archive",
        version="DJA-cutout-live",
        kind="observed_image_cutout",
    )
    response = cache.fetch(
        THROUGHPUT,
        "metadata/jwst_nircam_f090w_throughput.txt",
        survey="STScI",
        version="JDox NIRCam throughput v7.0",
        kind="instrument_reference",
    )
    scaling = OUT / "58441_flux_scale.json"
    rows, panels = imaging(cutout, source)
    controls = streak_controls(panels, rows)
    budget = line_budget(response, scaling)
    draw(panels, rows, source)
    write_json(
        OUT / "58441_optical_ir_context.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            input_images=receipt(cutout, url),
            throughput_curve=receipt(response, THROUGHPUT),
            line_fit_sha256=text_digest(scaling),
            bands=rows,
            streak_axis=controls,
            line_budget=budget,
            scope="The returned image bands are JWST NIRCam only. Flux errors use diagonal image weights. "
            "Centroids are positive-flux first moments with no cross-band astrometric calibration. "
            "Apertures do not match PSFs or deblend projected components.",
        ),
    )
    for row in rows:
        print(row["band"], row["centres"][1], row["apertures"][1]["aperture_flux_uJy"])
    print(controls)
    print(budget)


if __name__ == "__main__":
    main()

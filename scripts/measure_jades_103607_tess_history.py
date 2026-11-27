"""Measure historical optical variability on public TESS full-frame images."""

import json
from urllib.parse import urljoin
from zipfile import ZipFile

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from download_jades_103607_tess_history import TAG as INVENTORY_TAG
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.dates import AutoDateLocator, DateFormatter
from matplotlib.figure import Figure
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_tess_history"
NATIVE = OUT / "module_b_001_native_november_fade.json"
PRF = "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0004/"
DURATIONS = (1, 24)
OFFSETS = [(0, 0), *((x, y) for y in (-4, 0, 4) for x in (-4, 0, 4) if (x, y) != (0, 0))]
ZEROPOINT = 20.44


def read_cutout(product, coord):
    source = product["cutout"]
    path = ROOT / source["path"]
    if digest(path) != source["sha256"]:
        raise ValueError(f"TESS cutout does not match the recorded checksum: {path.name}")
    with ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.endswith(".fits")]
        if len(members) != 1:
            raise ValueError(f"Expected one TESSCut FITS file: {path.name}")
        with archive.open(members[0]) as stream, fits.open(stream, memmap=False) as hdu:
            table = hdu["PIXELS"].data
            header = hdu["PIXELS"].header.copy()
            x, y = (float(value) for value in WCS(hdu["APERTURE"].header).world_to_pixel(coord))
            if header["TIMESYS"] != "TDB" or header["TIMEREF"] != "SOLARSYSTEM":
                raise ValueError("Unexpected TESS FFI time convention")
            raw_time = np.asarray(table["TIME"], dtype=float)
            correction = np.asarray(table["TIMECORR"], dtype=float)
            flux = np.asarray(table["FLUX"], dtype=float)
            err = np.asarray(table["FLUX_ERR"], dtype=float)
            selected = (
                (table["QUALITY"] == 0) & np.isfinite(raw_time) & np.isfinite(correction)
                & np.all(np.isfinite(flux) & np.isfinite(err) & (err > 0), axis=(1, 2))
            )
    observed = Time(
        header["BJDREFI"] + header["BJDREFF"], (raw_time - correction)[selected],
        format="jd", scale="tdb",
    ).utc
    if not len(observed) or np.any(np.diff(observed.mjd) <= 0):
        raise ValueError("Empty or nonmonotonic TESS time series")
    if (
        abs(observed.mjd[0] - Time(header["DATE-OBS"]).mjd) > 1
        or abs(observed.mjd[-1] - Time(header["DATE-END"]).mjd) > 1
    ):
        raise ValueError("Valid TESS cadences disagree with the FITS observing date range")
    return header, members[0], table.size, (x, y), observed.mjd, flux[selected], err[selected]


def bracketing(nodes, value):
    index = int(np.clip(np.searchsorted(nodes, value), 1, len(nodes) - 1))
    low, high = nodes[index - 1:index + 1]
    fraction = float(np.clip((value - low) / (high - low), 0, 1))
    return [(low, 1 - fraction), (high, fraction)]


def prf_design(cache, product, header, xy, shape):
    # Both the FFI RAWX coordinates and the model filenames include overscan.
    x, y = xy
    column = header["1CRV5P"] + (x + 1 - header["1CRP5P"]) * header.get("1CDLT5P", 1)
    row = header["2CRV5P"] + (y + 1 - header["2CRP5P"]) * header.get("2CDLT5P", 1)
    camera, ccd = product["camera"], product["ccd"]
    grid_rows = (1, 513, 1025, 1536, 2048)
    grid_columns = (45, 557, 1069, 1580, 2092)
    models = {}
    receipts = []
    # The archived start_s0004 models have a common file timestamp per camera.
    stamp = "2019107181902" if camera == 4 else "2019107181901"
    for gy, _ in bracketing(grid_rows, row):
        for gx, _ in bracketing(grid_columns, column):
            name = f"tess{stamp}-prf-{camera}-{ccd}-row{gy:04d}-col{gx:04d}.fits"
            url = urljoin(PRF, f"cam{camera}_ccd{ccd}/{name}")
            path = cache.fetch(
                url, f"calibration/tess/{name}", survey="NASA TESS / MAST",
                version="start_s0004 UPDATED_2.0", kind="point_response_function",
            )
            with fits.open(path, memmap=False) as hdu:
                source = hdu[0]
                if (
                    source.data.shape != (117, 117) or source.header["NSAMP"] != 9
                    or source.header["CCD_RREF"] != gy or source.header["CCD_CREF"] != gx
                    or source.header["CAM"] != camera or source.header["CCD"] != ccd
                ):
                    raise ValueError(f"Inconsistent TESS PRF: {name}")
                models[gy, gx] = source.data.astype(float)
            receipts.append(receipt(path, url))
    yy, xx = np.indices(shape, dtype=float)
    profiles = []
    for dx, dy in OFFSETS:
        # The same four calibration nodes surround the small cutout. Clip at grid edges.
        model = sum(
            models[gy, gx] * wy * wx
            for gy, wy in bracketing(grid_rows, row + dy)
            for gx, wx in bracketing(grid_columns, column + dx)
        )
        native = map_coordinates(
            model, [58 + 9 * (yy - y - dy), 58 + 9 * (xx - x - dx)], order=1,
        )
        # Validate normalization against a complete 13 x 13 model at the same pixel phase.
        ty, tx = np.indices((13, 13), dtype=float)
        phase_norm = np.sum(map_coordinates(
            model, [58 + 9 * (ty - 6 - (y - round(y))), 58 + 9 * (tx - 6 - (x - round(x)))],
            order=1,
        ))
        if not 0.9 < phase_norm < 1.1:
            raise ValueError(f"Unexpected phase normalization: {phase_norm}")
        profiles.append(native / phase_norm)
    return np.asarray(profiles), [float(column), float(row)], receipts


def solve(design, differences, err):
    weights = 1 / err.reshape(len(err), -1) ** 2
    matrix = np.einsum("ni,tn,nj->tij", design, weights, design, optimize=True)
    rhs = (weights * differences.reshape(len(err), -1)) @ design
    covariance = np.linalg.inv(matrix)
    coefficients = np.einsum("tij,tj->ti", covariance, rhs)
    return coefficients, np.sqrt(np.diagonal(covariance, axis1=1, axis2=2))


def fit_frames(times, counts, err, profiles, column_adjustment=False):
    fitted = np.zeros((len(times), len(profiles)))
    errors = np.zeros_like(fitted)
    daily = np.floor(times).astype(int)
    yy, xx = np.indices(counts.shape[1:], dtype=float)
    xx -= np.mean(xx)
    yy -= np.mean(yy)
    if column_adjustment:
        # The constant is already in the span of the column indicators.
        background = [*(xx == col for col in np.unique(xx)), yy, yy**2]
    else:
        background = [np.ones_like(xx), xx, yy, xx**2 + yy**2]
    design_by_day = {}
    for day in np.unique(daily):
        group = daily == day
        reference = np.median(counts[~group], axis=0)
        grad_y, grad_x = np.gradient(reference)
        components = [*profiles, *background, reference, grad_x, grad_y]
        design = np.column_stack([component.ravel() for component in components])
        scale = np.linalg.norm(design, axis=0)
        design /= scale
        differences = counts[group] - reference
        coeff, errors_scaled = solve(design, differences, err[group])
        fitted[group] = (coeff / scale)[:, :len(profiles)]
        errors[group] = (errors_scaled / scale)[:, :len(profiles)]
        design_by_day[day] = (reference, design, scale)
    reference, design, scale = design_by_day[daily[0]]
    injected, _ = solve(design, counts[:1] - reference + 10 * profiles[0], err[:1])
    delta = float(injected[0, 0] / scale[0] - fitted[0, 0])
    if not np.isclose(delta, 10.0, atol=1e-7):
        raise ValueError(f"Matched-profile injection recovered {delta} instead of 10 e-/s")
    return fitted, errors, design_by_day, delta


def bin_frames(times, flux, err, hours, cadence_seconds):
    # UTC bins are aligned to whole clock hours or days. No interpolation through gaps.
    bins = np.floor(times * 24 / hours).astype(int)
    unique, start_indices, counts = np.unique(bins, return_index=True, return_counts=True)
    minimum = int(np.ceil(0.60 * hours * 3600 / cadence_seconds - 1e-7))
    keep = counts >= minimum
    weights = 1 / err**2
    total = np.add.reduceat(weights, start_indices, axis=0)[keep]
    fitted = np.add.reduceat(weights * flux, start_indices, axis=0)[keep] / total
    formal = 1 / np.sqrt(total)
    raw = fitted / formal
    center = np.median(raw, axis=0)
    inflation = np.maximum(1, 1.4826 * np.median(np.abs(raw - center), axis=0))
    empirical = formal * inflation
    return unique[keep] * hours / 24, counts[keep], minimum, fitted, formal, empirical, inflation


def color_limit(amplitude):
    return float(ZEROPOINT - 2.5 * np.log10(amplitude))


def quantify(times, flux, err, hours, cadence_seconds):
    start, counts, minimum, signal, formal, empirical, inflation = bin_frames(
        times, flux, err, hours, cadence_seconds,
    )
    score = signal / empirical
    peak = int(np.argmax(score[:, 0]))
    brightest_control = float(np.max(score[:, 1:]))
    same_window_controls = score[peak, 1:]
    required = np.maximum(0, 5 * empirical[:, 1:] - signal[:, 1:])
    # Transfer control residual scores onto the target's corresponding noise scale.
    target_amplitude = required / empirical[:, 1:] * empirical[:, :1]
    recovery = float(np.quantile(target_amplitude, 0.90, method="higher"))
    quiet = np.max(np.abs(score[:, 1:]), axis=1) < 5
    excursions = np.flatnonzero(quiet & (score[:, 0] >= 5))
    isolated = [
        dict(
            start_utc=mjd_iso(start[index]), score=float(score[index, 0]),
            delta_e_s=float(signal[index, 0]), uncertainty_e_s=float(empirical[index, 0]),
        )
        for index in excursions
    ]
    boundaries = Time([start[peak], start[peak] + hours / 24], format="mjd", scale="utc").isot
    result = dict(
        duration_hours=hours, windows=int(len(start)), minimum_cadences=minimum,
        control_windows=int(required.size), first_bin_utc=mjd_iso(start[0]),
        last_bin_end_utc=mjd_iso(start[-1] + hours / 24),
        target_median_formal_uncertainty_e_s=float(np.median(formal[:, 0])),
        target_median_empirical_uncertainty_e_s=float(np.median(empirical[:, 0])),
        error_inflation_by_site=inflation.tolist(),
        peak_target=dict(
            window_start_utc=boundaries[0], window_end_utc=boundaries[1],
            cadences=int(counts[peak]), delta_e_s=float(signal[peak, 0]),
            formal_uncertainty_e_s=float(formal[peak, 0]),
            empirical_uncertainty_e_s=float(empirical[peak, 0]),
            residual_score=float(score[peak, 0]),
            same_window_control_scores=same_window_controls.tolist(),
            same_window_control_max_abs_score=float(np.max(np.abs(same_window_controls))),
            same_window_control_count_abs_above_five=int(np.count_nonzero(np.abs(same_window_controls) >= 5)),
        ),
        max_control_score=brightest_control,
        control_scores_at_or_above_target_peak=int(np.count_nonzero(score[:, 1:] >= score[peak, 0])),
        target_windows_score_at_least_five=int(np.count_nonzero(score[:, 0] >= 5)),
        control_windows_score_at_least_five=int(np.count_nonzero(score[:, 1:] >= 5)),
        windows_with_all_controls_below_abs_score_five=int(np.count_nonzero(quiet)),
        target_excursions_in_those_windows=isolated,
        control_recovery_90pct_5sigma_e_s=recovery,
        approximate_incremental_tess_magnitude_at_recovery=color_limit(recovery),
    )
    return result, (start, signal, empirical, peak)


def mjd_iso(mjd):
    return Time(mjd, format="mjd", scale="utc").isot


def run_sector(cache, product, coord):
    header, member, count, xy, times, flux, err = read_cutout(product, coord)
    profiles, detector, receipts = prf_design(cache, product, header, xy, flux.shape[1:])
    fits_flux, fits_error, by_day, amplitude = fit_frames(times, flux, err, profiles)
    column_flux, column_error, _, column_amplitude = fit_frames(
        times, flux, err, profiles, column_adjustment=True,
    )
    cadence = header["TIMEDEL"] * 86400
    results, series, column_results = {}, {}, {}
    for hours in DURATIONS:
        results[str(hours)], series[hours] = quantify(times, fits_flux, fits_error, hours, cadence)
        column_results[str(hours)], _ = quantify(times, column_flux, column_error, hours, cadence)
    start, _, _, peak = series[1]
    use = (times >= start[peak]) & (times < start[peak] + 1 / 24)
    reference, design, _ = by_day[int(np.floor(start[peak]))]
    difference = flux[use] - reference
    coefficients, _ = solve(design, difference, err[use])
    nuisance = coefficients[:, len(profiles):] @ design[:, len(profiles):].T
    residual = difference - nuisance.reshape(difference.shape)
    weights = 1 / err[use] ** 2
    image = np.sum(residual * weights, axis=0) / np.sum(weights, axis=0)
    target_template = profiles[0] * results["1"]["peak_target"]["delta_e_s"]
    core = profiles[0] >= 0.02 * np.max(profiles[0])
    compare = results["1"]["peak_target"]
    compare["residual_image_rms_inside_target_profile_core_e_s_per_pixel"] = float(
        np.sqrt(np.mean(image[core] ** 2))
    )
    compare["residual_image_rms_outside_target_profile_core_e_s_per_pixel"] = float(
        np.sqrt(np.mean(image[~core] ** 2))
    )
    compare["target_model_cross_correlation_with_residual_image"] = float(
        np.corrcoef(target_template.ravel(), image.ravel())[0, 1]
    )
    compare["fraction_of_centered_image_variance_in_column_means"] = float(
        np.var(np.mean(image, axis=0)) / np.var(image)
    )
    summary = dict(
        **product, prf_receipts=receipts, file_member=member, cutout_target_xy=list(xy),
        ffi_target_column_row=detector, target_profile_fraction_in_cutout=float(np.sum(profiles[0])),
        cadence_seconds=float(cadence), fits_cadences=count, usable_cadences=int(len(times)),
        fits_start_utc=header["DATE-OBS"], fits_end_utc=header["DATE-END"],
        first_valid_midpoint_utc=mjd_iso(times[0]), last_valid_midpoint_utc=mjd_iso(times[-1]),
        valid_midpoint_dates_utc=sorted(set(Time(times, format="mjd", scale="utc").iso.astype("U10"))),
        matched_profile_injection_recovered_e_s=amplitude, windows=results,
        column_adjusted_matched_profile_injection_recovered_e_s=column_amplitude,
        column_adjusted_windows=column_results,
    )
    # The peak-window image retains the target and fitted control sources.
    return summary, series, (reference, image, xy)


def plot_result(rows):
    fig = Figure(figsize=(10.8, 9.4), layout="constrained")
    grid = fig.add_gridspec(3, 4, height_ratios=(1.10, 1, 1))
    _, _, (reference, difference, xy) = rows[-1]
    for col, (data, title) in enumerate((
        (reference, "S108 leave-day-out reference"),
        (difference, "S108 strongest positive hour"),
    )):
        ax = fig.add_subplot(grid[0, 2 * col:2 * col + 2])
        if col == 0:
            center = float(np.median(data))
            high = float(np.percentile(data, 97))
            rendered = ax.imshow(data, origin="lower", cmap="viridis", vmin=center - (high - center) / 4,
                                 vmax=high)
        else:
            bound = float(np.percentile(np.abs(data), 98))
            rendered = ax.imshow(data, origin="lower", cmap="RdBu_r", vmin=-bound, vmax=bound)
        ax.scatter(*xy, s=132, facecolors="none", edgecolors="white" if col == 0 else "#3d4324", lw=1.25)
        ax.set(title=title, xlabel="Cutout column (pixels)", ylabel="Cutout row (pixels)")
        fig.colorbar(rendered, ax=ax, label=r"e$^-$/s per pixel", shrink=0.83)
    # Three short panels per grid row, with an independent y scale in each sector.
    light_curves = [grid[index, :].subgridspec(1, 3) for index in (1, 2)]
    for index, (row, binned, _) in enumerate(rows):
        ax = fig.add_subplot(light_curves[index // 3][0, index % 3])
        utc, values, uncertainties, _ = binned[24]
        dates = Time(utc + 0.5, format="mjd", scale="utc").datetime
        ax.errorbar(dates, values[:, 0], uncertainties[:, 0], fmt="o", color="#226e7f",
                    ms=3, lw=0.8, capsize=1.6)
        ax.axhline(0, color="#666666", lw=0.8)
        ax.xaxis.set_major_formatter(DateFormatter("%d %b"))
        ax.xaxis.set_major_locator(AutoDateLocator(minticks=3, maxticks=5))
        ax.tick_params(axis="x", labelsize=8, rotation=26)
        ax.grid(axis="y", alpha=0.18)
        year = row["fits_start_utc"][:4]
        ax.set(title=f"Sector {row['sector']} · {year}", ylabel=r"Daily excess (e$^-$/s)")
    fig.suptitle("TESS observations at the JADES 103607 transient position", fontsize=16)
    save_figure(fig, FIGURES / TAG)


def plot_peak_images(rows):
    fig = Figure(figsize=(10.5, 7.2), layout="constrained")
    for ax, (row, _, (_, image, xy)) in zip(fig.subplots(2, 3).flat, rows, strict=True):
        peak = row["windows"]["1"]["peak_target"]
        bound = float(np.percentile(np.abs(image), 99))
        shown = ax.imshow(image, origin="lower", cmap="RdBu_r", vmin=-bound, vmax=bound)
        ax.scatter(*np.asarray(OFFSETS).T + np.asarray(xy)[:, None], marker="+", s=37,
                   color="#666666", alpha=0.66)
        ax.scatter(*xy, s=130, facecolors="none", edgecolors="#202020", lw=1.15)
        ax.set(
            title=f"S{row['sector']} · {peak['window_start_utc'][:16].replace('T', ' ')}",
            xlabel="Cutout column (pixels)", ylabel="Cutout row (pixels)",
        )
        fig.colorbar(shown, ax=ax, label=r"e$^-$/s per pixel", shrink=0.80)
    fig.suptitle("Observed TESS residuals at each sector’s highest-scoring hour (UTC)", fontsize=15)
    save_figure(fig, FIGURES / (TAG + "_peaks"))


def main():
    path = OUT / f"{INVENTORY_TAG}.json"
    inventory = json.loads(path.read_text(encoding="utf-8"))
    coord = SkyCoord(*inventory["position_icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    rows = []
    for product in inventory["tesscut_products"]:
        result, series, images = run_sector(cache, product, coord)
        rows.append((result, series, images))
        print(f"TESS {result['sector']}: {result['usable_cadences']} usable cadences", flush=True)
    newest = rows[-1][0]
    native = json.loads(NATIVE.read_text(encoding="utf-8"))
    jwst_mjd = min(
        exposure["expstart_mjd"]
        for band in native["bands"].values()
        for visit in band["inventory"]
        for exposure in visit["products"]
        if "2026-10" in mjd_iso(exposure["expstart_mjd"])
    )
    plot_result(rows)
    plot_peak_images(rows)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        tess_inventory_sha256=text_digest(path), jwst_native_fade_sha256=text_digest(NATIVE),
        position_icrs_deg=inventory["position_icrs_deg"],
        site_offsets_xy=OFFSETS, tess_zeropoint=ZEROPOINT,
        first_jwst_october_exposure_start_utc=mjd_iso(jwst_mjd),
        last_tess_fits_end_utc=newest["fits_end_utc"],
        hours_from_last_tess_fits_end_to_first_october_jwst=(
            (jwst_mjd - Time(newest["fits_end_utc"]).mjd) * 24
        ),
        sectors=[row for row, _, _ in rows],
        nearest_sectors_planned_coverage=inventory["sectors_during_2026_window"],
        method="Fit the MAST calibrated pixel response of a source at the fixed JWST position and "
        "eight spatial controls jointly, using QUALITY=0 and finite flux/errors at all pixels. "
        "Subtract a median reference excluding the UTC day under measurement. Use a local quadratic "
        "radial background, reference scene, and its two gradients as nuisance terms. "
        "Bin in fixed UTC hours and days, retaining windows with at least 60% cadence completeness. "
        "Subtract TIMECORR from BTJD to reconstruct cadence midpoint UTC. "
        "Use the larger of formal errors and formal errors scaled by the MAD of flux/formal error. "
        "A second fit uses one nuisance offset per detector column and quadratic row terms to "
        "test the effect of stripes in the residual images. Both fits list excursions that "
        "exceed score five while all eight same-window control scores have magnitude below five.",
        recovery_scope="Analytic matched-profile recovery on measured control residuals. Controls "
        "are rescaled to target uncertainty within corresponding windows. We quote the increment "
        "needed to exceed residual score 5 in at least 90% of such trials, conditional on observed "
        "windows and a fixed difference-image reference. Controls share images and are not "
        "independent global false-alarm trials. Magnitudes approximate the incremental TESS-band "
        "source flux using the instrument zeropoint, on the TESS Vega-based scale. Bright sources "
        "can saturate. Brief flares and signals in gaps or spanning much of a sector may be missed.",
    ))


if __name__ == "__main__":
    main()

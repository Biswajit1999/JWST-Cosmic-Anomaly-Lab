"""Measure calibrated historical DES z-band coadds at the JADES 103607 position."""

import argparse
import json
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from inventory_jades_103607_des_z import TAG as INVENTORY_TAG
from matplotlib.dates import DateFormatter, YearLocator
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from requests import RequestException
from scipy.ndimage import gaussian_filter, map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_des_z_history"
SCALE = 0.26
GRID = np.arange(-43, 44)
CENTER = len(GRID) // 2
FWHM = 1.5
RADIUS = 1.0
ANNULUS = (1.6, 2.5)
FLUX_TRIALS = (36, 100, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000, 9000, 14000)
MIN_NIGHT_EXPOSURES = 2
FGCM_FLAGS = {1, 2, 4}


def plan_nights(inventory, one_night_per_month):
    by_exposure = defaultdict(list)
    for row in inventory["footprints"]:
        by_exposure[row["expnum"], row["ccdnum"]].append(row)
    nights, rejected = defaultdict(list), Counter()
    exposures = set()
    for row in inventory["products"]:
        zp = row["fgcm"]
        if not zp or zp["calibration_flag"] not in FGCM_FLAGS:
            rejected["missing or rejected FGCM calibration"] += 1
            continue
        if zp["calibration_source"] != "FGCM" or zp["calibration_version"] != "y6a1_v2.1":
            raise ValueError("Unexpected DES zeropoint source/version")
        if set(row["cutouts"]) != {"science", "mask"}:
            rejected["missing science or mask cutout URL"] += 1
            continue
        footprints = by_exposure[row["expnum"], row["ccdnum"]]
        best = max(footprints, key=lambda x: x["approximate_signed_edge_arcsec"])
        if best["approximate_signed_edge_arcsec"] < 20:
            rejected["within 20 arcsec of archived CCD edge"] += 1
            continue
        if row["expnum"] in exposures:
            raise ValueError(f"Multiple eligible reductions or covering CCDs for exposure {row['expnum']}")
        exposures.add(row["expnum"])
        timestamp = datetime.fromisoformat(best["dateobs_center"])
        night = (timestamp - timedelta(hours=12)).date().isoformat()
        nights[night].append(dict(
            **row, archive_dateobs_center=best["dateobs_center"], field=best["object"],
            archive_exposure_s=best["exposure"], archive_edge_arcsec=best["approximate_signed_edge_arcsec"],
            archive_chip_fwhm_pixels=best["archive_chip_fwhm_pixels"],
        ))

    eligible = {night: exposures for night, exposures in nights.items() if len(exposures) >= MIN_NIGHT_EXPOSURES}
    by_month = defaultdict(list)
    for night in eligible:
        by_month[night[:7]].append(night)
    selected = set(eligible)
    if one_night_per_month:
        selected = {
            min(dates, key=lambda n: (-sum(float(row["archive_exposure_s"]) for row in eligible[n]), n))
            for dates in by_month.values()
        }
    return nights, eligible, selected, rejected


def angular_aperture():
    span = int(np.ceil(ANNULUS[1] / SCALE))
    dx, dy = np.meshgrid(np.arange(-span, span+1), np.arange(-span, span+1))
    distance = np.hypot(dx, dy) * SCALE
    core = distance <= RADIUS
    ring = (distance >= ANNULUS[0]) & (distance <= ANNULUS[1])
    use = core | ring
    return dx[use], dy[use], core[use]


DX, DY, CORE = angular_aperture()
CONTROL = [
    (x, y)
    for x in range(12, len(GRID)-11, 10)
    for y in range(12, len(GRID)-11, 10)
    if np.hypot(x-CENTER, y-CENTER) * SCALE >= 3.9
]


def aperture_flux(pixels):
    core = pixels[..., CORE]
    ring = pixels[..., ~CORE]
    good = np.isfinite(core).all(axis=-1) & (np.isfinite(ring).sum(axis=-1) >= 0.8 * ring.shape[-1])
    values = np.full(good.shape, np.nan)
    values[good] = core[good].sum(axis=-1) - CORE.sum() * np.nanmedian(ring[good], axis=-1)
    return values


def read_exposure(cache, row, grid):
    images, headers, receipts = {}, {}, {}
    for plane, hdu in (("science", "SCI"), ("mask", "MSK")):
        url = row["cutouts"][plane]
        rel = f"cross_detector/des_z_stamps/{row['filename']}.{plane}.fits"
        path = cache.fetch(
            url, rel, survey="NSF NOIRLab Astro Data Lab",
            version="DES DR2 single-epoch Final Cut", kind=f"observed_{plane}_cutout",
        )
        receipts[plane] = receipt(path, url)
        with fits.open(path, memmap=False) as fits_file:
            images[plane] = fits_file[hdu].data
            headers[plane] = fits_file[hdu].header
    header = headers["science"]
    if (header["EXPNUM"], header["CCDNUM"], header["BUNIT"], header["FILTER"].split()[0]) != (
        row["expnum"], row["ccdnum"], "electrons", "z"
    ):
        raise ValueError(f"Unexpected DES FITS exposure identity for {row['filename']}")
    row.update(files=receipts, date_obs_utc=header["DATE-OBS"], exposure_s=header["EXPTIME"])
    clock = datetime.fromisoformat(row["date_obs_utc"]).replace(tzinfo=timezone.utc)
    metadata_clock = datetime.fromisoformat(row["archive_dateobs_center"])
    if abs((clock - metadata_clock).total_seconds()) > 2 * header["EXPTIME"]:
        raise ValueError(f"DES FITS exposure time does not match the footprint {row['filename']}")
    wcs = WCS(header).celestial
    band = images["science"]
    mask = images["mask"]
    x, y = wcs.world_to_pixel(grid)
    mx, my = WCS(headers["mask"]).celestial.world_to_pixel(grid)
    if band.shape != mask.shape or not (np.allclose(mx, x, rtol=0, atol=1e-4)
                                       and np.allclose(my, y, rtol=0, atol=1e-4)):
        raise ValueError(f"DES science and mask stamps are misregistered {row['filename']}")
    if not (0 <= x[CENTER, CENTER] < band.shape[1] and 0 <= y[CENTER, CENTER] < band.shape[0]):
        row["status"] = "position outside returned science cutout"
        return None
    native_area = abs(np.linalg.det(wcs.pixel_scale_matrix)) * 3600**2
    native_scale = float(np.sqrt(native_area))
    seeing = float(header.get("FWHM", np.nan)) * native_scale
    zp = row["fgcm"]["ab_zp"]
    factor = float(3.631 * 10**((30-zp)/2.5) * (SCALE/native_scale)**2)
    row.update(fwhm_arcsec=seeing if np.isfinite(seeing) else None, resampled_njy_per_electron=factor)
    if not 0.5 <= seeing <= FWHM:
        row["status"] = "missing FWHM or outside Gaussian matching bounds"
        return None
    valid = np.isfinite(band) & (mask == 0)
    sigma = np.sqrt(FWHM**2 - seeing**2) / (2.354820045 * native_scale)
    kernel_support = gaussian_filter(valid.astype(float), sigma, mode="constant")
    smooth = gaussian_filter(np.where(valid, band, 0.0).astype(float), sigma, mode="constant")
    smooth = np.divide(smooth, kernel_support, out=np.full_like(smooth, np.nan), where=kernel_support >= 0.98)
    regridded = map_coordinates(smooth, (y, x), order=1, cval=np.nan, prefilter=False)
    if not np.isfinite(regridded).any():
        row["status"] = "no valid resampled sky pixels"
        return None
    image = (regridded - np.nanmedian(regridded)) * factor
    row["status"] = "usable"
    row["cutout_valid_fraction"] = float(np.isfinite(image).mean())
    return image


def median_image(images, minimum):
    data = np.asarray(images)
    stack = np.ma.median(np.ma.masked_invalid(data), axis=0).filled(np.nan)
    return np.where(np.isfinite(data).sum(axis=0) >= minimum, stack, np.nan)


def measure(residual):
    target = float(aperture_flux(residual[CENTER+DY, CENTER+DX]))
    x, y = np.asarray(CONTROL).T
    pixels = residual[y[:, None]+DY, x[:, None]+DX]
    raw = aperture_flux(pixels)
    valid = np.isfinite(raw)
    if valid.sum() < 20 or not np.isfinite(target):
        return None
    controls = raw[valid]
    bias = float(np.median(controls))
    scatter = float(1.4826 * np.median(np.abs(controls - bias)))
    if scatter <= 0:
        return None
    fractions = {}
    for profile in ("Gaussian", "Moffat beta=3"):
        r2 = (DX**2 + DY**2) * SCALE**2
        extent = GRID[:, None]**2 * SCALE**2 + GRID[None, :]**2 * SCALE**2
        if profile == "Gaussian":
            denominator = 2 * (FWHM / 2.354820045)**2
            template = np.exp(-r2 / denominator) / np.exp(-extent / denominator).sum()
        else:
            alpha2 = (FWHM / 2)**2 / (2**(1/3)-1)
            template = (1+r2/alpha2)**-3 / ((1+extent/alpha2)**-3).sum()
        model = pixels[None, valid] + np.asarray(FLUX_TRIALS)[:, None, None] * template
        recovered = (aperture_flux(model) - bias >= 5*scatter).mean(axis=1)
        fractions[profile] = dict(
            at_tested_fluxes=recovered.tolist(),
            first_tested_flux_with_90_percent_recovery_njy=next(
                (flux for flux, f in zip(FLUX_TRIALS, recovered, strict=True) if f >= 0.9), None
            ),
        )
    return dict(
        signed_aperture_excess_njy=target-bias, empirical_scatter_njy=scatter,
        residual_over_control_scatter=(target-bias)/scatter, control_median_njy=bias,
        control_indices=np.flatnonzero(valid).tolist(), control_residuals_njy=(controls-bias).tolist(),
        five_scatter_aperture_njy=5*scatter, injection_recovery=fractions,
    )


def season(night):
    year, month = map(int, night.split("-")[:2])
    return year if month >= 7 else year-1


def control_peaks(measured):
    peak = max(measured, key=lambda row: row["difference"]["residual_over_control_scatter"])
    site_scores = defaultdict(list)
    for row in measured:
        result = row["difference"]
        for index, flux in zip(result["control_indices"], result["control_residuals_njy"], strict=True):
            site_scores[index].append(flux / result["empirical_scatter_njy"])
    matched = [
        max(site_scores[index]) for index in peak["difference"]["control_indices"]
        if len(site_scores[index]) >= np.ceil(0.8 * len(measured))
    ]
    score = peak["difference"]["residual_over_control_scatter"]
    return dict(
        peak_night=peak["night"], residual_over_control_scatter=score,
        controls_covering_peak_and_80_percent_of_nights=len(matched),
        controls_reaching_peak_score=sum(value >= score for value in matched),
    )


def plot(measured, nightly, references, residuals):
    fig = Figure(figsize=(11.9, 7.5), layout="constrained")
    axes = fig.subplots(2, 3, gridspec_kw={"height_ratios": (1.32, 1)})
    axes[0, 1].remove()
    axes[0, 2].remove()
    time_ax = axes[0, 0]
    gs = time_ax.get_subplotspec().get_gridspec()
    time_ax.set_subplotspec(gs[0, :])
    dates = [datetime.fromisoformat(row["night"]) for row in measured]
    flux = [row["difference"]["signed_aperture_excess_njy"] for row in measured]
    scatters = [row["difference"]["empirical_scatter_njy"] for row in measured]
    time_ax.errorbar(dates, flux, yerr=scatters, fmt="o", ms=4, lw=1, color="#17687c", capsize=2)
    time_ax.axhline(0, linewidth=0.9, color="#555555")
    time_ax.set(ylabel="Signed z aperture excess (nJy)", title="DES nights at the JWST transient position")
    time_ax.xaxis.set_major_locator(YearLocator())
    time_ax.xaxis.set_major_formatter(DateFormatter("%Y"))
    picked = max(measured, key=lambda row: row["difference"]["residual_over_control_scatter"])["night"]
    half = (GRID[-1]+0.5) * SCALE
    for ax, label, img in zip(
        axes[1], (f"{picked} observing night", "Reference from other seasons", "Night minus reference"),
        (nightly[picked], references[picked], residuals[picked]), strict=True,
    ):
        centered = img - np.nanmedian(img)
        stretch = float(np.nanpercentile(np.abs(centered), 96))
        im = ax.imshow(centered, cmap="RdBu_r", extent=(-half, half, -half, half), origin="lower",
                       interpolation="nearest", vmin=-stretch, vmax=stretch)
        ax.add_patch(Circle((0, 0), RADIUS, edgecolor="black", facecolor="none"))
        ax.set(title=label, xlabel="East (arcsec)")
        fig.colorbar(im, ax=ax, shrink=0.76, label="nJy per grid pixel")
    axes[1, 0].set_ylabel("North (arcsec)")
    FIGURES.mkdir(parents=True, exist_ok=True)
    save_figure(fig, FIGURES / TAG)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-night-per-month", action="store_true",
                        help="Measure the largest eligible z exposure sum in each calendar month.")
    parser.add_argument("--max-exposures-per-night", type=int,
                        help="Cap each night using exposure time, then CCD seeing, as the ranking.")
    args = parser.parse_args()
    if args.max_exposures_per_night is not None and args.max_exposures_per_night < MIN_NIGHT_EXPOSURES:
        parser.error(f"At least {MIN_NIGHT_EXPOSURES} exposures per night are required.")
    path = OUT / f"{INVENTORY_TAG}.json"
    inventory = json.loads(path.read_text(encoding="utf-8"))
    sky = SkyCoord(*inventory["position_icrs_deg"], unit="deg")
    earth_offset = GRID * SCALE / 3600 * sky.ra.unit
    east, north = np.meshgrid(earth_offset, earth_offset)
    grid = sky.spherical_offsets_by(east, north)
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    nights, eligible, selected, rejected = plan_nights(inventory, args.one_night_per_month)

    nightly, rows = {}, []
    for night in sorted(selected):
        group = sorted(eligible[night], key=lambda row: (
            -float(row["archive_exposure_s"]),
            float(row["archive_chip_fwhm_pixels"] or np.inf),
            row["expnum"],
        ))[:args.max_exposures_per_night]
        images = []
        for row in group:
            try:
                image = read_exposure(cache, row, grid)
            except RequestException as exc:
                print(f"{row['filename']}: {exc}", flush=True)
                row["status"] = "failed HTTP retrieval"
                continue
            if image is not None:
                images.append(image)
        print(f"DES z night {night}: {len(images)}/{len(group)} usable images", flush=True)
        rows.append(dict(night=night, season=season(night), products=group, usable_images=len(images)))
        if len(images) >= MIN_NIGHT_EXPOSURES:
            nightly[night] = median_image(images, max(MIN_NIGHT_EXPOSURES, len(images)//2+1))

    references, residuals = {}, {}
    for row in rows:
        night = row["night"]
        reference_dates = [date for date in nightly if season(date) != season(night)]
        if night not in nightly or len(reference_dates) < 4:
            row["status"] = "too few good exposures or reference nights"
            continue
        refs = [nightly[date] for date in reference_dates]
        references[night] = median_image(refs, len(refs)//2+1)
        residuals[night] = nightly[night] - references[night]
        row["reference_nights"] = reference_dates
        row["difference"] = measure(residuals[night])
        row["status"] = "measured" if row["difference"] is not None else "insufficient valid aperture or controls"

    measured = [row for row in rows if row.get("difference") is not None]
    if not measured:
        raise ValueError("No measured DES z-band night stacks")
    plot(measured, nightly, references, residuals)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__), inventory_result_sha256=text_digest(path),
        position_icrs_deg=inventory["position_icrs_deg"], aperture_radius_arcsec=RADIUS,
        annulus_arcsec=ANNULUS, common_gaussian_fwhm_arcsec=FWHM, fgcm_flag_values=sorted(FGCM_FLAGS),
        one_night_per_month=args.one_night_per_month, excluded_products=dict(rejected),
        max_exposures_per_night=args.max_exposures_per_night,
        archived_nights_with_eligible_products={date: len(group) for date, group in sorted(nights.items())},
        eligible_nights=len(eligible), selected_nights=len(selected), measured_nights=len(measured),
        injection_flux_grid_total_njy=FLUX_TRIALS, injection_trigger_control_scatter=5,
        nightly_results=rows, peak_control_comparison=control_peaks(measured),
        scope="DES z-band SCI and MSK cutouts are matched to the FGCM calibrated image basename "
        "and reprojected after approximate Gaussian FWHM matching. Nightly medians have no WGT "
        "weighting and subtract a reference made from other observing seasons. UTC night boundaries "
        "are at noon. Reported scatter describes local aperture residuals including blended light. "
        "Gaussian and Moffat injections are performed on residual coadds to measure recovery of "
        "persistent light during the sampled exposures without modelling individual exposure durations.",
    ))


if __name__ == "__main__":
    main()

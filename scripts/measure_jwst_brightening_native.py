"""Measure signed aperture changes directly on the individual NIRCam exposures."""

import json
from types import SimpleNamespace
from warnings import catch_warnings, filterwarnings

import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning
from download_jwst_cross_detector import OUT, ROOT
from matplotlib import colormaps
from matplotlib.figure import Figure
from measure_jwst_brightening_offset import BANDS, patch
from measure_jwst_epoch_candidate import RADII, pixel_scale, show_stamp
from measure_jwst_module_b_brightening import ANNULUS, FILTERS, NEW, OLD
from search_jwst_cross_detector import window
from vet_jwst_field_b_candidate import FIGURES
from vet_jwst_module_b import native_alignment, native_products

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

STORAGE_MIB = 22_000
SUBPIXELS = 8
SUBPIXEL_CENTRES = (np.arange(SUBPIXELS)+.5)/SUBPIXELS-.5
VALID_FRACTION = .95
JUMP_DET = 4
VIEW_BANDS = ("F150W", "F200W", "F277W")


def fractional_apertures(image, x, y, scale):
    """Use only observed pixel values, weighting pixels intersected by the aperture."""
    half = int(np.ceil(ANNULUS[1]/scale))+1
    ix, iy = round(x), round(y)
    if min(ix, iy, image.shape[1]-ix-1, image.shape[0]-iy-1) < half:
        return np.full(len(RADII), np.nan), np.zeros(len(RADII))
    yy, xx = np.mgrid[iy-half:iy+half+1, ix-half:ix+half+1]
    dx, dy = xx-x, yy-y
    distance = np.hypot(dx, dy)
    patch = image[yy, xx]
    valid = np.isfinite(patch)
    sky = (distance >= ANNULUS[0]/scale) & (distance <= ANNULUS[1]/scale)
    if valid[sky].sum() < .8*sky.sum():
        return np.full(len(RADII), np.nan), np.zeros(len(RADII))
    background = np.median(patch[sky & valid])
    fluxes, fractions = [], []
    reach = np.sqrt(.5)
    for radius in RADII:
        r = radius/scale
        area = (distance <= r-reach).astype(float)
        edge = abs(distance-r) <= reach
        xx_sub = dx[edge, None, None]+SUBPIXEL_CENTRES[None, :, None]
        yy_sub = dy[edge, None, None]+SUBPIXEL_CENTRES[None, None, :]
        area[edge] = ((xx_sub**2+yy_sub**2) <= r**2).mean(axis=(1, 2))
        fraction = area[valid].sum()/area.sum()
        central_ok = valid[distance <= .07/scale].all()
        flux = np.sum((patch[valid]-background)*area[valid])
        fluxes.append(flux if fraction >= VALID_FRACTION and central_ok else np.nan)
        fractions.append(fraction)
    return np.asarray(fluxes), np.asarray(fractions)


def frame_photometry(frame, image, sky):
    x, y = frame["wcs"].world_to_pixel_values(*sky)
    xy = np.column_stack((np.atleast_1d(x), np.atleast_1d(y)))+frame["translation"]
    return [fractional_apertures(image, *point, frame["scale"]) for point in xy]


def positions_from_mosaic(record, target, references):
    _, wcs, _, scale = window(record, "sci", target, 16.5)
    x, y = wcs.world_to_pixel_values(*target)
    displacement = np.array([row["offset_arcsec"] for row in references if .85 <= row["new_to_old"] <= 1.15])
    return wcs.pixel_to_world_values(x+displacement[:, 0]/scale, y+displacement[:, 1]/scale)


def local_controls(wcs, target):
    scale = pixel_scale(wcs)
    x, y = wcs.world_to_pixel_values(*target)
    axis = np.arange(-6., 6.1, .60)
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    good = (np.hypot(xx, yy) >= 1.5) & (np.hypot(xx, yy) <= 6)
    return wcs.pixel_to_world_values(x+xx[good]/scale, y+yy[good]/scale)


def robust_sigma(values, minimum):
    values = values[np.isfinite(values)]
    if len(values) < minimum:
        return np.nan
    return float(1.4826*np.median(abs(values-np.median(values))))


def relative_gains(references):
    """Choose the largest reference aperture with enough coverage in every exposure."""
    trial_counts = {}
    for i in reversed(range(len(RADII))):
        old = references[:3, :, i]
        keep = (np.isfinite(old).sum(axis=0) >= 2) & (np.nansum(old, axis=0) > 0)
        anchors = np.nanmedian(old[:, keep], axis=0)
        ratios = references[:, keep, i]/anchors
        clean = [r[np.isfinite(r) & (r > .80) & (r < 1.20)] for r in ratios]
        trial_counts[f"{RADII[i]:.2f}"] = [len(r) for r in clean]
        if min(trial_counts[f"{RADII[i]:.2f}"]) >= 5:
            gains = [dict(aperture_gain=float(np.median(r)), reference_count=len(r),
                          reference_gain_mad=robust_sigma(r, 5)) for r in clean]
            return gains, dict(method="local_reference_sources", radius_arcsec=RADII[i], trial_counts=trial_counts)
    gains = [dict(aperture_gain=1., reference_count=None, reference_gain_mad=None) for _ in references]
    return gains, dict(method="pipeline_calibration", radius_arcsec=None, trial_counts=trial_counts)


def apertures(frames, refs, controls, sky):
    calibrated_images = [np.where((f["dq"] == 0) | (f["dq"] == JUMP_DET), f["image"], np.nan) for f in frames]
    paired = list(zip(frames, calibrated_images, strict=True))
    ref = np.asarray([[v[0] for v in frame_photometry(frame, image, refs)] for frame, image in paired])
    gains, audit = relative_gains(ref)
    target = [frame_photometry(frame, image, sky)[0] for frame, image in paired]
    zero_dq = [frame_photometry(frame, frame["clean"], sky)[0] for frame in frames]
    blank = np.asarray([[v[0] for v in frame_photometry(frame, image, controls)] for frame, image in paired])
    calibrated, background, records = [], [], []
    for i, (flux, fractions) in enumerate(target):
        gain = gains[i]["aperture_gain"]
        calibrated.append(flux/gain)
        background.append(blank[i]/gain)
        records.append(dict(**gains[i], usable_area_fraction=fractions.tolist(),
                            zero_dq_area_fraction=zero_dq[i][1].tolist(),
                            zero_dq_aperture_flux_njy=[float(v/gain) if np.isfinite(v) else None for v in zero_dq[i][0]],
                            raw_aperture_flux_njy=[float(v) if np.isfinite(v) else None for v in flux],
                            calibrated_aperture_flux_njy=[float(v/gain) if np.isfinite(v) else None for v in flux]))
    return np.asarray(calibrated), np.asarray(background), records, audit


def median_with_minimum(values, minimum, axis):
    count = np.isfinite(values).sum(axis=axis)
    valid = np.ma.masked_invalid(values)
    measured = np.ma.filled(np.ma.median(valid, axis=axis), np.nan)
    return np.where(count >= minimum, measured, np.nan)


def compare(target, background, records):
    report = []
    for i, radius in enumerate(RADII):
        epoch_fluxes = [median_with_minimum(e[:, i], 2, axis=0) for e in (target[:3], target[3:])]
        blank_medians = [median_with_minimum(e[:, :, i], 2, axis=0) for e in (background[:3], background[3:])]
        blank_difference = blank_medians[1]-blank_medians[0]
        uncertainty = robust_sigma(blank_difference, 40)
        measured = epoch_fluxes[1]-epoch_fluxes[0]
        all_pairs = []
        for older in range(3):
            for newer in range(3, 6):
                a, b = target[[older, newer], i]
                if np.isfinite(a+b):
                    all_pairs.append(dict(older_dither=older+1, newer_dither=newer-2, excess_njy=float(b-a)))
        audited = np.asarray([rec["zero_dq_aperture_flux_njy"][i] for rec in records], dtype=float)
        clean_old, clean_new = (median_with_minimum(e, 2, axis=0) for e in (audited[:3], audited[3:]))
        raw = np.asarray([rec["raw_aperture_flux_njy"][i] for rec in records], dtype=float)
        raw_old, raw_new = (median_with_minimum(e, 2, axis=0) for e in (raw[:3], raw[3:]))
        report.append(dict(
            radius_arcsec=radius, median_older_njy=float(epoch_fluxes[0]) if np.isfinite(epoch_fluxes[0]) else None,
            median_newer_njy=float(epoch_fluxes[1]) if np.isfinite(epoch_fluxes[1]) else None,
            difference_njy=float(measured) if np.isfinite(measured) else None,
            local_difference_mad_njy=uncertainty if np.isfinite(uncertainty) else None,
            matched_background_count=int(np.isfinite(blank_difference).sum()),
            paired_dither_differences=all_pairs,
            zero_dq_difference_njy=float(clean_new-clean_old) if np.isfinite(clean_old+clean_new) else None,
            unscaled_difference_njy=float(raw_new-raw_old) if np.isfinite(raw_old+raw_new) else None,
        ))
    return report


def stamps(frames, diagnostic):
    result = []
    for frame, row in zip(frames, diagnostic, strict=True):
        x, y = row["forced_center_pixels_xy"]
        scale = frame["scale"]
        half = int(np.ceil(.72/scale))+2
        x0, y0 = round(x)-half, round(y)-half
        image = frame["clean"][y0:y0+2*half+1, x0:x0+2*half+1].copy()
        result.append((image-np.nanmedian(image), x-x0, y-y0, scale))
    return result


def zero_dq_core(frames, diagnostics, sigma, radius, fluxes):
    """Fit a fixed Gaussian core and background only to DQ-zero pixels."""
    results = []
    for frame, row, calibrated in zip(frames, diagnostics, fluxes, strict=True):
        x, y = row["forced_center_pixels_xy"]
        dx, dy, image = patch(frame["clean"], x, y, frame["scale"], radius)
        radius_sq = dx**2+dy**2
        good = (radius_sq <= radius**2) & np.isfinite(image)
        if good.sum() < .65*(radius_sq <= radius**2).sum():
            raise ValueError(f"Core fit has too few unflagged pixels in {row['file']}")
        t = np.exp(-radius_sq[good]/(2*sigma**2))
        design = np.column_stack((t, np.ones_like(t)))
        pixels = image[good]
        amplitude, sky = np.linalg.lstsq(design, pixels, rcond=None)[0]
        brightest = int(np.argmax(pixels))
        others = np.arange(len(pixels)) != brightest
        without_peak = np.linalg.lstsq(design[others], pixels[others], rcond=None)[0][0]
        results.append(dict(fitted_core_amplitude_njy_per_pixel=float(amplitude/calibrated["aperture_gain"]),
                            without_brightest_pixel_njy_per_pixel=float(without_peak/calibrated["aperture_gain"]),
                            background_njy_per_pixel=float(sky), usable_pixels=int(good.sum())))
    amps = np.asarray([r["fitted_core_amplitude_njy_per_pixel"] for r in results])
    peak_removed = np.asarray([r["without_brightest_pixel_njy_per_pixel"] for r in results])
    return dict(fixed_sigma_arcsec=sigma, fit_radius_arcsec=radius, exposures=results,
                older_median_njy_per_pixel=float(np.median(amps[:3])),
                newer_median_njy_per_pixel=float(np.median(amps[3:])),
                all_newer_amplitudes_above_all_older=bool(amps[3:].min() > amps[:3].max()),
                all_newer_amplitudes_above_all_older_without_brightest_pixel=bool(
                    peak_removed[3:].min() > peak_removed[:3].max()))


def plot(records, native_stamps):
    figure = Figure(figsize=(12.8, 7.5), layout="constrained")
    grid = figure.add_gridspec(4, 6, height_ratios=(1.7, 1, 1, 1))
    ax = figure.add_subplot(grid[0, :])
    for i, band in enumerate(FILTERS):
        result = records[band]["epoch_changes"][-1]
        diffs = [v["excess_njy"] for v in result["paired_dither_differences"]]
        if diffs:
            ax.scatter(i+np.linspace(-.18, .18, len(diffs)), diffs, color="#809da6", s=14,
                       label="Individual dither pairs" if i == 0 else None, alpha=.85)
        value, error = result["difference_njy"], result["local_difference_mad_njy"]
        if value is not None:
            ax.errorbar(i, value, yerr=error, color="#993e24", marker="D", capsize=3, linestyle="none",
                        label="Difference of visit medians" if i == 0 else None)
    ax.set(ylabel="Signed 0.24″ flux change (nJy)", xticks=range(len(FILTERS)), xticklabels=FILTERS)
    ax.axhline(0, color=".55", lw=.7)
    ax.grid(axis="y", alpha=.2)
    ax.legend(frameon=False, loc="upper left")
    masked_gray = colormaps["gray"].with_extremes(bad="#8e64a6")
    for j, band in enumerate(VIEW_BANDS, 1):
        values = native_stamps[band]
        floor, top = np.nanpercentile(np.concatenate([d[0].ravel() for d in values]), [8, 98.5])
        for i, (image, x, y, scale) in enumerate(values):
            ax = figure.add_subplot(grid[j, i])
            show_stamp(ax, image, x, y, scale, (floor, top), masked_gray, .24, half_arcsec=.61)
            ax.set_title(f"{band} · {2025 if i<3 else 2026} · {i%3+1}", fontsize=8)
    figure.suptitle("JADES 103607 · JWST/NIRCam brightness change before image resampling\n"
                    "Native calibrated pixels · purple marks masked pixels · common display scale within each filter",
                    fontsize=12)
    save_figure(figure, FIGURES/"module_b_001_native_epoch_flux")


def native_metadata(receipt, visit, band, detector):
    metadata = []
    for product in receipt["products"]:
        h = fits.getheader(ROOT/product["path"])
        if (h["VISIT_ID"], h["FILTER"], h["DETECTOR"], h["PUPIL"]) != (visit[2:], band, detector, "CLEAR"):
            raise ValueError(f"Exposure metadata mismatch for {product['file']}")
        if h["S_PHOTOM"] != "COMPLETE" or h.get("S_RESAMP") == "COMPLETE":
            raise ValueError(f"Expected unresampled calibrated exposure: {product['file']}")
        metadata.append(dict(filename=product["file"], calibration_pipeline=h["CAL_VER"],
                             crds_context=h["CRDS_CTX"], detector=h["DETECTOR"],
                             exposure_start_mjd=h["EXPSTART"], start_utc=h["DATE-BEG"]))
    return metadata


def measure_band(cache, band, selections, comparison, fit_data):
    sky = comparison["icrs_deg"]
    base = selections[band, OLD]
    prior = comparison["difference_photometry"][band]["local_flux_scale"]["references"]
    refs = positions_from_mosaic(base, sky, prior)
    controls = local_controls(window(base, "sci", sky, 16.5)[1], sky)
    frames, audit, sources = [], [], []
    for visit in (OLD, NEW):
        mosaic = selections[band, visit]
        receipt = native_products(cache, SimpleNamespace(pointing="001"), band, mosaic, sky)
        metadata = native_metadata(receipt, visit, band, mosaic["detector"])
        with catch_warnings():
            filterwarnings("ignore", message="'(datfix|obsfix)' made the change", category=FITSFixedWarning)
            visit_frames, diagnostics = native_alignment(base, receipt, sky)
        if len(visit_frames) != 3:
            raise ValueError(f"Expected three science exposures per visit in {band}")
        frames.extend(visit_frames)
        audit.extend(diagnostics)
        sources.append(dict(visit=visit, fits_headers=metadata, **receipt))
    target, background, fluxes, calibration = apertures(frames, refs, controls, sky)
    report = compare(target, background, fluxes)
    record = dict(sources=sources, epoch_changes=report, relative_throughput=calibration,
                  native_exposures=[dict(**a, fractional_photometry=b) for a, b in zip(audit, fluxes, strict=True)])
    views = None
    if band in VIEW_BANDS:
        views = stamps(frames, audit)
        sigma = fit_data[band]["astrometry"]["host_core"]["core_sigma_arcsec"]
        record["zero_dq_core_fits"] = [zero_dq_core(frames, audit, sigma*factor, BANDS[band], fluxes)
                                       for factor in (.85, 1., 1.15)]
    return record, views


def main():
    photo_path = OUT/"module_b_001_six_filter_change.json"
    prior_path = OUT/"module_b_products.json"
    data, search = (json.loads(path.read_text()) for path in (photo_path, prior_path))
    fit_path = OUT/"module_b_001_mast_comparison.json"
    fit_data = json.loads(fit_path.read_text())["bands"]
    sky = data["icrs_deg"]
    mosaics = [r for r in search["images"] if r["pointing"] == "001"
               and r["detector"] == ("NRCB4" if r["band"]=="F200W" else "NRCBLONG")]
    mosaics.extend(r for r in data["additional_products"] if r["visit"] in {OLD, NEW})
    selections = {(m["band"], m["visit"]): m for m in mosaics}
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    records, views = {}, {}
    for band in FILTERS:
        records[band], views[band] = measure_band(cache, band, selections, data, fit_data)
        print(band, records[band]["epoch_changes"][-1], flush=True)
    plot(records, views)
    write_json(OUT/"module_b_001_native_epoch_flux.json", dict(
        executed_utc=utcnow(), candidate=data["candidate"], icrs_deg=sky, native_band_measurements=records,
        photometry_subpixels_per_axis=SUBPIXELS, minimum_usable_aperture_fraction=VALID_FRACTION,
        primary_allowed_pixel_flags=[0, JUMP_DET], fixed_core_widths_from_sha256=text_digest(fit_path),
        storage_cap_mib=STORAGE_MIB, six_filter_sha256=text_digest(photo_path), search_products_sha256=text_digest(prior_path),
        script_sha256=text_digest(__file__),
        scope="Calibrated native SCI pixels. Primary photometry allows valid ramp fits marked only JUMP_DET."
              " The DQ-zero aperture and core measurements test stricter masking. No missing-pixel imputation."
              " DAWN positions of local references anchor the translations and relative throughput."
              " Dither-pair combinations share exposures.",
    ))


if __name__ == "__main__":
    main()

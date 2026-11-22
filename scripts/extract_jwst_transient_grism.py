"""Force a detector-space extraction at the transient and neighboring traces."""

import json
from contextlib import ExitStack

import asdf
import numpy as np
from astropy.io import fits
from download_jwst_cross_detector import OUT
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import verified_path
from measure_jwst_transient_grism_geometry import (
    dispersed_trace,
    imaging_transform,
    pointing_header,
)
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

BINS = np.linspace(3.18, 3.94, 15)
CONTROL_OFFSETS = (-34, -26, -18, 18, 26, 34)
BACKGROUND_RADIUS = 12
CORE_RADIUS = 2.2


def trace_for_column(record, calibrations, target, neighbour):
    _, sci = pointing_header(verified_path(record))
    with ExitStack() as stack:
        refs = {k: stack.enter_context(asdf.open(verified_path(calibrations[v])))
                for k, v in record["references"].items()}
        (offset,) = [item for item in refs["R_FILOFF"]["filters"]
                     if (item["filter"], item["pupil"]) == ("F356W", "GRISMR")]
        wcs = imaging_transform(sci, refs["R_DISTOR"]["model"], offset)
        order = refs["R_SPCWCS"]
        x, y, wave = dispersed_trace(order, 1, *wcs.inverse(*target), sci["VELOSYS"])
        hx, hy, _ = dispersed_trace(order, 1, *wcs.inverse(*neighbour), sci["VELOSYS"])
    cols = np.arange(max(4, int(np.ceil(np.interp(BINS[0], wave, x)))),
                     min(2044, int(np.floor(np.interp(BINS[-1], wave, x)))+1))
    return cols, np.interp(cols, x, wave), np.interp(cols, x, y), np.interp(cols, hx, hy)


def flat_strip(rate, corrections, columns, center):
    margin = BACKGROUND_RADIUS + max(map(abs, CONTROL_OFFSETS)) + 2
    low, high = int(np.floor(center.min()))-margin, int(np.ceil(center.max()))+margin+1
    if not (4 <= low < high <= 2044):
        raise ValueError("Diagnostic cross-dispersion strip exceeds the science detector")
    region = np.s_[low:high, columns[0]:columns[-1]+1]
    with fits.open(verified_path(rate), memmap=False) as h, fits.open(
            verified_path(corrections[rate["references"]["R_FLAT"]]), memmap=False) as flat:
        if h["SCI"].header["BUNIT"] != "DN/s":
            raise ValueError("Unexpected unit for the rate image")
        if h[0].header["EXP_TYPE"] != "NRC_WFSS" or h[0].header["PUPIL"] != "GRISMR":
            raise ValueError("Unexpected mode for the detector image")
        divisor = flat["SCI"].data[region]
        good = (h["DQ"].data[region] == 0) & (flat["DQ"].data[region] == 0) & (divisor > 0)
        sci = np.where(good, h["SCI"].data[region], np.nan)
        result = np.full(sci.shape, np.nan)
        np.divide(sci, divisor, out=result, where=good)
    return result, np.arange(low, high)


def extract_column(values, rows, target_y, galaxy_y, degree):
    """Sum the fixed core after a side-window fit, leaving contaminated pixels missing."""
    dy = rows-target_y
    usable = np.isfinite(values) & (abs(rows-galaxy_y) >= 3.5)
    keep = usable & (abs(dy) >= 4) & (abs(dy) <= BACKGROUND_RADIUS)
    core = abs(dy) <= CORE_RADIUS
    if not usable[core].all() or min(np.sum(keep & (dy < 0)), np.sum(keep & (dy > 0))) < degree+2:
        return np.nan
    for _ in range(2):
        fit = np.polynomial.polynomial.polyfit(dy[keep], values[keep], degree)
        residual = values-np.polynomial.polynomial.polyval(dy, fit)
        base = residual[keep]
        center = np.median(base)
        scatter = 1.4826*np.median(abs(base-center))
        if scatter <= 0:
            break
        limited = keep & (abs(residual-center) <= 3.5*scatter)
        if min(np.sum(limited & (dy < 0)), np.sum(limited & (dy > 0))) < degree+2:
            break
        keep = limited
    return float(residual[core].sum())


def bin_rates(wavelength, rate):
    result, samples = [], []
    for a, b in zip(BINS[:-1], BINS[1:], strict=True):
        d = rate[(wavelength >= a) & (wavelength < b) & np.isfinite(rate)]
        samples.append(len(d))
        result.append(float(np.median(d)*1000) if len(d) >= 15 else np.nan)
    return np.array(result), samples


def check_one(frame, headers, models, refs, target, galaxy):
    cols, wave, center, nearest = trace_for_column(headers, models, target, galaxy)
    science, rows = flat_strip(frame, refs, cols, center)
    offsets = (0,)+CONTROL_OFFSETS
    by_method, sample_count = [], []
    for degree in (1, 2):
        fit_vectors = []
        for offset in offsets:
            raw = np.array([extract_column(values, rows, y+offset, h, degree)
                            for values, y, h in zip(science.T, center, nearest, strict=True)])
            b, n = bin_rates(wave, raw)
            fit_vectors.append(b)
            if offset == 0:
                sample_count.append(n)
        by_method.append(fit_vectors)
    # Select observed pixels with the nearest integer row for each trace offset.
    display_y = np.arange(-12, 13)[:, None]
    sampled = science[np.rint(center).astype(int)+display_y-rows[0], np.arange(science.shape[1])]
    sampled = sampled-np.nanmedian(sampled, axis=0)
    grid_wave = np.linspace(BINS[0], BINS[-1], 749)
    strips = [np.interp(grid_wave, wave, d, left=np.nan, right=np.nan) for d in sampled]
    print(frame["filename"], "host relative trace offset", round(float(np.median(nearest-center)), 2),
          "target samples", sample_count[0], flush=True)
    return np.array(by_method), sample_count, np.array(strips), float(np.median(nearest-center))


def to_values(array):
    return np.where(np.isfinite(array), array, None).tolist()


def plot_diagnostic(strips, result, offsets):
    fig = Figure(figsize=(9.4, 6.6), layout="constrained")
    axes = fig.subplots(3, 1, gridspec_kw={"height_ratios": [1.45, 1, 1]}, sharex=True)
    stack = np.nanmedian(strips, axis=0)
    extent = (BINS[0], BINS[-1], -12.5, 12.5)
    lim = np.nanpercentile(abs(stack), 97)
    axes[0].imshow(stack, origin="lower", aspect="auto", extent=extent, cmap="RdBu_r", vmin=-lim, vmax=lim)
    axes[0].axhline(0, color="#247474", linewidth=1, label="Predicted transient trace")
    axes[0].axhline(np.median(offsets), color="#d6a13c", linewidth=1, label="Closest sky match in archive")
    axes[0].set(ylabel="Cross-dispersion pixel", title="Observed rate residuals aligned on the predicted trace")
    axes[0].legend(loc="lower right", fontsize=8)
    centers = (BINS[:-1]+BINS[1:])/2
    for i, ax in enumerate(axes[1:]):
        raw = np.nanmedian(result[:, i], axis=0)
        for ctrl in raw[1:]:
            ax.plot(centers, ctrl, color="#a7a9ae", alpha=.48, lw=.7)
        for j in range(2):
            half = np.nanmedian(result[j*3:j*3+3, i, 0], axis=0)
            ax.plot(centers, half, ".-", color=("#a5574c", "#476c9c")[j], lw=.8,
                    label=("Exposures 1–3", "Exposures 4–6")[j])
        ax.axhline(0, color="#878c8b", linewidth=.65)
        ax.plot(centers, raw[0], "o-", color="#183f42", ms=3, lw=1.2, label="All exposures")
        ax.set_ylabel("mDN/s per column")
        ax.set_title(("Linear" if i == 0 else "Quadratic")+" local background · grey traces show spatial controls",
                     fontsize=9)
        ax.legend(loc="upper right", ncols=3, fontsize=7)
    axes[-1].set_xlabel("Wavelength along transient trace (µm) · 54 nm median bins")
    fig.suptitle("GOODS-S transient · forced detector extraction\n"
                 "Six GO 4549 GRISMR exposures · flux calibration and contamination unresolved", fontsize=13)
    save_figure(fig, FIGURES / "field_b_transient_forced_grism")


def main():
    rate_path = OUT / "field_b_transient_grism_rates.json"
    receipt = json.loads(rate_path.read_text())
    pointings_path = OUT / "field_b_transient_grism_products.json"
    pointings = json.loads(pointings_path.read_text())
    audit_path = OUT / "field_b_transient_wfss.json"
    audit = json.loads(audit_path.read_text())
    if len(receipt["native_rate_frames"]) != 6 or pointings["audit_sha256"] != text_digest(audit_path):
        raise ValueError("Require six matching public rate exposures and consistent source provenance")
    headers = {r["rate_filename"]: r for r in pointings["rate_headers"]}
    nearest = {r["filename"].replace("_x1d.fits", "_rate.fits"): r["nearest"]["icrs_deg"]
               for r in audit["x1d_exposures"]}
    frames, all_counts, strips, offsets = [], [], [], []
    for rate in receipt["native_rate_frames"]:
        name = rate["filename"]
        result, n, pixels, offset = check_one(rate, headers[name], pointings["calibrations"],
                                              receipt["calibrations"], receipt["position_degrees"], nearest[name])
        frames.append(result)
        strips.append(pixels)
        offsets.append(offset)
        all_counts.append(dict(filename=name, host_relative_trace_pixels=offset, usable_target_columns=n))
    frames = np.array(frames)
    result = np.nanmedian(frames, axis=0)
    control = np.nanmedian(result[:, 1:], axis=1)
    scatter = 1.4826*np.nanmedian(abs(result[:, 1:]-control[:, None]), axis=1)
    plot_diagnostic(strips, frames, offsets)
    write_json(OUT / "field_b_transient_forced_grism.json", {
        "executed_utc": utcnow(), "candidate": receipt["candidate"], "wavelength_bin_edges_um": BINS.tolist(),
        "model_order": [1, 2], "trace_offsets_pix": [0, *CONTROL_OFFSETS],
        "detector_rate_mdn_per_s": to_values(frames), "stack_median_mdn_per_s": to_values(result),
        "spatial_control_mad_mdn_per_s": to_values(scatter), "frame_coverage": all_counts,
        "scope": "Flat-corrected detector rates; no point-spread aperture correction or contamination model",
        "rate_products_sha256": text_digest(rate_path), "pointings_sha256": text_digest(pointings_path),
        "wfss_audit_sha256": text_digest(audit_path), "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

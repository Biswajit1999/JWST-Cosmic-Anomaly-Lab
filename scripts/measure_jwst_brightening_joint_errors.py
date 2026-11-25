"""Estimate shared-epoch photometric uncertainties from matched local apertures."""

import json

import numpy as np
from astropy.coordinates import SkyCoord
from compare_jwst_brightening_mast import REGION_ARCSEC, common_grid
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from measure_jwst_brightening_november import EPOCHS, FILTERS, STORAGE_MIB
from measure_jwst_brightening_offset import BANDS, convolve, local_alignment
from measure_jwst_epoch_candidate import aperture, pixel_scale, show_stamp, verified_path
from measure_jwst_module_b_brightening import ANNULUS
from read_jwst_i2d_cutout import MastImageSections
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

PREFIX = "module_b_001_joint_epoch_errors"
RADII = (.18, .24, .30)
SMOOTHING = (.014, .034)
CONTRASTS = np.array([[-1, 1, 0], [-1, 0, 1], [0, -1, 1]])
NAMES = ("Oct 2026 - Oct 2025", "Nov 2026 - Oct 2025", "Nov 2026 - Oct 2026")
SEED = 103607


class PinnedSections(MastImageSections):
    """Use the exact previously inventoried byte ranges and verify their checksums."""

    def __init__(self, cache, receipt):
        self.cache = cache
        self.filename = receipt["filename"]
        self.size = receipt["object_bytes"]
        self.url = receipt["url"]
        self.etag = receipt["etag"]
        self.ranges = []
        self.expected = {(r["first_byte"], r["last_byte"]): r for r in receipt["ranges"]}

    def read(self, start, stop):
        expected = self.expected[start, stop-1]
        if not (ROOT/expected["path"]).is_file():
            super().read(start, stop)
        return verified_path(expected)


def joint_images(stamps, band, smoothing):
    reference = stamps[0]
    _, wcs, x, y, _, _ = reference
    scale = pixel_scale(wcs)
    images = [common_grid(reference, reference)[0]]
    shifts = [np.zeros(2)]
    gains = [1.]
    relative = [0.]
    registration = []
    for stamp in stamps[1:]:
        a, b = common_grid(reference, stamp)
        reg = local_alignment(a, b, x, y, scale, BANDS.get(band, .22))
        images.append(b)
        shifts.append(-np.array(reg["old_to_new_xy_arcsec"][::-1])/scale)
        gains.append(reg["new_to_old_scale"])
        relative.append(reg["matched_variance_arcsec2"])
        registration.append(reg)
    blur = np.sqrt(max(relative) - np.array(relative) + smoothing**2)
    matched = [convolve(im, shift, sigma/scale)/gain
               for im, shift, sigma, gain in zip(images, shifts, blur, gains, strict=True)]
    return matched, (x, y, scale), registration, blur.tolist()


def control_apertures(images, x, y, scale, radius):
    step = round(.61/scale)
    shifts = range(-8*step, 9*step, step)
    # Use a fixed annulus of control centers, including real residual sky structure.
    offsets = [(dx, dy) for dx in shifts for dy in shifts if 1.5 <= np.hypot(dx, dy)*scale <= 5.]
    target = np.array([aperture(im, x, y, scale, radius, ANNULUS) for im in images])
    rows = np.array([[aperture(im, x+dx, y+dy, scale, radius, ANNULUS) for im in images]
                     for dx, dy in offsets])
    good = np.isfinite(rows).all(axis=1)
    if not np.isfinite(target).all() or good.sum() < 65:
        raise ValueError("Target apertures or matched local controls have insufficient coverage")
    return target, rows[good], np.asarray(offsets)[good]*scale


def robust_scale(rows, axis=0):
    return 1.4826*np.median(np.abs(rows-np.median(rows, axis=axis, keepdims=True)), axis=axis)


def summarize(target, control, offsets, rng):
    signal = CONTRASTS@target
    samples = control@CONTRASTS.T
    scatter = robust_scale(samples)
    median = np.median(samples, axis=0)
    # The winsorization is applied in epoch space so all derived contrasts
    # share a consistent covariance, including exact algebraic dependence.
    centre = np.median(control, axis=0)
    spread = robust_scale(control)
    joint = np.cov(np.clip(control, centre-4*spread, centre+4*spread), rowvar=False)
    covariance = CONTRASTS@joint@CONTRASTS.T
    block_labels = np.floor(offsets/1.5).astype(int)
    _, membership = np.unique(block_labels, axis=0, return_inverse=True)
    blocks = [np.flatnonzero(membership == i) for i in range(membership.max()+1)]
    replicates = []
    for draw in rng.integers(len(blocks), size=(800, len(blocks))):
        chosen = np.concatenate([blocks[i] for i in draw])
        replicates.append(robust_scale(samples[chosen]))
    errors = np.quantile(replicates, [.16, .84], axis=0)
    # Compare the rise and fade to the same empirical positions without
    # interpreting ranks from a selected source as a blind detection probability.
    joint_score = np.minimum((samples[:, 0]-median[0])/scatter[0], -(samples[:, 2]-median[2])/scatter[2])
    score = min((signal[0]-median[0])/scatter[0], -(signal[2]-median[2])/scatter[2])
    return dict(
        signed_contrasts_njy=signal.tolist(), local_mad_njy=scatter.tolist(),
        local_median_njy=median.tolist(), scatter_block_bootstrap_16_84_njy=errors.tolist(),
        winsorized_contrast_covariance_njy2=covariance.tolist(),
        target_aperture_njy=target.tolist(), control_count=len(control), spatial_block_count=len(blocks),
        controls_at_least_target_rise_and_fade_score=int((joint_score >= score).sum()),
        target_rise_and_fade_score=float(score), control_offsets_arcsec=offsets.tolist(),
        control_epoch_fluxes_njy=control.tolist(),
    )


def figure(bands):
    palette = ("#245c82", "#af5029", "#76588b", "#24734b")
    fig = Figure(figsize=(7.1, 3.85), layout="constrained")
    ax = fig.add_subplot()
    for i, (band, colour) in enumerate(zip(FILTERS, palette, strict=True)):
        row = bands[band]["primary"]
        ax.errorbar(np.arange(3)+(i-1.5)*.13, row["signed_contrasts_njy"], yerr=row["local_mad_njy"],
                    fmt="o", markersize=4, capsize=2.5, color=colour, label=band)
    ax.axhline(0, color="#526875", lw=.7)
    ax.grid(axis="y", alpha=.2)
    ax.set(xticks=range(3), xticklabels=["Oct 2026\nminus Oct 2025", "Nov 2026\nminus Oct 2025",
                                      "Nov 2026\nminus Oct 2026"],
           ylabel="Signed excess in 0.24 arcsec aperture (nJy)")
    ax.legend(loc="upper right", frameon=False, ncols=2)
    save_figure(fig, FIGURES/PREFIX)


def f200w_images(images, centre):
    fig = Figure(figsize=(7.2, 1.85), layout="constrained")
    old, bright, late = images
    x, y, scale = centre
    half = .49
    n = round(half/scale)
    ix, iy = round(x), round(y)
    crop = np.s_[iy-n:iy+n+1, ix-n:ix+n+1]
    image_range = np.nanpercentile(np.concatenate([im[crop].ravel() for im in images]), [4, 98.5])
    rise, fade = bright-old, late-bright
    residual_limit = np.nanpercentile(np.abs(np.r_[rise[crop].ravel(), fade[crop].ravel()]), 98)
    scenes = [
        (old, "Oct 2025", "gray", image_range),
        (bright, "Oct 2026", "gray", image_range),
        (late, "Nov 2026", "gray", image_range),
        (rise, "Oct 2026 - 2025", "RdBu_r", (-residual_limit, residual_limit)),
        (fade, "Nov - Oct 2026", "RdBu_r", (-residual_limit, residual_limit)),
    ]
    for ax, (image, label, cmap, bounds) in zip(fig.subplots(1, 5), scenes, strict=True):
        show_stamp(ax, image, x, y, scale, bounds, cmap, .24, half_arcsec=half)
        ax.set_title(label, fontsize=8)
    save_figure(fig, FIGURES/"module_b_001_preprint_f200w")


def main():
    prior = OUT/"module_b_001_november_photometry.json"
    report = json.loads(prior.read_text(encoding="utf-8"))
    sky = SkyCoord(*report["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    rng = np.random.default_rng(SEED)
    bands = {}
    for band in FILTERS:
        records = {(r["epoch"]): r for r in report["official_mosaics"] if r["filter"] == band}
        stamps = [PinnedSections(cache, records[epoch]).cutout(sky, REGION_ARCSEC)[0] for epoch in EPOCHS]
        rows = {}
        for smoothing in SMOOTHING:
            images, centre, registration, blur = joint_images(stamps, band, smoothing)
            for radius in RADII:
                target, controls, offsets = control_apertures(images, *centre, radius)
                rows[f"{smoothing:.3f}/{radius:.2f}"] = summarize(target, controls, offsets, rng)
            if smoothing == SMOOTHING[0]:
                if band == "F200W":
                    f200w_images(images, centre)
                primary = rows["0.014/0.24"]
                bands[band] = dict(primary=primary, registration=registration, Gaussian_sigma_arcsec=blur)
                print(band, primary["signed_contrasts_njy"], primary["local_mad_njy"],
                      primary["control_count"], flush=True)
        bands[band]["alternative_apertures"] = {k: v for k, v in rows.items() if k != "0.014/0.24"}
    figure(bands)
    write_json(OUT/f"{PREFIX}.json", dict(
        executed_utc=utcnow(), icrs_deg=report["icrs_deg"], prior_photometry_sha256=text_digest(prior),
        script_sha256=text_digest(__file__), epoch_labels=EPOCHS, contrast_labels=NAMES,
        aperture_radius_arcsec=.24, sky_annulus_arcsec=ANNULUS, bootstrap_seed=SEED,
        bands=bands,
        scope="Common spatial resolution and relative throughput for all three epochs. Local aperture MAD bars."
              " A block bootstrap describes the precision of the scatter estimate."
              " Covariances are empirical and contain shared epochs."
              " No search-wide significance or calibration completeness is implied.",
    ))


if __name__ == "__main__":
    main()

"""Test recovery of synthetic point sources added to observed CANDELS epoch images."""

import json

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from measure_jades_103607_hst_history import (
    ANNULUS,
    GRID,
    RADIUS,
    SCALE,
    SOURCE,
    control_centers,
    difference,
    epoch_start,
    load_series,
)
from scipy.ndimage import map_coordinates
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_hst_injection_recovery"
PSF_ROOT = (
    "https://www.stsci.edu/files/live/sites/www/files/home/hst/"
    "instrumentation/wfc3/data-analysis/psf/_documents"
)
BANDS = ("F125W", "F160W")
FLUXES_NJY = (0, 15, 30, 45, 60, 80, 100, 125, 150, 200, 250, 300)
PSF_SAMPLE_ARCSEC = 0.128 / 4
ROTATIONS_DEG = (0, 45, 90)
PHASES_PIXEL = ((0, 0), (0.35, -0.35), (-0.35, 0.35))
THRESHOLD = 5.0
PSF_HALF_PIXELS = 39
APERTURE_HALF = int(np.ceil(ANNULUS[1] / SCALE)) + 1
LOCAL = np.arange(-APERTURE_HALF, APERTURE_HALF + 1)
DISTANCE = np.hypot(LOCAL[:, None], LOCAL[None, :]) * SCALE
SOURCE_MASK = DISTANCE <= RADIUS
SKY_MASK = (DISTANCE >= ANNULUS[0]) & (DISTANCE < ANNULUS[1])


def psf_templates(cache, band):
    name = f"PSFSTD_WFC3IR_{band}.fits"
    url = f"{PSF_ROOT}/{name}"
    path = cache.fetch(
        url,
        f"instrument/{name}",
        survey="HST WFC3 calibration",
        version="Anderson ISR 2016-12",
        kind="observed_empirical_psf_reference",
    )
    with fits.open(path, memmap=False) as handle:
        cube = handle[0].data.astype(float)
    if cube.shape != (9, 101, 101) or not np.isfinite(cube).all():
        raise ValueError("Unexpected WFC3/IR PSF cube")
    coords = np.arange(-PSF_HALF_PIXELS, PSF_HALF_PIXELS + 1)
    yy, xx = np.meshgrid(coords, coords, indexing="ij")
    scale_ratio = SCALE / PSF_SAMPLE_ARCSEC
    bound = np.s_[PSF_HALF_PIXELS - APERTURE_HALF : PSF_HALF_PIXELS + APERTURE_HALF + 1]
    templates = []
    for angle in ROTATIONS_DEG:
        theta = np.deg2rad(angle)
        cosine, sine = np.cos(theta), np.sin(theta)
        for dx, dy in PHASES_PIXEL:
            xi, yi = (xx - dx) * scale_ratio, (yy - dy) * scale_ratio
            px = 50 + xi * cosine + yi * sine
            py = 50 - xi * sine + yi * cosine
            for model in cube:
                sampled = np.maximum(map_coordinates(model, (py, px), order=1, prefilter=False), 0)
                sampled /= sampled.sum()
                templates.append(sampled[bound, bound])
    return np.asarray(templates), receipt(path, url)


def signed_aperture(stamps):
    counts = stamps[..., SOURCE_MASK]
    sky = stamps[..., SKY_MASK]
    usable = np.isfinite(counts).all(axis=-1) & (
        np.isfinite(sky).sum(axis=-1) >= 0.8 * sky.shape[-1]
    )
    if not usable.all():
        raise ValueError("Insufficient coverage in an injection aperture")
    return counts.sum(axis=-1) - counts.shape[-1] * np.nanmedian(sky, axis=-1)


def measure_epoch(image, baseline, templates):
    n = len(GRID) // 2
    stamps = np.asarray(
        [
            image[
                y - APERTURE_HALF : y + APERTURE_HALF + 1, x - APERTURE_HALF : x + APERTURE_HALF + 1
            ]
            for x, y in [(n, n), *control_centers()]
        ]
    )
    valid = np.isfinite(stamps[..., SOURCE_MASK]).all(axis=-1) & (
        np.isfinite(stamps[..., SKY_MASK]).sum(axis=-1) >= 0.8 * SKY_MASK.sum()
    )
    if not valid[0]:
        raise ValueError("Injection target has insufficient image coverage")
    stamps = stamps[valid]
    comparison = signed_aperture(stamps)
    background = baseline["control_median_njy"]
    scatter = baseline["empirical_scatter_njy"]
    if len(stamps) != len(baseline["controls_njy"]) + 1 or not np.allclose(
        comparison,
        [baseline["target_minus_leave_one_out_njy"], *baseline["controls_njy"]],
        rtol=1e-6,
        atol=1e-5,
    ):
        raise ValueError("Original CANDELS aperture photometry did not reproduce")
    measurements = []
    for flux in FLUXES_NJY:
        perturbed = stamps[:, None] + flux * templates[None, :]
        possible = signed_aperture(perturbed) - background
        detected = possible >= THRESHOLD * scatter
        increment = (possible[0] + background) - comparison[0]
        measurements.append(
            dict(
                injected_model_total_njy=flux,
                center_recovered_increment_njy_range=[
                    float(increment.min()),
                    float(increment.max()),
                ],
                center_recovered_increment_njy_median=float(np.median(increment)),
                center_detections=int(detected[0].sum()),
                center_trials=int(detected.shape[1]),
                control_detections=int(detected[1:].sum()),
                control_trials=int(detected[1:].size),
            )
        )
    return dict(
        epoch=baseline["epoch"],
        date_start_utc=epoch_start(baseline),
        baseline_aperture_njy=baseline["differential_njy"],
        baseline_scatter_njy=scatter,
        control_sites=len(stamps) - 1,
        control_grid_indices=np.flatnonzero(valid[1:]).tolist(),
        synthetic_recovery=measurements,
    )


def analyze(cache, sky, band):
    baseline_path = OUT / f"module_b_001_hst_{band.lower()}_history.json"
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    rows, images = load_series(cache, sky, band)
    if [r["ranges"] for r in rows] != [r["ranges"] for r in baseline["coverage"]]:
        raise ValueError("Source FITS byte-range provenance differs from the baseline")
    models, psf = psf_templates(cache, band)
    observed = [r for r in baseline["coverage"] if r["center_weight"] > 0]
    epochs = [measure_epoch(difference(images, i), row, models) for i, row in enumerate(observed)]
    summary = []
    for i, flux in enumerate(FLUXES_NJY):
        all_results = [e["synthetic_recovery"][i] for e in epochs]
        controls_detected = sum(row["control_detections"] for row in all_results)
        controls_tested = sum(row["control_trials"] for row in all_results)
        fit_by_epoch = [row["control_detections"] / row["control_trials"] for row in all_results]
        summary.append(
            dict(
                injected_model_total_njy=flux,
                control_detection_fraction=controls_detected / controls_tested,
                control_fraction_by_epoch_range=[min(fit_by_epoch), max(fit_by_epoch)],
                center_template_detection_fraction=sum(r["center_detections"] for r in all_results)
                / sum(r["center_trials"] for r in all_results),
            )
        )
    return dict(
        baseline_file=baseline_path.name,
        baseline_sha256=text_digest(baseline_path),
        psf_reference=psf,
        model_template_count=len(models),
        epochs=epochs,
        synthetic_summary=summary,
        first_tested_flux_at_90_percent_recovery_njy=next(
            (
                row["injected_model_total_njy"]
                for row in summary
                if row["control_detection_fraction"] >= 0.90
            ),
            None,
        ),
    )


def paired_observed_history(bands):
    """Test for excess in both bands at the same positions and nominal epochs."""
    rows_by_band = {}
    grids_by_band = {}
    for band in BANDS:
        path = OUT / f"module_b_001_hst_{band.lower()}_history.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        grids_by_band[band] = {
            row["epoch"]: row["control_grid_indices"] for row in bands[band]["epochs"]
        }
        rows_by_band[band] = {
            row["epoch"]: row for row in data["coverage"] if row["center_weight"] > 0
        }
    paired = sorted(set(rows_by_band[BANDS[0]]) & set(rows_by_band[BANDS[1]]))
    if not paired:
        raise ValueError("No matched CANDELS epochs")
    common_sites = sorted(
        set.intersection(*[set(grids_by_band[band][epoch]) for band in BANDS for epoch in paired])
    )
    if len(common_sites) < 20:
        raise ValueError("Insufficient matched control sites across paired epochs")
    target, null_by_epoch = [], []
    for epoch in paired:
        scores, null = {}, []
        for band in BANDS:
            row = rows_by_band[band][epoch]
            scale = row["empirical_scatter_njy"]
            grid = dict(zip(grids_by_band[band][epoch], row["controls_njy"], strict=True))
            controls = np.array([grid[index] for index in common_sites])
            scores[band] = row["differential_njy"] / scale
            null.append((controls - row["control_median_njy"]) / scale)
        target.append(
            dict(epoch=epoch, residual_over_scatter=scores, min_score=min(scores.values()))
        )
        null_by_epoch.append(np.minimum(*null))
    strongest = max(target, key=lambda row: row["min_score"])
    null = np.max(null_by_epoch, axis=0)
    return dict(
        statistic="Maximum across paired epochs of the smaller same-epoch signed residual "
        "over its own band's control scatter.",
        epochs=target,
        strongest_target_pair=strongest,
        control_site_count=len(null),
        common_grid_indices=common_sites,
        controls_meeting_or_exceeding_target=int(np.count_nonzero(null >= strongest["min_score"])),
        control_max_scores=null.tolist(),
        scope="Identical spatial controls in both bands and every tested epoch. Their rank "
        "provides context, with shared reference data and no independent-trial probability.",
    )


def plot(bands):
    fig = Figure(figsize=(9.2, 4.7), layout="constrained")
    axes = fig.subplots(1, 2, sharex=True, sharey=True)
    for axis, (band, result) in zip(axes, bands.items(), strict=True):
        rows = result["synthetic_summary"]
        limits = np.asarray([r["control_fraction_by_epoch_range"] for r in rows]).T
        axis.fill_between(FLUXES_NJY, *limits, color="#256982", alpha=0.18, label="Epoch range")
        axis.plot(
            FLUXES_NJY,
            [r["control_detection_fraction"] for r in rows],
            ".-",
            color="#256982",
            label="Nearby control sites",
        )
        axis.plot(
            FLUXES_NJY,
            [r["center_template_detection_fraction"] for r in rows],
            "s--",
            ms=3,
            color="#b4572d",
            label="Actual target residuals",
        )
        axis.axvline(30, color=".45", linestyle=":", linewidth=1)
        axis.set(title=band, xlim=(0, 300), ylim=(-0.025, 1.025))
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False, loc="lower right", fontsize=8)
    axes[0].set_ylabel("Fraction above the five-scatter threshold")
    fig.supxlabel("Artificial source model flux density in each HST band (nJy)")
    fig.suptitle("JADES 103607 | artificial sources added to observed HST differences")
    save_figure(fig, FIGURES / TAG)


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    bands = {band: analyze(cache, sky, band) for band in BANDS}
    plot(bands)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            position_source_sha256=text_digest(SOURCE),
            historical_reduction_script_sha256=text_digest(
                ROOT / "scripts/measure_jades_103607_hst_history.py"
            ),
            observation_type="Observed CANDELS mosaics and simulated source increments",
            rotations_deg=ROTATIONS_DEG,
            subpixel_offsets_output_pixels=PHASES_PIXEL,
            threshold_in_local_scatter=THRESHOLD,
            bands=bands,
            paired_observed_history=paired_observed_history(bands),
            scope="A source in one already drizzled epoch adds linearly to its leave-one-out "
            "difference. Orientation, detector position and centroid are varied using a grid "
            "of resampled WFC3 empirical PSFs. The flux is normalized over the sampled model "
            "and applied one site at a time. Recovery fractions are descriptive paired trials "
            "on observed residuals. The trials share data. They omit source shot noise, native "
            "exposure injection, color conversion, light-curve dilution and unobserved dates.",
        ),
    )


if __name__ == "__main__":
    main()

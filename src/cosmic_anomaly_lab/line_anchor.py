"""Observed narrow-line widths and conditional instrument bounds, not unique LSFs."""

import numpy as np
from scipy.optimize import least_squares

from .line_profiles import C_KMS, profile_bins
from .spectroscopy import pixel_edges

OIII = (0.4960295, 0.5008240)  # Vacuum wavelengths, microns.


def doublet_model(wave, p):
    """Two independent nonnegative areas, common observed Gaussian width/shift."""
    continuum, slope, area1, area2, width, velocity = p
    result = continuum + slope * (wave - np.mean(OIII)) / 0.01
    for center, area in zip(OIII, [area1, area2], strict=True):
        result = result + area * profile_bins(wave, center * (1 + velocity / C_KMS), width)
    return result


def profile_crossing(grid, delta, best, level, side):
    """Interpolate the first level crossing outwards from the fitted minimum."""
    chosen = np.flatnonzero(grid >= best) if side == "upper" else np.flatnonzero(grid <= best)[::-1]
    for a, b in zip(chosen[:-1], chosen[1:], strict=False):
        if delta[a] <= level <= delta[b]:
            fraction = (level - delta[a]) / (delta[b] - delta[a])
            return float(grid[a] + fraction * (grid[b] - grid[a]))
    return None


def fit_oiii(spec, redshift, half_window=900, profile=True):
    wave = spec.wavelength.value / (1 + redshift)
    region = (wave > OIII[0] * (1 - (half_window + 200) / C_KMS)) & (
        wave < OIII[1] * (1 + (half_window + 200) / C_KMS)
    )
    wave, flux, error, good = (
        wave[region],
        spec.flux.value[region],
        spec.error.value[region],
        spec.good[region],
    )
    if len(wave) < 3:
        return dict(status="NO_WAVELENGTH_COVERAGE")
    edges = pixel_edges(wave)
    windows = [(w * (1 - half_window / C_KMS), w * (1 + half_window / C_KMS)) for w in OIII]
    fit_mask = np.zeros(len(wave), bool)
    coverage, pixels = [], []
    for lo, hi in windows:
        overlap = np.maximum(0, np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo))
        mask = good & (overlap > 0)
        coverage.append(float(np.sum(overlap[good]) / (hi - lo)))
        pixels.append(int(mask.sum()))
        fit_mask |= mask
    base = dict(half_window_kms=half_window, coverage=coverage, valid_pixels=pixels)
    if min(coverage) < 0.8 or min(pixels) < 3 or fit_mask.sum() < 12:
        return dict(**base, status="INSUFFICIENT_COVERAGE")
    lower, upper = [-10, -30, 0, 0, 2, -400], [30, 30, 300, 300, 500, 400]
    scale = max(float(np.percentile(abs(flux[fit_mask]), 90)), float(np.median(error[fit_mask])))
    y, err = flux / scale, error / scale
    attempts = []
    for width in [60, 120, 220]:
        p = [float(np.clip(np.median(y[fit_mask]), -9, 29)), 0, 0.2, 0.6, width, 0]
        result = least_squares(
            lambda p: ((doublet_model(wave, p) - y) / err)[fit_mask],
            p,
            bounds=(lower, upper),
            x_scale="jac",
            max_nfev=600,
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-10,
        )
        if result.success:
            attempts.append(result)
    if not attempts:
        return dict(**base, status="NO_CONVERGED_FIT")
    best = min(attempts, key=lambda r: np.sum(r.fun**2))
    chi2 = float(np.sum(best.fun**2))
    p = best.x.copy()
    p[:4] *= scale
    n = int(fit_mask.sum())
    result = dict(
        **base,
        status="CONDITIONAL_GAUSSIAN_WIDTH",
        observed_sigma_kms=float(p[4]),
        chi2=chi2,
        n=n,
        k=6,
        reduced_chi2=chi2 / (n - 6),
        parameter_vector=p.tolist(),
        rest_wave_um=wave.tolist(),
        flux_uJy=[float(v) if np.isfinite(v) else None for v in flux],
        error_uJy=[float(v) if np.isfinite(v) else None for v in error],
        mask=fit_mask.tolist(),
        model_uJy=doublet_model(wave, p).tolist(),
    )
    if profile:
        grid = np.unique(np.r_[np.arange(2, 351, 1), best.x[4]])
        values = []
        nuisance = [0, 1, 2, 3, 5]
        for width in grid:

            def prediction(q):
                pars = best.x.copy()
                pars[4] = width
                pars[nuisance] = q
                return doublet_model(wave, pars)

            fit = least_squares(
                lambda q: ((prediction(q) - y) / err)[fit_mask],
                best.x[nuisance],
                bounds=(np.array(lower)[nuisance], np.array(upper)[nuisance]),
                x_scale="jac",
                max_nfev=400,
            )
            if not fit.success:
                raise ValueError("Unconverged width-profile point")
            values.append(float(np.sum(fit.fun**2)) - chi2)
        result.update(profile_sigma_kms=grid.tolist(), profile_delta_chi2=values)
        for level, tag in [(1, "delta1"), (3.84, "delta3p84")]:
            for side in ["lower", "upper"]:
                result[f"sigma_{tag}_{side}"] = profile_crossing(
                    grid, np.array(values), best.x[4], level, side
                )
    return result


def resolution_lower_bound(observed_sigma_upper):
    """R >= c/(FWHM_factor*sigma_obs_upper), conditional on Gaussian broadening.

    Intrinsic variance cannot be negative. No finite upper bound on resolving
    power follows from a single observed width with unknown intrinsic width.
    """
    if observed_sigma_upper is None:
        return None
    if not np.isfinite(observed_sigma_upper) or observed_sigma_upper <= 0:
        raise ValueError("Positive finite observed width required")
    return float(C_KMS / (2.354820045 * observed_sigma_upper))

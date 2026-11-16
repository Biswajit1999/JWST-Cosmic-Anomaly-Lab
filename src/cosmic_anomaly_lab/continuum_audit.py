"""Exploratory line-core measurements requiring continuum on both sides."""

import numpy as np

from .spectroscopy import pixel_edges

SIDEBANDS = {
    "legacy": ((-3000, -1500), (2000, 3000)),
    "symmetric_far": ((-3000, -2000), (2000, 3000)),
    "near": ((-1500, -800), (800, 1300)),
}


def bilateral_core(velocity, flux, error, good, center_angstrom, half_width=600,
                   sidebands=SIDEBANDS["legacy"]):
    """Signed core with >=3 valid continuum pixels on EACH side.

    Formal errors include the linear continuum covariance but not released-pixel
    covariance or systematic calibration errors. These are finite-window cores,
    not deblended total fluxes. Near windows may contain line wings.
    """
    v, f, e = (np.asarray(a, dtype=float) for a in (velocity, flux, error))
    g = np.asarray(good, dtype=bool).copy()
    if v.ndim != 1 or any(a.shape != v.shape for a in (f, e, g)):
        raise ValueError("Equal one-dimensional arrays are required")
    if len(v) < 2 or not np.isfinite(v).all() or np.any(np.diff(v) <= 0):
        raise ValueError("Velocity must be finite and strictly increasing")
    (lo, hi), (lo2, hi2) = sidebands
    if not (lo < hi < -half_width < half_width < lo2 < hi2 and center_angstrom > 0):
        raise ValueError("Continuum windows must bracket and exclude the positive-width core")
    g &= np.isfinite(f) & np.isfinite(e) & (e > 0)
    edges = pixel_edges(v)
    overlap = np.maximum(0, np.minimum(edges[1:], half_width)
                         - np.maximum(edges[:-1], -half_width))
    core = g & (overlap > 0)
    blue, red = g & (v >= lo) & (v <= hi), g & (v >= lo2) & (v <= hi2)
    coverage = float(overlap[g].sum() / (2 * half_width))
    result = dict(good_coverage=coverage, core_pixels=int(core.sum()),
                  blue_pixels=int(blue.sum()), red_pixels=int(red.sum()),
                  flux=None, error=None, continuum_intercept=None,
                  continuum_slope_per_1000kms=None, continuum_variance_fraction=None)
    if coverage < .8 or core.sum() < 4:
        return dict(result, status="INSUFFICIENT_CORE")
    if min(blue.sum(), red.sum()) < 3:
        return dict(result, status="INSUFFICIENT_BILATERAL_CONTINUUM")
    side = blue | red
    design = np.column_stack([np.ones_like(v), v / 1000])
    info = (design[side].T / e[side]**2) @ design[side]
    if np.linalg.cond(info) > 1e10:
        return dict(result, status="SINGULAR_CONTINUUM")
    cov = np.linalg.inv(info)
    beta = cov @ ((design[side].T / e[side]**2) @ f[side])
    w = overlap[core] * center_angstrom / 299792.458
    cw = w @ design[core]
    cv = cw @ cov @ cw
    variance = np.sum((w * e[core])**2) + cv
    return dict(result, status="MEASURED_CORE", flux=float(w @ f[core] - cw @ beta),
                error=float(np.sqrt(variance)), continuum_intercept=float(beta[0]),
                continuum_slope_per_1000kms=float(beta[1]),
                continuum_variance_fraction=float(cv / variance))

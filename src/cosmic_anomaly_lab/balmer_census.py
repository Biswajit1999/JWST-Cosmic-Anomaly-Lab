"""Catalogue diagnostics with explicit covariance and missing-data semantics."""

import re

import numpy as np

from .spectroscopy import pixel_edges

LINES = {"ha": "HBaA_6563", "hb": "HBaB_4861", "hg": "HBaG_4340"}


def normalize_filter(value):
    """Accept only exact released disperser names; never repair arbitrary bytes."""
    if isinstance(value, bytes):
        try:
            value = value.decode("ascii")
        except UnicodeDecodeError:
            return None
    value = str(value).strip()
    if re.fullmatch(r"g(?:140m_f(?:070|100)lp|235m_f170lp|395m_f290lp)(?:_ext)?", value):
        return value
    return None


def contrast_score(numerator, denominator, numerator_error, denominator_error,
                   reference, fractional_error=0.0, correlation=0.0):
    """Signed (N-rD)/error, a sensitivity diagnostic, NOT a calibrated z-score.

    Statistical cross-line correlation is explicit. Fractional systematic terms
    are independent per line for this sensitivity calculation, not measured errors.
    Negative input fluxes remain signed; missing/invalid errors return NaN.
    """
    if not (reference > 0 and fractional_error >= 0 and -1 <= correlation <= 1):
        raise ValueError("Invalid contrast assumptions")
    n, d, en, ed = np.broadcast_arrays(
        *map(lambda a: np.asarray(a, dtype=float),
             [numerator, denominator, numerator_error, denominator_error])
    )
    valid = np.isfinite(n) & np.isfinite(d) & np.isfinite(en) & np.isfinite(ed)
    valid &= (en > 0) & (ed > 0)
    variance = (en**2 + (reference * ed)**2 - 2 * reference * correlation * en * ed
                + fractional_error**2 * (n**2 + (reference * d)**2))
    valid &= variance > 0
    return np.divide(n - reference * d, np.sqrt(np.maximum(variance, 0)),
                     out=np.full(n.shape, np.nan), where=valid)


def aperture_minimum(values):
    """Require every extraction; missing aperture must not count as replication."""
    a = np.asarray(values, dtype=float)
    return np.where(np.isfinite(a).all(axis=0), np.min(a, axis=0), np.nan)


def signed_core_integral(velocity, flux_density, error, good, center_angstrom, half_width=600):
    """Signed finite-window line-core integral after a sideband continuum fit.

    Inputs are F_lambda in a common per-Angstrom unit. The returned integral
    keeps that unit times Angstrom. Released pixel covariance is unavailable;
    propagation includes the fitted continuum's diagonal-noise uncertainty.
    This is NOT a deblended total line flux or a broad/narrow decomposition.
    """
    v, f, e = map(lambda a: np.asarray(a, dtype=float), [velocity, flux_density, error])
    valid = np.asarray(good, dtype=bool) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    edges = pixel_edges(v)
    overlap = np.maximum(0, np.minimum(edges[1:], half_width) - np.maximum(edges[:-1], -half_width))
    core = valid & (overlap > 0)
    side = valid & (((v >= -3000) & (v <= -1500)) | ((v >= 2000) & (v <= 3000)))
    coverage = float(overlap[valid].sum() / (2 * half_width))
    result = dict(half_width_kms=half_width, good_coverage=coverage,
                  core_pixels=int(core.sum()), sideband_pixels=int(side.sum()))
    if coverage < .8 or core.sum() < 4 or side.sum() < 6:
        return dict(**result, status="INSUFFICIENT_COVERAGE", flux=None, error=None)
    design = np.column_stack([np.ones_like(v), v / 1000])
    information = (design[side].T / e[side]**2) @ design[side]
    if np.linalg.cond(information) > 1e10:
        return dict(**result, status="SINGULAR_CONTINUUM", flux=None, error=None)
    covariance = np.linalg.inv(information)
    beta = covariance @ ((design[side].T / e[side]**2) @ f[side])
    weights = overlap[core] * center_angstrom / 299792.458
    continuum_weights = weights @ design[core]
    integral = weights @ f[core] - continuum_weights @ beta
    variance = np.sum((weights * e[core])**2) + continuum_weights @ covariance @ continuum_weights
    return dict(**result, status="MEASURED_CORE", flux=float(integral), error=float(np.sqrt(variance)),
                continuum_intercept=float(beta[0]), continuum_slope_per_1000kms=float(beta[1]))

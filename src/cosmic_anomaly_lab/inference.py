"""Explicit Gaussian likelihoods, censoring, and predictive model separation."""

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.special import log_ndtr


def gaussian_loglike(data, model, *, sigma=None, covariance=None, calibration_fraction=0):
    data, model = np.asarray(data, float), np.asarray(model, float)
    if data.shape != model.shape or data.ndim != 1 or not np.all(np.isfinite(data + model)):
        raise ValueError("Finite matching one-dimensional arrays required")
    if calibration_fraction < 0:
        raise ValueError("Calibration uncertainty must be nonnegative")
    if (sigma is None) == (covariance is None):
        raise ValueError("Provide exactly one of sigma or covariance")
    if covariance is None:
        sigma = np.asarray(sigma, float)
        if sigma.shape != data.shape or np.any(sigma <= 0) or not np.all(np.isfinite(sigma)):
            raise ValueError("Positive finite errors required")
        cov = np.diag(sigma**2)
    else:
        cov = np.asarray(covariance, float).copy()
    if (
        cov.shape != (len(data), len(data))
        or not np.all(np.isfinite(cov))
        or not np.allclose(cov, cov.T)
    ):
        raise ValueError("Finite symmetric covariance required")
    cov += calibration_fraction**2 * np.outer(model, model)
    factor = cho_factor(cov, lower=True)
    residual = data - model
    return float(
        -0.5
        * (
            residual @ cho_solve(factor, residual)
            + 2 * np.log(np.diag(factor[0])).sum()
            + len(data) * np.log(2 * np.pi)
        )
    )


def upper_limit_loglike(limit, model, sigma):
    """Independent censored Gaussian P(measurement < limit), not a zero detection."""
    limit, model, sigma = np.broadcast_arrays(limit, model, sigma)
    if np.any(sigma <= 0) or not np.all(np.isfinite(limit + model + sigma)):
        raise ValueError("Finite values and positive sigma required")
    return float(log_ndtr((limit - model) / sigma).sum())


def discriminants(prediction_a, prediction_b, measurement_sigma, model_sigma=0):
    """Expected per-observable separation; assumes independent Gaussian errors."""
    a, b, err, merr = np.broadcast_arrays(
        prediction_a, prediction_b, measurement_sigma, model_sigma
    )
    if np.any(err <= 0) or np.any(merr < 0) or not np.all(np.isfinite(a + b + err + merr)):
        raise ValueError("Finite predictions and positive measurement uncertainties required")
    return np.abs(a - b) / np.sqrt(err**2 + merr**2)

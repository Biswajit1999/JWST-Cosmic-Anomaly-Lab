"""Instrument-convolved empirical line shapes, not a dense-envelope solver."""

import numpy as np
from scipy.optimize import least_squares
from scipy.special import ndtr
from scipy.stats import exponnorm

from .spectroscopy import pixel_edges

C_KMS = 299792.458


def profile_bins(wave, center, sigma_kms, tail_kms=0):
    """Pixel-averaged Gaussian or symmetric normal-Laplace density in velocity.

    Return 1000 * density(km/s)^-1. Gaussian sigma includes intrinsic+LSF
    quadrature; tail is a symmetric exponential scale. This conserves integrated
    profile area under pixel integration, including nonuniform wavelength bins.
    """
    if center <= 0 or sigma_kms <= 0 or tail_kms < 0:
        raise ValueError("Invalid line-profile parameters")
    velocity_edges = C_KMS * (pixel_edges(wave) / center - 1)
    if tail_kms == 0:
        cdf = ndtr(velocity_edges / sigma_kms)
    else:
        shape = tail_kms / sigma_kms
        cdf = 0.5 * (
            exponnorm.cdf(velocity_edges, shape, scale=sigma_kms)
            + exponnorm.sf(-velocity_edges, shape, scale=sigma_kms)
        )
    return 1000 * np.diff(cdf) / np.diff(velocity_edges)


def fit_complex(
    wave,
    flux,
    error,
    good,
    *,
    main_line,
    components,
    resolution,
    family="gaussian",
    broad=True,
    starts=(400, 1000, 1800),
):
    """Local linear continuum + tied narrow lines + optional main broad line.

    components is a list of [(rest wavelength, relative amplitude), ...] per
    nuisance amplitude, allowing a tied doublet. All wavelengths are in Angstrom
    in the same rest frame. Return conditional fits and residuals, no class.
    """
    wave, flux, error, good = map(np.asarray, (wave, flux, error, good))
    if family not in ["gaussian", "exponential"]:
        raise ValueError("Unknown profile family")
    if wave.ndim != 1 or any(x.shape != wave.shape for x in [flux, error, good]):
        raise ValueError("Matching one-dimensional spectra required")
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("Positive resolving power required")
    mask = good & np.isfinite(flux) & np.isfinite(error) & (error > 0)
    if mask.sum() < 25:
        raise ValueError("Insufficient line-window coverage")
    scale = max(float(np.median(np.abs(flux[mask]))), float(np.median(error[mask])))
    y, err = flux / scale, error / scale
    x = (wave - main_line) / (wave.max() - wave.min())
    sigma_lsf = C_KMS / (2.35482 * resolution)
    count = len(components)

    def predict(p):
        model = p[0] + p[1] * x
        for amp, group in zip(p[4 : 4 + count], components, strict=True):
            for center, relative in group:
                model = model + amp * relative * profile_bins(
                    wave, center * (1 + p[3] / C_KMS), np.hypot(p[2], sigma_lsf)
                )
        if broad:
            amp, sigma, shift = p[4 + count : 7 + count]
            tail = p[-1] if family == "exponential" else 0
            model = model + amp * profile_bins(
                wave, main_line * (1 + shift / C_KMS), np.hypot(sigma, sigma_lsf), tail
            )
        return model

    lower = [-50, -100, 30, -600] + [0] * count
    upper = [100, 100, 500, 600] + [10000] * count
    if broad:
        lower += [0, 100, -900]
        upper += [10000, 3000, 900]
        if family == "exponential":
            lower += [20]
            upper += [4000]
    solutions = []
    for width in starts if broad else (100,):
        initial = [float(np.median(y[mask])), 0, 120, 0] + [0.1] * count
        initial[4] = max(0.1, float(np.max(y[mask])) * 0.1)
        if broad:
            initial += [max(0.1, float(np.max(y[mask]))), width, 0]
            if family == "exponential":
                initial += [800]
        solution = least_squares(
            lambda p: ((predict(p) - y) / err)[mask],
            initial,
            bounds=(lower, upper),
            x_scale="jac",
            max_nfev=1800,
            ftol=1e-9,
            xtol=1e-9,
            gtol=1e-9,
        )
        if solution.success and np.all(np.isfinite(solution.fun)):
            solutions.append(solution)
    if not solutions:
        raise ValueError("No converged multi-start fit")
    best = min(solutions, key=lambda s: np.sum(s.fun**2))
    chi2 = float(np.sum(best.fun**2))
    n, k = int(mask.sum()), len(best.x)
    near_bound = ((best.x - np.array(lower)) < 1e-7 * (np.array(upper) - lower)) | (
        (np.array(upper) - best.x) < 1e-7 * (np.array(upper) - lower)
    )
    return dict(
        family=family if broad else "narrow_only",
        chi2=chi2,
        n=n,
        k=k,
        bic=float(chi2 + k * np.log(n)),
        reduced_chi2=chi2 / (n - k),
        parameter_values=best.x.tolist(),
        parameter_indices_near_bound=np.flatnonzero(near_bound).tolist(),
        broad_sigma_kms=float(best.x[5 + count]) if broad else None,
        exponential_scale_kms=float(best.x[-1]) if broad and family == "exponential" else None,
        resolution=float(resolution),
        model=predict(best.x) * scale,
        residual=np.where(mask, (flux - predict(best.x) * scale) / error, np.nan),
        wavelength=wave,
        mask=mask,
        observed_flux=flux,
        observed_error=error,
    )

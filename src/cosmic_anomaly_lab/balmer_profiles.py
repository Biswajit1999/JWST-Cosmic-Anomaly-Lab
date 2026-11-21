"""Empirical pixel-integrated Balmer emission profiles, not physical retrievals."""

import numpy as np
from scipy.optimize import least_squares
from scipy.special import ndtr

from .spectroscopy import pixel_edges

C_KMS = 299792.458
REST = np.array([4341.65, 4862.64])  # Vacuum Angstrom, Hgamma and Hbeta.
OIII = 4364.44  # Vacuum auroral-line wavelength; empirical nuisance amplitude.


def bin_gaussian(edges, center, sigma):
    """Unit integral Gaussian averaged over the actual, unmasked input pixels."""
    return np.diff(ndtr((edges-center)/sigma))/np.diff(edges)


def pair_model(waves, redshift, parameters, *, continuum_order=1, kinematics="free"):
    """Evaluate observed profiles with no LSF deconvolution or intrinsic width claim.

    Parameters: Hg/Hb integrated fluxes, Hg/Hb continuum intercepts, shared or
    independent velocity centers and observed sigmas, [OIII]4363 nuisance flux,
    then optional Hg/Hb continuum slopes (per 1000 km/s).
    """
    p = np.asarray(parameters)
    free = kinematics == "free"
    shifts = p[4:6] if free else np.repeat(p[4], 2)
    widths = p[6:8] if free else np.repeat(p[5], 2)
    oi = 8 if free else 6
    slope = p[oi+1:oi+3] if continuum_order else np.zeros(2)
    models = []
    for i, w in enumerate(waves):
        center0 = REST[i]*(1+redshift)
        center = center0*(1+shifts[i]/C_KMS)
        sigma = center*widths[i]/C_KMS
        edges = pixel_edges(w)
        velocity = (w/center0-1)*C_KMS
        y = p[2+i] + slope[i]*velocity/1000 + p[i]*bin_gaussian(edges, center, sigma)
        if i == 0:
            ocenter = OIII*(1+redshift)*(1+shifts[i]/C_KMS)
            y = y + p[oi]*bin_gaussian(edges, ocenter, ocenter*widths[i]/C_KMS)
        models.append(y)
    return models


def fit_pair(wavelength, flux, error, good, redshift, *, window_kms=3000,
             continuum_order=1, kinematics="free"):
    """Fit both lines including independent continua and a signed 4363 nuisance.

    All amplitudes remain signed. Noise is diagonal and errors are the unscaled
    local-Jacobian approximation, not profile intervals or detection probabilities.
    Pixel integration precedes masking, so gaps never become artificial large bins.
    """
    if continuum_order not in (0, 1) or kinematics not in ("free", "shared"):
        raise ValueError("Unsupported continuum or kinematics")
    if redshift <= -1 or not np.isfinite(redshift):
        raise ValueError("Invalid redshift")
    w, y, e = [np.asarray(a, float) for a in (wavelength, flux, error)]
    g = np.asarray(good, bool)
    if w.ndim != 1 or not (w.shape == y.shape == e.shape == g.shape) or not np.isfinite(w).all():
        raise ValueError("Expected matching one-dimensional spectra with finite wavelengths")
    g = g & np.isfinite(y+e) & (e > 0)
    edges = pixel_edges(w)
    waves, ys, es, masks, coverage = [], [], [], [], []
    for rest in REST:
        v = (w/(rest*(1+redshift))-1)*C_KMS
        region = abs(v) <= window_kms
        valid = region & g
        blue = valid & (v < -800)
        red = valid & (v > 800)
        center = rest * (1 + redshift)
        lo, hi = center * (1 + np.array([-350, 350]) / C_KMS)
        overlap = np.maximum(0, np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo))
        core_fraction = float(overlap[g].sum() / (hi - lo))
        coverage.append(dict(valid_pixels=int(valid.sum()), expected_pixels=int(region.sum()),
                             blue_pixels=int(blue.sum()), red_pixels=int(red.sum()),
                             core_fraction=core_fraction))
        if (region.sum() < 12 or valid.sum() < .8*region.sum() or min(blue.sum(), red.sum()) < 3
                or core_fraction < .75):
            return dict(status="INSUFFICIENT_COVERAGE", coverage=coverage)
        waves.append(w[region])
        ys.append(y[region])
        es.append(e[region])
        masks.append(g[region])
    free = kinematics == "free"
    nkin = 2 if free else 1
    initial = [100, 200, *[float(np.median(a[b])) for a, b in zip(ys, masks, strict=True)]]
    initial += [0]*nkin + [150]*nkin + [0] + ([0, 0] if continuum_order else [])
    lower = [-np.inf]*4 + [-300]*nkin + [30]*nkin + [-np.inf]*(1+2*continuum_order)
    upper = [np.inf]*4 + [300]*nkin + [600]*nkin + [np.inf]*(1+2*continuum_order)

    def residual(p):
        model = pair_model(waves, redshift, p, continuum_order=continuum_order, kinematics=kinematics)
        return np.concatenate([((yy-mm)/ee)[mask] for yy, mm, ee, mask in
                               zip(ys, model, es, masks, strict=True)])

    fits = []
    for width in [80, 150, 300]:
        guess = np.array(initial, float)
        guess[4+nkin:4+2*nkin] = width
        fits.append(least_squares(residual, guess, bounds=(lower, upper), x_scale="jac", max_nfev=1500))
    fit = min(fits, key=lambda x: np.dot(x.fun, x.fun))
    info = fit.jac.T @ fit.jac
    rank = np.linalg.matrix_rank(fit.jac)
    dof = len(fit.fun)-len(fit.x)
    result = dict(status="FIT", converged=bool(fit.success), parameters=fit.x.tolist(),
                  parameter_errors=None, covariance=None, chi2=float(fit.fun@fit.fun),
                  dof=dof, chi2_dof=float(fit.fun@fit.fun/dof), max_abs_residual=float(abs(fit.fun).max()),
                  jacobian_rank=int(rank), n_parameters=len(fit.x), coverage=coverage,
                  at_boundary=bool(np.any(fit.active_mask != 0)), ratio=None, ratio_error=None,
                  hg_flux=float(fit.x[0]), hb_flux=float(fit.x[1]), hg_error=None, hb_error=None,
                  eligible=False, fit_quality_eligible=False,
                  observed_sigma_hg_kms=float(fit.x[4+nkin]),
                  observed_sigma_hb_kms=float(fit.x[4+nkin+(1 if free else 0)]),
                  velocity_hg_kms=float(fit.x[4]), velocity_hb_kms=float(fit.x[4+(1 if free else 0)]))
    if not fit.success or rank != len(fit.x):
        result["status"] = "NONCONVERGED_OR_SINGULAR"
        return result
    covariance = np.linalg.inv(info)
    errors = np.sqrt(np.diag(covariance))
    result.update(covariance=covariance.tolist(), parameter_errors=errors.tolist(),
                  hg_error=float(errors[0]), hb_error=float(errors[1]))
    if fit.x[0] >= 5*errors[0] and fit.x[1] >= 5*errors[1]:
        gradient = np.array([1/fit.x[1], -fit.x[0]/fit.x[1]**2])
        result.update(eligible=True, ratio=float(fit.x[0]/fit.x[1]),
                      ratio_error=float(np.sqrt(gradient @ covariance[:2, :2] @ gradient)),
                      fit_quality_eligible=not result["at_boundary"] and result["chi2_dof"] <= 2)
    return result

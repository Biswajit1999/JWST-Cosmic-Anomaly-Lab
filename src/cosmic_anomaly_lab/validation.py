"""Seeded simulation experiments with known truth, separate from observations."""

import numpy as np
from astropy import units as u
from scipy.stats import beta

from .anomaly_detection import ErrorAwareDistance, calibrated_pvalues
from .physics import blackbody_lnu, eddington_luminosity, observed_fnu
from .spectroscopy import gaussian_template, linear_fit


def binomial_interval(success, count, confidence=0.95):
    if count <= 0 or not 0 <= success <= count:
        raise ValueError("Invalid binomial counts")
    alpha = 1 - confidence
    lo = 0.0 if success == 0 else beta.ppf(alpha / 2, success, count - success + 1)
    hi = 1.0 if success == count else beta.ppf(1 - alpha / 2, success + 1, count - success)
    return [float(lo), float(hi)]


def simulate(seed=20261110):
    rng = np.random.default_rng(seed)
    wave = np.linspace(3.5, 3.7, 121)
    profile = gaussian_template(wave, 3.6, 0.009)
    design = np.column_stack([np.ones(len(wave)), wave - 3.6, profile])
    truth = np.array([1.5, 0.2, 0.15])
    errors = np.full(len(wave), 0.15)
    fits = []
    for _ in range(500):
        coef, cov, chi2, residual = linear_fit(
            design, design @ truth + rng.normal(0, errors), errors
        )
        fits.append([coef[2], np.sqrt(cov[2, 2])])
    fits = np.asarray(fits)
    pulls = (fits[:, 0] - truth[2]) / fits[:, 1]
    covered = int(np.sum(np.abs(pulls) <= 1.959964))
    # Synthetic inliers share a distribution across independent data splits.
    intrinsic = np.array([[0.09, 0.035], [0.035, 0.16]])

    def sample(n, shift=0):
        sigma = rng.uniform(0.03, 0.15, (n, 2))
        cov = np.array([np.diag(s * s) for s in sigma])
        x = rng.multivariate_normal([0.5, 0.3], intrinsic, n) + rng.normal(size=(n, 2)) * sigma
        x[:, 0] += shift
        return x, cov

    train, tc = sample(1500)
    cal, cc = sample(1000)
    controls, cov = sample(1000)
    injected, ic = sample(500, 1.5)
    model = ErrorAwareDistance().fit(train, tc)
    cs = model.score(cal, cc)
    control_p = calibrated_pvalues(cs, model.score(controls, cov))
    injection_p = calibrated_pvalues(cs, model.score(injected, ic))
    fp, recovered = int(np.sum(control_p <= 0.01)), int(np.sum(injection_p <= 0.01))
    wavelengths = np.geomspace(0.1, 5, 300) * u.um
    blackbody = observed_fnu(blackbody_lnu(wavelengths, 5000 * u.K, 1000 * u.R_sun), 4)
    results = dict(
        data_kind="SIMULATED",
        seed=seed,
        line_recovery=dict(
            trials=500,
            truth_integrated_flux_arbitrary_units=0.15,
            mean_pull=float(pulls.mean()),
            std_pull=float(pulls.std(ddof=1)),
            coverage_95=covered / 500,
            coverage_interval=binomial_interval(covered, 500),
        ),
        anomaly=dict(
            training=1500,
            calibration=1000,
            controls=1000,
            injections=500,
            threshold_p=0.01,
            false_positives=fp,
            false_positive_rate=fp / 1000,
            false_positive_interval=binomial_interval(fp, 1000),
            recovered=recovered,
            recovery_rate=recovered / 500,
            recovery_interval=binomial_interval(recovered, 500),
            limitation="Gaussian 2D population and artificial 1.5 dex shift; not astrophysical discovery performance",
        ),
        physics=dict(
            eddington_luminosity_1e7_Msun_erg_s=float(eddington_luminosity(1e7 * u.M_sun).value),
            blackbody_temperature_K=5000,
            blackbody_radius_Rsun=1000,
            redshift=4,
            flux_min_uJy=float(blackbody.min().value),
            flux_max_uJy=float(blackbody.max().value),
        ),
    )
    return results, pulls, control_p, injection_p, wavelengths, blackbody

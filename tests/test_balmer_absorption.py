import numpy as np
import pytest
from scipy.special import ndtr

from cosmic_anomaly_lab.balmer_absorption import BalmerModel, fit_balmer


def test_zero_absorption_matches_analytic_blurred_gaussian_on_uneven_pixels():
    v = np.linspace(-3000, 3000, 101) + 2 * np.sin(np.arange(101))
    m = BalmerModel(v, 2700, "absorbed", step=2)
    p = [0, 0, 0, 50, 0, 0, 2, 450, 0, 1200, 30, 0, 100, 0]
    sigma = np.hypot(450, m.sigma_lsf)
    expected = 2000 * np.diff(ndtr((m.edges - 30) / sigma)) / np.diff(m.edges)
    assert m.predict(p) == pytest.approx(expected, abs=2e-5)


def test_absorption_precedes_convolution_and_retains_narrow_peak():
    v = np.linspace(-1000, 1000, 81)
    m = BalmerModel(v, 2700, "absorbed")
    p = [0.1, 0, 1.5, 40, 0, 0, 4, 500, 4, 1500, 0, 4, 120, 0]
    absorbed = m.predict(p)
    p[11] = 0
    emission = m.predict(p)
    assert np.all(absorbed <= emission + 1e-12)
    assert absorbed[40] > absorbed[34]  # Central emission survives broad-only absorption.
    p[11] = 4
    assert np.max(absorbed - BalmerModel(v, 1000, "absorbed").predict(p)) > 1


def test_noiseless_absorbed_profile_recovery():
    v = np.linspace(-4000, 4000, 161)
    m = BalmerModel(v, 2700, "absorbed")
    p = [0.2, 0, 1.3, 70, 20, 0.05, 5, 550, 4, 1400, -20, 2, 110, 0]
    flux = m.predict(p)
    fitted = fit_balmer(
        v, flux, np.full(len(v), 0.05), np.ones(len(v), bool), 2700, "absorbed", starts=2
    )
    assert fitted["chi2"] < 0.01
    assert fitted["parameters"]["abs_sigma"] == pytest.approx(110, abs=1)
    assert fitted["parameters"]["tau"] == pytest.approx(2, abs=0.02)

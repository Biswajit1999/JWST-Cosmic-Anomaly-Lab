import numpy as np
import pytest
from astropy import units as u

from cosmic_anomaly_lab.line_anchor import (
    doublet_model,
    fit_oiii,
    profile_crossing,
    resolution_lower_bound,
)
from cosmic_anomaly_lab.spectroscopy import Spectrum


def test_width_recovery_and_missing_pixels():
    rest = np.linspace(0.493, 0.504, 400)
    p = [0.2, 0.01, 0.7, 2, 85, 30]
    flux = doublet_model(rest, p)
    s = Spectrum(
        rest * 6 * u.um,
        flux * u.uJy,
        np.full(len(rest), 0.02) * u.uJy,
        np.ones(len(rest), bool),
        {},
    )
    fit = fit_oiii(s, 5, profile=False)
    assert fit["observed_sigma_kms"] == pytest.approx(85, abs=0.01)
    assert fit["chi2"] < 1e-4
    s.good[:] = False
    assert fit_oiii(s, 5)["status"] == "INSUFFICIENT_COVERAGE"


def test_intrinsic_instrument_degeneracy_and_resolution_bound():
    observed = 100
    for intrinsic in [0, 60, 80]:
        instrument = np.sqrt(observed**2 - intrinsic**2)
        assert resolution_lower_bound(instrument) >= resolution_lower_bound(observed)
    assert resolution_lower_bound(None) is None
    with pytest.raises(ValueError):
        resolution_lower_bound(-1)


def test_profile_crossing_and_unbounded_endpoint():
    x = np.arange(0, 11, dtype=float)
    delta = (x - 5) ** 2
    assert profile_crossing(x, delta, 5, 1, "upper") == 6
    assert profile_crossing(x, delta, 5, 1, "lower") == 4
    assert profile_crossing(x, delta, 5, 100, "upper") is None

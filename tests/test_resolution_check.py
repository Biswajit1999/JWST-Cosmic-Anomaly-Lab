import numpy as np
import pytest
from astropy import units as u

from cosmic_anomaly_lab.resolution_check import halpha_core_contrast
from cosmic_anomaly_lab.spectroscopy import Spectrum


def test_signed_core_contrast_distinguishes_peak_and_dip():
    velocity = np.linspace(-1000, 1000, 501)
    wave = 0.656461 * (1 + velocity / 299792.458) * 6
    spec = Spectrum(
        wave * u.um, np.ones(501) * u.uJy, np.full(501, 0.1) * u.uJy, np.ones(501, dtype=bool), {}
    )
    flat = halpha_core_contrast(spec, 5)
    assert flat["contrast_uJy"] == pytest.approx(0, abs=1e-12)
    spec.flux += np.exp(-0.5 * (velocity / 100) ** 2) * u.uJy
    assert halpha_core_contrast(spec, 5)["contrast_uJy"] < 0
    spec.flux -= 2 * np.exp(-0.5 * (velocity / 100) ** 2) * u.uJy
    dip = halpha_core_contrast(spec, 5)
    assert dip["contrast_uJy"] > 0
    assert dip["statistical_error_uJy"] > 0
    spec.good[:] = False
    assert halpha_core_contrast(spec, 5)["status"] == "INSUFFICIENT_COVERAGE"

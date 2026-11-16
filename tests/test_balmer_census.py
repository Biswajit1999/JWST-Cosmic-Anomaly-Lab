import numpy as np
import pytest

from cosmic_anomaly_lab.balmer_census import (
    aperture_minimum,
    contrast_score,
    normalize_filter,
    signed_core_integral,
)


def test_missing_filter_is_not_repaired_into_instrument_association():
    assert normalize_filter(b"g235m_f170lp") == "g235m_f170lp"
    assert normalize_filter("g395m_f290lp_ext") == "g395m_f290lp_ext"
    assert normalize_filter(b"g235m_f170lp\x00junk") is None
    assert normalize_filter(b"\xffg235m_f170lp") is None
    assert normalize_filter("") is None


def test_signed_contrast_covariance_and_systematic_sensitivity():
    n, d, en, ed, ref = -2, 3, 1, 2, 2
    assert contrast_score(n, d, en, ed, ref) == pytest.approx(-8 / np.sqrt(17))
    assert abs(contrast_score(n, d, en, ed, ref, correlation=-0.5)) < abs(
        contrast_score(n, d, en, ed, ref))
    assert abs(contrast_score(n, d, en, ed, ref, fractional_error=0.2)) < abs(
        contrast_score(n, d, en, ed, ref))
    assert np.isnan(contrast_score(n, d, 0, ed, ref))
    assert np.isnan(contrast_score(np.nan, d, en, ed, ref))
    with pytest.raises(ValueError):
        contrast_score(n, d, en, ed, ref, correlation=1.1)


def test_aperture_agreement_requires_both_and_keeps_opposite_signs():
    result = aperture_minimum([[4, 4, 4], [3, np.nan, -2]])
    np.testing.assert_allclose(result, [3, np.nan, -2], equal_nan=True)


def test_core_integration_keeps_absorption_and_propagates_continuum_error():
    v = np.arange(-4000, 4001, 50.)
    core = abs(v) < 300
    f = 5 + v / 1000 - 2 * core
    e = np.ones_like(v) * .1
    r = signed_core_integral(v, f, e, np.ones_like(v, bool), 10000)
    expected = -2 * core.sum() * 50 * 10000 / 299792.458
    assert r["flux"] == pytest.approx(expected)
    assert r["error"] > .1 * np.sqrt(core.sum()) * 50 * 10000 / 299792.458
    missing = signed_core_integral(v, f, e, abs(v) > 500, 10000)
    assert missing["flux"] is None and missing["status"] == "INSUFFICIENT_COVERAGE"

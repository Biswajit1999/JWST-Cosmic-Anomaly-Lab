import numpy as np
import pytest

from cosmic_anomaly_lab.balmer_census import signed_core_integral
from cosmic_anomaly_lab.continuum_audit import bilateral_core


def inputs():
    v = np.arange(-4000, 4001, 50.)
    return v, 5 + v / 1000 - 2 * (abs(v) < 300), np.full(v.shape, .1), np.ones(v.shape, bool)


def test_signed_flux_and_formal_error_match_legacy_with_full_coverage():
    v, f, e, g = inputs()
    r, old = bilateral_core(v, f, e, g, 10000), signed_core_integral(v, f, e, g, 10000)
    assert r['flux'] < 0
    assert r['flux'] == pytest.approx(old['flux'])
    assert r['error'] == pytest.approx(old['error'])


def test_one_sided_continuum_is_not_accepted_as_bracketing():
    v, f, e, g = inputs()
    g[v < -1400] = False
    assert signed_core_integral(v, f, e, g, 10000)['status'] == 'MEASURED_CORE'
    result = bilateral_core(v, f, e, g, 10000)
    assert result['status'] == 'INSUFFICIENT_BILATERAL_CONTINUUM'
    assert result['blue_pixels'] == 0 and result['flux'] is None


def test_invalid_errors_cannot_supply_missing_continuum():
    v, f, e, g = inputs()
    e[v > 1500] = 0
    assert bilateral_core(v, f, e, g, 10000)['red_pixels'] == 0


def test_invalid_windows_and_grid_are_rejected():
    v, f, e, g = inputs()
    with pytest.raises(ValueError):
        bilateral_core(v, f, e, g, 10000, sidebands=((-100, 100), (800, 1300)))
    with pytest.raises(ValueError):
        bilateral_core(v[::-1], f, e, g, 10000)

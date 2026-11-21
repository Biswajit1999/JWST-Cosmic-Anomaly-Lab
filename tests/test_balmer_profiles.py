import numpy as np
import pytest

from cosmic_anomaly_lab.balmer_profiles import fit_pair, pair_model


def test_observed_profile_recovery_with_masked_pixels():
    redshift = 6.7
    wave = np.arange(3.2e4, 3.9e4, 10.)
    p = [100, 200, .1, .2, 25, -15, 120, 100, 8, .02, -.01]
    # Build distant-line windows separately on the same regularly sampled grid.
    f1, f2 = pair_model([wave, wave], redshift, p)
    flux = np.where(wave < 36000, f1, f2)
    error = np.full(len(wave), .01)
    good = np.ones(len(wave), bool)
    good[np.argmin(abs(wave-4341.65*(1+redshift)))] = False
    result = fit_pair(wave, flux, error, good, redshift)
    assert result['ratio'] == pytest.approx(.5, abs=1e-5)
    assert result['observed_sigma_hg_kms'] == pytest.approx(120, abs=.01)
    assert result['velocity_hb_kms'] == pytest.approx(-15, abs=.01)
    assert result['chi2'] < 1e-6


def test_one_sided_coverage_is_rejected():
    z = 6.7
    wave = np.arange(3.2e4, 3.9e4, 10.)
    good = wave > 4341.65*(1+z)
    result = fit_pair(wave, np.ones(len(wave)), np.ones(len(wave)), good, z)
    assert result['status'] == 'INSUFFICIENT_COVERAGE'

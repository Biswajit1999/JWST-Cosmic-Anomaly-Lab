import numpy as np
import pytest
from astropy import units as u

from cosmic_anomaly_lab.open_search import (
    band_measurements,
    growth_efold_myr,
    normalized_shape,
    peer_mask,
    shape_rarity,
)
from cosmic_anomaly_lab.spectroscopy import Spectrum


def test_band_flux_covariance_and_mask():
    spec = Spectrum(
        np.arange(1, 6) * u.um,
        np.full(5, -2.0) * u.uJy,
        np.ones(5) * u.uJy,
        np.ones(5, dtype=bool),
        {},
    )
    measured = band_measurements(spec, 0, [0.5, 2, 3.5], minimum_pixels=2)
    assert measured["status"] == "MEASURED"
    assert measured["flux_uJy"] == pytest.approx([-2, -2])
    cov = np.asarray(measured["covariance_uJy2"])
    assert cov[0, 0] == pytest.approx(5 / 9)
    assert cov[0, 1] == pytest.approx(1 / 9)
    spec.good[1] = False
    assert band_measurements(spec, 0, [0.5, 2, 3.5])["status"] == "INSUFFICIENT_COVERAGE"


def test_shape_is_amplitude_invariant_and_signed():
    flux = np.array([-1.0, 2, 3])
    assert normalized_shape(flux) == pytest.approx(normalized_shape(flux * 123))
    assert normalized_shape(flux)[0] < 0
    with pytest.raises(ValueError):
        normalized_shape(np.zeros(4))


def test_peers_exclude_self_other_surveys_and_mismatched_depth():
    mask = peer_mask(
        np.array(["A", "A", "B", "A", "A"]),
        np.array([3.0, 3.1, 3.1, 6.0, 3.1]),
        np.array([10.0, 11.0, 11.0, 11.0, 40.0]),
        0,
    )
    assert mask.tolist() == [False, True, False, False, False]


def test_rare_shape_scores_above_known_peer_shape_and_scales_consistently():
    peers = normalized_shape(np.random.default_rng(1).normal(1, 0.05, (80, 12)))
    covariance = np.eye(12) * 0.01**2
    ordinary = shape_rarity(np.ones(12), covariance, peers, np.random.default_rng(2))
    rare_flux = np.r_[np.ones(6), np.full(6, 3.0)]
    rare = shape_rarity(rare_flux, covariance, peers, np.random.default_rng(2))
    scaled = shape_rarity(rare_flux * 5, covariance * 25, peers, np.random.default_rng(2))
    assert rare["perturbation_q10"] > ordinary["perturbation_q90"]
    assert rare["score"] == pytest.approx(scaled["score"])
    assert rare["perturbation_q10"] == pytest.approx(scaled["perturbation_q10"])


def test_accretion_clock_dimensions_and_assumptions():
    assert 49 < growth_efold_myr() < 51
    assert growth_efold_myr(duty_cycle=0.5) == pytest.approx(2 * growth_efold_myr())
    assert growth_efold_myr(eddington_ratio=2) == pytest.approx(growth_efold_myr() / 2)
    with pytest.raises(ValueError):
        growth_efold_myr(efficiency=1)

import numpy as np
import pytest

from cosmic_anomaly_lab.transit_followup import combine_event_depths, event_depth


def test_partial_exposure_and_signed_linear_baseline_fit():
    time = np.arange(-30, 31) * .015 + .003
    duration, cadence = .12, .015
    fraction = np.clip((duration / 2 + cadence / 2 - abs(time)) / cadence, 0, 1)
    flux = 1 + .006 * time + .0025 * fraction
    result = event_depth(time, flux, np.full_like(time, .0001), 0, duration, cadence)
    assert result["depth"] == pytest.approx(-.0025)
    assert result["error"] > 0


def test_missing_bracket_and_empty_combination():
    time = np.arange(-20, 5) * .02
    result = event_depth(time, np.ones_like(time), np.ones_like(time), 0, .12, .02)
    assert result["status"] == "INSUFFICIENT_SAMPLES"
    assert combine_event_depths([result])["depth"] is None

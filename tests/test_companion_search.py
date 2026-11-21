import numpy as np
import pytest

from cosmic_anomaly_lab.companion_search import (
    duration_resolved_period_grid,
    known_transit_mask,
    prepare_companion_quarter,
)
from cosmic_anomaly_lab.transits import phase_distance


def test_known_deep_transits_are_excluded_before_detrending():
    time = np.arange(0, 14, .02)
    known = dict(koi_period=.9, koi_time0bk=.2, koi_duration=4.8)
    in_known = abs(phase_distance(time, .9, .2)) < .09
    outside_known = ~known_transit_mask(time, [known], .02)
    hidden = (abs(time - 7.85) < .06) & outside_known
    flux = 1 - .2*in_known - .001*hidden
    t, f, _, removed, valid = prepare_companion_quarter(
        time, flux, np.full_like(time, 1e-4), np.zeros(len(time)), [known], .02
    )
    assert valid == len(time) and removed == len(time) - outside_known.sum()
    assert not known_transit_mask(t, [known], .02).any()
    assert f.min() == pytest.approx(.999)
    assert np.median(f) == pytest.approx(1)


def test_period_grid_bounds_end_to_end_transit_drift():
    baseline, shortest = 680, .08
    periods = duration_resolved_period_grid(2, 50, baseline, shortest)
    assert periods[[0, -1]].tolist() == [2, 50]
    assert np.max(baseline*np.diff(periods)/periods[:-1]) <= shortest/2
    assert len(periods) < 60000

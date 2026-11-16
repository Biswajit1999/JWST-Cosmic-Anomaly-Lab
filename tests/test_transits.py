import numpy as np
import pytest

from cosmic_anomaly_lab.transits import flatten_segments, local_depth, phase_distance


def test_transit_depth_survives_gap_aware_detrending_and_bad_rows():
    t = np.r_[np.arange(0, 10, .02), np.arange(12, 22, .02)]
    transit = abs(phase_distance(t, 2.5, 2.)) < .05
    f = (1+.0005*t)*(1-.01*transit)
    q = np.zeros(len(t), int)
    q[210:215] = 1
    e = np.full(len(t), .0001)
    ix, y, err = flatten_segments(t, f, e, q)
    assert not np.isin(np.arange(210, 215), ix).any()
    assert not ((t[ix] > 9) & (t[ix] < 13)).any()
    depth = local_depth(t[ix], y, err, 2.5, 2., .10)
    assert depth['depth'] == pytest.approx(.01, abs=.0001)


def test_local_depth_keeps_brightenings_signed_and_missing_unmeasured():
    t = np.arange(0, 10, .02)
    f = 1+.01*(abs(phase_distance(t, 2., 1.)) < .05)
    e = np.full(len(t), .001)
    r = local_depth(t, f, e, 2., 1., .1)
    assert r['depth'] == pytest.approx(-.01)
    r = local_depth(t[:2], f[:2], e[:2], 2., 1., .1)
    assert r['depth'] is None

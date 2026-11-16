import numpy as np
import pytest

from cosmic_anomaly_lab.exposure_audit import combine_extractions


def test_signed_weighted_mean_missingness_and_omission():
    f = [[-2, 10, 1], [2, np.nan, 2]]
    e = [[1, 0, -1], [2, 1, 0]]
    mean, error, count = combine_extractions(f, e)
    assert mean[0] == pytest.approx(-1.2)
    assert error[0] == pytest.approx(np.sqrt(.8))
    assert count.tolist() == [2, 0, 0]
    assert np.isnan(mean[1:]).all() and np.isnan(error[1:]).all()
    mean, error, count = combine_extractions(f, e, omit=0)
    assert mean[0] == 2 and error[0] == 2 and count[0] == 1


def test_bad_shapes_and_exclusion_rejected():
    with pytest.raises(ValueError):
        combine_extractions([1], [1])
    with pytest.raises(ValueError):
        combine_extractions([[1]], [[1]], omit=1)

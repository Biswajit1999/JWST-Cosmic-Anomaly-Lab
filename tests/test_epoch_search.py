import numpy as np

from cosmic_anomaly_lab.epoch_search import aperture_score, fit_registration_residual, local_peaks


def test_masked_hole_and_signed_point_perturbations():
    rng = np.random.default_rng(29)
    old = rng.normal(0, 1, size=(360, 360)).astype(np.float32)
    new = old + rng.normal(0, .35, size=old.shape).astype(np.float32) + 1
    new[105:108, 120:123] += 10
    new[230:233, 250:253] -= 10
    new[150:153, 160:163] += 50
    valid = np.ones(old.shape, bool)
    valid[145:158, 155:168] = False
    residual, interior, _ = fit_registration_residual(old, new, valid)
    weights = np.ones_like(old)
    signal, score, _ = aperture_score(residual, weights, weights, interior, 3)
    y, x = local_peaks(score, 1, 5)[0]
    assert abs(y-106) < 3 and abs(x-121) < 3 and signal[y, x] > 0
    y, x = local_peaks(score, -1, 5)[0]
    assert abs(y-231) < 3 and abs(x-251) < 3 and signal[y, x] < 0
    assert np.isnan(score[151, 161])

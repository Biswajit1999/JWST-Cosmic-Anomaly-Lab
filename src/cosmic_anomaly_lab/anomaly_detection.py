"""Interpretable error-aware distances with held-out calibration, no class assignment."""

import numpy as np
from sklearn.covariance import MinCovDet


class ErrorAwareDistance:
    """Robust intrinsic covariance plus each observation's measurement covariance.

    Complete measured features only. Missing values are rejected, never scored
    as anomalous. Requires exchangeability for calibration; survey shift can
    invalidate that assumption and is reported rather than explained away.
    """

    def fit(self, features, measurement_covariance):
        x, cov = self._validate(features, measurement_covariance)
        if len(x) < max(20, 5 * x.shape[1]):
            raise ValueError("Insufficient training sample")
        model = MinCovDet(random_state=42).fit(x)
        self.location = model.location_
        intrinsic = model.covariance_ - np.median(cov[model.support_], axis=0)
        eigenvalues, eigenvectors = np.linalg.eigh(intrinsic)
        self.intrinsic = (eigenvectors * np.maximum(eigenvalues, 1e-4)) @ eigenvectors.T
        return self

    @staticmethod
    def _validate(features, measurement_covariance):
        x, cov = np.asarray(features, float), np.asarray(measurement_covariance, float)
        if x.ndim != 2 or cov.shape != (len(x), x.shape[1], x.shape[1]):
            raise ValueError("Feature/covariance shapes do not match")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(cov)):
            raise ValueError("Missing features cannot be assigned an anomaly score")
        if not np.allclose(cov, cov.transpose(0, 2, 1)) or np.any(np.linalg.eigvalsh(cov) < -1e-12):
            raise ValueError("Measurement covariance must be positive semidefinite")
        return x, cov

    def score(self, features, measurement_covariance):
        x, cov = self._validate(features, measurement_covariance)
        delta = x - self.location
        solved = np.linalg.solve(self.intrinsic[None, :, :] + cov, delta[:, :, None])[:, :, 0]
        return np.sqrt(np.maximum(0, np.sum(delta * solved, axis=1)))


def calibrated_pvalues(calibration_scores, test_scores):
    calibration = np.asarray(calibration_scores)
    tests = np.asarray(test_scores)
    if (
        len(calibration) < 20
        or not np.all(np.isfinite(calibration))
        or not np.all(np.isfinite(tests))
    ):
        raise ValueError("At least 20 finite held-out calibration scores required")
    return (1 + np.sum(calibration[:, None] >= tests[None, :], axis=0)) / (len(calibration) + 1)


def line_ratio_features(frame, min_sn=5):
    """log10([OIII]/Hb), log10([OIII]/[OII]); retain shared-numerator covariance."""
    cols = ["line_oiii", "line_hb", "line_oii"]
    flux = frame[cols].to_numpy(float)
    err = frame[[c + "_err" for c in cols]].to_numpy(float)
    good = np.all(np.isfinite(flux + err) & (err > 0) & (flux > 0), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        good &= np.all(flux / err >= min_sn, axis=1)
    f, e = flux[good], err[good]
    feature = np.log10(f[:, [0, 0]] / f[:, [1, 2]])
    variance = (e / f / np.log(10)) ** 2
    covariance = np.zeros((len(f), 2, 2))
    covariance[:, 0, 0] = variance[:, 0] + variance[:, 1]
    covariance[:, 1, 1] = variance[:, 0] + variance[:, 2]
    covariance[:, 0, 1] = covariance[:, 1, 0] = variance[:, 0]
    return good, feature, covariance

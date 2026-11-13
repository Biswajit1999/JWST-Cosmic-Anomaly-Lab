"""Synthetic stress tests: descriptive rarity must not be mistaken for novelty."""

from pathlib import Path

import numpy as np
from scipy.stats import beta

from cosmic_anomaly_lab.open_search import normalized_shape, shape_rarity
from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(20261113)
x = np.linspace(-1, 1, 12)
noise = 0.04
cov = np.eye(12) * noise**2


def ordinary(n):
    slope = rng.normal(0, 0.15, (n, 1))
    curve = rng.normal(0, 0.04, (n, 1))
    return np.exp(slope * x + curve * x**2) + rng.normal(0, noise, (n, 12))


peers = normalized_shape(ordinary(300))


def scores(fluxes):
    return np.array(
        [shape_rarity(f, cov, peers, rng, draws=128)["perturbation_q10"] for f in fluxes]
    )


threshold = float(np.quantile(scores(ordinary(500)), 0.99))
outcomes = []
for name in [
    "heldout_ordinary",
    "continuum_step",
    "smooth_bump",
    "single_band_spike",
    "negative_band_artifact",
]:
    fluxes = ordinary(500)
    if name == "continuum_step":
        fluxes[:, 6:] *= 2
    elif name == "smooth_bump":
        fluxes += np.exp(-0.5 * (x / 0.3) ** 2)
    elif name == "single_band_spike":
        fluxes[:, 5] += 2
    elif name == "negative_band_artifact":
        fluxes[:, 5] -= 2
    values = scores(fluxes)
    k, n = int((values > threshold).sum()), len(values)
    outcomes.append(
        dict(
            scenario=name,
            trials=n,
            selected=k,
            fraction=k / n,
            binomial_95_interval=[
                float(beta.ppf(0.025, k, n - k + 1)) if k else 0.0,
                float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0,
            ],
            median_score=float(np.median(values)),
        )
    )
write_json(
    root / "results/open_search/synthetic_validation.json",
    dict(
        executed_utc=utcnow(),
        seed=20261113,
        peers=300,
        null_calibration_trials=500,
        threshold_q99=float(threshold),
        outcomes=outcomes,
        interpretation="Synthetic shape sensitivity only. Artifacts can score as outliers. "
        "Not a real-source false-positive rate or discovery threshold.",
        assumptions="12 bands, log-quadratic continuum family, independent 0.04 errors, "
        "fixed empirical peers. No instrument/extraction/population realism.",
        analysis_sha256=text_digest(root / "src/cosmic_anomaly_lab/open_search.py"),
        script_sha256=text_digest(Path(__file__)),
    ),
)
print(outcomes)

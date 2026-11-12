"""Bounded synthetic shape-recovery study; not real-source discovery calibration."""

from pathlib import Path

import numpy as np

from cosmic_anomaly_lab.line_profiles import C_KMS, fit_complex, profile_bins
from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.validation import binomial_interval

root = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(20261111)
wave = np.linspace(6400, 6750, 151)
resolution = 1300
sigma_lsf = C_KMS / (2.35482 * resolution)
results = []
for truth_family in ["gaussian", "exponential"]:
    for index in range(50):
        sigma = rng.uniform(400, 1100)
        tail = rng.uniform(500, 1300) if truth_family == "exponential" else 0
        error = np.full(len(wave), rng.uniform(0.02, 0.06))
        truth = (
            1
            + 0.4 * profile_bins(wave, 6564.61, np.hypot(120, sigma_lsf))
            + 3 * profile_bins(wave, 6564.61, np.hypot(sigma, sigma_lsf), tail)
        )
        flux = truth + rng.normal(0, error)
        fits = [
            fit_complex(
                wave,
                flux,
                error,
                np.ones(len(wave), bool),
                main_line=6564.61,
                components=[[(6564.61, 1)]],
                resolution=resolution,
                family=f,
            )
            for f in ["gaussian", "exponential"]
        ]
        results.append(
            dict(
                truth=truth_family,
                trial=index,
                sigma_truth=float(sigma),
                tail_truth=float(tail),
                error_sigma=float(error[0]),
                delta_bic_gaussian_minus_exponential=fits[0]["bic"] - fits[1]["bic"],
            )
        )
    print(truth_family, "50 completed", flush=True)
summary = dict(
    data_kind="SIMULATED",
    seed=20261111,
    executed_utc=utcnow(),
    code_sha256=text_digest(root / "src/cosmic_anomaly_lab/line_profiles.py"),
    threshold_delta_bic=10,
    model="Isolated Gaussian narrow + Gaussian or normal-Laplace broad line; known LSF; independent normal errors",
    limitation="Correctly specified synthetic family and fixed resolution; does not include blends, absorption, extraction errors, wrong continua or physical-class degeneracy",
)
for family in ["gaussian", "exponential"]:
    selected = [r for r in results if r["truth"] == family]
    count = sum(r["delta_bic_gaussian_minus_exponential"] > 10 for r in selected)
    summary[family] = dict(
        trials=len(selected),
        exponential_selected=count,
        selection_fraction=count / len(selected),
        binomial_interval_95=binomial_interval(count, len(selected)),
    )
write_json(root / "results/profile_simulation_validation.json", dict(**summary, trials=results))
print(summary)

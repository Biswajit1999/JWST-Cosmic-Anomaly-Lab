import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from cosmic_anomaly_lab.batch_download import download_batch
from cosmic_anomaly_lab.expanded_search import triage_flags
from cosmic_anomaly_lab.line_profiles import C_KMS, fit_complex, profile_bins
from cosmic_anomaly_lab.spectroscopy import pixel_edges


def job(name):
    return dict(
        filename=name,
        url="https://example.org/" + name,
        survey="test",
        version="1",
        kind="synthetic_test",
    )


def reader(url, cap):
    return b"abc", dict(final_url=url, etag=None, last_modified=None)


def test_batch_budget_cache_and_failed_denominator(tmp_path):
    outcome = download_batch(
        tmp_path, [job("a"), job("b")], budget_bytes=5, workers=2, reader=reader
    )
    assert [o["status"] for o in outcome] == ["DOWNLOADED", "FAILED"]
    records = json.loads((tmp_path / "data/manifests/downloads.json").read_text())
    assert len(records) == 1 and records[0]["bytes"] == 3
    assert not (tmp_path / "data/raw/b").exists()
    outcome = download_batch(tmp_path, [job("a")], budget_bytes=5, reader=reader)
    assert outcome[0]["status"] == "CACHED"
    (tmp_path / "data/raw/a").write_bytes(b"bad")
    with pytest.raises(ValueError, match="changed cache"):
        download_batch(tmp_path, [job("a")], budget_bytes=5, reader=reader)


def test_batch_reject_escape_and_oversized_product(tmp_path):
    with pytest.raises(ValueError, match="inside"):
        download_batch(tmp_path, [job("../../escape")], budget_bytes=10, reader=reader)
    outcome = download_batch(tmp_path, [job("a")], budget_bytes=10, maximum_bytes=2, reader=reader)
    assert outcome[0]["status"] == "FAILED"


@pytest.mark.parametrize("tail", [0, 500, 1500])
def test_integrated_line_profiles_conserve_flux(tail):
    # Nonuniform bins; broad wings well inside the integration domain.
    velocity = np.sinh(np.linspace(-4, 4, 1001)) * 1000
    wave = 6564.61 * (1 + velocity / C_KMS)
    profile = profile_bins(wave, 6564.61, 400, tail)
    dv = np.diff(pixel_edges(wave)) / 6564.61 * C_KMS
    assert np.all(profile >= -1e-12)
    assert np.sum(profile * dv) == pytest.approx(1000, rel=1e-6)


def test_recover_gaussian_width_with_instrument_convolution():
    rng = np.random.default_rng(914)
    wave = np.linspace(6400, 6750, 151)
    resolution = 1300
    sigma_lsf = C_KMS / (2.35482 * resolution)
    truth = (
        1
        + 0.4 * profile_bins(wave, 6564.61, np.hypot(120, sigma_lsf))
        + 3 * profile_bins(wave, 6564.61, np.hypot(800, sigma_lsf))
    )
    error = np.full(len(wave), 0.025)
    fit = fit_complex(
        wave,
        truth + rng.normal(0, error),
        error,
        np.ones(len(wave), bool),
        main_line=6564.61,
        components=[[(6564.61, 1)]],
        resolution=resolution,
    )
    assert fit["broad_sigma_kms"] == pytest.approx(800, rel=0.1)
    assert 0.5 < fit["reduced_chi2"] < 1.5


def test_screen_requires_both_diagnostics_and_retains_missing():
    config = yaml.safe_load(
        (Path(__file__).parents[1] / "configs/expanded_search.yaml").read_text()
    )["triage"]
    assert not triage_flags({}, config)["joint_followup_priority"]
    r = dict(
        continuum_ratio=3,
        continuum_lower_3sigma=2.5,
        alternate_continuum_lower_3sigma=2.4,
        minimum_broad_delta_bic=20,
        minimum_broad_sn=8,
    )
    assert triage_flags(r, config)["joint_followup_priority"]
    r["minimum_broad_delta_bic"] = -2
    assert not triage_flags(r, config)["joint_followup_priority"]

import numpy as np
import pandas as pd
import pytest
import yaml
from astropy import units as u

from cosmic_anomaly_lab.catalogues import select_pilot
from cosmic_anomaly_lab.classification import HypothesisRegistry
from cosmic_anomaly_lab.spectroscopy import Spectrum, fit_halpha_blend, gaussian_template


def test_config_changes_selection_and_empty_is_safe():
    frame = pd.DataFrame(
        {
            "file": ["a", "b"],
            "source_id": ["one", "two"],
            "grating": ["PRISM"] * 2,
            "grade": [3, 3],
            "sn50": [10, 1],
            "z_best": [2.5, 2.6],
            "survey": ["JADES"] * 2,
            "root": ["jades"] * 2,
            "srcid": [1, 2],
        }
    )
    config = {"spectroscopy": {"redshift_bins": [[2, 3]], "per_survey_per_bin": [5]}}
    assert len(select_pilot(frame, config=config)) == 1
    config["spectroscopy"]["min_archive_median_pixel_sn"] = 100
    assert len(select_pilot(frame, config=config)) == 0


def test_halpha_recovery_uses_instrument_lsf():
    z = 3.0
    wave = np.linspace(2.45, 2.8, 240) * u.um
    center = 6564.61 * (1 + z)
    resolution = 100 * 1.3
    width = center * np.hypot(1 / (2.35482 * resolution), 100 / 299792.458)
    profile = gaussian_template(wave.to_value(u.AA), center, width)
    unit = u.erg / u.s / u.cm**2 / u.AA
    injected_flux = 2e-17
    flam = (1e-20 + injected_flux * profile) * unit
    error = np.ones(len(wave)) * 1e-21 * unit
    spec = Spectrum(
        wave,
        flam.to(u.uJy, equivalencies=u.spectral_density(wave)),
        error.to(u.uJy, equivalencies=u.spectral_density(wave)),
        np.ones(len(wave), bool),
        {},
    )
    out = fit_halpha_blend(spec, z, np.array([0.6, 5.3]), np.array([100, 100]))
    narrow = next(x for x in out["fits"] if x["model"] == "narrow_blend")
    assert narrow["flux_erg_s_cm2"] == pytest.approx(injected_flux, rel=1e-7)
    assert narrow["rest_ew_A"] == pytest.approx(injected_flux / 1e-20 / (1 + z), rel=1e-7)
    assert narrow["bic"] == min(x["bic"] for x in out["fits"])


def test_hypotheses_do_not_turn_toy_models_into_physical_classes(tmp_path):
    p = tmp_path / "h.yaml"
    p.write_text(
        yaml.safe_dump(
            {
                "populations": {"AGN": {}, "COMPOSITE": {}},
                "hypotheses": {
                    "obscured_AGN": {
                        "requires": ["xray", "sed"],
                        "implemented_forward_model": False,
                    }
                },
            }
        )
    )
    registry = HypothesisRegistry(p)
    assert len(registry.assess_labels({"AGN": 0.8, "COMPOSITE": 0.7})["labels"]) == 2
    assert not registry.applicability("obscured_AGN", ["xray", "sed"])["testable"]
    with pytest.raises(ValueError):
        registry.assess_labels({"MAGIC": 1})

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from scipy.stats import norm

from cosmic_anomaly_lab.anomaly_detection import (
    ErrorAwareDistance,
    calibrated_pvalues,
    line_ratio_features,
)
from cosmic_anomaly_lab.crossmatching import candidates
from cosmic_anomaly_lab.imaging import aperture_measurement
from cosmic_anomaly_lab.inference import discriminants, gaussian_loglike, upper_limit_loglike
from cosmic_anomaly_lab.physics import blackbody_lnu, eddington_luminosity, observed_fnu
from cosmic_anomaly_lab.spectroscopy import (
    Spectrum,
    continuum_break,
    gaussian_template,
    linear_fit,
    pixel_edges,
    read_spectrum,
)
from cosmic_anomaly_lab.validation import simulate


def test_wraparound_and_ambiguity():
    left = pd.DataFrame(dict(ra_deg=[359.99999, 30], dec_deg=[0, 0]))
    right = pd.DataFrame(dict(ra_deg=[0.00001, 359.99998, 31], dec_deg=[0, 0, 0]))
    result = candidates(left, right, 0.2)
    assert len(result) == 2 and result.ambiguous.all()
    assert result.separation_arcsec.max() < 0.1


def test_polar_coordinates():
    a = pd.DataFrame(dict(ra_deg=[0], dec_deg=[89.99999]))
    b = pd.DataFrame(dict(ra_deg=[180], dec_deg=[89.99999]))
    assert len(candidates(a, b, 0.1)) == 1


def test_gaussian_normalization_and_covariance():
    y = np.array([1.0, 2.0])
    m = np.array([0.5, 2.1])
    s = np.array([0.3, 0.4])
    assert gaussian_loglike(y, m, sigma=s) == pytest.approx(norm.logpdf(y, m, s).sum())
    assert gaussian_loglike(y, m, covariance=np.diag(s * s)) == pytest.approx(
        gaussian_loglike(y, m, sigma=s)
    )
    cov = np.diag(s * s) + 0.1**2 * np.outer(m, m)
    assert gaussian_loglike(y, m, sigma=s, calibration_fraction=0.1) == pytest.approx(
        gaussian_loglike(y, m, covariance=cov)
    )


def test_upper_limits_are_censored():
    assert upper_limit_loglike(3.0, 0.0, 1.0) == pytest.approx(np.log(norm.cdf(3)))
    assert upper_limit_loglike(3.0, 6.0, 1.0) < upper_limit_loglike(3.0, 0.0, 1.0)


@pytest.mark.parametrize("sigma", [-1, 0, np.nan])
def test_invalid_errors(sigma):
    with pytest.raises(ValueError):
        gaussian_loglike([1], [1], sigma=[sigma])


def test_physical_units_and_flux_transform():
    assert eddington_luminosity(1 * u.M_sun).to_value(u.erg / u.s) == pytest.approx(
        1.257e38, rel=0.003
    )
    from astropy.cosmology import Planck18

    lnu = 1e28 * u.erg / u.s / u.Hz
    fnu = observed_fnu(lnu, 4)
    expected = 5 * lnu / (4 * np.pi * Planck18.luminosity_distance(4) ** 2)
    assert fnu.to_value(u.uJy) == pytest.approx(expected.to_value(u.uJy))
    first = blackbody_lnu(1 * u.um, 5000 * u.K, 1 * u.R_sun)
    second = blackbody_lnu(1 * u.um, 5000 * u.K, 2 * u.R_sun)
    assert (second / first).value == pytest.approx(4)


def test_pixel_integrated_line_and_negative_flux_recovery():
    wave = np.linspace(-1, 1, 100)
    template = gaussian_template(wave, 0, 0.02)
    assert np.sum(template * np.diff(pixel_edges(wave))) == pytest.approx(1, abs=1e-10)
    design = np.column_stack([np.ones(len(wave)), template])
    coef, cov, chi2, _ = linear_fit(design, design @ np.array([2.0, -0.7]), np.ones(len(wave)))
    assert coef[1] == pytest.approx(-0.7)
    assert chi2 < 1e-18 and cov[1, 1] > 0


def test_fits_units_mask_and_sort_validation(tmp_path):
    wave = np.linspace(0.7, 5.2, 20)

    def write(wave, unit="mJy"):
        cols = [
            fits.Column(name="wave", format="D", unit="um", array=wave),
            fits.Column(name="flux", format="D", unit=unit, array=np.ones(20) * 2),
            fits.Column(name="err", format="D", unit=unit, array=np.r_[0, np.ones(19) * 0.1]),
        ]
        fits.HDUList(
            [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols, name="SPEC1D")]
        ).writeto(tmp_path / "s.fits", overwrite=True)

    write(wave)
    spec = read_spectrum(tmp_path / "s.fits")
    assert spec.flux[1].to_value(u.uJy) == 2000 and spec.good.sum() == 19
    write(wave[::-1])
    with pytest.raises(ValueError):
        read_spectrum(tmp_path / "s.fits")
    write(wave, unit=None)
    with pytest.raises(ValueError):
        read_spectrum(tmp_path / "s.fits")


def test_break_does_not_divide_by_noise():
    wave = np.linspace(0.33, 0.43, 101) * u.um
    spec = Spectrum(wave, np.zeros(101) * u.uJy, np.ones(101) * u.uJy, np.ones(101, dtype=bool), {})
    assert continuum_break(spec, 0)["ratio"] is None
    spec.flux = np.where(wave.value < 0.38, 10, 20) * u.uJy
    assert continuum_break(spec, 0)["ratio"] == pytest.approx(2)


def test_missing_features_refused_and_uncertainty_downweights():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(100, 2))
    cov = np.tile(np.eye(2) * 0.01, (100, 1, 1))
    model = ErrorAwareDistance().fit(x, cov)
    with pytest.raises(ValueError):
        model.score([[np.nan, 1]], np.eye(2)[None])
    assert (
        model.score([[5, 5]], np.eye(2)[None] * 10)[0]
        < model.score([[5, 5]], np.eye(2)[None] * 0.01)[0]
    )


def test_shared_ratio_covariance_and_limits():
    df = pd.DataFrame(
        {
            "line_oiii": [10, -1],
            "line_hb": [5, 5],
            "line_oii": [2, 2],
            "line_oiii_err": [1, 1],
            "line_hb_err": [0.5, 0.5],
            "line_oii_err": [0.2, 0.2],
        }
    )
    good, x, cov = line_ratio_features(df)
    assert good.tolist() == [True, False] and cov[0, 0, 1] > 0
    assert x[0, 0] == pytest.approx(np.log10(2))
    assert calibrated_pvalues(np.arange(100), [100])[0] == pytest.approx(1 / 101)


def test_aperture_calibration_and_background():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [51, 51]
    wcs.wcs.crval = [10, 20]
    wcs.wcs.cdelt = [-0.05 / 3600, 0.05 / 3600]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    image = np.ones((101, 101)) * 3
    image[50, 50] += 100
    result = aperture_measurement(image, np.ones_like(image) * 4, wcs, 10, 20, jy_per_unit=1e-8)
    assert result["aperture_flux_uJy"] == pytest.approx(1)
    assert result["diagonal_error_uJy"] > 0


def test_discriminant_prediction():
    assert discriminants([1, 2], [4, 6], [1, 2]).tolist() == [3, 2]


def test_injection_recovery_is_calibrated():
    result, *_ = simulate()
    assert 0.90 < result["line_recovery"]["coverage_95"] < 0.99
    assert abs(result["line_recovery"]["mean_pull"]) < 0.15
    assert result["anomaly"]["false_positive_rate"] < 0.03
    assert result["anomaly"]["recovery_rate"] > 0.5

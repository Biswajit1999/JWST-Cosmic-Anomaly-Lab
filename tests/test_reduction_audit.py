import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits

from cosmic_anomaly_lab.reduction_audit import (
    collapse_spatial,
    read_dr4_spectrum,
    signed_spatial_moment,
)


def test_dr4_units_variants_and_missing_values(tmp_path):
    wave = np.linspace(1, 2, 20)
    flux = (np.full(20, -2.0) * u.uJy).to_value(
        u.erg / u.s / u.cm**2 / u.AA, equivalencies=u.spectral_density(wave * u.um)
    )
    error = abs(flux) / 4
    flux[5] = np.nan
    cols = [fits.Column(name="WAVELENGTH", array=wave, format="D", unit="um")]
    for name, values in [
        ("FLUX", flux),
        ("FLUX_ERR", error),
        ("FLUX_NOD2", 2 * flux),
        ("FLUX_ERR_NOD2", 2 * error),
    ]:
        cols.append(fits.Column(name=name, array=values, format="D", unit="erg s-1 cm-2 AA-1"))
    path = tmp_path / "dr4.fits"
    fits.HDUList(
        [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(cols, name="EXTRACT5PIX1D")]
    ).writeto(path)
    spec = read_dr4_spectrum(path)
    assert not spec.good[5]
    assert spec.flux.value[0] == pytest.approx(-2)
    assert read_dr4_spectrum(path, nods=2).flux.value[0] == pytest.approx(-4)
    with pytest.raises(ValueError):
        read_dr4_spectrum(path, aperture=4)


def test_spatial_collapse_propagates_weights_and_missing_coverage():
    wave = np.arange(1, 6, dtype=float)
    science = np.full((3, 5), 2.0)
    weight = np.full((3, 5), 4.0)
    p, v = collapse_spatial(wave, science, weight, 0.5, 2)
    assert p == pytest.approx([2, 2, 2])
    assert v == pytest.approx(np.full(3, 5 / 36))
    weight[0, :] = 0
    p, v = collapse_spatial(wave, science, weight, 0.5, 2)
    assert np.isnan(p[0]) and np.isnan(v[0])


def test_signed_moment_retains_negative_pixels_and_rejects_low_sn():
    y = np.arange(-2, 3, dtype=float)
    flux = np.array([-0.2, 1, 3, 2, -0.1])
    result = signed_spatial_moment(y, flux, np.full(5, 0.001), half_width=2)
    assert result["centroid_pixels"] == pytest.approx(y @ flux / flux.sum())
    assert signed_spatial_moment(y, flux, np.ones(5) * 100)["status"] == "LOW_SIGNED_SUM_SN"

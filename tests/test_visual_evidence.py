import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from cosmic_anomaly_lab.visual_evidence import (
    BANDS,
    check_registration,
    colour_composite,
    read_images,
)


def image_wcs():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [10, -30]
    wcs.wcs.crpix = [61, 61]
    wcs.wcs.cdelt = [-0.05 / 3600, 0.05 / 3600]
    return wcs


def test_rgb_preserves_calibrated_positive_ratios_at_saturation():
    channels = np.array([[[100.0, 20, 10], [1, 0.2, 0.1]]])
    rgb = colour_composite(channels, np.ones((1, 2), dtype=bool), scale=0.1)
    assert rgb.max() <= 1
    np.testing.assert_allclose(rgb[..., 0] / rgb[..., 1], 5)
    np.testing.assert_allclose(rgb[..., 1] / rgb[..., 2], 2)


def test_rgb_negative_and_missing_pixels_are_display_only():
    channels = np.array([[[-1.0, 1, 0], [np.nan, np.nan, np.nan]]])
    original = channels.copy()
    rgb = colour_composite(channels, np.array([[True, False]]), scale=1)
    assert rgb[0, 0, 0] == 0
    np.testing.assert_allclose(rgb[0, 1], 0.22)
    np.testing.assert_equal(channels, original)
    with pytest.raises(ValueError, match="Nonfinite"):
        colour_composite(channels, np.ones((1, 2), dtype=bool), scale=1)


def test_reject_misregistered_channels():
    reference = image_wcs()
    assert check_registration(reference, reference, (120, 120)) < 1e-5
    shifted = image_wcs()
    shifted.wcs.crpix[0] += 0.2
    with pytest.raises(ValueError, match="reprojection"):
        check_registration(reference, shifted, (120, 120))


@pytest.mark.parametrize(
    "bands,blank_band,expect_rgb",
    [(BANDS, None, True), (BANDS[2:], None, False), (BANDS, BANDS[1], False)],
)
def test_fits_calibration_coverage_and_missing_band(tmp_path, bands, blank_band, expect_rgb):
    rng = np.random.default_rng(71)
    hdus = [fits.PrimaryHDU()]
    for band in bands:
        data = 10 + rng.normal(0, 0.1, (120, 120))
        data[57:64, 57:64] += 20
        weight = np.ones_like(data) * 100
        weight[0, 0] = 0
        if band == blank_band:
            weight[:] = 0
        for version, values in [("SCI", data), ("WHT", weight)]:
            header = image_wcs().to_header()
            header["FILTER"] = band
            header["EXTVER"] = version
            header["PHOTFNU"] = 1e-8
            hdus.append(fits.ImageHDU(values, header=header))
    path = tmp_path / "images.fits"
    fits.HDUList(hdus).writeto(path)
    images, rgb, detail = read_images(path, 10, -30)
    assert (rgb is not None) == expect_rgb
    assert set(images) == set(bands) - {blank_band}
    for band, im in images.items():
        assert np.isnan(im["data"][0, 0])
        assert detail["bands"][band]["display_background_uJy_pixel"] == pytest.approx(0.1, abs=1e-4)
        assert im["data"][60, 60] == pytest.approx(0.2, abs=0.005)
    for aperture in detail["apertures"]:
        assert aperture["status"] == (
            "UNAVAILABLE" if aperture["filter"] == blank_band else "MEASURED"
        )

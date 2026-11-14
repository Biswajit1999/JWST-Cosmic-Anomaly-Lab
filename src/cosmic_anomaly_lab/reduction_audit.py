"""Cross-reduction readers and signed spatial diagnostics with explicit limits."""

import numpy as np
from astropy import units as u
from astropy.io import fits

from .spectroscopy import Spectrum, pixel_edges


def read_dr4_spectrum(path, aperture=5, nods=3):
    """Read the documented DR4 boxcar product; no undocumented error rescaling."""
    if aperture not in (3, 5) or nods not in (2, 3):
        raise ValueError("DR4 supports 3/5-pixel, 2/3-nod variants")
    with fits.open(path, checksum=True) as hdus:
        hdus.verify("exception")
        table = hdus[f"EXTRACT{aperture}PIX1D"]
        flux_name, error_name = (
            ("FLUX", "FLUX_ERR") if nods == 3 else ("FLUX_NOD2", "FLUX_ERR_NOD2")
        )
        for name in ["WAVELENGTH", flux_name, error_name]:
            if not table.columns[name].unit:
                raise ValueError(f"Missing DR4 column unit: {name}")
        wave = (np.array(table.data["WAVELENGTH"]) * u.Unit(table.columns["WAVELENGTH"].unit)).to(
            u.um
        )
        if not np.all(np.isfinite(wave.value)) or np.any(np.diff(wave.value) <= 0):
            raise ValueError("Increasing finite wavelength grid required")
        flux = (np.array(table.data[flux_name]) * u.Unit(table.columns[flux_name].unit)).to(
            u.uJy, equivalencies=u.spectral_density(wave)
        )
        error = (np.array(table.data[error_name]) * u.Unit(table.columns[error_name].unit)).to(
            u.uJy, equivalencies=u.spectral_density(wave)
        )
        good = np.isfinite(flux.value + error.value) & (error.value > 0)
        good &= (wave.value >= 0.7) & (wave.value <= 5.2)
        meta = {
            key: hdus[0].header.get(key)
            for key in ["RA", "DEC", "TARGET", "PROPOSID", "MJD-BEG", "MJD-END", "HLSPVER", "DOI"]
        }
        meta.update(
            aperture_pixels=aperture,
            nods=nods,
            quality_mask="Finite flux/error and positive error; no DQ column supplied",
        )
    if good.sum() < 10:
        raise ValueError("Fewer than ten usable DR4 spectral bins")
    return Spectrum(wave, flux, error, good, meta)


def collapse_spatial(wave, science, inverse_variance, low, high):
    """Wavelength-overlap mean per row; WHT is treated as inverse variance.

    This is a remeasurement of the released rectified image, not a detector-level
    independent extraction. Resampling covariance is unavailable and omitted.
    """
    edges = pixel_edges(wave)
    weights = np.maximum(0, np.minimum(edges[1:], high) - np.maximum(edges[:-1], low))
    valid = np.isfinite(science) & np.isfinite(inverse_variance) & (inverse_variance > 0)
    overlap = valid * weights[None, :]
    total = overlap.sum(axis=1)
    usable = total >= 0.8 * (high - low)
    a = np.divide(overlap, total[:, None], out=np.zeros_like(overlap), where=total[:, None] > 0)
    profile = np.sum(a * np.where(valid, science, 0), axis=1)
    variance = np.sum(np.divide(a * a, inverse_variance, out=np.zeros_like(a), where=valid), axis=1)
    return np.where(usable, profile, np.nan), np.where(usable, variance, np.nan)


def signed_spatial_moment(y, profile, variance, half_width=3, minimum_sn=5):
    """Signed centroid inside a stated aperture, conditional on positive total S/N.

    Do not zero negative background-subtracted pixels. An unstable signed moment
    can lie outside the aperture; return unavailable instead of inventing a center.
    """
    mask = (abs(y) <= half_width) & np.isfinite(profile + variance) & (variance > 0)
    if mask.sum() < 3:
        return dict(status="INSUFFICIENT_SPATIAL_COVERAGE")
    yy, ff, vv = y[mask], profile[mask], variance[mask]
    total = float(ff.sum())
    error = float(np.sqrt(vv.sum()))
    result = dict(total_uJy=total, statistical_error_uJy=error, half_width_pixels=half_width)
    if total <= minimum_sn * error:
        return dict(**result, status="LOW_SIGNED_SUM_SN", centroid_pixels=None)
    center = float(yy @ ff / total)
    if abs(center) > half_width:
        return dict(**result, status="UNSTABLE_SIGNED_MOMENT", centroid_pixels=None)
    uncertainty = float(np.sqrt(np.sum((yy - center) ** 2 * vv)) / total)
    return dict(
        **result,
        status="CONDITIONAL_MOMENT",
        centroid_pixels=center,
        centroid_statistical_error_pixels=uncertainty,
    )

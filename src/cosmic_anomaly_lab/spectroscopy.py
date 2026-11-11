"""Conservative spectral diagnostics with explicit units and conditional errors."""

from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from scipy.special import ndtr


@dataclass
class Spectrum:
    wavelength: u.Quantity
    flux: u.Quantity
    error: u.Quantity
    good: np.ndarray
    metadata: dict

    def as_specutils(self):
        from specutils import Spectrum as SpecutilsSpectrum

        return SpecutilsSpectrum(
            spectral_axis=self.wavelength,
            flux=self.flux,
            uncertainty=StdDevUncertainty(self.error),
            mask=~self.good,
        )

    def rest_wavelength(self, redshift):
        if not np.isfinite(redshift) or redshift <= -1:
            raise ValueError("Invalid redshift")
        return self.wavelength / (1 + redshift)


def read_spectrum(path, wavelength_range=(0.7, 5.2)):
    """DJA SPEC1D and standard JWST EXTRACT1D; units must be in FITS columns.

    Read released flux as-is; no extra path-loss or photometric rescaling. Flux
    calibration is a product-level assumption, not independently established here.
    """
    with fits.open(path, checksum=True) as hdus:
        hdus.verify("exception")
        extension = "SPEC1D" if "SPEC1D" in hdus else "EXTRACT1D"
        hdu = hdus[extension]
        names = {name.lower(): name for name in hdu.columns.names}
        wave_name = names.get("wave", names.get("wavelength"))
        flux_name = names.get("flux")
        error_name = names.get("err", names.get("flux_error"))
        if any(name is None for name in [wave_name, flux_name, error_name]):
            raise ValueError("Required spectrum columns missing")
        for name in [wave_name, flux_name, error_name]:
            if not hdu.columns[name].unit:
                raise ValueError(f"Missing FITS unit for {name}")
        wave = (np.array(hdu.data[wave_name], float) * u.Unit(hdu.columns[wave_name].unit)).to(u.um)
        if (
            not np.all(np.isfinite(wave))
            or np.any(wave.value <= 0)
            or np.any(np.diff(wave.value) <= 0)
        ):
            raise ValueError("Wavelengths must be finite, positive and strictly increasing")
        flux = (np.array(hdu.data[flux_name], float) * u.Unit(hdu.columns[flux_name].unit)).to(
            u.uJy, equivalencies=u.spectral_density(wave)
        )
        error = (np.array(hdu.data[error_name], float) * u.Unit(hdu.columns[error_name].unit)).to(
            u.uJy, equivalencies=u.spectral_density(wave)
        )
        good = np.isfinite(flux.value) & np.isfinite(error.value) & (error.value > 0)
        good &= (wave.value >= wavelength_range[0]) & (wave.value <= wavelength_range[1])
        if "dq" in names:
            good &= np.asarray(hdu.data[names["dq"]]) == 0
        if "npix" in names:
            good &= np.asarray(hdu.data[names["npix"]]) > 0
        meta = {
            key: hdu.header.get(key, hdus[0].header.get(key))
            for key in ["GRATING", "FILTER", "VERSION", "CRDS_CTX", "CAL_VER", "SRCRA", "SRCDEC"]
        }
    if good.sum() < 10:
        raise ValueError("Fewer than 10 usable spectral bins")
    return Spectrum(wave, flux, error, good, meta)


def continuum_break(spec, redshift, blue=(0.34, 0.36), red=(0.40, 0.42)):
    """Rest-frame fnu red/blue continuum proxy, not a unique Balmer-break interpretation."""
    wave = spec.rest_wavelength(redshift).to_value(u.um)
    measurements = []
    for lo, hi in [blue, red]:
        in_band = (wave >= lo) & (wave <= hi)
        mask = spec.good & in_band
        if mask.sum() < 3 or mask.sum() < 0.8 * in_band.sum():
            return dict(status="INSUFFICIENT_COVERAGE", ratio=None, ratio_error=None)
        weights = 1 / spec.error.value[mask] ** 2
        mean = np.sum(weights * spec.flux.value[mask]) / weights.sum()
        measurements.append((mean, np.sqrt(1 / weights.sum())))
    (b, be), (r, re) = measurements
    result = dict(
        blue_uJy=float(b),
        blue_error_uJy=float(be),
        red_uJy=float(r),
        red_error_uJy=float(re),
        ratio=None,
        ratio_error=None,
        status="LOW_CONTINUUM_SN",
    )
    if b > 3 * be and r > 3 * re:
        result.update(
            ratio=float(r / b),
            ratio_error=float(np.sqrt((re / b) ** 2 + (r * be / b**2) ** 2)),
            status="MEASURED_CONDITIONAL_ON_REDSHIFT",
        )
    return result


def pixel_edges(wave):
    wave = np.asarray(wave, float)
    if len(wave) < 2 or np.any(np.diff(wave) <= 0):
        raise ValueError("Increasing wavelength grid required")
    mid = (wave[1:] + wave[:-1]) / 2
    return np.r_[wave[0] - (mid[0] - wave[0]), mid, wave[-1] + (wave[-1] - mid[-1])]


def gaussian_template(wave, center, sigma):
    """Unit integrated flux, pixel-averaged Gaussian; wavelength and sigma in same units."""
    if not np.isfinite(center + sigma) or sigma <= 0:
        raise ValueError("Positive finite sigma required")
    edges = pixel_edges(wave)
    return np.diff(ndtr((edges - center) / sigma)) / np.diff(edges)


def linear_fit(design, y, error):
    design, y, error = np.asarray(design), np.asarray(y), np.asarray(error)
    if np.any(error <= 0) or not np.all(np.isfinite(design)) or not np.all(np.isfinite(y + error)):
        raise ValueError("Finite inputs with positive errors required")
    a, b = design / error[:, None], y / error
    coef, _, rank, _ = np.linalg.lstsq(a, b, rcond=None)
    if rank < design.shape[1]:
        raise ValueError("Unidentifiable model")
    covariance = np.linalg.inv(a.T @ a)
    residual = (y - design @ coef) / error
    chi2 = float(residual @ residual)
    return coef, covariance, chi2, residual


def fit_halpha_blend(
    spec, redshift, resolution_wave, resolution_r, *, lsf_scale=1.3, error_scale=1
):
    """Local empirical Halpha+[NII] blend diagnostic; NOT a black-hole estimator.

    Linear continuum plus one Gaussian at vacuum Halpha 6564.61 A. The blend
    cannot separate [NII], broad Halpha, absorption or velocity shifts. Compare
    conditional intrinsic sigma=100 and 1500 km/s at a fixed redshift. Gaussian
    LSF interpolates the instrument reference R(lambda). BIC is descriptive.
    """
    if not np.isfinite(redshift) or redshift <= -1 or lsf_scale <= 0 or error_scale <= 0:
        raise ValueError("Invalid redshift or scale")
    wave = spec.wavelength.to_value(u.AA)
    center = 6564.61 * (1 + redshift)
    rest = wave / (1 + redshift)
    window = (rest >= 6200) & (rest <= 6900)
    mask = spec.good & window
    # Require continuum on both sides and no substantial coverage gaps.
    if (
        mask.sum() < 10
        or np.count_nonzero(mask & (rest < 6400)) < 2
        or np.count_nonzero(mask & (rest > 6750)) < 2
        or mask.sum() < 0.8 * window.sum()
    ):
        return dict(status="INSUFFICIENT_COVERAGE", fits=[])
    if not resolution_wave.min() <= center / 1e4 <= resolution_wave.max():
        return dict(status="NO_LSF_COVERAGE", fits=[])
    resolution = np.interp(center / 1e4, resolution_wave, resolution_r) * lsf_scale
    fluxunit = u.erg / u.s / u.cm**2 / u.AA
    conversion = (np.ones(len(wave)) * u.uJy).to_value(
        fluxunit, equivalencies=u.spectral_density(spec.wavelength)
    )
    y = spec.flux.value * conversion / 1e-20
    err = spec.error.value * conversion / 1e-20 * error_scale
    x = (wave - center) / (350 * (1 + redshift))
    output = []
    for name, velocity in [("continuum", None), ("narrow_blend", 100), ("broad_blend", 1500)]:
        design = np.column_stack([np.ones(len(wave)), x])
        if velocity is not None:
            sigma = center * np.hypot(1 / (2.35482 * resolution), velocity / 299792.458)
            profile = gaussian_template(wave, center, sigma)
            design = np.column_stack([design, profile])
        coef, cov, chi2, residual = linear_fit(design[mask], y[mask], err[mask])
        k, n = len(coef), int(mask.sum())
        result = dict(
            model=name,
            chi2=chi2,
            dof=n - k,
            bic=float(chi2 + k * np.log(n)),
            n_pixels=n,
            intrinsic_sigma_km_s=velocity,
            effective_resolving_power=float(resolution),
            continuum_erg_s_cm2_A=float(coef[0] * 1e-20),
            continuum_error_erg_s_cm2_A=float(np.sqrt(cov[0, 0]) * 1e-20),
            flux_erg_s_cm2=None,
            flux_error_erg_s_cm2=None,
            rest_ew_A=None,
            rest_ew_error_A=None,
            sn=None,
            upper_limit_3sigma_erg_s_cm2=None,
            normalized_residual_rms=float(np.sqrt(np.mean(residual**2))),
        )
        if velocity is not None:
            amp, amp_err = coef[2] * 1e-20, np.sqrt(cov[2, 2]) * 1e-20
            result.update(
                flux_erg_s_cm2=float(amp),
                flux_error_erg_s_cm2=float(amp_err),
                sn=float(amp / amp_err),
                upper_limit_3sigma_erg_s_cm2=float(max(0, amp) + 3 * amp_err),
            )
            if coef[0] > 3 * np.sqrt(cov[0, 0]):
                grad = np.array([-coef[2] / coef[0] ** 2, 0, 1 / coef[0]]) / (1 + redshift)
                result.update(
                    rest_ew_A=float(coef[2] / coef[0] / (1 + redshift)),
                    rest_ew_error_A=float(np.sqrt(grad @ cov @ grad)),
                )
        output.append(result)
    return dict(
        status="EMPIRICAL_BLEND_DIAGNOSTIC",
        fits=output,
        caveat="Conditional diagonal errors; fixed z, LSF, continuum and profile. Halpha and [NII] unresolved.",
    )

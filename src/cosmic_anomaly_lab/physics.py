"""Unit-aware analytic baselines, not stellar atmospheres or radiative transfer."""

import numpy as np
from astropy import constants as c
from astropy import units as u
from astropy.cosmology import Planck18


def eddington_luminosity(mass):
    """Fully ionized pure hydrogen, Thomson scattering, spherical force balance."""
    mass = u.Quantity(mass).to(u.kg)
    if np.any(mass <= 0 * u.kg):
        raise ValueError("Mass must be positive")
    return (4 * np.pi * c.G * mass * c.m_p * c.c / c.sigma_T).to(u.erg / u.s)


def accretion_luminosity(mdot, efficiency=0.1):
    if not 0 < efficiency <= 1:
        raise ValueError("Efficiency must be in (0,1]")
    mdot = u.Quantity(mdot).to(u.kg / u.s)
    if np.any(mdot < 0 * u.kg / u.s):
        raise ValueError("Negative accretion rate")
    return (efficiency * mdot * c.c**2).to(u.erg / u.s)


def blackbody_lnu(wavelength, temperature, radius):
    """Spherical opaque isotropic blackbody; no lines, dust, envelope, or IGM."""
    wave = u.Quantity(wavelength).to(u.m)
    temp, rad = u.Quantity(temperature).to(u.K), u.Quantity(radius).to(u.m)
    if np.any(wave <= 0 * u.m) or np.any(temp <= 0 * u.K) or np.any(rad <= 0 * u.m):
        raise ValueError("Positive wavelength, temperature and radius required")
    nu = (c.c / wave).to(u.Hz)
    bnu = 2 * c.h * nu**3 / c.c**2 / np.expm1((c.h * nu / (c.k_B * temp)).decompose().value)
    return (4 * np.pi**2 * rad**2 * bnu).to(u.erg / u.s / u.Hz)


def observed_fnu(rest_lnu, redshift, magnification=1):
    """Fnu at nu_rest/(1+z) = (1+z)*Lnu/(4 pi DL^2), Hogg 1999 Eq. 22.

    No IGM/dust correction. Magnification is an input, never inferred here.
    """
    if redshift <= 0 or magnification <= 0:
        raise ValueError("Positive redshift and magnification required")
    dl = Planck18.luminosity_distance(redshift)
    return (magnification * (1 + redshift) * rest_lnu / (4 * np.pi * dl**2)).to(u.uJy)

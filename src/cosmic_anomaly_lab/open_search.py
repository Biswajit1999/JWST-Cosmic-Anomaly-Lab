"""Class-agnostic broad spectral shapes and matched-peer descriptive rarity."""

import numpy as np
from astropy import units as u
from astropy.constants import G, c, m_p, sigma_T

from .spectroscopy import pixel_edges


def band_measurements(spec, redshift, edges, minimum_coverage=0.8, minimum_pixels=2):
    """Integrate signed fnu with pixel-overlap weights; never fill missing bands.

    Covariance includes pixels split across adjacent bands, conditional on the
    released diagonal pixel errors. The wavelength average is not an AB filter.
    """
    rest = spec.rest_wavelength(redshift).to_value(u.um)
    pixels = pixel_edges(rest)
    edges = np.asarray(edges, dtype=float)
    if np.any(np.diff(edges) <= 0):
        raise ValueError("Band edges must increase")
    overlap = np.maximum(
        0,
        np.minimum(edges[1:, None], pixels[None, 1:])
        - np.maximum(edges[:-1, None], pixels[None, :-1]),
    )
    overlap[:, ~spec.good] = 0
    coverage = overlap.sum(axis=1) / np.diff(edges)
    count = (overlap > 0).sum(axis=1)
    result = dict(coverage=coverage.tolist(), pixel_count=count.tolist())
    if np.any(coverage < minimum_coverage) or np.any(count < minimum_pixels):
        return dict(**result, status="INSUFFICIENT_COVERAGE")
    weights = overlap / overlap.sum(axis=1)[:, None]
    flux = weights @ np.where(spec.good, spec.flux.to_value(u.uJy), 0)
    error = np.where(spec.good, spec.error.to_value(u.uJy), 0)
    covariance = (weights * error**2) @ weights.T
    return dict(
        **result,
        status="MEASURED",
        flux_uJy=flux.tolist(),
        covariance_uJy2=covariance.tolist(),
    )


def normalized_shape(flux):
    """Normalize amplitude without removing signed continuum or spectral slope."""
    flux = np.asarray(flux, dtype=float)
    norm = np.sqrt(np.mean(flux**2, axis=-1, keepdims=True))
    if np.any(~np.isfinite(norm)) or np.any(norm <= 0):
        raise ValueError("Finite nonzero RMS flux required")
    return flux / norm


def peer_mask(survey, redshift, sn, index, maximum_z_gap=0.15, maximum_sn_ratio=2):
    """Strictly exclude the target; no angular matching or class labels."""
    mask = (
        (survey == survey[index])
        & (np.abs(redshift - redshift[index]) / (1 + redshift[index]) < maximum_z_gap)
        & (sn >= sn[index] / maximum_sn_ratio)
        & (sn <= sn[index] * maximum_sn_ratio)
    )
    mask[index] = False
    return mask


def shape_rarity(flux, covariance, peer_shapes, rng, *, floor=0.05, draws=128):
    """RMS robust peer residual, with conditional noise perturbation quantiles.

    A score is NOT chi-square, a p-value, a posterior or physical-model failure.
    The peer MAD contains astrophysical spread AND measurement noise. Shared
    normalization couples bins; perturbations propagate this, but the scalar
    distance deliberately remains descriptive rather than a likelihood.
    """
    flux, covariance = np.asarray(flux), np.asarray(covariance)
    median = np.median(peer_shapes, axis=0)
    spread = np.maximum(1.4826 * np.median(np.abs(peer_shapes - median), axis=0), floor)
    norm = np.sqrt(np.mean(flux**2))
    denominator = np.sqrt(spread**2 + np.diag(covariance) / norm**2)
    shape = normalized_shape(flux)
    residual = (shape - median) / denominator
    simulated = rng.multivariate_normal(flux, covariance, size=draws, check_valid="raise")
    scores = np.sqrt(np.mean(((normalized_shape(simulated) - median) / denominator) ** 2, axis=1))
    return dict(
        score=float(np.sqrt(np.mean(residual**2))),
        perturbation_q10=float(np.quantile(scores, 0.1)),
        perturbation_q90=float(np.quantile(scores, 0.9)),
        dominant_bin=int(np.argmax(np.abs(residual))),
        peer_median=median.tolist(),
        peer_scatter=spread.tolist(),
        normalized_flux=shape.tolist(),
        residual=residual.tolist(),
    )


def growth_efold_myr(efficiency=0.1, eddington_ratio=1, duty_cycle=1):
    """Idealized radiatively efficient accretion clock, not an evolutionary solver."""
    if not 0 < efficiency < 1 or eddington_ratio <= 0 or not 0 < duty_cycle <= 1:
        raise ValueError("Invalid growth assumptions")
    eddington_time = (c * sigma_T / (4 * np.pi * G * m_p)).to_value(u.Myr)
    return float(eddington_time * efficiency / (1 - efficiency) / eddington_ratio / duty_cycle)

"""A signed local core contrast; not an absorption or physical-class classifier."""

import numpy as np

from .open_search import band_measurements


def halpha_core_contrast(spec, redshift):
    """Mean of two shoulders minus central flux, conditional on a fixed redshift.

    Windows in km/s: [-450,-250], [-120,120], [250,450]. Positive contrast
    describes a local dip, not its cause. Instrument LSF, emission components,
    neighbouring lines and correlated extraction errors can alter the result.
    """
    velocity_edges = np.array([-450, -250, -120, 120, 250, 450])
    edges = 0.656461 * (1 + velocity_edges / 299792.458)
    bands = band_measurements(spec, redshift, edges, minimum_pixels=1)
    if bands["status"] != "MEASURED":
        return bands
    if min(np.array(bands["pixel_count"])[[0, 2, 4]]) < 2:
        return dict(status="INSUFFICIENT_SAMPLING", pixel_count=bands["pixel_count"])
    weights = np.array([0.5, 0, -1, 0, 0.5])
    flux = np.asarray(bands["flux_uJy"])
    covariance = np.asarray(bands["covariance_uJy2"])
    contrast = float(weights @ flux)
    error = float(np.sqrt(weights @ covariance @ weights))
    return dict(
        status="MEASURED_CONDITIONAL_CONTRAST",
        redshift=float(redshift),
        contrast_uJy=contrast,
        statistical_error_uJy=error,
        conditional_contrast_over_error=contrast / error,
        velocity_edges_kms=velocity_edges.tolist(),
        bands=bands,
    )

"""Diagnostic combinations of released extractions, preserving signed flux."""

import numpy as np


def combine_extractions(flux, error, omit=None):
    """Inverse-variance mean, NOT a reconstruction of the survey stack.

    Inputs share a wavelength grid. Formal errors assume independent diagonal
    noise; nodded backgrounds can violate this assumption. No clipping is done.
    """
    f, e = np.asarray(flux, float), np.asarray(error, float)
    if f.ndim != 2 or f.shape != e.shape:
        raise ValueError("Matched exposure-by-wavelength matrices required")
    valid = np.isfinite(f) & np.isfinite(e) & (e > 0)
    if omit is not None:
        if not 0 <= omit < len(f):
            raise ValueError("Omitted exposure index out of range")
        valid[omit] = False
    weights = np.zeros_like(e)
    np.divide(1, e**2, out=weights, where=valid)
    total = weights.sum(axis=0)
    mean = np.divide((weights * np.where(valid, f, 0)).sum(axis=0), total,
                     out=np.full(f.shape[1], np.nan), where=total > 0)
    err = np.sqrt(np.divide(1, total, out=np.full_like(total, np.nan), where=total > 0))
    return mean, err, valid.sum(axis=0)

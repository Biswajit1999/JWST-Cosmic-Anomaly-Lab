"""Pixel-domain diagnostics for two calibrated, registered epoch mosaics."""

import numpy as np
from scipy import ndimage as ndi


def disk(radius):
    extent = int(np.ceil(radius))
    y, x = np.mgrid[-extent:extent+1, -extent:extent+1]
    return (x*x+y*y <= radius*radius).astype(np.float32)


def fit_registration_residual(old, new, valid):
    """Robust small-shift and small symmetric PSF-change model on reference pixels."""
    gy, gx = np.gradient(old)
    laplace = ndi.laplace(old)
    models = [old, gx, gy, laplace]
    interior = ndi.binary_erosion(valid, iterations=3)
    sample = np.s_[::5, ::5]
    good = interior[sample]
    values = new[sample][good].astype(float)
    design = np.column_stack([np.ones(len(values)), *[a[sample][good] for a in models]])
    if len(values) < 2000:
        raise ValueError("Too little common sky to model differences")
    keep = np.ones(len(values), bool)
    for _ in range(6):
        beta, _, rank, _ = np.linalg.lstsq(design[keep], values[keep], rcond=None)
        if rank != design.shape[1]:
            raise ValueError("Insufficient variation to fit the epoch model")
        beta[1:] = np.clip(beta[1:], [.8, -.45, -.45, -.65], [1.2, .45, .45, .65])
        residual = values - design @ beta
        center = np.median(residual)
        beta[0] += center
        scatter = 1.4826*np.median(abs(residual-center))
        if scatter <= 0:
            break
        keep = abs(residual-center) <= 4*scatter
    predicted = np.full_like(old, beta[0])
    for coefficient, term in zip(beta[1:], models, strict=True):
        predicted += np.float32(coefficient)*term
    residual = new-predicted
    metadata = dict(background=float(beta[0]), flux_scale=float(beta[1]),
                    derivative_x=float(beta[2]), derivative_y=float(beta[3]),
                    laplacian=float(beta[4]), fit_pixels=int(keep.sum()),
                    trial_pixels=len(keep), displaced_or_psf_residuals_possible=True)
    return residual, interior, metadata


def aperture_score(diff, w_old, w_new, valid, radius, block=180):
    """Locally normalize circular aperture sums using empirical sky scatter."""
    aperture = disk(radius)
    good = ndi.convolve(valid.astype(np.float32), aperture, mode="constant") == aperture.sum()
    signal = ndi.convolve(np.where(valid, diff, 0), aperture, mode="constant")
    v = np.zeros_like(diff, dtype=np.float32)
    v[valid] = 1/w_old[valid]+1/w_new[valid]
    formal = np.sqrt(ndi.convolve(v, aperture, mode="constant"))
    ratio = np.zeros_like(diff, dtype=np.float32)
    np.divide(signal, formal, where=good & (formal > 0), out=ratio)
    z = np.full_like(diff, np.nan, dtype=np.float32)
    noise = np.full_like(diff, np.nan, dtype=np.float32)
    for y in range(0, diff.shape[0], block):
        for x in range(0, diff.shape[1], block):
            box = np.s_[y:y+block, x:x+block]
            sub = ratio[box][good[box]]
            if len(sub) < block**2*.12:
                continue
            median = np.median(sub)
            scale = 1.4826*np.median(abs(sub-median))
            if np.isfinite(scale) and scale > 0:
                z[box] = np.where(good[box], (ratio[box]-median)/scale, np.nan)
                noise[box] = np.where(good[box], formal[box]*scale, np.nan)
    return signal, z, noise


def local_peaks(score, sign, minimum_score, separation_pixels=6):
    """Find separated trial locations without reversing the sign of stored fluxes."""
    q = np.where(np.isfinite(score), score*sign, -np.inf)
    maxima = ndi.maximum_filter(q, size=2*separation_pixels+1, mode="constant", cval=-np.inf)
    y, x = np.nonzero((q >= minimum_score) & (q == maxima))
    return sorted(zip(y, x, strict=True), key=lambda pair: -q[pair])

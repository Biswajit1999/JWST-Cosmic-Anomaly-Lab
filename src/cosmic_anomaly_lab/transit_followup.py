"""Measure signed depths at a locked ephemeris with local baselines."""

import numpy as np


def event_depth(time, flux, error, center, duration, cadence):
    """Fit a cadence-averaged box and linear baseline, retaining signed depth.

    Exposure is approximated by a contiguous box with width TIMEDEL. Errors
    are conditional on the timing and shape and omit between-cadence covariance.
    """
    x = np.asarray(time, float) - center
    flux, error = np.asarray(flux, float), np.asarray(error, float)
    box = np.maximum(
        0, np.minimum(x + cadence / 2, duration / 2) - np.maximum(x - cadence / 2, -duration / 2)
    ) / cadence
    in_event = box > 0
    left = (x < -duration) & (x >= -3 * duration)
    right = (x > duration) & (x <= 3 * duration)
    result = dict(midpoint=float(center), in_points=int(in_event.sum()),
                  left_points=int(left.sum()), right_points=int(right.sum()),
                  depth=None, error=None, chi2_dof=None)
    if in_event.sum() < 3 or min(left.sum(), right.sum()) < 3:
        return dict(result, status="INSUFFICIENT_SAMPLES")
    use = in_event | left | right
    design = np.column_stack([np.ones(use.sum()), x[use] / duration, -box[use]])
    weighted = design / error[use, None]
    rhs = flux[use] / error[use]
    beta, _, rank, _ = np.linalg.lstsq(weighted, rhs, rcond=None)
    if rank < len(beta):
        return dict(result, status="SINGULAR")
    residual = rhs - weighted @ beta
    chi2_dof = float(residual @ residual / (len(rhs) - len(beta)))
    covariance = np.linalg.inv(weighted.T @ weighted)
    return dict(result, status="MEASURED", depth=float(beta[2]),
                error=float(np.sqrt(covariance[2, 2] * max(1, chi2_dof))), chi2_dof=chi2_dof)


def measure_events(time, flux, error, period, epoch, duration, cadence):
    """Record all predicted events within the observation span, including gaps."""
    start = int(np.ceil((time[0] - epoch) / period))
    stop = int(np.floor((time[-1] - epoch) / period))
    for cycle in range(start, stop + 1):
        center = epoch + cycle * period
        lo, hi = np.searchsorted(time, [center - 3 * duration, center + 3 * duration])
        yield dict(cycle=cycle, **event_depth(time[lo:hi], flux[lo:hi], error[lo:hi],
                                             center, duration, cadence))


def combine_event_depths(events):
    """Inverse-variance mean with an event scatter floor, under independence."""
    measured = [r for r in events if r["status"] == "MEASURED"]
    if not measured:
        return dict(events=0, positive=0, depth=None, error=None, contrast=None, event_scatter=None)
    depths = np.array([r["depth"] for r in measured])
    errors = np.array([r["error"] for r in measured])
    weights = 1 / errors**2
    mean = float(np.average(depths, weights=weights))
    scatter = float(np.sum(weights * (depths - mean)**2) / (len(depths) - 1)) if len(depths) > 1 else None
    error = float(np.sqrt(max(1, scatter or 0) / weights.sum()))
    return dict(events=len(depths), positive=int((depths > 0).sum()), depth=mean, error=error,
                contrast=mean/error, event_scatter=scatter)

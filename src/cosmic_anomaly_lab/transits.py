"""Exploratory transit screening: quality masks and diagnostics, not validation."""

import numpy as np
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter


def flatten_segments(time, flux, error, quality, window_days=2.):
    """Divide by a running median separately across >0.5-day gaps.

    Retain downward outliers (possible transits); reject flagged/nonfinite rows.
    Trim one day at each segment edge, avoiding incomplete median windows.
    Diagonal errors omit detrending covariance. Returns original row indices.
    """
    t, f, e = (np.asarray(a, float) for a in (time, flux, error))
    good = np.isfinite(t+f+e) & (f > 0) & (e > 0) & (np.asarray(quality) == 0)
    ix = np.where(good)[0]
    ix = ix[np.argsort(t[ix])]
    if len(ix) < 100:
        return np.array([], int), np.array([]), np.array([])
    if np.any(np.diff(t[ix]) <= 0):
        raise ValueError("Duplicate or non-increasing retained times")
    selected, flattened, errors = [], [], []
    for group in np.split(ix, np.where(np.diff(t[ix]) > .5)[0]+1):
        if len(group) < 100 or np.ptp(t[group]) <= 2*window_days:
            continue
        dt = np.median(np.diff(t[group]))
        size = max(3, int(round(window_days/dt)) | 1)
        trend = median_filter(f[group], size=size, mode="reflect")
        interior = (t[group] >= t[group[0]]+window_days/2) & (t[group] <= t[group[-1]]-window_days/2)
        selected.extend(group[interior])
        flattened.extend(f[group][interior]/trend[interior])
        errors.extend(e[group][interior]/trend[interior])
    return np.array(selected, int), np.array(flattened), np.array(errors)


def phase_distance(time, period, epoch):
    return (np.asarray(time)-epoch+period/2) % period-period/2


def local_depth(time, flux, error, period, epoch, duration):
    """Local out-minus-in weighted depth, signed and with formal error."""
    phase = abs(phase_distance(time, period, epoch))
    inside = phase < duration/2
    outside = (phase > duration) & (phase < min(3*duration, period/2))
    if inside.sum() < 3 or outside.sum() < 6:
        return dict(depth=None, error=None, in_points=int(inside.sum()), status="INSUFFICIENT_SAMPLES")
    def mean(mask):
        w = 1/np.asarray(error)[mask]**2
        return float(np.sum(np.asarray(flux)[mask]*w)/w.sum()), float(1/w.sum())
    inn, vi = mean(inside)
    off, vo = mean(outside)
    return dict(depth=off-inn, error=float(np.sqrt(vi+vo)), in_points=int(inside.sum()), status="MEASURED")


def scan_transits(time, flux, error):
    """Fixed 0.8–20 day BLS scan and inverted-flux comparison; no FAP claim."""
    t, y, e = map(lambda x: np.asarray(x, float), (time, flux, error))
    scatter = 1.4826*np.median(abs(y-np.median(y)))
    effective_error = np.maximum(e, scatter)
    model = BoxLeastSquares(t, y, effective_error)
    result = model.autopower([.06, .10, .16, .24], minimum_period=.8, maximum_period=20,
                            minimum_n_transit=3, frequency_factor=5, objective="snr")
    index = int(np.argmax(result.power))
    period, epoch, duration = (float(result[k][index]) for k in ["period", "transit_time", "duration"])
    inverted = BoxLeastSquares(t, 2-y, effective_error).power(result.period, [.06, .10, .16, .24], objective="snr")
    stats = model.compute_stats(period, duration, epoch)
    odd, even = stats["depth_odd"], stats["depth_even"]
    secondary = stats["depth_phased"]
    odd_even = float(abs(odd[0]-even[0])/np.hypot(odd[1], even[1]))
    events = []
    indices = np.rint((t-epoch)/period).astype(int)
    for n in np.unique(indices):
        near = abs(t-(epoch+n*period)) < min(3*duration, period/2)
        d = local_depth(t[near], y[near], effective_error[near], period, epoch, duration)
        events.append(dict(event_number=int(n), midpoint_bkjd=epoch+n*period, **d))
    answer = dict(period_days=period, epoch_bkjd=epoch, duration_days=duration,
        depth=float(result.depth[index]), depth_error=float(result.depth_err[index]),
        depth_snr=float(result.depth_snr[index]), odd_even_contrast=odd_even,
        secondary_depth=float(secondary[0]), secondary_error=float(secondary[1]),
        secondary_snr=float(secondary[0]/secondary[1]),
        inverted_max_snr=float(np.max(inverted.power)),
        empirical_point_scatter=float(scatter), periods_searched=len(result.period),
        observed_events=sum(e["status"] == "MEASURED" for e in events),
        positive_events=sum(e["depth"] is not None and e["depth"] > 0 for e in events))
    return answer, events, result, effective_error

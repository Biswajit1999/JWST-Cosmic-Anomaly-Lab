"""Known-transit masking and local detrending for companion searches."""

import numpy as np
import pandas as pd

from .transits import phase_distance


def duration_resolved_period_grid(minimum, maximum, baseline, shortest_duration):
    """Bound maximum end-to-end drift to half the shortest trial duration."""
    if not 0 < minimum < maximum or baseline <= 0 or shortest_duration <= 0:
        raise ValueError("Positive baseline, duration and ordered period bounds required")
    max_log_step = np.log1p(shortest_duration / (2 * baseline))
    points = int(np.ceil(np.log(maximum / minimum) / max_log_step)) + 1
    return np.geomspace(minimum, maximum, points)


def known_transit_mask(time, planets, cadence, multiplier=1.0):
    """Flag exposure centers within one full catalogue duration of transit center."""
    time = np.asarray(time, float)
    masked = np.zeros(len(time), dtype=bool)
    for planet in planets:
        half_width = multiplier * float(planet["koi_duration"]) / 24 + cadence / 2
        masked |= abs(phase_distance(time, planet["koi_period"], planet["koi_time0bk"])) <= half_width
    return masked


def prepare_companion_quarter(time, flux, error, quality, planets, cadence, *, window_days=3.0,
                               mask_multiplier=1.0):
    """Estimate a running median from unmasked data in contiguous time segments."""
    time, flux, error = (np.asarray(v, float) for v in [time, flux, error])
    valid = np.isfinite(time+flux+error) & (flux > 0) & (error > 0) & (np.asarray(quality) == 0)
    t, f, e = time[valid], flux[valid], error[valid]
    if np.any(np.diff(t) <= 0):
        raise ValueError("Expected an increasing time grid")
    masked = known_transit_mask(t, planets, cadence, multiplier=mask_multiplier)
    selected = []
    for idx in np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > .5) + 1):
        if len(idx) < 100 or np.ptp(t[idx]) < window_days:
            continue
        values = pd.Series(np.where(masked[idx], np.nan, f[idx]),
                           index=pd.to_timedelta(t[idx]-t[idx[0]], unit="D"))
        trend = values.rolling(pd.Timedelta(days=window_days), center=True, min_periods=30).median().to_numpy()
        keep = (~masked[idx] & np.isfinite(trend) & (trend > 0)
                & (t[idx] >= t[idx[0]] + window_days/2) & (t[idx] <= t[idx[-1]] - window_days/2))
        selected.append((t[idx][keep], f[idx][keep]/trend[keep], e[idx][keep]/trend[keep]))
    if not selected:
        return np.array([]), np.array([]), np.array([]), int(masked.sum()), int(valid.sum())
    t, f, e = (np.concatenate(v) for v in zip(*selected, strict=True))
    return t, f, e, int(masked.sum()), int(valid.sum())

"""Observed host-star pixels and transit difference image for the known control."""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import colormaps
from matplotlib.colors import AsinhNorm
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transits import phase_distance
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/exoplanet_pilot"
p = json.loads((out / "pixel_product.json").read_text())
assert digest(root / p["path"]) == p["sha256"]
r = pd.read_csv(out / "screening.csv")
r = r[r.kic == p["kic"]].iloc[0]
with fits.open(root / p["path"]) as h:
    assert h[0].header["KEPLERID"] == p["kic"] and h[0].header["QUARTER"] == p["quarter"]
    d = h[1].data
    t = np.array(d["TIME"])
    cube = np.array(d["FLUX"], float)
    good = np.isfinite(t) & (d["QUALITY"] == 0)
    aperture = (h[2].data.astype(int) & 2) != 0
    wcs = WCS(h[2].header).celestial
    x, y = wcs.world_to_pixel_values(h[0].header["RA_OBJ"], h[0].header["DEC_OBJ"])
    flux_unit = h[1].columns["FLUX"].unit
    assert flux_unit == "e-/s"
phase = abs(phase_distance(t, r.period_days, r.epoch_bkjd))
inside = good & (phase < r.duration_days/2)
outside = good & (phase > r.duration_days) & (phase < 3*r.duration_days)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    mean = np.nanmedian(cube[good], axis=0)
    difference = np.nanmean(cube[outside], axis=0)-np.nanmean(cube[inside], axis=0)
fig = Figure(figsize=(10, 4.8), layout="constrained")
axes = fig.subplots(1, 2)
high = float(np.nanpercentile(mean, 99))
im = axes[0].imshow(mean, origin="lower", cmap=colormaps["magma"].with_extremes(bad="#909090"), norm=AsinhNorm(linear_width=high/20, vmin=0, vmax=high))
fig.colorbar(im, ax=axes[0], label="Median observed flux (electrons/s)", shrink=.7)
lim = float(np.nanmax(abs(difference)))
im = axes[1].imshow(difference, origin="lower", cmap=colormaps["RdBu_r"].with_extremes(bad="#909090"), vmin=-lim, vmax=lim)
fig.colorbar(im, ax=axes[1], label="Out − in transit (electrons/s)", shrink=.7)
for ax in axes:
    ax.contour(aperture.astype(float), levels=[.5], colors="cyan", linewidths=.8)
    ax.plot(x, y, "+", color="#74ee6f", ms=12)
    ax.set(xlabel="Cutout column (≈4 arcsec/pixel)", ylabel="Cutout row")
axes[0].set_title("Actual Kepler Q2 pixels: host star")
axes[1].set_title("Observed transit difference image")
fig.suptitle("TrES-2 / KIC 11446443 · known-planet control\nCross: catalogue coordinate; cyan: pipeline aperture; planet is not spatially resolved", fontsize=11)
save_figure(fig, root / "figures/exoplanet_pilot/tres2_observed_pixels")
write_json(out / "pixel_summary.json", dict(executed_utc=utcnow(), input_sha256=p["sha256"],
    script_sha256=text_digest(Path(__file__)), in_transit_cadences=int(inside.sum()),
    comparison_cadences=int(outside.sum()), unit=flux_unit, pixel_scale_arcsec=3.988966,
    method="Median quality-zero image; mean local-phase out minus in image",
    centroid_significance_measured=False, psf_fit_performed=False,
    limitation="Difference image can retain pointing, background, and crowding effects; not validation"))

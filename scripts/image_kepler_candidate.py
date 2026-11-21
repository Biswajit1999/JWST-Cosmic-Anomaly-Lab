"""Inspect actual Kepler pixel data around the strongest later-quarter events."""

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

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_followup"
products = json.loads((out / "products.json").read_text())
duration = json.loads((out / "summary.json").read_text())["duration_days"]
events = pd.read_csv(out / "events.csv").query("segment >= 4 and status == 'MEASURED' and reduction == 'PDCSAP_FLUX'")
quarters = ((events.depth / events.error)**2).groupby(events.segment).sum().nlargest(2).index.tolist()
cache = DownloadCache(root, 4 * 1024**3)
fig = Figure(figsize=(11, 9), layout="constrained")
axes = fig.subplots(2, 2)
records = []
for row, quarter in enumerate(quarters):
    matches = []
    for p in products["products"]:
        with fits.open(root / p["path"]) as h:
            if h[0].header["QUARTER"] == quarter:
                matches.append(p)
    assert len(matches) == 1
    name = Path(matches[0]["path"]).name.replace("_llc.fits", "_lpd-targ.fits.gz")
    assert name in products["ancillary"]["listed_target_pixel_files"]
    url = products["ancillary"]["pixel_directory"]["url"] + name
    path = cache.fetch(url, "exoplanets/" + name, survey="Kepler", version="public archive",
                       kind="observed_target_pixel_cube")
    with fits.open(path) as h:
        assert h[0].header["KEPLERID"] == 11446424 and h[0].header["QUARTER"] == quarter
        assert h[1].header["BJDREFI"] == 2454833 and h[1].header["TIMESYS"] == "TDB"
        assert h[1].columns["FLUX"].unit == "e-/s"
        data = h[1].data
        time, cube = np.asarray(data["TIME"]), np.asarray(data["FLUX"], float)
        good = np.isfinite(time) & (data["QUALITY"] == 0)
        aperture = (h[2].data.astype(int) & 2) != 0
        wcs = WCS(h[2].header).celestial
        x, y = wcs.world_to_pixel_values(h[0].header["RA_OBJ"], h[0].header["DEC_OBJ"])
    differences, used = [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        median = np.nanmedian(cube[good], axis=0)
        for center in events.loc[events.segment == quarter, "midpoint"]:
            phase = time - center
            inside = good & (abs(phase) <= duration/2)
            left = good & (phase >= -3*duration) & (phase < -duration)
            right = good & (phase <= 3*duration) & (phase > duration)
            if min(inside.sum(), left.sum(), right.sum()) < 3:
                continue
            differences.append((np.nanmean(cube[left], axis=0) + np.nanmean(cube[right], axis=0))/2
                               - np.nanmean(cube[inside], axis=0))
            used.append(dict(midpoint_bkjd=center, in_points=int(inside.sum()),
                             left_points=int(left.sum()), right_points=int(right.sum())))
        assert differences
        difference = np.nanmean(differences, axis=0)
    high = float(np.nanpercentile(median, 99))
    limit = float(np.nanmax(abs(difference)))
    im = axes[row, 0].imshow(median, origin="lower", cmap=colormaps["magma"].with_extremes(bad="#999999"),
                            norm=AsinhNorm(linear_width=high/20, vmin=0, vmax=high))
    fig.colorbar(im, ax=axes[row, 0], label="Median flux (electrons/s)", shrink=.8)
    im = axes[row, 1].imshow(difference, origin="lower",
                            cmap=colormaps["RdBu_r"].with_extremes(bad="#999999"), vmin=-limit, vmax=limit)
    fig.colorbar(im, ax=axes[row, 1], label="Out − in (electrons/s)", shrink=.8)
    for ax in axes[row]:
        ax.contour(aperture.astype(float), levels=[.5], colors="cyan", linewidths=.8)
        ax.plot(x, y, "+", color="#29cb64", ms=10)
        ax.set(xlabel="Cutout column", ylabel="Cutout row")
    axes[row, 0].set_title(f"Q{quarter} · median observed pixels")
    axes[row, 1].set_title(f"Q{quarter} · difference across {len(used)} events")
    peak_y, peak_x = np.unravel_index(np.nanargmax(difference), difference.shape)
    records.append(dict(quarter=int(quarter), path=path.relative_to(root).as_posix(), url=url, sha256=digest(path),
                        target_column=float(x), target_row=float(y), brightest_difference_column=int(peak_x),
                        brightest_difference_row=int(peak_y), events=used))
fig.suptitle("K08223.01 · observed Kepler pixels at the catalogue ephemeris\nGreen = target coordinate, cyan = pipeline aperture, gray = missing", fontsize=12)
save_figure(fig, root / "figures/kepler_followup/k08223_observed_pixels")
write_json(out / "pixel_audit.json", dict(executed_utc=utcnow(), products=records,
    selection="Two later quarters ranked by sum of squared PDCSAP depth-to-error ratios",
    script_sha256=text_digest(Path(__file__)), centroid_significance_measured=False,
    localization="Uncalibrated difference images with no PSF fit or neighbour exclusion"))
print(records)

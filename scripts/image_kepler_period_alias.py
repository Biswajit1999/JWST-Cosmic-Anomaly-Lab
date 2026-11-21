"""Compare actual pixel images of selected and other transits of Kepler-1437 b."""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import colormaps
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transits import phase_distance
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_long_period"
audit = json.loads((out / "alias_audit.json").read_text())
kic = audit["kic"]
inventory = json.loads((out / "products.json").read_text())["targets"]
target = next(row for row in inventory if row["kic"] == kic)
overlap = pd.read_csv(out / "alias_overlap.csv").query("split == 'replication'")
quarter = int(overlap.groupby("quarter").size().sort_index().idxmax())
for p in target["products"]:
    with fits.open(root / p["path"]) as h:
        if h[0].header["QUARTER"] == quarter:
            corresponding = p
            break
else:
    raise ValueError("No light curve found for selected quarter")
base = target["directory"]["url"].replace("lightcurves", "target_pixel_files")
cache = DownloadCache(root, 4096*1024**2)
listing = cache.fetch(base, f"exoplanets/kic{kic}_tpf_directory.html", survey="Kepler",
                      version="public archive", kind="directory")
name = Path(corresponding["path"]).name.replace("_llc.fits", "_lpd-targ.fits.gz")
assert name in re.findall(r'href="([^"]+_lpd-targ.fits.gz)"', listing.read_text())
path = cache.fetch(base + name, "exoplanets/" + name, survey="Kepler",
                   version="public archive", kind="observed_target_pixel_cube")
with fits.open(path) as h:
    assert h[0].header["KEPLERID"] == kic and h[0].header["QUARTER"] == quarter
    assert h[1].header["BJDREFI"] == 2454833 and h[1].header["TIMESYS"] == "TDB"
    assert h[1].columns["FLUX"].unit == "e-/s"
    d = h[1].data
    time, cube = np.array(d["TIME"]), np.array(d["FLUX"], float)
    good = np.isfinite(time) & (d["QUALITY"] == 0)
    sample_pixels = np.isfinite(cube[good]).all(axis=0)
    aperture = (h[2].data.astype(int) & 2) != 0
    wcs = WCS(h[2].header).celestial
    x, y = wcs.world_to_pixel_values(h[0].header["RA_OBJ"], h[0].header["DEC_OBJ"])
reference = audit["reference"]
period, epoch, duration = reference["koi_period"], reference["koi_time0bk"], reference["koi_duration"]/24
alias = audit["variants"]["original_peak"]
time, cube = time[good], cube[good][:, sample_pixels]
median = np.full(sample_pixels.shape, np.nan)
median[sample_pixels] = np.median(cube, axis=0)
grouped = {"Alias-selected transits": [], "Other published transits": []}
counts = []
for cycle in range(int(np.ceil((time.min()-epoch)/period)), int(np.floor((time.max()-epoch)/period))+1):
    center = epoch + period * cycle
    d = time - center
    inside = abs(d) <= duration/2
    left = (d >= -3*duration) & (d < -duration)
    right = (d <= 3*duration) & (d > duration)
    if min(inside.sum(), left.sum(), right.sum()) < 3:
        continue
    belongs = abs(float(phase_distance(center, alias["period_days"], alias["epoch_bkjd"]))) < duration
    label = "Alias-selected transits" if belongs else "Other published transits"
    grouped[label].append((cube[left].mean(axis=0) + cube[right].mean(axis=0))/2 - cube[inside].mean(axis=0))
    counts.append(dict(cycle=cycle, group=label, midpoint_bkjd=center, in_cadences=int(inside.sum())))
diff = {}
for label, copies in grouped.items():
    diff[label] = np.full(sample_pixels.shape, np.nan)
    if copies:
        diff[label][sample_pixels] = np.mean(copies, axis=0)
limit = max(float(np.nanmax(abs(image))) for image in diff.values() if np.isfinite(image).any())
fig = Figure(figsize=(12, 4.8), layout="constrained")
axes = fig.subplots(1, 3)
for index, (ax, (label, image)) in enumerate(zip(axes, [("Median observed field", median), *diff.items()], strict=True)):
    cmap = colormaps["magma" if index == 0 else "RdBu_r"].with_extremes(bad="#9f9f9f")
    im = ax.imshow(image, origin="lower", cmap=cmap,
                   vmin=None if not index else -limit, vmax=None if not index else limit)
    ax.contour(aperture.astype(float), levels=[.5], colors="cyan", linewidths=.8)
    ax.plot(x, y, "+", color="#48db69", ms=10)
    ax.set(title=label + (f" ({len(grouped[label])})" if index else ""),
           xlabel="Detector column in cutout", ylabel="Detector row in cutout")
    fig.colorbar(im, ax=ax, label="Flux (e⁻/s)" if not index else "Out − in (e⁻/s)", shrink=.72)
fig.suptitle(f"Kepler-1437 · actual Kepler Q{quarter} target pixels\nCyan = photometric aperture, green cross = catalogued target", fontsize=11)
save_figure(fig, root / "figures/kepler_long_period/kepler1437_pixels")
write_json(out / "alias_pixel_audit.json", dict(executed_utc=utcnow(), quarter=quarter,
    product=dict(url=base+name, path=path.relative_to(root).as_posix(), sha256=digest(path)),
    script_sha256=text_digest(__file__), comparison_cadences=counts,
    selection="Earliest later quarter with highest count of usable alias events",
    quantitative_centroid_or_blend_exclusion=False))
print(quarter, {k: len(v) for k, v in grouped.items()})

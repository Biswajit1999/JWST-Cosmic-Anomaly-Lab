"""Render an actual Kepler target and a locked-ephemeris difference image."""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import colormaps
from matplotlib.figure import Figure

from cosmic_anomaly_lab.companion_search import known_transit_mask
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_companions"
vetting = pd.read_csv(out / "vetting.csv")
leads = vetting.query("status == 'UNMATCHED_REPLICATION_REQUIRES_VETTING' and role == 'search_host'")
if leads.empty:
    leads = vetting.query("kic == 8644288 and rank == 1 and status == 'CATALOGUE_MATCH_OR_ALIAS'")
if leads.empty:
    leads = vetting.query("kic == 11904151 and rank == 1 and status == 'CATALOGUE_MATCH_OR_ALIAS'")
target = leads.iloc[0]
kic, rank = int(target.kic), int(target["rank"])
rows = pd.read_csv(out / "events.csv").query(
    "kic == @kic and rank == @rank and sign == 'observed' and split == 'later' and status == 'MEASURED'"
)
quarter = int(rows.groupby("quarter").size().sort_index().idxmax())
fit = pd.read_csv(out / "candidates.csv").query("kic == @kic and rank == @rank and sign == 'observed'").iloc[0]
star = next(r for r in json.loads((out / "selection.json").read_text())["targets"] if r["kic"] == kic)
products = next(r for r in json.loads((out / "products.json").read_text())["targets"] if r["kic"] == kic)
for p in products["products"]:
    with fits.open(root / p["path"]) as h:
        if h[0].header["QUARTER"] == quarter:
            break
else:
    raise ValueError(f"No light-curve record for quarter {quarter}")
name = Path(p["path"]).name.replace("_llc.fits", "_lpd-targ.fits.gz")
url = products["listing"]["url"].replace("lightcurves", "target_pixel_files") + name
cache = DownloadCache(root, 4096*1024**2)
file = cache.fetch(url, "exoplanets/"+name, survey="Kepler", version="public archive",
                   kind="observed_target_pixel_cube")
with fits.open(file) as h:
    assert h[0].header["KEPLERID"] == kic and h[0].header["QUARTER"] == quarter
    assert h[1].header["BJDREFI"] == 2454833 and h[1].header["TIMESYS"] == "TDB"
    assert h[1].columns["FLUX"].unit == "e-/s"
    t, flux, flag = (np.asarray(h[1].data[key]) for key in ["TIME", "FLUX", "QUALITY"])
    cadence = h[1].header["TIMEDEL"]
    aperture = (h[2].data.astype(int) & 2) != 0
    x, y = WCS(h[2].header).celestial.world_to_pixel_values(h[0].header["RA_OBJ"], h[0].header["DEC_OBJ"])
good = np.isfinite(t) & (flag == 0) & ~known_transit_mask(t, star["masks"], cadence)
comparisons, used = [], []
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    median = np.nanmedian(flux[good], axis=0)
    for center in rows[rows.quarter == quarter].midpoint:
        delta = t-center
        duration = fit.duration_days
        inside = good & (abs(delta) <= duration/2)
        before = good & (delta >= -3*duration) & (delta < -duration)
        after = good & (delta > duration) & (delta <= 3*duration)
        if min(inside.sum(), before.sum(), after.sum()) < 3:
            continue
        comparisons.append((np.nanmean(flux[before], axis=0) + np.nanmean(flux[after], axis=0))/2
                           - np.nanmean(flux[inside], axis=0))
        used.append(float(center))
    assert comparisons
    change = np.nanmean(comparisons, axis=0)
limit = float(np.nanmax(abs(change)))
fig = Figure(figsize=(9.8, 4.7), layout="constrained")
for ax, image, cmap, title, low, high in zip(fig.subplots(1, 2), [median, change], ["magma", "RdBu_r"],
        ["Median observed field", f"Out minus in transit, {len(comparisons)} events"],
        [None, -limit], [None, limit], strict=True):
    im = ax.imshow(image, origin="lower", cmap=colormaps[cmap].with_extremes(bad="#999999"),
                   vmin=low, vmax=high)
    ax.contour(aperture.astype(float), levels=[.5], colors="cyan", linewidths=.8)
    ax.plot(x, y, "+", color="#44e076", ms=9)
    ax.set(title=title, xlabel="Detector column in cutout", ylabel="Detector row in cutout")
    fig.colorbar(im, ax=ax, label="Flux (e⁻/s)", shrink=.75)
fig.suptitle(f"KIC {kic} · actual Kepler Q{quarter} pixels · {fit.period_days:.5f} days\nCyan = aperture, green cross = target coordinate",
             fontsize=11)
save_figure(fig, root / f"figures/kepler_companions/kic{kic}_target_pixels")
write_json(out / "pixel_audit.json", dict(executed_utc=utcnow(), kic=kic, rank=rank, quarter=quarter,
    product=dict(url=url, path=file.relative_to(root).as_posix(), sha256=digest(file)),
    centers_bkjd=used, assessment=target.status, script_sha256=text_digest(__file__),
    measured_centroid_uncertainty=False))
print(kic, quarter, len(used))

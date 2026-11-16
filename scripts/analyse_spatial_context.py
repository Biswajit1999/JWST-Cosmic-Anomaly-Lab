"""Render actual NIRCam data and audit nearby spectroscopic catalogue entries."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astropy.table import Table
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import BANDS, image_panel, read_images, save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/spatial_context"
products = json.loads((out / "image_products.json").read_text())["products"]
catalogue = Table.read(root / "data/raw/jades_dr4.fits", hdu="Obs_info")
coords = SkyCoord(catalogue["RA_TARG"]*u.deg, catalogue["Dec_TARG"]*u.deg)
neighbors, records = [], []
fig = Figure(figsize=(15, 15), layout="constrained")
axes = fig.subplots(len(products), 4)
for row, p in enumerate(products):
    if p["status"] != "DOWNLOADED":
        for ax in axes[row]:
            ax.text(.5, .5, "Image unavailable", ha="center")
        continue
    path = root / p["path"]
    assert digest(path) == p["sha256"]
    images, rgb, detail = read_images(path, p["ra"], p["dec"])
    center = SkyCoord(p["ra"]*u.deg, p["dec"]*u.deg)
    sep = coords.separation(center).arcsec
    nearby = np.where((sep > .1) & (sep <= 3))[0]
    for i in nearby:
        q = catalogue[i]
        z = float(q["z_Spec"])
        dv = 299792.458 * (z-p["z"]) / (1+p["z"]) if np.isfinite(z) and z > 0 and abs(z-p["z"])/(1+p["z"]) < .01 else None
        neighbors.append(dict(target=p["dr4_id"], neighbor=str(q["Unique_ID"]), ra=float(q["RA_TARG"]),
            dec=float(q["Dec_TARG"]), separation_arcsec=float(sep[i]), neighbor_z=z if z>0 else None,
            neighbor_z_flag=str(q["z_Spec_flag"]), target_z=p["z"], small_dz_velocity_proxy_kms=dv,
            neighbor_nircam_dr5_id=int(q["NIRCam_DR5_ID"]),
            projected_kpc_if_at_target_z=float(sep[i]*Planck18.kpc_proper_per_arcmin(p["z"]).value/60),
            physical_association="NOT_ESTABLISHED"))
    for col, band in enumerate([None, *BANDS]):
        ax = axes[row, col]
        image_panel(ax, images, rgb, SimpleNamespace(ra=p["ra"], dec=p["dec"]), band=band)
        ref = images.get(band, next(iter(images.values())))
        # Labels denote catalogue positions, not an automated merge or classification.
        for i in nearby:
            if sep[i] > 1.2:
                continue
            x, y = ref["wcs"].world_to_pixel_values(float(catalogue["RA_TARG"][i]), float(catalogue["Dec_TARG"][i]))
            ax.plot(x, y, "+", color="#ffd166", ms=9)
            ax.annotate(str(catalogue["NIRSpec_ID"][i]), (x, y), xytext=(8, -14),
                        textcoords="offset points", color="#ffd166", fontsize=8)
        name = "F444W / F277W / F150W RGB" if band is None else band.replace("-CLEAR", "")
        ax.set_title(p["dr4_id"] + "\n" + name, fontsize=9)
    record = dict(**p, image_details=detail, catalog_neighbors_within_3arcsec=len(nearby))
    records.append(record)
fig.suptitle("ACTUAL JWST NIRCam IMAGES · four Balmer follow-up targets\nCyan circle: target coordinate; gold crosses: nearby catalogue targets; 1″ bars; PSFs unmatched", fontsize=14)
save_figure(fig, root / "figures/spatial_context/target_images")
pd.DataFrame(neighbors).to_csv(out / "catalogue_neighbors.csv", index=False)
write_json(out / "summary.json", dict(executed_utc=utcnow(), products=records,
    catalog_sha256=digest(root / "data/raw/jades_dr4.fits"), script_sha256=text_digest(Path(__file__)),
    image_reader_sha256=text_digest(root / "src/cosmic_anomaly_lab/visual_evidence.py"),
    separation_method="Spherical coordinates, 0.1 < separation <= 3 arcsec; no automatic merging",
    cosmology="Astropy Planck18; projected distances conditional on target redshift",
    morphology_or_blending_measured=False, physical_classification="NOT_ESTABLISHED"))
print(pd.DataFrame(neighbors).to_string(index=False))

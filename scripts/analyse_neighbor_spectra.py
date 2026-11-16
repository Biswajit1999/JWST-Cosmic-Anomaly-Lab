"""Spatially distinct nearby targets: released slit polygons and signed cores."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib.figure import Figure
from matplotlib.path import Path as PolygonPath

from cosmic_anomaly_lab.continuum_audit import SIDEBANDS, bilateral_core
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import image_panel, read_images, save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/spatial_context"
image = next(p for p in json.loads((out / "image_products.json").read_text())["products"] if p["dr4_id"].endswith("_89464"))
near = json.loads((out / "neighbor_product.json").read_text())
dja = json.loads((root / "results/multiwavelength/dja_products.json").read_text())
target = next(p for p in dja if p["dr4_id"].endswith("_89464"))
coordinates = {89464: (image["ra"], image["dec"], image["z"]), 38432: (189.18617000457152, 62.27086361455187, 6.7170105817648125)}
center = SkyCoord(image["ra"]*u.deg, image["dec"]*u.deg)
records, geometry, hashes = [], [], {}
images, rgb, _ = read_images(root / image["path"], image["ra"], image["dec"])
fig = Figure(figsize=(13, 5), layout="constrained")
axes = fig.subplots(1, 2)
from types import SimpleNamespace

for ax, band in zip(axes, [None, "F444W-CLEAR"], strict=True):
    image_panel(ax, images, rgb, SimpleNamespace(ra=image["ra"], dec=image["dec"]), band=band)
    ax.set_title("NIRCam F444W/F277W/F150W" if band is None else "F444W with released S_REGION polygons")
for sid, product, color in [(89464, target, "#50d3e0"), (38432, near, "#ffd166")]:
    assert digest(root / product["path"]) == product["sha256"]
    hashes[product["path"]] = product["sha256"]
    ra, dec, z = coordinates[sid]
    with fits.open(root / product["path"]) as h:
        head = h["SPEC1D"].header
        assert int(head["SOURCEID"]) == sid and int(head["PROGRAM"]) == 1181
        assert SkyCoord(head["SRCRA"]*u.deg, head["SRCDEC"]*u.deg).separation(SkyCoord(ra*u.deg, dec*u.deg)).arcsec < .1
        region = h["SCI"].header["S_REGION"]
    assert region.startswith("POLYGON ICRS")
    vertices = np.array(region.split()[2:], float).reshape(-1, 2)
    sky = SkyCoord(vertices[:, 0]*u.deg, vertices[:, 1]*u.deg).transform_to(center.skyoffset_frame())
    xy = np.column_stack([sky.lon.arcsec, sky.lat.arcsec])
    for other, (ora, odec, _) in coordinates.items():
        sk = SkyCoord(ora*u.deg, odec*u.deg).transform_to(center.skyoffset_frame())
        point = np.array([sk.lon.arcsec, sk.lat.arcsec])
        distances = []
        for a, b in zip(xy, np.roll(xy, -1, axis=0), strict=True):
            v = b-a
            fraction = np.clip(np.dot(point-a, v)/np.dot(v, v), 0, 1)
            distances.append(float(np.linalg.norm(point-(a+fraction*v))))
        geometry.append(dict(slit_target=sid, tested_target=other, centroid_inside=bool(PolygonPath(xy).contains_point(point)),
                             minimum_edge_distance_arcsec=min(distances), region=region,
                             scope="Released SCI header polygon, not every nod or PSF throughput"))
    for ax in axes:
        wcs = images["F444W-CLEAR"]["wcs"]
        x, y = wcs.world_to_pixel_values(ra, dec)
        ax.plot(x, y, "+", color=color, ms=9)
        ax.annotate(str(sid), (x, y), xytext=(9, -20 if sid == 38432 else 17), textcoords="offset points", color=color)
    px, py = images["F444W-CLEAR"]["wcs"].world_to_pixel_values(vertices[:, 0], vertices[:, 1])
    axes[1].plot(np.r_[px, px[0]], np.r_[py, py[0]], color=color, lw=1.2, label=f"{sid} G395M")
    spec = read_spectrum(root / product["path"], wavelength_range=(.6, 5.5))
    unit = u.erg/u.s/u.cm**2/u.AA
    f = spec.flux.to_value(unit, equivalencies=u.spectral_density(spec.wavelength))/1e-20
    e = spec.error.to_value(unit, equivalencies=u.spectral_density(spec.wavelength))/1e-20
    for line, rest in [("hg", .434165), ("hb", .486264)]:
        wavelength = rest*(1+z)
        velocity = (spec.wavelength.value/wavelength-1)*299792.458
        for width in [400, 600]:
            for name, windows in SIDEBANDS.items():
                records.append(dict(source_id=sid, line=line, fixed_z=z, half_width_kms=width, sidebands=name,
                    **bilateral_core(velocity, f, e, spec.good, wavelength*1e4, width, windows)))
axes[1].legend(fontsize=8, loc="lower right")
fig.suptitle("GN-89464 and GN-38432: observed spatial context\nDistinct catalogue entries; 0.707″ apart; proximity does not establish contamination or a bound pair", fontsize=12)
save_figure(fig, root / "figures/spatial_context/gn89464_neighbor")
pd.DataFrame(geometry).to_csv(out / "slit_geometry.csv", index=False)
data = pd.DataFrame(records)
data.to_csv(out / "neighbor_signed_cores.csv", index=False)
ratios = []
for (sid, width, bands), group in data.groupby(["source_id", "half_width_kms", "sidebands"]):
    group = group.set_index("line")
    n, d = group.loc["hg"], group.loc["hb"]
    eligible = n.status == d.status == "MEASURED_CORE" and n.flux/n.error >= 5 and d.flux/d.error >= 5
    ratios.append(dict(source_id=int(sid), half_width_kms=int(width), sidebands=bands, eligible=bool(eligible),
                       ratio=float(n.flux/d.flux) if eligible else None))
pd.DataFrame(ratios).to_csv(out / "neighbor_ratios.csv", index=False)
write_json(out / "neighbor_summary.json", dict(executed_utc=utcnow(), input_sha256=hashes,
    script_sha256=text_digest(Path(__file__)), physical_association="NOT_ESTABLISHED",
    psf_contamination_quantified=False, raw_detector_reduction=False))
print(pd.DataFrame(geometry).drop(columns=["region", "scope"]).to_string(index=False))
print(pd.DataFrame(ratios).to_string(index=False))

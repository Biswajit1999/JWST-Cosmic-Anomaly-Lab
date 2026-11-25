"""Identify the origin of the F090W diffraction spike at GOODS-N 58441."""

import json
from io import StringIO
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from download_jades_reionization_sample import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.colors import AsinhNorm
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from read_jwst_i2d_cutout import MastImageSections
from scipy.ndimage import gaussian_filter, map_coordinates

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

MOSAIC = "jw06434-c1034_t000_nircam_clear-f090w_i2d.fits"
STEM = "jades_58441_diffraction_spike"
GAIA_RADIUS = 3
HALF_ARCSEC = 18


def gaia_star(cache, position, epoch):
    params = {
        "-source": "I/355/gaiadr3",
        "-c": f"{position.ra.deg:.9f} {position.dec.deg:.9f}",
        "-c.rs": GAIA_RADIUS,
        "-out.max": 30,
        "-out": "Source,RA_ICRS,DE_ICRS,Gmag,Plx,pmRA,pmDE",
    }
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode(params)
    path = cache.fetch(
        url,
        "metadata/jades_58441_spike_gaia.tsv",
        survey="CDS VizieR",
        version="Gaia DR3 I/355/gaiadr3",
        kind="public_foreground_catalogue",
    )
    content = path.read_text(encoding="utf-8")
    if "#Name: I/355/gaiadr3" not in content or "#Truncated" in content:
        raise ValueError("Incomplete Gaia catalogue response")
    tab = pd.read_csv(StringIO(content), comment="#", sep="\t", dtype=str)
    columns = ["RA_ICRS", "DE_ICRS", "Gmag", "Plx", "pmRA", "pmDE"]
    tab[columns] = tab[columns].apply(pd.to_numeric, errors="coerce")
    rows = []
    for row in tab.dropna(subset=columns[:2]).to_dict(orient="records"):
        sky = SkyCoord(row["RA_ICRS"], row["DE_ICRS"], unit="deg")
        moved = sky
        if np.isfinite([row["pmRA"], row["pmDE"]]).all():
            elapsed = (epoch - Time(2016.0, format="jyear")).to_value(u.yr)
            moved = sky.spherical_offsets_by(
                row["pmRA"] * elapsed * u.mas, row["pmDE"] * elapsed * u.mas
            )
        data = {
            key: value if np.isfinite(value) else None
            for key, value in row.items()
            if key in columns
        }
        rows.append(
            dict(
                data,
                source_id=row["Source"],
                separation_j2016_arcsec=float(sky.separation(position).arcsec),
                separation_propagated_arcsec=float(moved.separation(position).arcsec),
                propagated_icrs_deg=[float(moved.ra.deg), float(moved.dec.deg)],
            )
        )
    return dict(
        receipt=receipt(path, url),
        radius_arcsec=GAIA_RADIUS,
        propagated_epoch_utc=epoch.isot,
        rows=rows,
    )


def find_ridge(stamp, star_xy):
    image, wcs, x, y, x0, y0 = stamp
    centre = np.array([x, y])
    delta = centre - np.array(star_xy)
    direction = delta / np.linalg.norm(delta)
    normal = np.array([-direction[1], direction[0]])
    scale = float(wcs.proj_plane_pixel_scales()[0].to_value(u.arcsec))
    off = np.linspace(-0.50, 0.50, 101)
    profiles = []
    for displacement in (-2.0, -1.0, 0.0, 1.0, 2.0):
        pos = centre[:, None] + (displacement * direction[:, None] + off * normal[:, None]) / scale
        values = map_coordinates(image, pos[::-1], order=1, mode="constant", cval=np.nan)
        near = np.flatnonzero((abs(off) <= 0.22) & np.isfinite(values))
        if not len(near):
            raise ValueError("Missing spike profile near the target")
        i = near[np.argmax(values[near])]
        background = float(np.nanmedian(values[abs(off) > 0.33]))
        profiles.append(
            dict(
                along_ridge_arcsec=displacement,
                offset_from_star_target_axis_arcsec=float(off[i]),
                peak_minus_sideband_MJy_sr=float(values[i] - background),
                offset_arcsec=off.tolist(),
                sky_subtracted_MJy_sr=(values - background).tolist(),
            )
        )
    star = wcs.pixel_to_world(star_xy[0] + x0, star_xy[1] + y0)
    target = wcs.pixel_to_world(x + x0, y + y0)
    return dict(
        bright_source_icrs_deg=[float(star.ra.deg), float(star.dec.deg)],
        target_separation_arcsec=float(star.separation(target).arcsec),
        position_angle_at_star_deg=float(star.position_angle(target).deg),
        profiles=profiles,
    ), star


def plot(stamp, star_xy, geometry):
    image, wcs, x, y, *_ = stamp
    scale = float(wcs.proj_plane_pixel_scales()[0].to_value(u.arcsec))
    cy, cx = image.shape
    extent = np.array([-x, cx - x, -y, cy - y]) * scale
    finite = image[np.isfinite(image)]
    lo = float(np.median(finite))
    vmax = float(np.percentile(finite, 99.8))
    fig = Figure(figsize=(10, 4.9), layout="constrained")
    image_ax, profile_ax = fig.subplots(1, 2)
    image_ax.imshow(
        image,
        origin="lower",
        extent=extent,
        cmap="gray",
        norm=AsinhNorm(linear_width=0.016, vmin=lo - 0.045, vmax=vmax),
    )
    image_ax.add_patch(Circle((0, 0), 0.4, fill=False, edgecolor="#47dcf0", lw=1.4))
    mark = (np.array(star_xy) - [x, y]) * scale
    image_ax.annotate(
        "Gaia foreground star",
        xy=mark,
        xytext=(mark[0] - 5, mark[1] - 6),
        color="white",
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="white", lw=0.8),
    )
    image_ax.set(
        xlabel="Mosaic x offset (arcsec)",
        ylabel="Mosaic y offset (arcsec)",
        title="SAPPHIRES F090W · 8 January 2025",
    )
    for p in geometry["profiles"]:
        t = p["along_ridge_arcsec"]
        profile_ax.plot(
            p["offset_arcsec"],
            p["sky_subtracted_MJy_sr"],
            label=f"{t:+.0f}″",
            lw=2 if t == 0 else 1,
        )
    profile_ax.axvline(0, color="0.4", lw=0.75, ls="--")
    profile_ax.set(
        xlabel="Offset across star–target axis (arcsec)",
        ylabel="Signal above sidebands (MJy/sr)",
        title="Cross-sections near GOODS-N 58441",
    )
    profile_ax.legend(title="Along axis from target", frameon=False, fontsize=8)
    fig.suptitle("A foreground star’s diffraction spike crosses GOODS-N 58441", fontsize=13)
    save_figure(fig, ROOT / f"figures/{STEM}")


def main():
    path = OUT / "jades_58441_f090w_exposures.json"
    native = json.loads(path.read_text())
    target = SkyCoord(*native["position_icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, 28_000 * 1024**2)
    stamp, mosaic = MastImageSections(cache, MOSAIC).cutout(target, HALF_ARCSEC)
    image = stamp[0]
    centred = np.where(np.isfinite(image), image - np.nanmedian(image), 0.0)
    y, x = np.unravel_index(gaussian_filter(centred, 2.5).argmax(), image.shape)
    geometry, star = find_ridge(stamp, (x, y))
    match = gaia_star(cache, star, Time(mosaic["date_utc"], scale="utc"))
    # PA_V3 is recorded in each native science extension.
    orientations = []
    for exposure in native["exposures"]:
        if exposure["status"] != "MEASURED":
            continue
        rec = exposure["ranges"][1]
        native_header = ROOT / rec["path"]
        if digest(native_header) != rec["sha256"]:
            raise ValueError("Native science header checksum mismatch")
        hdr = fits.Header.fromstring(native_header.read_bytes().decode("ascii", errors="replace"))
        if hdr["EXTNAME"] != "SCI":
            raise ValueError("Native spectrum receipt did not identify the science header")
        orientations.append(dict(exposure=exposure["exposure"], pa_v3_deg=float(hdr["PA_V3"])))
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            native_exposure_summary_sha256=text_digest(path),
            mosaic=mosaic,
            geometry=geometry,
            gaia=match,
            orientations=orientations,
            scope="Optical origin check. Bright-star centre is a smoothed peak near a saturated core. "
            "Gaia positions were linearly propagated with proper motion alone. "
            "Coadd and native images of SAPPHIRES share exposures.",
        ),
    )
    plot(stamp, (x, y), geometry)
    print(
        geometry["target_separation_arcsec"], geometry["position_angle_at_star_deg"], match["rows"]
    )
    print(
        [
            (
                row["along_ridge_arcsec"],
                row["offset_from_star_target_axis_arcsec"],
                row["peak_minus_sideband_MJy_sr"],
            )
            for row in geometry["profiles"]
        ]
    )


if __name__ == "__main__":
    main()

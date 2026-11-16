"""Observed 2D spectral sheets for three stable-core priorities and a nod control."""

import json
import re
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib.colors import AsinhNorm
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/jades_balmer_census"
queue = pd.read_csv(out / "followup_queue.csv").set_index("dr4_id")
measurements = pd.read_csv(out / "signed_core_measurements.csv")
products = json.loads((out / "spectral_products.json").read_text())
table = pd.read_csv(out / "aperture_measurements.csv.gz")
table = table[table.aperture_pixels == 5].set_index("dr4_id")
rest = dict(ha=.656452, hb=.486264, hg=.434165)
selected = []
for sid, target in queue.iterrows():
    line = "ha" if target.direction == "low_ha_over_hb" else "hg"
    m = measurements[measurements.dr4_id == sid]
    n = m[m.line == line].set_index(["aperture_pixels", "nods", "half_width_kms"])
    d = m[m.line == "hb"].set_index(["aperture_pixels", "nods", "half_width_kms"])
    good = (n.flux / n.error >= 5) & (d.flux / d.error >= 5)
    if len(good) == 8 and good.all():
        selected.append(sid)
control = "goods-n-mediumhst_25734"
selected.append(control)
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
records = []
for sid in selected:
    target = queue.loc[sid]
    images = {}
    for product in products:
        if product["dr4_id"] != sid or product["status"] != "DOWNLOADED":
            continue
        filename = Path(product["path"]).name.replace("_x1d.fits", "_s2d.fits")
        tier = sid.rsplit("_", 1)[0]
        grating, filt = product["grating_filter"].split("_")
        listing = root / f"data/raw/metadata/dr4_{tier}_{grating}_{filt}.html"
        assert filename in re.findall(r'href="([^\"]+)"', listing.read_text())
        url = urljoin(product["directory_url"], filename)
        path = cache.fetch(url, "independent/" + filename, survey="JADES", version="DR4",
                           kind="observed_spectral_image")
        with fits.open(path) as hdus:
            h = hdus[0].header
            sep = SkyCoord(h["RA"] * u.deg, h["DEC"] * u.deg).separation(
                SkyCoord(target.ra * u.deg, target.dec * u.deg)).arcsec
            assert sep <= .2 and int(h["PROPOSID"]) == target.pid
            assert int(h["TARGET"]) == int(sid.rsplit("_", 1)[1])
            wave = (hdus["WAVELENGTH"].data * u.Unit(hdus["WAVELENGTH"].header["BUNIT"])).to_value(u.um)
            assert np.all(np.diff(wave) > 0)
            data = {}
            for nods, extension in [(3, "FLUX"), (2, "FLUX_2NOD")]:
                image = hdus[extension]
                error = hdus["FLUX_ERR" if nods == 3 else "FLUX_ERR_2NOD"].data
                one_d = read_dr4_spectrum(root / product["path"], 5, nods, wavelength_range=(.6, 5.5))
                assert np.allclose(wave, one_d.wavelength.value, rtol=0, atol=1e-8)
                scale = (1 * u.Unit(image.header["BUNIT"])).to_value(u.erg / u.s / u.cm**2 / u.AA) / 1e-20
                y = image.header["CRVAL2"] + (np.arange(image.data.shape[0]) + 1 - image.header["CRPIX2"]) * image.header["CDELT2"]
                valid = np.isfinite(error) & (error > 0) & one_d.good[None, :]
                data[nods] = (y, np.where(valid, image.data * scale, np.nan))
            images[product["grating_filter"]] = (wave.copy(), data)
        records.append(dict(dr4_id=sid, path=path.relative_to(root).as_posix(), sha256=digest(path),
                            url=url, header_separation_arcsec=float(sep),
                            selection="nod_sensitivity_control" if sid == control else "both_cores_SN5_in_all_8_variants"))
    fig = Figure(figsize=(15, 7), layout="constrained")
    axes = fig.subplots(2, 3)
    for col, (line, wavelength) in enumerate(rest.items()):
        label = str(table.loc[sid, line + "_filter"]).removesuffix("_ext")
        if label not in images:
            for ax in axes[:, col]:
                ax.text(.5, .5, "No accepted image", ha="center", transform=ax.transAxes)
            continue
        wave, data = images[label]
        v = (wave / (wavelength * (1 + target.z)) - 1) * 299792.458
        keep = abs(v) <= 2500
        values = np.concatenate([data[n][1][:, keep].ravel() for n in [2, 3]])
        finite = values[np.isfinite(values)]
        if finite.size == 0 or keep.sum() < 2:
            for ax in axes[:, col]:
                ax.text(.5, .5, "No image coverage supported\nby usable 1D measurements", ha="center", transform=ax.transAxes)
                ax.set_title(line.upper())
            continue
        scale = max(float(np.percentile(abs(finite), 99)), 1e-8)
        norm = AsinhNorm(linear_width=scale / 10, vmin=-scale, vmax=scale)
        for row, nods in enumerate([3, 2]):
            ax = axes[row, col]
            y, flux = data[nods]
            mesh = ax.pcolormesh(v[keep], y, flux[:, keep], shading="nearest", cmap="RdBu_r", norm=norm)
            ax.axvline(0, color="black", ls=":", lw=.7)
            ax.axhline(0, color="black", ls=":", lw=.7)
            ax.set(xlabel="Velocity at catalogue z (km/s)", ylabel="Released spatial coordinate\n(fractional shutter)",
                   title=f"{line.upper()} · {nods} nods · {label}")
        fig.colorbar(mesh, ax=axes[:, col].tolist(), label="Fλ per image pixel (10⁻²⁰ cgs/Å)",
                     shrink=.8, extend="both", ticks=[-scale, 0, scale], format="%.2g")
    fig.suptitle(f"TELESCOPE SPECTRAL IMAGES · {sid}\nSigned flux; same stretch between nods within each line; not sky-registered images", fontsize=12)
    save_figure(fig, root / "figures/jades_balmer_images" / sid)
    print(sid, "spectral image sheet saved", flush=True)
write_json(out / "spectral_images.json", dict(executed_utc=utcnow(), products=records,
    selection_is_exploratory=True, targets=selected,
    input_sha256={"scripts/image_jades_balmer_queue.py": text_digest(Path(__file__)),
                  "results/jades_balmer_census/signed_core_measurements.csv": text_digest(out / "signed_core_measurements.csv"),
                  "results/jades_balmer_census/followup_queue.csv": text_digest(out / "followup_queue.csv"),
                  "results/jades_balmer_census/spectral_products.json": text_digest(out / "spectral_products.json"),
                  "results/jades_balmer_census/aperture_measurements.csv.gz": digest(out / "aperture_measurements.csv.gz")},
    mask="Finite positive 2D error and valid released 5px 1D wavelength bin for the same nod variant; signed flux retained",
    interpretation="Visual contamination/nod inspection; not deblending, detector reprocessing or independent confirmation"))

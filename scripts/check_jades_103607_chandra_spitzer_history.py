"""Measure historical Chandra count limits and Spitzer confusion at JADES 103607."""

import gzip
import json
import tarfile
from io import BytesIO

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from scipy.stats import binom, chi2
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

SOURCE = OUT / "module_b_001_native_november_fade.json"
TAG = "module_b_001_chandra_spitzer_history"
CDFS_NAME = "luo2017-cdfs-data-products03.tar.gz"
CDFS_URL = f"https://zenodo.org/records/8021940/files/{CDFS_NAME}?download=1"
IRSA = "https://irsa.ipac.caltech.edu/data/SPITZER/GOODS/images/DR2/South"


def centered_image(hdu, position, half):
    x, y = map(float, WCS(hdu.header, fix=False).celestial.world_to_pixel(position))
    ix, iy = round(x), round(y)
    ny, nx = hdu.data.shape
    if not (half <= ix < nx - half and half <= iy < ny - half):
        raise ValueError("Target cutout outside archived image")
    data = np.array(hdu.data[iy - half : iy + half + 1, ix - half : ix + half + 1], dtype=float)
    return data, (x - ix + half, y - iy + half)


def chandra(cache, sky):
    path = cache.fetch(
        CDFS_URL,
        f"cross_detector/{CDFS_NAME}",
        survey="Chandra Deep Field-South",
        version="7 Ms Luo et al. 2017 Zenodo 8021940",
        kind="observed_survey_images",
    )
    summaries, display = {}, {}
    with tarfile.open(path, "r:gz") as bundle:
        for band in ("0p5to2", "2to7", "0p5to7"):
            label = "CDFS-7Ms-" + band
            names = [
                f"{label}-asca-im-bin1-01.fits.gz",
                f"{label}-bin1-02.emap.fits.gz",
            ]
            image_data, centers, headers = [], [], []
            for name in names:
                with bundle.extractfile("./" + name) as raw:
                    data = gzip.decompress(raw.read())
                with fits.open(BytesIO(data)) as hdul:
                    image, center = centered_image(hdul[0], sky, 45)
                    image_data.append(image)
                    centers.append(center)
                    headers.append(hdul[0].header.copy())
            counts, eff_exposure = image_data
            if (
                not np.allclose(centers, centers[0], atol=0.02, rtol=0)
                or headers[1]["BUNIT"] != "cm**2 s"
            ):
                raise ValueError("Chandra products are not on the same count/exposure grid")
            arcsec_per_pix = abs(headers[0]["CDELT2"]) * 3600
            yy, xx = np.indices(counts.shape)
            x, y = centers[0]
            r = np.hypot(xx - x, yy - y) * arcsec_per_pix
            annulus = (r >= 8) & (r <= 20)
            annulus_count = int(counts[annulus].sum())
            annulus_pixels = int(annulus.sum())
            options = {}
            for radius in (1.0, 1.5, 2.0, 3.0):
                use = r <= radius
                total = int(counts[use].sum())
                pixels = int(use.sum())
                bg = annulus_count * pixels / annulus_pixels
                p = pixels / (pixels + annulus_pixels)
                upper = float(chi2.ppf(0.95, 2 * (total + 1)) / 2)
                options[str(radius)] = dict(
                    observed_counts=total,
                    expected_background_counts=bg,
                    aperture_pixels=pixels,
                    signed_net_counts=total - bg,
                    local_background_binomial_tail=float(
                        binom.sf(total - 1, total + annulus_count, p)
                    ),
                    counts_upper95_no_background_subtraction=upper,
                )
            summaries[band] = dict(
                archive_members=names,
                background_annulus_arcsec=[8, 20],
                background_annulus_counts=annulus_count,
                background_annulus_pixels=annulus_pixels,
                center_exposure_cm2_s=float(eff_exposure[round(y), round(x)]),
                radii_arcsec=options,
                observation_start=headers[0]["DATE-OBS"],
                observation_end=headers[0]["DATE-END"],
            )
            display[band] = counts, x, y, arcsec_per_pix
    return dict(product=receipt(path, CDFS_URL), bands=summaries), display


def irac_aperture(data, exp, x, y):
    yy, xx = np.ogrid[: data.shape[0], : data.shape[1]]
    r = np.hypot(xx - x, yy - y) * 0.6
    ap, sky = r <= 1.8, (r >= 2.4) & (r <= 4.2)
    if not np.all(exp[ap | sky] > 0) or not np.isfinite(data[ap | sky]).all():
        return np.nan
    return float(data[ap].sum() - ap.sum() * np.median(data[sky]))


def spitzer(cache, sky):
    results, display = {}, {}
    for channel in (1, 2):
        names = [f"s_irac_{channel}_s2_v0.30_{kind}.fits" for kind in ("sci", "exp")]
        images, centroids, headers, used = [], [], [], []
        for name in names:
            url = f"{IRSA}/{name}"
            path = cache.fetch(
                url,
                f"cross_detector/{name}",
                survey="Spitzer GOODS-S",
                version="IRAC epoch 2 DR2 v0.30",
                kind="observed_survey_images",
            )
            used.append(receipt(path, url))
            with fits.open(path, memmap=False) as hdul:
                image, center = centered_image(hdul[0], sky, 52)
                images.append(image)
                centroids.append(center)
                headers.append(hdul[0].header.copy())
        data, exp = images
        x, y = centroids[0]
        w = WCS(headers[0], fix=False).celestial
        scale = abs(np.linalg.det(w.pixel_scale_matrix)) ** 0.5 * 3600
        if not np.allclose(centroids, [x, y], rtol=0, atol=0.02) or abs(scale - 0.6) > 0.001:
            raise ValueError("Weight and science mosaics must share the IRAC pixel grid")
        units = headers[0]["BUNIT"].strip()
        if units != "DN/sec" or headers[1]["BUNIT"].strip() != "sec":
            raise ValueError(f"Unexpected GOODS IRAC science/exposure units {units}")
        factor = float(headers[0]["FLUXCONV"]) * 1000
        data *= factor
        local = float(exp[round(y), round(x)])
        target = irac_aperture(data, exp, x, y)
        controls = np.asarray(
            [
                irac_aperture(data, exp, i, j)
                for i in range(9, 97, 9)
                for j in range(9, 97, 9)
                if np.hypot(i - x, j - y) > 10
            ]
        )
        usable = controls[np.isfinite(controls)]
        scatter = None
        if len(usable) >= 20:
            scatter = float(1.4826 * np.median(abs(usable - np.median(usable))))
        results[str(channel)] = dict(
            products=used,
            calibration_njy_per_dn_second=factor,
            center_exposure_map_value=local,
            exposure_header_unit=headers[1].get("BUNIT"),
            target_blended_aperture_njy=target if np.isfinite(target) else None,
            control_aperture_count=len(usable),
            control_aperture_njy=usable.tolist(),
            control_scatter_njy=scatter,
            five_times_control_scatter_njy=5 * scatter if scatter else None,
        )
        display[str(channel)] = data, exp, x, y
    return dict(observation_window_utc=["2004-08-12", "2004-08-18"], channels=results), display


def figure(chandra_images, irac_images):
    fig = Figure(figsize=(10, 4.2), layout="constrained")
    axes = fig.subplots(1, 3)
    for ax, band, label in zip(axes[:2], ("0p5to2", "2to7"), ("0.5–2 keV", "2–7 keV"), strict=True):
        image, x, y, scale = chandra_images[band]
        ax.imshow(np.sqrt(image), origin="lower", cmap="Blues", vmin=0, vmax=2.5)
        ax.add_patch(Circle((x, y), 1.5 / scale, fill=False, color="#b94b34", lw=1.1))
        ax.set(xlim=(x - 19, x + 19), ylim=(y - 19, y + 19), title=f"Chandra · {label}")
    data, exp, x, y = irac_images["1"]
    if exp[round(y), round(x)] == 0:
        data, exp, x, y = irac_images["2"]
        label = "4.5"
    else:
        label = "3.6"
    image = np.where(exp > 0, data, np.nan)
    scale = np.nanpercentile(image, 96)
    axes[2].imshow(image, origin="lower", cmap="magma", vmin=-scale * 0.1, vmax=scale)
    axes[2].add_patch(Circle((x, y), 1.8 / 0.6, fill=False, color="#86e0ed", lw=1.1))
    axes[2].set(xlim=(x - 15, x + 15), ylim=(y - 15, y + 15), title=f"Spitzer IRAC · {label} µm")
    for ax in axes:
        ax.set(xticks=[], yticks=[])
    fig.suptitle("JADES 103607 · pre-JWST archive images and target apertures", fontsize=13)
    save_figure(fig, FIGURES / TAG)


def main():
    coords = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*coords, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    xr, ximages = chandra(cache, sky)
    ir, iimages = spitzer(cache, sky)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=coords,
            position_source_sha256=text_digest(SOURCE),
            chandra=xr,
            spitzer=ir,
            scope="Chandra counts use the full 7 Ms stack and have no PSF or spectral flux correction. "
            "IRAC uses blended apertures and empirical spatial scatter from one legacy observing epoch.",
        ),
    )
    figure(ximages, iimages)
    for band, r in xr["bands"].items():
        print(band, r["center_exposure_cm2_s"], r["radii_arcsec"]["1.5"])
    for key, row in ir["channels"].items():
        print(
            "IRAC",
            key,
            row["center_exposure_map_value"],
            row["target_blended_aperture_njy"],
            row["five_times_control_scatter_njy"],
        )


if __name__ == "__main__":
    main()

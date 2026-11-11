"""WCS-registered diagnostic apertures for DJA cutouts; no intrinsic morphology claims."""

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


def aperture_measurement(
    image,
    weight,
    wcs,
    ra_deg,
    dec_deg,
    *,
    jy_per_unit,
    radius_arcsec=0.3,
    sky_annulus_arcsec=(0.8, 1.2),
):
    """Pixel-centre aperture and clipped-median sky, diagonal WHT error only.

    Covariance from drizzling, aperture correction, PSF mismatch, and blending
    are unmodelled. Returned fluxes are aperture diagnostics, not total SED fluxes.
    """
    if image.shape != weight.shape or image.ndim != 2:
        raise ValueError("Image and weight must be matching two-dimensional arrays")
    if jy_per_unit <= 0 or not 0 < radius_arcsec < sky_annulus_arcsec[0] < sky_annulus_arcsec[1]:
        raise ValueError("Invalid aperture or calibration")
    x0, y0 = wcs.world_to_pixel_values(ra_deg, dec_deg)
    scales = proj_plane_pixel_scales(wcs.celestial) * 3600
    yy, xx = np.indices(image.shape)
    radius = np.hypot((xx - x0) * scales[0], (yy - y0) * scales[1])
    ap = radius <= radius_arcsec
    sky = (radius >= sky_annulus_arcsec[0]) & (radius < sky_annulus_arcsec[1])
    valid = np.isfinite(image) & np.isfinite(weight) & (weight > 0)
    # Reject clipped apertures instead of reporting artificially low fluxes.
    margin = sky_annulus_arcsec[1] / np.min(scales)
    if min(x0, y0, image.shape[1] - 1 - x0, image.shape[0] - 1 - y0) < margin:
        raise ValueError("Sky annulus extends beyond cutout")
    if ap.sum() == 0 or (ap & valid).sum() < 0.95 * ap.sum() or (sky & valid).sum() < 20:
        raise ValueError("Insufficient valid aperture/sky pixels")
    sky_values = sigma_clip(image[sky & valid], sigma=3, maxiters=5).compressed()
    bg = np.median(sky_values)
    bgvar = np.pi / 2 * np.var(sky_values, ddof=1) / len(sky_values)
    mask = ap & valid
    flux = np.sum(image[mask] - bg) * jy_per_unit * 1e6
    error = np.sqrt(np.sum(1 / weight[mask]) + mask.sum() ** 2 * bgvar) * jy_per_unit * 1e6
    return dict(
        aperture_flux_uJy=float(flux),
        diagonal_error_uJy=float(error),
        sn=float(flux / error),
        upper_limit_3sigma_uJy=float(max(0, flux) + 3 * error),
        aperture_radius_arcsec=radius_arcsec,
        background_uJy_pixel=float(bg * jy_per_unit * 1e6),
        n_aperture_pixels=int(mask.sum()),
        pixel_scale_arcsec=float(np.mean(scales)),
        flags=[
            "PSF_UNMATCHED",
            "NO_APERTURE_CORRECTION",
            "BLENDING_UNASSESSED",
            "DRIZZLE_COVARIANCE_UNMODELLED",
        ],
    )


def measure_cutout(path, ra_deg, dec_deg):
    output, panels = [], []
    with fits.open(path) as hdus:
        hdus.verify("exception")
        for hdu in hdus:
            if hdu.header.get("EXTVER") != "SCI":
                continue
            band = hdu.header["FILTER"]
            weights = [
                h
                for h in hdus
                if h.header.get("FILTER") == band and h.header.get("EXTVER") == "WHT"
            ]
            if len(weights) != 1 or "PHOTFNU" not in hdu.header:
                raise ValueError("Missing or ambiguous inverse variance/calibration")
            wcs = WCS(hdu.header).celestial
            result = aperture_measurement(
                hdu.data, weights[0].data, wcs, ra_deg, dec_deg, jy_per_unit=hdu.header["PHOTFNU"]
            )
            output.append(dict(filter=band, **result))
            panels.append((band, np.array(hdu.data), wcs))
    if not output:
        raise ValueError("No SCI image extensions")
    return output, panels

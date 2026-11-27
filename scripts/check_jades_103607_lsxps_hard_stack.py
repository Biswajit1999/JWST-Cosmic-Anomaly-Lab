"""Inspect the archived image behind the lone hard-band LSXPS flag at JADES 103607."""

import json
import warnings
from urllib.parse import urlsplit

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from measure_jades_103607_swift_xrt_history import PAIR_ARCSEC, SXPS_API
from measure_jades_103607_swift_xrt_history import TAG as HISTORY
from scipy.stats import poisson
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_lsxps_hard_stack"
NATIVE_APERTURE_PIXELS = 12


def archive_images(cache, obs):
    query = dict(
        APIFunc="getSXPSDatasetImages", APIVersion="1.0.5", whichCat="LSXPS",
        ObsID=obs, bands=["hard"], types="all", getRegions=True,
    )
    path = cache.fetch(
        SXPS_API, f"metadata/{TAG}_{obs}_products.json", survey="UK Swift Science Data Centre",
        version="LSXPS hard-band stack product locator", kind="observed_stack_products", json_payload=query,
    )
    products = json.loads(path.read_text(encoding="utf-8"))
    if products.get("OK") != 1 or "SupersededBy" in products:
        raise ValueError(f"Expected an active LSXPS stack {obs}")
    paths, sources = {}, {"locator": receipt(path, SXPS_API)}
    for key, field in (
        ("counts", "Image"), ("background", "BackgroundMap"),
        ("background_with_catalogue_sources", "BackgroundMapWithSources"),
    ):
        url = products["Hard"][field]
        u = urlsplit(url)
        if (u.scheme, u.hostname) != ("https", "www.swift.ac.uk") or not u.path.startswith(f"/LSXPS/obs/{obs}/"):
            raise ValueError(f"Unexpected LSXPS archive image URL {url}")
        file = cache.fetch(
            url, f"cross_detector/lsxps/{u.path.rsplit('/', 1)[-1]}", survey="UK Swift Science Data Centre",
            version="LSXPS archived hard-band stack", kind="observed_image_or_calibration",
        )
        paths[key] = file
        sources[key] = receipt(file, url)
    return paths, sources, query


def cutout(paths, point):
    arrays, transforms = {}, {}
    xmin = ymin = 0
    for label, path in paths.items():
        with fits.open(path) as hdus:
            wcs = WCS(hdus[0].header, fix=False).celestial
            transforms[label] = wcs
            image = hdus[0].data
            if label == "counts":
                x, y = map(float, wcs.world_to_pixel(point))
                xmin, ymin = round(x)-30, round(y)-30
            stamp = image[ymin:ymin+61, xmin:xmin+61]
            if stamp.shape != (61, 61) or not np.isfinite(stamp).all() or (stamp < 0).any():
                raise ValueError(f"Invalid local LSXPS {label} image")
            arrays[label] = np.array(stamp, dtype=float)
    xx, yy = np.meshgrid(np.arange(xmin, xmin+61), np.arange(ymin, ymin+61))
    sky = transforms["counts"].pixel_to_world(xx, yy)
    for label, other in transforms.items():
        xo, yo = other.world_to_pixel(sky)
        if not (np.allclose(xx, xo, atol=.01, rtol=0) and np.allclose(yy, yo, atol=.01, rtol=0)):
            raise ValueError(f"LSXPS background map is not on the science grid: {label}")
    scale = float(np.sqrt(abs(np.linalg.det(transforms["counts"].pixel_scale_matrix))) * 3600)
    icrs = SkyCoord(sky.icrs.ra, sky.icrs.dec, frame="icrs")
    east, north = point.spherical_offsets_to(icrs)
    return arrays, xx-x, yy-y, east.arcsec, north.arcsec, scale


def counts_in(arrays, region):
    counts = float(arrays["counts"][region].sum())
    if counts != int(counts):
        raise ValueError("LSXPS count map contains non-integral aperture counts")
    bg = float(arrays["background_with_catalogue_sources"][region].sum())
    return dict(
        count_pixels=int(np.count_nonzero(region)), observed_counts=int(counts),
        background_with_catalogue_sources_counts=bg,
        background_without_catalogue_sources_counts=float(arrays["background"][region].sum()),
        signed_excess_counts=counts-bg, poisson_tail_at_fixed_map_background=float(poisson.sf(counts-1, bg)),
    )


def plot(obs, arrays, east, north, neighbors, radius):
    fig = Figure(figsize=(9.8, 5.1), layout="constrained")
    axes = fig.subplots(1, 2)
    net = arrays["counts"] - arrays["background_with_catalogue_sources"]
    limit = float(np.percentile(np.abs(net), 99.5))
    for ax, title, pixels, cmap, norm in (
        (axes[0], "Observed hard-band counts", arrays["counts"], "Blues", None),
        (axes[1], "Counts minus catalogued-source background model", net, "RdBu_r",
         TwoSlopeNorm(vcenter=0, vmin=-limit, vmax=limit)),
    ):
        mesh = ax.pcolormesh(east, north, pixels, shading="nearest", cmap=cmap, norm=norm)
        ax.add_patch(Circle((0, 0), radius, fill=False, color="#ad552c", lw=1.35))
        ax.scatter(0, 0, marker="+", color="#ad552c", s=55, label="JWST position")
        ax.scatter(*np.asarray([row["offset_east_north_arcsec"] for row in neighbors]).T,
                   marker="x", color="#527043", s=33, label="Nearby Chandra sources")
        ax.set(xlim=(-66, 66), ylim=(-66, 66), aspect="equal", xlabel="East (arcsec)",
               ylabel="North (arcsec)", title=title)
        fig.colorbar(mesh, ax=ax, shrink=.72, label="Counts per pixel")
    axes[0].legend(loc="upper right", fontsize=8)
    fig.suptitle(f"LSXPS 2–10 keV stack {obs} at JADES 103607")
    FIGURES.mkdir(parents=True, exist_ok=True)
    save_figure(fig, FIGURES / TAG)


def main():
    source = OUT / f"{HISTORY}.json"
    archive = json.loads(source.read_text(encoding="utf-8"))
    flags = [flag for flag in archive["lsxps"]["detection_flags"] if flag["stack"] and flag["band"] == "Hard"]
    if len(flags) != 1:
        raise ValueError("Expected one flagged hard-band stack in the archived LSXPS search")
    flag = flags[0]
    obs = flag["obsid"]
    dataset = next(row for row in archive["lsxps"]["datasets"] if row["obsid"] == obs)
    source_point = SkyCoord(*archive["position_icrs_deg"], unit="deg", frame="icrs")
    neighbors = archive["nearest_chandra_sources"]
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    paths, receipts, query = archive_images(cache, obs)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*RADECSYS.*", category=FITSFixedWarning)
        arrays, dx, dy, east, north, scale = cutout(paths, source_point)
    fixed = counts_in(arrays, np.hypot(dx, dy) <= NATIVE_APERTURE_PIXELS)
    if fixed["observed_counts"] != flag["counts"] or not np.isclose(
        fixed["background_with_catalogue_sources_counts"], flag["background_counts"], atol=2e-4, rtol=0
    ):
        raise ValueError("Image aperture did not reproduce the LSXPS forced count measurement")
    local = {"jwst_position": counts_in(arrays, np.hypot(east, north) <= PAIR_ARCSEC)}
    for index, row in enumerate(neighbors, start=1):
        x, y = row["offset_east_north_arcsec"]
        local[f"chandra_neighbor_{index}"] = counts_in(arrays, np.hypot(east-x, north-y) <= PAIR_ARCSEC)
    plot(obs, arrays, east, north, neighbors, NATIVE_APERTURE_PIXELS * scale)
    print(f"Hard stack {obs}, aperture {fixed}, local 6 arcsec apertures {local}", flush=True)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__), input_result_sha256=text_digest(source),
        position_icrs_deg=archive["position_icrs_deg"], instruments=["Swift XRT", "Chandra catalogue"],
        obsid=obs, energy_kev=[2, 10],
        stack_time_span_utc=Time([dataset["start_mjd"], dataset["stop_mjd"]], format="mjd", scale="utc").isot.tolist(),
        locator_request=query, receipts=receipts, aperture_radius_pixels=NATIVE_APERTURE_PIXELS,
        aperture_radius_arcsec=NATIVE_APERTURE_PIXELS*scale, fixed_aperture=fixed,
        paired_aperture_radius_arcsec=PAIR_ARCSEC, local_apertures=local,
        scope="Direct reproduction of the pipeline aperture from its archived count and background maps. "
        "The conditional-on-model-background Poisson tail omits background-model uncertainty and survey trials. "
        "The small apertures sample blended flux in the Swift PSF. The coadd cannot date an excess.",
    ))


if __name__ == "__main__":
    main()

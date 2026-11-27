"""Check Swift coverage and the closest sampled epoch to the 2024 Euclid excess."""

import json
from urllib.parse import urlencode
from xml.etree import ElementTree

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from check_jades_103607_independent_2026 import END, HEASARC, SOURCE, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from scipy.ndimage import map_coordinates
from scipy.stats import chi2, poisson
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_swift_coverage"
SWIFT = "https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs"
EUCLID = "2024-08-06T09:27:39"
CONE = 0.5
NARROW_CENTER_ARCMIN = 8.5
XRT_SOURCE_ARCSEC = 20
XRT_BACKGROUND_ARCSEC = (32, 60)


def swift_log(cache, point, tag, window=None):
    query = (
        "SELECT obsid,name,ra,dec,start_time,stop_time,xrt_exposure,uvot_exposure,archive_date "
        "FROM swiftmastr WHERE CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{point.ra.deg:.12f},{point.dec.deg:.12f},{CONE}))=1"
    )
    if window is not None:
        lo, hi = Time(window, scale="utc").mjd
        query += f" AND start_time < {hi:.8f} AND stop_time >= {lo:.8f}"
    query += " ORDER BY start_time"
    url = f"{HEASARC}?{urlencode(dict(REQUEST='doQuery', FORMAT='votable', LANG='ADQL', QUERY=query))}"
    path = cache.fetch(
        url, f"metadata/{TAG}_{tag}.vot", survey="NASA HEASARC Swift",
        version="SWIFTMASTR positional pointing log", kind="observed_pointing_metadata",
    )
    status = [node.get("value") for node in ElementTree.parse(path).findall(".//{*}INFO")
              if node.get("name") == "QUERY_STATUS"]
    if status != ["OK"]:
        raise ValueError(f"Unexpected Swift archive status {status}")
    table = Table.read(path, format="votable", use_names_over_ids=True)
    rows = json.loads(table.to_pandas().to_json(orient="records", double_precision=12))
    if len(table):
        separations = point.separation(SkyCoord(table["ra"], table["dec"], unit="deg")).arcmin
        for row, separation in zip(rows, separations, strict=True):
            row["pointing_separation_arcmin"] = float(separation)
            row["start_utc"] = Time(row["start_time"], format="mjd", scale="utc").isot
    return rows, receipt(path, url)


def product(cache, observation, subdir, suffix):
    obs_id = observation["obsid"]
    month = observation["start_utc"][:7].replace("-", "_")
    name = f"sw{obs_id}{suffix}"
    url = f"{SWIFT}/{month}/{obs_id}/{subdir}/{name}"
    path = cache.fetch(
        url, f"cross_detector/swift_coverage/{name}", survey="NASA HEASARC Swift",
        version="Public pipeline FITS products", kind="observed_swift_fits",
    )
    return path, receipt(path, url)


def position_on_image(hdu, point):
    wcs = WCS(hdu.header, fix=False).celestial
    x, y = map(float, wcs.world_to_pixel(point))
    scale = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600)
    if not (0 <= round(x) < hdu.data.shape[1] and 0 <= round(y) < hdu.data.shape[0]):
        raise ValueError(f"Point does not fall within {hdu.name}")
    return x, y, scale


def measure_xrt(cache, observation, point):
    events, event_receipt = product(cache, observation, "xrt/event", "xpcw3po_cl.evt.gz")
    expmap, exposure_receipt = product(cache, observation, "xrt/products", "xpc_ex.img.gz")
    with fits.open(events) as event_hdus, fits.open(expmap) as exposure_hdus:
        event = event_hdus["EVENTS"]
        exposure = exposure_hdus[0]
        eh = event.header
        if (eh["OBS_ID"], eh["DATAMODE"], exposure.header["OBS_ID"], exposure.header["BUNIT"]) != (
            observation["obsid"], "PHOTON", observation["obsid"], "s"
        ):
            raise ValueError("Unexpected Swift XRT event or exposure map")
        pairs = (("CRPIX", "TCRPX"), ("CRVAL", "TCRVL"), ("CDELT", "TCDLT"))
        for map_key, event_key in pairs:
            if not np.allclose(
                [exposure.header[f"{map_key}{i}"] for i in (1, 2)],
                [eh[f"{event_key}{i}"] for i in (2, 3)], atol=1e-5, rtol=0,
            ):
                raise ValueError(f"Swift XRT event and exposure map WCS mismatch in {map_key}")
        x, y, scale = position_on_image(exposure, point)
        dx = (np.asarray(event.data["X"], dtype=float) - 1 - x) * scale
        dy = (np.asarray(event.data["Y"], dtype=float) - 1 - y) * scale
        er = np.hypot(dx, dy)
        pi = event.data["PI"]
        py, px = np.ogrid[:exposure.data.shape[0], :exposure.data.shape[1]]
        pr = np.hypot(px-x, py-y) * scale
        weights = np.asarray(exposure.data, dtype=float)
        core = pr <= XRT_SOURCE_ARCSEC
        ring = (pr >= XRT_BACKGROUND_ARCSEC[0]) & (pr <= XRT_BACKGROUND_ARCSEC[1])
        mean_core = float(weights[core].mean())
        annulus_sum = float(weights[ring].sum())
        if np.min(weights[core]) <= 0 or annulus_sum <= 0:
            raise ValueError("Swift position lacks full XRT aperture and usable annulus exposure")
        factor = float(weights[core].sum() / annulus_sum)
        results = []
        for low, high in ((0.3, 10), (0.5, 2), (2, 7)):
            band = (pi >= round(low*100)) & (pi < round(high*100))
            n = int(np.count_nonzero(band & (er <= XRT_SOURCE_ARCSEC)))
            bg = int(np.count_nonzero(band & (er >= XRT_BACKGROUND_ARCSEC[0]) & (er <= XRT_BACKGROUND_ARCSEC[1])))
            expectation = bg * factor
            upper = float(chi2.ppf(.95, 2*(n+1)) / 2)
            results.append(dict(
                energy_kev_lower_inclusive_upper_exclusive=[low, high], source_counts=n,
                annulus_counts=bg, estimated_aperture_background_counts=expectation,
                signed_excess_counts=n-expectation,
                poisson_tail_at_plugin_background=float(poisson.sf(n-1, expectation)),
                gross_95_percent_aperture_counts_upper=upper,
                gross_95_percent_aperture_rate_per_s_upper=upper/mean_core,
            ))
        band = (pi >= 30) & (pi < 1000) & (er <= 130)
        wcs = WCS(exposure.header, fix=False).celestial
        positions = wcs.pixel_to_world(event.data["X"][band]-1, event.data["Y"][band]-1).icrs
        east, north = point.spherical_offsets_to(SkyCoord(positions.ra, positions.dec, frame="icrs"))
        datum = dict(
            events=event_receipt, exposure_map=exposure_receipt, target_image_pixel_xy=[x, y],
            point_exposure_s=float(weights[round(y), round(x)]),
            mean_source_aperture_exposure_s=mean_core, aperture_radius_arcsec=XRT_SOURCE_ARCSEC,
            annulus_arcsec=XRT_BACKGROUND_ARCSEC, background_exposure_area_scale=factor,
            date_obs_utc=eh["DATE-OBS"], date_end_utc=eh["DATE-END"], bands=results,
        )
        return datum, (east.arcsec, north.arcsec)


def uvw2_image(cache, observation, point):
    exp, exp_receipt = product(cache, observation, "uvot/image", "uw2_ex.img.gz")
    sci, sci_receipt = product(cache, observation, "uvot/image", "uw2_sk.img.gz")
    with fits.open(exp) as exp_hdus, fits.open(sci) as sci_hdus:
        pairs = [(a, sci_hdus[a.name]) for a in exp_hdus if a.data is not None and a.data.ndim == 2]
        if len(pairs) != 1:
            raise ValueError("Expected exactly one public Swift UVW2 exposure in this visit")
        exposure, science = pairs[0]
        x, y, scale = position_on_image(exposure, point)
        positions = position_on_image(science, point)
        if (exposure.header["FILTER"], science.header["FILTER"], exposure.header["OBS_ID"]) != (
            "UVW2", "UVW2", observation["obsid"]
        ) or exposure.data.shape != science.data.shape or not np.allclose((x, y, scale), positions, atol=1e-4):
            raise ValueError("Unmatched Swift UVW2 science and exposure image")
        axis = np.arange(-42, 43) / 3600 * point.ra.unit
        east, north = np.meshgrid(axis, axis)
        grid = point.spherical_offsets_by(east, north)
        cx, cy = WCS(science.header, fix=False).celestial.world_to_pixel(grid)
        valid = np.where(exposure.data > 0, science.data, np.nan)
        cutout = map_coordinates(valid.astype(float), (cy, cx), order=0, cval=np.nan)
        extent = (-42.5, 42.5, -42.5, 42.5)
        coverage = dict(
            science_image=sci_receipt, exposure_map=exp_receipt, target_pixel_xy=[x, y],
            point_exposure_s=float(exposure.data[round(y), round(x)]),
            date_obs_utc=exposure.header["DATE-OBS"], date_end_utc=exposure.header["DATE-END"],
            scope="Observed UVW2 image and exposure coverage. Photometric sensitivity remains unmeasured.",
        )
        return coverage, cutout, extent


def plot(observation, xrt, offsets, uv, extent):
    fig = Figure(figsize=(8.7, 4.6), layout="constrained")
    xray, uvot = fig.subplots(1, 2)
    xray.scatter(*offsets, s=12, color="#215877", alpha=.75)
    xray.add_patch(Circle((0, 0), XRT_SOURCE_ARCSEC, color="#92532b", fill=False))
    for radius in XRT_BACKGROUND_ARCSEC:
        xray.add_patch(Circle((0, 0), radius, color="#838383", fill=False, linestyle=":"))
    counts = xrt["bands"][0]["source_counts"]
    exposure = xrt["mean_source_aperture_exposure_s"]
    xray.set(aspect="equal", xlim=(-125, 125), ylim=(-125, 125),
             title=f"XRT events in 0.3–10 keV\n{counts} aperture counts in {exposure:.0f} s",
             xlabel="East (arcsec)", ylabel="North (arcsec)")
    im = uvot.imshow(uv, extent=extent, origin="lower", cmap="bone",
                    vmin=0, vmax=max(1, np.nanpercentile(uv, 99)))
    uvot.scatter(0, 0, marker="+", color="#d17943", s=58)
    uvot.set(title="UVOT UVW2 science image", xlabel="East (arcsec)", ylabel="North (arcsec)")
    fig.colorbar(im, ax=uvot, label="Observed image counts", shrink=.78)
    fig.suptitle(f"Swift {observation['obsid']} | {observation['start_utc'][:10]} | JWST transient position")
    FIGURES.mkdir(parents=True, exist_ok=True)
    save_figure(fig, FIGURES / TAG)


def main():
    coords = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    point = SkyCoord(*coords, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    history, history_source = swift_log(cache, point, "history")
    recent, recent_source = swift_log(cache, point, "2026", (START, END))
    narrow = [row for row in history if row["pointing_separation_arcmin"] <= NARROW_CENTER_ARCMIN]
    epoch = Time(EUCLID, scale="utc")
    closest = min(narrow, key=lambda row: abs(row["start_time"] - epoch.mjd))
    xrt, offsets = measure_xrt(cache, closest, point)
    uvot, image, extent = uvw2_image(cache, closest, point)
    plot(closest, xrt, offsets, image, extent)
    hours_before_euclid = float((epoch - Time(xrt["date_obs_utc"], scale="utc")).to_value("hour"))
    print(f"{len(history)} Swift pointings, {len(narrow)} within 8.5 arcmin, {len(recent)} in 2026 window. "
          f"Closest usable XRT snapshot to Euclid starts {hours_before_euclid:.2f} hours earlier.", flush=True)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        position_source_sha256=text_digest(SOURCE), position_icrs_deg=coords,
        archive_audit_cone_deg=CONE, history_source=history_source, recent_source=recent_source,
        history_rows=history, recent_window_utc=[START, END], recent_rows=recent,
        narrow_center_radius_arcmin=NARROW_CENTER_ARCMIN, narrow_center_count=len(narrow),
        nearest_to_euclid_from_narrow_centers=closest, euclid_reference_utc=EUCLID,
        xrt_hours_before_euclid=hours_before_euclid, xrt=xrt, uvot=uvot,
        scope="Pointing cones describe archive metadata. Image and event coverage are verified only "
        "for the sampled August 2024 visit. XRT count-rate bounds retain gross aperture counts "
        "and use the aperture-mean exposure. No PSF, energy-response or spectral flux conversion "
        "is applied. The count-rate limits constrain the Swift snapshot time only.",
    ))


if __name__ == "__main__":
    main()

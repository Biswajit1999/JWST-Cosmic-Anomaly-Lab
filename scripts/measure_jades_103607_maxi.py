"""Compare public MAXI scan photometry with individual JWST integrations."""

import gzip
import json
import re
import tarfile
from io import BytesIO

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from check_jades_103607_einstein_probe import POSITION, STORAGE_MIB
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from measure_jwst_epoch_candidate import verified_path
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

TAG = "module_b_001_maxi_2026"
ARCHIVE = "https://maxi.riken.jp/mxondem/api/results"
JOBS = {
    "October": ("20261127092354_q0acckfu2rm08l1r", 61322, 61324, "jw12577001001"),
    "November": ("20261127092373_g90xwlullvjbdll2", 61368, 61370, "jw12577004001"),
}
BANDS = ("2.0-6.0", "6.0-20.0")
COLORS = ("#176b7a", "#a84b43")
SOURCE_RADIUS = 1.5
BACKGROUND_INNER_RADIUS = 1.6
BACKGROUND_RADIUS = 3.0
SECONDS_PER_DAY = 86400


def mjd_utc(values):
    return Time(values, format="mjd", scale="utc").isot


def pin(cache, label, work_id, filename):
    url = f"{ARCHIVE}/{work_id}/{filename}"
    path = cache.fetch(url, f"maxi_2026/{label}/{filename}",
                       survey="MAXI GSC", version="MAXI on-demand 3.95", kind="public_extraction")
    return path, receipt(path, url)


def integrations(position_result, visit):
    files = {
        product["file"]: product
        for band in position_result["bands"].values()
        for tier in band["inventory"] if tier["visit"] == visit
        for product in tier["products"]
    }
    windows = {}
    for filename, product in sorted(files.items()):
        header = fits.getheader(verified_path(product))
        lo, hi = header["EXPSTART"], header["EXPEND"]
        if header["TIMESYS"] != "UTC" or not 0 < hi-lo < 0.1:
            raise ValueError(f"Unexpected JWST exposure time in {filename}")
        # The two NIRCam channels give the same interval for simultaneous integrations.
        windows.setdefault((lo, hi), []).append(filename)
    if not windows:
        raise ValueError(f"No native integrations for {visit}")
    merged = []
    for (lo, hi), names in sorted(windows.items()):
        products = [files[name] for name in names]
        if merged and lo < merged[-1]["stop_mjd"]:
            merged[-1]["stop_mjd"] = max(hi, merged[-1]["stop_mjd"])
            merged[-1]["products"].extend(products)
        else:
            merged.append(dict(start_mjd=lo, stop_mjd=hi, products=products))
    for row in merged:
        row["utc"] = mjd_utc([row["start_mjd"], row["stop_mjd"]]).tolist()
    return merged


def tar_fits(tar, member):
    with tar.extractfile(member) as stream:
        return fits.open(BytesIO(gzip.decompress(stream.read())))


def check_times(primary, published, native_rows, start_mjd, stop_mjd):
    if primary["TIMESYS"] != "TT" or primary["TIMEUNIT"] != "s":
        raise ValueError("Unexpected MAXI raw time definition")
    ref = Time(primary["MJDREFI"], primary["MJDREFF"], format="mjd", scale="tt")
    for key, query_mjd in (("TSTART", start_mjd), ("TSTOP", stop_mjd)):
        calculated = (ref + TimeDelta(primary[key], format="sec")).utc.mjd
        if abs(calculated-query_mjd)*SECONDS_PER_DAY > 1:
            raise ValueError("MAXI header time and query interval disagree")
    for key in ("START", "STOP"):
        raw_mjd = (ref + TimeDelta(native_rows[key], format="sec")).utc.mjd
        # A published scan may combine detectors and use the extrema of their times.
        seconds_difference = np.min(np.abs(published[key][:, None]-raw_mjd[None, :]), axis=1)
        if np.max(seconds_difference)*SECONDS_PER_DAY > 1:
            raise ValueError("MAXI light-curve times do not reproduce raw detector times")


def statistics(rate, error):
    weights = 1/error**2
    mean = np.sum(rate*weights)/weights.sum()
    uncertainty = 1/np.sqrt(weights.sum())
    return dict(
        inverse_variance_mean=float(mean), statistical_error=float(uncertainty),
        nominal_constant_rate_upper_3sigma=float(max(0, mean)+3*uncertainty),
        maximum_single_scan_z=float(np.max(rate/error)),
    )


def light_curves(path, tag, start_mjd, stop_mjd, windows):
    with tarfile.open(path, mode="r:gz") as tar, tar_fits(tar, tag+"_gsc.flc.gz") as hdus:
        primary = hdus[0].header
        raw = hdus["LCRAWDATA"].data
        intervals = hdus["BINLC0"].data
        check_times(primary, intervals, raw, start_mjd, stop_mjd)
        lo = np.asarray(intervals["START"])
        hi = np.asarray(intervals["STOP"])
        if not np.all((lo >= start_mjd) & (hi <= stop_mjd) & (hi > lo)):
            raise ValueError("A MAXI scan falls outside its provenance interval")
        overlaps = [
            [j for j, w in enumerate(windows) if min(b, w["stop_mjd"]) > max(a, w["start_mjd"])]
            for a, b in zip(lo, hi, strict=True)
        ]
        scans = [
            dict(start_mjd=float(a), stop_mjd=float(b), scan_utc=mjd_utc([a, b]).tolist(),
                 scan_wall_span_s=float((b-a)*SECONDS_PER_DAY), overlapping_jwst_windows=overlap,
                 jwst_overlap_s=float(sum(
                     min(b, windows[j]["stop_mjd"])-max(a, windows[j]["start_mjd"]) for j in overlap
                 )*SECONDS_PER_DAY), bands={})
            for a, b, overlap in zip(lo, hi, overlaps, strict=True)
        ]
        paired = np.asarray([bool(item) for item in overlaps])
        band_stats = {}
        for i, band in enumerate(BANDS):
            data = hdus[f"BINLC{i}"].data
            if not (np.array_equal(lo, data["START"]) and np.array_equal(hi, data["STOP"])):
                raise ValueError("MAXI scan intervals differ across the two bands")
            rate, error = data["RATE"], data["ERROR"]
            if not (np.isfinite(rate).all() and np.isfinite(error).all() and (error > 0).all()):
                raise ValueError("Expected finite MAXI scan rates and positive supplied errors")
            dat = np.loadtxt(BytesIO(tar.extractfile(f"{tag}_{band}keV_gsclc.dat").read()), ndmin=2)
            if not np.allclose(dat, np.column_stack((lo, hi, rate, error)), rtol=0, atol=5.1e-7):
                raise ValueError("MAXI FITS and ASCII photometry disagree")
            for row, net, sigma in zip(scans, rate, error, strict=True):
                row["bands"][band] = dict(
                    signed_rate=float(net), statistical_error=float(sigma),
                    nominal_upper_3sigma=float(max(0, net)+3*sigma),
                )
            band_stats[band] = dict(
                all_sampled_scans=statistics(rate, error),
                overlapping_scans=statistics(rate[paired], error[paired]) if paired.any() else None,
            )
        return scans, band_stats, len(raw), primary["MXONDEM"]


def images(path, tag, point):
    products = {}
    with tarfile.open(path, mode="r:gz") as tar:
        for band in BANDS:
            with tar_fits(tar, f"{tag}_{band}keV_image.img.gz") as hdus:
                hdu = hdus["IMG"]
                frame = WCS(hdu.header).celestial
                x, y = frame.world_to_pixel(point)
                if not 0 <= x < hdu.data.shape[1] or not 0 <= y < hdu.data.shape[0]:
                    raise ValueError("MAXI target falls outside archived count map")
                products[band] = (np.array(hdu.data), frame)
        regions = {}
        for kind, expected_radii in [
            ("src", [SOURCE_RADIUS]), ("bgd", [BACKGROUND_RADIUS, BACKGROUND_INNER_RADIUS]),
        ]:
            name = f"{tag}_{kind}.reg"
            with tar.extractfile(name) as stream:
                regions[kind] = stream.read().decode()
            entries = re.findall(r'circle\(([\d.]+),(-?[\d.]+),([\d.]+)"\)', regions[kind])
            if len(entries) != len(expected_radii):
                raise ValueError(f"Unexpected aperture definition in {name}")
            for (ra, dec, radius), expected in zip(entries, expected_radii, strict=True):
                center = SkyCoord(float(ra), float(dec), frame="fk5", unit="deg")
                if abs(float(radius)/3600-expected) > 1e-7 or point.separation(center).arcsec > 3:
                    raise ValueError(f"Unexpected aperture geometry in {name}")
    return products, regions


def plot_image(fig, ax, grid, point, band):
    image, frame = grid
    mappable = ax.imshow(image, origin="lower", cmap="bone", interpolation="nearest")
    x, y = frame.world_to_pixel(point)
    degrees_per_pixel = frame.proj_plane_pixel_scales()[0].to_value("deg")
    for radius, style in [
        (SOURCE_RADIUS, "-"), (BACKGROUND_INNER_RADIUS, "--"), (BACKGROUND_RADIUS, "--"),
    ]:
        ax.add_patch(Circle((x, y), radius/degrees_per_pixel, facecolor="none",
                            edgecolor="#ffb85c", linestyle=style, lw=.8))
    ax.plot(x, y, "+", color="#ffb85c", ms=5)
    ax.set_title(f"{band} keV counts", fontsize=10)
    ax.coords[0].set_axislabel("Right ascension (J2000)", fontsize=8)
    ax.coords[1].set_axislabel("Declination (J2000)", fontsize=8)
    for coord in ax.coords:
        coord.set_ticks(number=3)
        coord.set_ticklabel(size=7)
    bar = fig.colorbar(mappable, ax=ax, shrink=.7, ticks=range(int(image.max())+1))
    bar.ax.tick_params(labelsize=7)


def render(point, results, maps):
    fig = Figure(figsize=(13.7, 7.05))
    grid = fig.add_gridspec(2, 3, width_ratios=(1.05, 1.05, 1.65), left=.075, right=.978,
                           bottom=.105, top=.855, wspace=.56, hspace=.46)
    for i, (label, row) in enumerate(results.items()):
        for j, band in enumerate(BANDS):
            picture = maps[label][band]
            plot_image(fig, fig.add_subplot(grid[i, j], projection=picture[1]), picture,
                       point, band)
        ax = fig.add_subplot(grid[i, 2])
        origin = row["queried_mjd"][0]
        for j, exposure in enumerate(row["jwst_integrations"]):
            ax.axvspan(24*(exposure["start_mjd"]-origin), 24*(exposure["stop_mjd"]-origin),
                       color="#e4b25f", alpha=.3, lw=0, label="JWST exposures" if j == 0 else None)
        for k, band in enumerate(BANDS):
            rows = row["scans"]
            centers = [24*((s["start_mjd"]+s["stop_mjd"])/2-origin) for s in rows]
            ax.errorbar(centers, [s["bands"][band]["signed_rate"] for s in rows],
                        yerr=[s["bands"][band]["statistical_error"] for s in rows],
                        xerr=[12*(s["stop_mjd"]-s["start_mjd"]) for s in rows],
                        fmt=("o" if k == 0 else "s"), ms=3, capsize=1.8, color=COLORS[k],
                        label=f"{band} keV", alpha=.91)
        ax.axhline(0, color=".5", linewidth=.6)
        ax.set(xlim=(0, 48), xlabel=f"Hours since {mjd_utc(origin)[:10]} 00:00 UTC",
               ylabel="Photon rate (cm⁻² s⁻¹)", title=f"{label} 2026 · {len(row['scans'])} scans")
        ax.yaxis.label.set_size(9)
        ax.xaxis.label.set_size(9)
        if i == 0:
            ax.legend(loc="upper left", fontsize=7.5, frameon=False, ncol=2)
    fig.suptitle("MAXI/GSC observations during the JADES 103607 JWST visits", fontsize=14, y=.967)
    save_figure(fig, FIGURES / TAG)


def main():
    source = json.loads(POSITION.read_text(encoding="utf-8"))
    point = SkyCoord(*source["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB * 1024**2)
    results, maps = {}, {}
    for label, (work_id, start, stop, visit) in JOBS.items():
        name = f"JADES103607_2026{label[:3]}"
        index, page = pin(cache, label, work_id, "index.html")
        html = index.read_text(encoding="latin-1")
        if not all(
            re.search(r'href\s*=\s*["\']?'+re.escape(f"{name}_{kind}.tgz"), html, re.I)
            for kind in ("lc", "img")
        ):
            raise ValueError(f"MAXI {label} page does not list both required data products")
        lc, light = pin(cache, label, work_id, name+"_lc.tgz")
        img, counts = pin(cache, label, work_id, name+"_img.tgz")
        windows = integrations(source, visit)
        scans, bands, native_rows, version = light_curves(lc, name, start, stop, windows)
        maps[label], ds9_regions = images(img, name, point)
        results[label] = dict(
            queried_mjd=[start, stop], work_id=work_id, on_demand_version=version,
            regions_ds9=ds9_regions, products=dict(result_page=page, light_curve=light, images=counts),
            jwst_visit=visit, jwst_integrations=windows, raw_detector_scan_rows=native_rows,
            fractional_day_scan_bin_request=0.00005, scans=scans, bands=bands,
        )
    render(point, results, maps)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), icrs_deg=source["icrs_deg"], position_sha256=text_digest(POSITION),
        script_sha256=text_digest(__file__), archive_query_icrs_deg=[53.18866983, -27.81253940],
        source_radius_deg=SOURCE_RADIUS, background_radii_deg=[BACKGROUND_INNER_RADIUS, BACKGROUND_RADIUS],
        rate_unit="ph cm^-2 s^-1", windows=results,
        limitations=[
            "The scan interval is wall time. Effective detector area and usable exposure vary "
            "throughout the scan.",
            "The 3-sigma levels use the supplied statistical uncertainties and a zero floor on "
            "the point estimate. They are nominal Gaussian rate bounds, without a calibration "
            "or systematic-background allowance.",
            "Combined rates and upper levels assume constant intensity over the sampled scans.",
            "The source aperture blends positions across 1.5 degrees of radius. A counterpart "
            "fainter than the limits, or active between scans, remains unconstrained.",
        ],
    ))


if __name__ == "__main__":
    main()

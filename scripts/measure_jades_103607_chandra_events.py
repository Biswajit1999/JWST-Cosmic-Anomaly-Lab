"""Search calibrated Chandra photon times in discrete GOODS-S observations."""

import argparse
import json
import re
import tarfile
from warnings import catch_warnings, simplefilter

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS, FITSFixedWarning
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from matplotlib.path import Path as Polygon
from scipy.stats import binom, chi2
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

CXC = "https://cxc.cfa.harvard.edu/cdaftp/byobsid"
SOURCE = OUT / "module_b_001_native_november_fade.json"
JOURNAL = ROOT / "data/raw/cross_detector/luo2017-cdfs-data-products03.tar.gz"
EARLIER = OUT / "module_b_001_chandra_spitzer_history.json"
RADIUS, BG = 1.5, (8, 20)
EDGE_GUARD = 8
WINDOW = 10000
MINIMUM_GOOD_TIME = 5000
ENERGY_EV = (500, 6000)


def search_observations(year):
    prior = json.loads(EARLIER.read_text())["chandra"]["product"]
    if digest(JOURNAL) != prior["sha256"]:
        raise ValueError("Survey observation journal archive checksum changed")
    with tarfile.open(JOURNAL, "r:gz") as bundle:
        text = bundle.extractfile("./tbl01.txt").read().decode("ascii")
    matches = re.findall(
        r"^\s*(\d+)\s+(\d{4} \w{3} \d{2}), (\d{2}:\d{2})\s+([\d.]+)", text, re.MULTILINE
    )
    return [
        dict(obsid=int(number), journal_start=f"{date} {start}", survey_clean_exposure_ks=float(ks))
        for number, date, start, ks in matches
        if date[:4] == str(year)
    ], prior


def download_paths(cache, obsid):
    folder = f"{CXC}/{obsid % 10}/{obsid}/primary/"
    key = f"chandra_evt_{obsid}"
    listing = cache.fetch(
        folder,
        f"metadata/{key}_listing.html",
        survey="Chandra Data Archive",
        version="Public primary product listing",
        kind="archive_metadata",
    )
    html = listing.read_text(encoding="utf-8", errors="replace")
    names = {}
    for product in ("fov1", "evt2"):
        pattern = rf"acisf{obsid:05d}(?:_\d+)?N\d+_{product}\.fits\.gz"
        found = sorted(set(re.findall(pattern, html)))
        if len(found) != 1:
            raise ValueError(f"Expected one {product} for observation {obsid}, got {found}")
        names[product] = folder + found[0]
    return names, receipt(listing, folder)


def fetch_product(cache, url, product):
    path = cache.fetch(
        url,
        f"cross_detector/{url.rsplit('/', 1)[-1]}",
        survey="Chandra Data Archive",
        version="Public calibrated ACIS",
        kind=product,
    )
    return path, receipt(path, url)


def sky_pixel(h, columns, sky):
    header = h.copy()
    for col in columns:
        key = f"TCUNI{col}"
        if header.get(key) in {"degree", "degrees"}:
            header[key] = "deg"
    with catch_warnings():
        simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header, keysel=["pixel"], colsel=columns, fix=False)
    coord = np.array(wcs.wcs_world2pix(sky.ra.deg, sky.dec.deg, 1), dtype=float)
    scale = np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600
    if not np.isclose(scale, 0.492, rtol=0.001) or not np.all(np.isfinite(coord)):
        raise ValueError("Unexpected ACIS sky-pixel solution")
    return coord, scale


def signed_margin(vertices, center, scale):
    corners = vertices - center
    ends = np.roll(corners, -1, axis=0) - corners
    lengths = np.sum(ends**2, axis=1)
    keep = lengths > 0
    frac = np.clip(-np.sum(corners[keep] * ends[keep], axis=1) / lengths[keep], 0, 1)
    gap = float(np.min(np.linalg.norm(corners[keep] + ends[keep] * frac[:, None], axis=1)) * scale)
    return gap if Polygon(vertices).contains_point(center) else -gap


def footprint(path, sky, obsid):
    with fits.open(path, memmap=False) as info:
        h = info["FOV"]
        if int(h.header["OBS_ID"]) != obsid:
            raise ValueError("ACIS footprint ObsID mismatch")
        pos, scale = sky_pixel(h.header, [1, 2], sky)
        ccds = {}
        for row in h.data:
            if str(row["SHAPE"]).casefold() != "polygon" or int(row["CCD_ID"]) not in range(4):
                continue
            xy = np.column_stack((row["X"], row["Y"]))
            xy = xy[np.isfinite(xy).all(axis=1) & np.any(xy != 0, axis=1)]
            ccds[int(row["CCD_ID"])] = signed_margin(xy, pos, scale)
    if not ccds:
        raise ValueError("No ACIS-I chip polygons in footprint")
    chip = max(ccds, key=ccds.get)
    return chip, ccds[chip], pos


def statistic(target, background, good_seconds):
    n = len(target)
    b = len(background)
    ratio = RADIUS**2 / (BG[1] ** 2 - BG[0] ** 2)
    return dict(
        source_counts=n,
        annulus_counts=b,
        background_expected=b * ratio,
        signed_excess=n - b * ratio,
        binomial_background_tail=float(binom.sf(n - 1, n + b, ratio / (1 + ratio))),
        upper95_aperture_counts=float(chi2.ppf(0.95, 2 * (n + 1)) / 2),
        good_time_s=good_seconds,
    )


def windows(path, sky, chip, obsid, expected_pos):
    with fits.open(path, memmap=False) as h:
        events = h["EVENTS"]
        hdr = events.header
        if (
            int(hdr["OBS_ID"]) != obsid
            or hdr["TIMEUNIT"].strip() != "s"
            or hdr["TIMESYS"].strip() != "TT"
        ):
            raise ValueError("Unexpected ACIS event identity or time coordinates")
        xcol = events.columns.names.index("x") + 1
        ycol = events.columns.names.index("y") + 1
        pos, scale = sky_pixel(hdr, [xcol, ycol], sky)
        if not np.allclose(pos, expected_pos, rtol=0, atol=0.05):
            raise ValueError("Event and FOV sky-pixel solutions disagree")
        (gti,) = [
            item for item in h[2:] if item.name == "GTI" and item.header.get("CCD_ID") == chip
        ]
        starts = np.asarray(gti.data["START"])
        ends = np.asarray(gti.data["STOP"])
        if np.any(ends < starts) or np.any(starts[1:] < ends[:-1]):
            raise ValueError("Bad ordering of ACIS good time intervals")
        r = np.hypot(events.data["x"] - pos[0], events.data["y"] - pos[1]) * scale
        use = (
            (events.data["ccd_id"] == chip)
            & (events.data["energy"] >= ENERGY_EV[0])
            & (events.data["energy"] < ENERGY_EV[1])
        )
        time = np.asarray(events.data["time"])
        target = time[use & (r <= RADIUS)]
        bg = time[use & (r >= BG[0]) & (r <= BG[1])]
        mjdref = Time(hdr["MJDREF"], format="mjd", scale="tt")

        def utc(seconds):
            return (mjdref + TimeDelta(float(seconds), format="sec")).utc.isot

        cuts = []
        empty = []
        for lo in np.arange(starts[0], ends[-1], WINDOW):
            hi = min(lo + WINDOW, ends[-1])
            exposure = float(np.maximum(0, np.minimum(ends, hi) - np.maximum(starts, lo)).sum())
            if exposure < MINIMUM_GOOD_TIME:
                empty.append(dict(start_utc=utc(lo), end_utc=utc(hi), good_time_s=exposure))
                continue
            a = target[(target >= lo) & (target < hi)]
            b = bg[(bg >= lo) & (bg < hi)]
            cuts.append(dict(start_utc=utc(lo), end_utc=utc(hi), **statistic(a, b, exposure)))
        full = statistic(target, bg, float(np.sum(ends - starts)))
        return dict(
            start_utc=utc(starts[0]),
            end_utc=utc(ends[-1]),
            observation=full,
            windows=cuts,
            excluded_short_or_gap_windows=empty,
            source_events_utc=[utc(t) for t in target],
        )


def analyze(cache, row, sky):
    obsid = row["obsid"]
    urls, listing = download_paths(cache, obsid)
    fov, meta = fetch_product(cache, urls["fov1"], "observed_detector_footprint")
    chip, margin, coords = footprint(fov, sky, obsid)
    row.update(
        primary_archive_listing=listing,
        footprint_product=meta,
        chip=chip,
        edge_margin_arcsec=margin,
    )
    if margin < BG[1] + EDGE_GUARD:
        row["status"] = "insufficient distance from detector edge"
        return
    event, event_meta = fetch_product(cache, urls["evt2"], "observed_calibrated_photon_events")
    row.update(
        event_product=event_meta, status="measured", **windows(event, sky, chip, obsid, coords)
    )


def plot(rows, year, path):
    fig = Figure(figsize=(10, 4.75), layout="constrained")
    axes = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": (1.65, 1)})
    for i, row in enumerate(rows):
        d = row["observation"]
        axes[0].scatter(
            i - 0.08, d["source_counts"], color="#2f6c83", label="Observed" if not i else None
        )
        axes[0].scatter(
            i + 0.08,
            d["background_expected"],
            facecolors="none",
            edgecolors="#a8552a",
            label="Local background" if not i else None,
        )
        axes[0].plot(
            [i - 0.24, i + 0.24],
            [d["upper95_aperture_counts"]] * 2,
            color=".38",
            lw=1,
            ls="--",
            label="95% upper bound" if not i else None,
        )
        for j, w in enumerate(row["windows"]):
            x = i + (0.48 * (j + 0.5) / max(1, len(row["windows"])) - 0.24)
            axes[1].scatter(
                x,
                w["source_counts"],
                c=w["good_time_s"] * 0.001,
                cmap="viridis",
                vmin=5,
                vmax=10,
                s=23,
            )
    axes[0].set(ylabel="Full observation counts")
    axes[0].legend(ncols=3, frameon=False, fontsize=8)
    if len(rows) <= 17:
        axes[1].set_xticks(range(len(rows)), [r["obsid"] for r in rows], rotation=35)
    axes[1].set(
        ylabel="10 ks counts",
        xlabel="Chandra observation ID"
        if len(rows) <= 17
        else "Chandra observations in time order",
    )
    fig.suptitle(f"JADES 103607 · {year} Chandra photons in a 1.5″ aperture", fontsize=13)
    save_figure(fig, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--year", type=int, choices=(2000, 2007, 2010, 2014, 2015, 2016), default=2007
    )
    year = parser.parse_args().year
    source = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*source, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    rows, prior = search_observations(year)
    if not rows:
        raise ValueError(f"No {year} pointings in the CDF-S input journal")
    for row in rows:
        analyze(cache, row, sky)
        print(
            row["obsid"], row["status"], row.get("observation", {}).get("source_counts"), flush=True
        )
    good = [r for r in rows if r["status"] == "measured"]
    matches = []
    for visit in good:
        matches.extend(dict(obsid=visit["obsid"], kind="10 ks", **w) for w in visit["windows"])
        matches.append(dict(obsid=visit["obsid"], kind="whole pointing", **visit["observation"]))
    if not matches:
        raise ValueError("No useful exposure at the target")
    strongest = min(matches, key=lambda w: w["binomial_background_tail"])
    tag = f"module_b_001_chandra_{year}_events"
    plot(good, year, FIGURES / tag)
    write_json(
        OUT / f"{tag}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=source,
            position_sha256=text_digest(SOURCE),
            input_observing_journal=prior,
            campaign_year=year,
            aperture_arcsec=RADIUS,
            background_annulus_arcsec=BG,
            energy_ev=ENERGY_EV,
            duration_s=WINDOW,
            minimum_gti_s=MINIMUM_GOOD_TIME,
            chip_edge_guard_arcsec=EDGE_GUARD,
            observations=rows,
            number_of_tests=len(matches),
            best_trial=strongest,
            best_bonferroni_p=min(1.0, len(matches) * strongest["binomial_background_tail"]),
            scope="Fixed unbarycentred windows in ACIS GTIs and a spatial background from the same chip. "
            "Count bounds retain background. A time-varying PSF, detector response and full effective "
            "area model are needed for flux limits. Other years and gaps are outside a campaign result.",
        ),
    )


if __name__ == "__main__":
    main()

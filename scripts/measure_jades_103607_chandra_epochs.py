"""Measure X-ray counts at JADES 103607 in independent archival Chandra epochs."""

import json
import tarfile
from io import BytesIO
from pathlib import PurePosixPath

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from read_jwst_i2d_cutout import MastImageSections
from scipy.ndimage import map_coordinates
from scipy.stats import binom, chi2
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

BASE = "https://cxc.cfa.harvard.edu/cdaftp/arcftp"
SOURCE = OUT / "module_b_001_native_november_fade.json"
STEM = "module_b_001_chandra_epochs"
FILES = {
    2000: ("CDFS2000merge.tar", "2000"),
    2007: ("CDFS2007merge.tar", "2007"),
    2010: ("CDFS2010N002merge.tar", "2010N002"),
}
APERTURES = (1.5, 2.0, 3.0)
BG = (8, 20)
HALF = 44


class CxcSections(MastImageSections):
    """Read pinned byte ranges from the uncompressed CXC archive bundles."""

    def __init__(self, cache, url):
        self.cache = cache
        self.url = url
        self.filename = PurePosixPath(url).name
        with cache.session.head(url, timeout=(15, 40)) as response:
            response.raise_for_status()
            self.size = int(response.headers["Content-Length"])
            self.etag = response.headers["ETag"]
        self.ranges = []

    def members(self, expected):
        cursor = 0
        result = {}
        while cursor + 512 <= self.size and len(result) < len(expected):
            raw = self.read(cursor, cursor + 512).read_bytes()
            if raw == bytes(512):
                break
            member = tarfile.TarInfo.frombuf(raw, "utf-8", "strict")
            start = cursor + 512
            if member.name in expected and member.isfile():
                result[member.name] = (start, member.size)
            cursor = start + ((member.size + 511) // 512) * 512
        if set(result) != set(expected):
            raise ValueError(f"Missing requested FITS members in {self.filename}")
        return result

    def image_stamp(self, entry, sky):
        start, size = entry
        with BytesIO(self.read(start, start + min(172800, size)).read_bytes()) as source:
            h = fits.Header.fromfile(source)
            loc = start + source.tell()
        if h["NAXIS"] != 2 or h.get("BSCALE", 1) != 1 or h.get("BZERO", 0) != 0:
            raise ValueError("Unrecognised encoding in the CXC archive")
        ncols, nrows = h["NAXIS1"], h["NAXIS2"]
        dtype = {16: ">i2", 32: ">i4", -32: ">f4", -64: ">f8"}[h["BITPIX"]]
        bpr = ncols * np.dtype(dtype).itemsize
        wcs = WCS(h, fix=False).celestial
        x, y = map(float, wcs.world_to_pixel(sky))
        ix, iy = round(x), round(y)
        if not (HALF <= ix < ncols - HALF and HALF <= iy < nrows - HALF):
            raise ValueError("Target is outside the supported CXC image stamp")
        beginning = loc + bpr * (iy - HALF)
        end = loc + bpr * (iy + HALF + 1)
        if end > start + size:
            raise ValueError("CXC FITS stripe overruns the archive member")
        block = self.read(beginning, end)
        mapped = np.memmap(block, dtype=dtype, mode="r", shape=(2 * HALF + 1, ncols))
        return (
            np.array(mapped[:, ix - HALF : ix + HALF + 1], dtype=float),
            (x - ix + HALF, y - iy + HALF),
            h,
        )


def measurement(cache, year, sky):
    archive, suffix = FILES[year]
    source = CxcSections(cache, f"{BASE}/cdfs{year}/{archive}")
    names = [f"acismCDFS_{suffix}_{kind}.fits" for kind in ("full_img2", "expmap")]
    members = source.members(names)
    counts, (x, y), hdr = source.image_stamp(members[names[0]], sky)
    exposure, (ex, ey), exp_hdr = source.image_stamp(members[names[1]], sky)
    scale = abs(float(hdr["CDELT2"])) * 3600
    if not (
        np.allclose(
            (hdr["CRVAL1"], hdr["CRVAL2"]),
            (exp_hdr["CRVAL1"], exp_hdr["CRVAL2"]),
            atol=1e-8,
            rtol=0,
        )
        and np.allclose(
            (hdr["CDELT1"], hdr["CDELT2"]),
            (exp_hdr["CDELT1"], exp_hdr["CDELT2"]),
            atol=1e-10,
            rtol=0,
        )
        and max(abs(x - ex), abs(y - ey)) <= 1
        and np.all(counts >= 0)
        and np.all(counts == np.rint(counts))
    ):
        raise ValueError("Counts and exposure products are not aligned native maps")
    yy, xx = np.indices(counts.shape)
    exposure = map_coordinates(exposure, [yy + ey - y, xx + ex - x], order=1, prefilter=False)
    r = np.hypot(xx - x, yy - y) * scale
    annulus = (r >= BG[0]) & (r <= BG[1]) & (exposure > 0)
    bg_counts = int(counts[annulus].sum())
    bg_pixels = int(annulus.sum())
    fits_time = float(exp_hdr["ONTIME"])
    apertures = {}
    for radius in APERTURES:
        use = r <= radius
        if not np.all(exposure[use] > 0):
            raise ValueError("The claimed X-ray aperture lacks exposure")
        n = int(counts[use].sum())
        n_pix = int(use.sum())
        expected = bg_counts * n_pix / bg_pixels
        prob = float(binom.sf(n - 1, bg_counts + n, n_pix / (bg_pixels + n_pix)))
        apertures[str(radius)] = dict(
            counts=n,
            pixels=n_pix,
            expected_background=expected,
            signed_excess=n - expected,
            local_background_p=prob,
            upper95_counts=float(chi2.ppf(0.95, 2 * (n + 1)) / 2),
        )
    record = dict(
        group=year,
        header_utc=[exp_hdr["DATE-OBS"], exp_hdr["DATE-END"]],
        full_band_kev=[0.5, 6.0],
        pointing_count=round(float(exp_hdr["EXPOSURE"])),
        nominal_ontime_s=fits_time,
        central_exposure_map_value=float(exposure[round(y), round(x)]),
        exposure_bunit=exp_hdr.get("BUNIT"),
        exposure_to_counts_sampling_shift_dy_dx_pix=[ey - y, ex - x],
        background_annulus_arcsec=list(BG),
        background_pixels=bg_pixels,
        background_counts=bg_counts,
        apertures=apertures,
        archive_url=source.url,
        archive_etag=source.etag,
        tar_members={k: dict(offset=v[0], bytes=v[1]) for k, v in members.items()},
        ranges=source.ranges,
    )
    return record, (counts, x, y, scale)


def plot(records, images):
    fig = Figure(figsize=(9.5, 5.4), layout="constrained")
    layout = fig.add_gridspec(2, 3, height_ratios=(1.4, 1))
    graph = fig.add_subplot(layout[0, :])
    for i, (record, (data, x, y, scale)) in enumerate(zip(records, images, strict=True)):
        m = record["apertures"]["1.5"]
        graph.scatter(
            i - 0.1, m["counts"], color="#345f81", label="Target counts" if i == 0 else None
        )
        graph.scatter(
            i + 0.1,
            m["expected_background"],
            facecolors="none",
            edgecolors="#87583b",
            label="Expected background" if i == 0 else None,
        )
        graph.plot(
            [i - 0.24, i + 0.24],
            [m["upper95_counts"]] * 2,
            color=".31",
            ls="--",
            label="95% upper bound" if i == 0 else None,
        )
        ax = fig.add_subplot(layout[1, i])
        ax.imshow(data, cmap="Blues", origin="lower", vmin=0, vmax=4)
        ax.add_patch(Circle((x, y), 1.5 / scale, fill=False, ec="#bb5039", lw=1))
        ax.set(
            xlim=(x - 18, x + 18),
            ylim=(y - 18, y + 18),
            xticks=[],
            yticks=[],
            title=f"{record['group']} count image",
        )
    graph.set(
        xticks=range(len(records)),
        xticklabels=[r["group"] for r in records],
        ylabel="Counts in 1.5″ aperture",
        ylim=(-0.5, None),
    )
    graph.legend(frameon=False, ncols=3, fontsize=9)
    fig.suptitle("JADES 103607 · Chandra Deep Field South by observing group", fontsize=13)
    save_figure(fig, FIGURES / STEM)


def main():
    coords = json.loads(SOURCE.read_text())["icrs_deg"]
    sky = SkyCoord(*coords, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    records, images = [], []
    for year in FILES:
        record, display = measurement(cache, year, sky)
        records.append(record)
        images.append(display)
        print(year, record["central_exposure_map_value"], record["apertures"]["1.5"], flush=True)
    plot(records, images)
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=coords,
            position_source_sha256=text_digest(SOURCE),
            epochs=records,
            scope="Counts and local annular background in public, disjoint 0.5–6 keV count images. "
            "The upper bound is on aperture counts without background subtraction or PSF correction. "
            "Exposure-map values verify coverage in the released maps and are not flux calibration.",
        ),
    )


if __name__ == "__main__":
    main()

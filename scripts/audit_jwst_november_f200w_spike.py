"""Check a November single-filter excess directly on the three detector frames."""

import json

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.wcs import WCS
from compare_jwst_brightening_mast import REGION_ARCSEC
from download_jwst_cross_detector import OUT, ROOT
from matplotlib.figure import Figure
from measure_jwst_brightening_joint_errors import PinnedSections
from measure_jwst_epoch_candidate import (
    aperture,
    photometry,
    pixel_scale,
    show_stamp,
    verified_path,
)
from read_jwst_i2d_cutout import S3, MastImageSections, end_of_block
from search_jwst_november_catalogues import STORAGE_MIB
from vet_jwst_field_b_candidate import FIGURES
from vet_jwst_november_catalogues import fit_difference

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

STEM = "november_goods_s_f200w_spike"
SATURATED = 2
EXPECTED_SOURCE_LABEL = 77


class NativeSections(MastImageSections):
    """Read exact byte stripes of fully calibrated detector SCI and DQ."""

    def __init__(self, cache, filename):
        self.cache = cache
        self.filename = filename
        self.url = f"{S3}/{filename[:7]}/{filename[:13]}/{filename}"
        response = cache.session.head(self.url, timeout=(15, 30))
        response.raise_for_status()
        self.etag = response.headers["ETag"]
        self.size = int(response.headers["Content-Length"])
        self.ranges = []

    def native(self, position):
        header, offset = self.header(0)
        if header["FILENAME"] != self.filename or header["INSTRUME"] != "NIRCAM" or header["FILTER"] != "F200W":
            raise ValueError(f"Unexpected detector header for {self.filename}")
        extensions = {}
        for _ in range(10):
            h, data = self.header(offset)
            if h["NAXIS"] != 2 or h["XTENSION"] != "IMAGE":
                raise ValueError("Unexpected extension before both detector images")
            key = h["EXTNAME"]
            if key in ("SCI", "DQ"):
                extensions[key] = h, data
            if len(extensions) == 2:
                break
            offset = data+end_of_block(h["NAXIS1"]*h["NAXIS2"]*abs(h["BITPIX"])//8)
        sci, base = extensions["SCI"]
        dq, dq_base = extensions["DQ"]
        if sci["BUNIT"] != "MJy/sr" or sci["BITPIX"] != -32 or dq["BITPIX"] != 32:
            raise ValueError("Unexpected SCI or DQ image encoding")
        wcs = WCS(sci, fix=False).celestial
        x, y = map(float, wcs.world_to_pixel(position))
        scale = pixel_scale(wcs)
        half = int(np.ceil(.84/scale))+1
        ix, iy = round(x), round(y)
        width, height = sci["NAXIS1"], sci["NAXIS2"]
        if min(ix, iy, width-1-ix, height-1-iy) < half:
            raise ValueError("Sky position has insufficient detector coverage")
        start, stop = iy-half, iy+half+1
        stride = 4*width
        cols = np.s_[:, ix-half:ix+half+1]
        image = np.memmap(self.read(base+start*stride, base+stop*stride), dtype=">f4", mode="r",
                          shape=(stop-start, width))[cols].astype(float)
        raw = np.memmap(self.read(dq_base+start*stride, dq_base+stop*stride), dtype=">i4", mode="r",
                        shape=(stop-start, width))[cols].astype(np.int64)
        flags = (raw+dq.get("BZERO", 0)).astype(np.uint32)
        if sci.get("BZERO", 0) or sci.get("BSCALE", 1) != 1 or dq.get("BSCALE", 1) != 1:
            raise ValueError("Unexpected detector image scaling")
        image *= float(sci["PIXAR_SR"])*1e15
        px, py = x-ix+half, y-iy+half
        yy, xx = np.indices(image.shape)
        radial = np.hypot(xx-px, yy-py)*scale
        maximum = np.unravel_index(np.nanargmax(np.where(radial <= .18, image, np.nan)), image.shape)
        iy_peak, ix_peak = maximum
        row = dict(file=self.filename, upstream_etag=self.etag, date_utc=header["DATE-BEG"],
                   predicted_detector_xy=[x, y], brightest_pixel_xy=[int(ix+ix_peak-half), int(iy+iy_peak-half)],
                   brightest_pixel_dq=int(flags[maximum]), brightest_pixel_njy=float(image[maximum]),
                   saturation_flagged_pixels_within_core=int(((flags & SATURATED != 0) & (radial <= .18)).sum()),
                   ranges=self.ranges)
        for label, plane in (("unmasked", image), ("dq_zero", np.where(flags == 0, image, np.nan))):
            value = aperture(plane, px, py, scale, .12, (.30, .48))
            row[f"{label}_aperture_njy"] = value if np.isfinite(value) else None
        return row, (image, flags, px, py, scale)


def plot(frames):
    fig = Figure(figsize=(8.1, 3.3), layout="constrained")
    for i, (ax, (image, flags, x, y, scale)) in enumerate(zip(fig.subplots(1, 3), frames, strict=True), 1):
        show_stamp(ax, image, x, y, scale, (-2., 10.), "gray", .12, half_arcsec=.36)
        yy, xx = np.where(flags & SATURATED)
        ax.scatter((xx-x)*scale, (yy-y)*scale, marker="x", s=45, linewidths=1.3, color="#ce4156")
        ax.set(xlim=(-.36, .36), ylim=(-.36, .36), title=f"NRCA3 exposure {i}")
    fig.suptitle("F200W residual near GO 8060 target 2510042\n"
                 "Native science pixels · red × saturation flag · cyan 0.12″ aperture")
    save_figure(fig, FIGURES/STEM)


def main():
    catalogue_path = OUT/"november_goods_s_catalogue.json"
    catalogue = json.loads(catalogue_path.read_text())
    receipt = next(r["product"] for r in catalogue["products"] if r["visit"] == "12577-003" and r["band"] == "F200W")
    table = QTable.read(verified_path(receipt), format="ascii.ecsv")
    source, = table[table["label"] == EXPECTED_SOURCE_LABEL]
    sky = source["sky_centroid"].icrs
    prior = OUT/"november_goods_s_catalogue_vetting.json"
    data = json.loads(prior.read_text())
    nearby = next(x for x in data["candidates"] if x["candidate"] == "GSN-53.195383-27.779348")
    center = SkyCoord(*nearby["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    contrasts = {}
    for band, pair in nearby["bands"].items():
        stamps = [PinnedSections(cache, r).cutout(center, REGION_ARCSEC)[0] for r in pair["products"]]
        _result, view = fit_difference(*stamps, band, 1)
        _, wcs, _, _, x0, y0 = stamps[0]
        x, y = np.array(wcs.world_to_pixel(sky))-[x0, y0]
        contrasts[band] = photometry(view[2], x, y, view[-1], (.12, .18, .24), (.30, .48))
    rows, frames = [], []
    for i in (1, 2, 3):
        filename = f"jw12577003001_03101_{i:05d}_nrca3_cal.fits"
        row, frame = NativeSections(cache, filename).native(sky)
        rows.append(row)
        frames.append(frame)
        print(filename, row["brightest_pixel_dq"], row["brightest_pixel_njy"], flush=True)
    plot(frames)
    write_json(OUT/f"{STEM}.json", dict(
        executed_utc=utcnow(), icrs_deg=[sky.ra.deg, sky.dec.deg], separation_from_fading_source_arcsec=center.separation(sky).arcsec,
        st_sci_catalogue=receipt, source_label=EXPECTED_SOURCE_LABEL, signed_mosaic_photometry=contrasts, exposures=rows,
        script_sha256=text_digest(__file__), vetting_sha256=text_digest(prior),
        scope="Byte-pinned Science/DQ observations on NRCA3. Native positions use calibrated WCS without local translation."
              " Aperture values null when the aperture has missing or DQ-masked pixels.",
    ))


if __name__ == "__main__":
    main()

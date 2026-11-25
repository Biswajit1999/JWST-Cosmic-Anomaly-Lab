"""Read exact SCI and WHT stripes from public STScI JWST FITS mosaics."""

import re

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from measure_jwst_epoch_candidate import pixel_scale

from cosmic_anomaly_lab.provenance import digest

FITS_BLOCK = 2880
S3 = "https://stpubdata.s3.us-east-1.amazonaws.com/jwst/public"


def end_of_block(size):
    return (size+FITS_BLOCK-1)//FITS_BLOCK*FITS_BLOCK


class MastImageSections:
    """Pin a public S3 object and cache each exact byte range used by a cutout."""

    def __init__(self, cache, filename):
        pattern = re.fullmatch(r"(jw\d{5})-(o\d{3}|c\d{4})_t\d+_nircam_clear-f\d+w_i2d\.fits", filename)
        if not pattern:
            raise ValueError("Expected a NIRCam imaging stage-3 association file")
        program, observation = pattern.groups()
        self.filename = filename
        self.cache = cache
        self.url = f"{S3}/{program}/L3/t/{observation}/{filename}"
        response = cache.session.head(self.url, timeout=(15, 30))
        response.raise_for_status()
        self.etag = response.headers["ETag"]
        if not re.fullmatch(r'"[a-f0-9]+(?:-\d+)?"', self.etag):
            raise ValueError("Missing stable S3 version validator")
        self.size = int(response.headers["Content-Length"])
        self.ranges = []

    def read(self, start, stop):
        end = stop-1
        extent = f"bytes {start}-{end}/{self.size}"
        validator = self.etag[1:-1]
        name = f"cross_detector/mast_stripes/{self.filename}/{validator}-{start}-{end}.bytes"
        session = self.cache.session
        old_headers = session.headers.copy()
        old_hooks = session.hooks["response"]

        def require_range(response, **_kwargs):
            response.raise_for_status()
            if response.status_code != 206 or response.headers.get("Content-Range") != extent:
                raise ValueError(f"Expected the exact byte range {extent} from {self.filename}")
            if response.headers.get("ETag") != self.etag:
                raise ValueError(f"Upstream ETag changed while reading {self.filename}")

        session.headers.update({"Range": f"bytes={start}-{end}", "If-Match": self.etag,
                                "Accept-Encoding": "identity"})
        session.hooks["response"] = [*old_hooks, require_range]
        try:
            path = self.cache.fetch(self.url, name, survey="MAST JWST public S3",
                                    version=f"{self.etag} {extent}", kind="observed_fits_byte_range")
        finally:
            session.headers = old_headers
            session.hooks["response"] = old_hooks
        if path.stat().st_size != stop-start:
            raise ValueError(f"Unexpected cached byte count for {self.filename}")
        self.ranges.append(dict(url=self.url, path=path.relative_to(self.cache.root).as_posix(),
                                sha256=digest(path), first_byte=start, last_byte=end))
        return path

    def header(self, start):
        data = b""
        for block in range(3):
            offset = start+block*8*FITS_BLOCK
            data += self.read(offset, offset+8*FITS_BLOCK).read_bytes()
            for i in range(len(data)-8*FITS_BLOCK, len(data), 80):
                if data[i:i+8] == b"END     ":
                    length = end_of_block(i+80)
                    return fits.Header.fromstring(data[:length].decode("ascii")), start+length
        raise ValueError(f"Missing FITS END card at byte {start} in {self.filename}")

    def image_extensions(self):
        primary, start = self.header(0)
        if primary["NAXIS"] != 0 or primary["INSTRUME"] != "NIRCAM" or primary["S_RESAMP"] != "COMPLETE":
            raise ValueError(f"Unexpected primary header in {self.filename}")
        if primary["FILENAME"] != self.filename:
            raise ValueError(f"Upstream mosaic filename disagreement for {self.filename}")
        images = {}
        for _ in range(10):
            hdr, offset = self.header(start)
            if hdr["XTENSION"] != "IMAGE":
                raise ValueError(f"Expected image extension at byte {start} in {self.filename}")
            if hdr["EXTNAME"] in {"SCI", "WHT"}:
                images[hdr["EXTNAME"]] = hdr, offset
            if set(images) == {"SCI", "WHT"}:
                return primary, images
            shape = [hdr[f"NAXIS{i}"] for i in range(1, hdr["NAXIS"]+1)]
            start = offset+end_of_block(int(np.prod(shape))*abs(hdr["BITPIX"])//8)
        raise ValueError(f"Missing SCI or WHT in {self.filename}")

    def stripe(self, entry, bounds):
        header, data_start = entry
        x0, x1, y0, y1 = bounds
        if header["NAXIS"] != 2 or header["BITPIX"] != -32 or header.get("BSCALE", 1) != 1:
            raise ValueError(f"Unexpected image encoding in {self.filename}")
        if header.get("BZERO", 0) != 0 or x1 > header["NAXIS1"] or y1 > header["NAXIS2"]:
            raise ValueError(f"Invalid image range in {self.filename}")
        row_bytes = 4*header["NAXIS1"]
        path = self.read(data_start+y0*row_bytes, data_start+y1*row_bytes)
        data = np.memmap(path, dtype=">f4", mode="r", shape=(y1-y0, header["NAXIS1"]))
        return np.array(data[:, x0:x1], dtype=np.float32)

    def cutout(self, sky, half_arcsec):
        primary, images = self.image_extensions()
        h = images["SCI"][0]
        if h["BUNIT"] != "MJy/sr":
            raise ValueError(f"Unexpected science units in {self.filename}")
        wcs = WCS(h, fix=False).celestial
        target_x, target_y = map(float, wcs.world_to_pixel(sky))
        scale = pixel_scale(wcs)
        pixels = int(np.ceil(half_arcsec/scale))
        if not (0 <= target_x < h["NAXIS1"] and 0 <= target_y < h["NAXIS2"]):
            raise ValueError(f"Target lies outside the mosaic in {self.filename}")
        x0, x1 = max(0, round(target_x)-pixels), min(h["NAXIS1"], round(target_x)+pixels+1)
        y0, y1 = max(0, round(target_y)-pixels), min(h["NAXIS2"], round(target_y)+pixels+1)
        bounds = (x0, x1, y0, y1)
        data, weight = (self.stripe(images[name], bounds) for name in ("SCI", "WHT"))
        data[~np.isfinite(weight) | (weight <= 0)] = np.nan
        receipt = dict(filename=self.filename, url=self.url, etag=self.etag, object_bytes=self.size,
                       science_array_shape=[h["NAXIS2"], h["NAXIS1"]], cutout_bounds_x0_x1_y0_y1=list(bounds),
                       date_utc=primary["DATE-BEG"], exposure_start_mjd=primary["EXPSTART"], filter=primary["FILTER"],
                       pipeline_version=primary["CAL_VER"], crds_context=primary["CRDS_CTX"],
                       native_pixel_scale_arcsec=scale, header_pixel_area_sr=h["PIXAR_SR"], ranges=self.ranges)
        return (data, wcs, target_x-x0, target_y-y0, x0, y0), receipt

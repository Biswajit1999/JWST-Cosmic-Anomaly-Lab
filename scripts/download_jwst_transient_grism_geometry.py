"""Fetch bright-epoch pointing headers, calibration models and direct F356W images."""

import json

from astropy.io import fits
from audit_jwst_transient_wfss import MAST
from download_jwst_cross_detector import OUT, ROOT

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

STORAGE_MIB = 11500
REFERENCE_KEYS = ("R_SPCWCS", "R_DISTOR", "R_FILOFF", "R_WAVRAN")
RANGE = "bytes=0-57599"
CRDS = "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"


def fetch(cache, url, name, kind, **details):
    path = cache.fetch(url, "cross_detector/"+name, survey=details.pop("survey", "MAST JWST"),
                       version=details.pop("version", "GO 4549 public archive"), kind=kind)
    return {
        "url": url, "path": path.relative_to(ROOT).as_posix(),
        "sha256": digest(path), "bytes": path.stat().st_size, **details,
    }


def main():
    audit_path = OUT / "field_b_transient_wfss.json"
    audit = json.loads(audit_path.read_text())
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    output = OUT / "field_b_transient_grism_products.json"
    result = {
        "candidate": audit["candidate"], "position_degrees": audit["position_degrees"],
        "calibrations": {}, "rate_headers": [], "direct_images": [], "audit_sha256": text_digest(audit_path),
        "storage_cap_mib": STORAGE_MIB, "script_sha256": text_digest(__file__),
    }

    def save():
        write_json(output, {"executed_utc": utcnow(), **result})

    seen = set()
    for exposure in audit["x1d_exposures"]:
        spec = ROOT / exposure["path"]
        if digest(spec) != exposure["sha256"]:
            raise ValueError(f"Archived extraction checksum failed: {spec.name}")
        head = fits.getheader(spec)
        config = head["CRDS_CTX"]
        references = {k: head[k].removeprefix("crds://") for k in REFERENCE_KEYS}
        for name in set(references.values())-seen:
            result["calibrations"][name] = fetch(cache, CRDS+name, name, "calibration_reference",
                                                 version=config, survey="JWST CRDS")
            seen.add(name)
        name = exposure["filename"].replace("_x1d.fits", "_rate.fits")
        url = MAST + "uri=mast%3AJWST%2Fproduct%2F"+name
        cache.session.headers["Range"] = RANGE
        try:
            receipt = fetch(cache, url, name+".header.bin", "observation_header_range",
                            version="bytes 0-57599", byte_range=RANGE)
        finally:
            del cache.session.headers["Range"]
        if receipt["bytes"] != 57600:
            raise ValueError(f"Range response has incorrect length for {name}")
        result["rate_headers"].append(dict(**receipt, rate_filename=name, x1d_path=exposure["path"],
                                           x1d_sha256=exposure["sha256"], references=references))
        print(name, flush=True)
        save()

    for direct in audit["direct_image_catalogues"]:
        name = direct["filename"].replace("_cat.ecsv", "_i2d.fits")
        url = MAST + "uri=mast%3AJWST%2Fproduct%2F"+name
        result["direct_images"].append(fetch(cache, url, name, "observed_mast_direct_mosaic"))
        print(name, flush=True)
        save()


if __name__ == "__main__":
    main()

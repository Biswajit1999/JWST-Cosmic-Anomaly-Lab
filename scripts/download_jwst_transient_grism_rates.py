"""Retrieve native rate frames for the six exposures with full grism coverage."""

import json

from astropy.io import fits
from audit_jwst_transient_wfss import MAST
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import CRDS, STORAGE_MIB, fetch

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json


def main():
    geometry_path = OUT / "field_b_transient_grism_geometry.json"
    geometry = json.loads(geometry_path.read_text())
    products_path = OUT / "field_b_transient_grism_products.json"
    products = json.loads(products_path.read_text())
    if geometry["products_sha256"] != text_digest(products_path):
        raise ValueError("Geometry refers to a different pointing ledger")
    indices = {r["rate_filename"]: r for r in products["rate_headers"]}
    selected = [r for r in geometry["exposure_coverage"] if r["nominal"]["fraction_of_band"] > .99]
    if len(selected) != 6:
        raise ValueError("Expected six full-coverage observations for this fixed execution")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    receipts, references = [], {}
    path = OUT / "field_b_transient_grism_rates.json"
    for row in selected:
        name = row["filename"]
        header = fits.getheader(ROOT / indices[name]["x1d_path"])
        needed = {key: header[key].removeprefix("crds://") for key in ("R_FLAT", "R_PHOTOM")}
        for n in sorted(set(needed.values())-set(references)):
            references[n] = fetch(cache, CRDS+n, n, "calibration_reference",
                                  version=header["CRDS_CTX"], survey="JWST CRDS")
        raw = fetch(cache, MAST+"uri=mast%3AJWST%2Fproduct%2F"+name, name, "observed_native_grism_rate")
        receipts.append({**raw, "filename": name, "references": needed})
        write_json(path, {
            "executed_utc": utcnow(), "candidate": geometry["candidate"], "position_degrees": geometry["position_degrees"],
            "geometry_sha256": text_digest(geometry_path), "script_sha256": text_digest(__file__),
            "native_rate_frames": receipts, "calibrations": references,
        })
        print(name, flush=True)


if __name__ == "__main__":
    main()

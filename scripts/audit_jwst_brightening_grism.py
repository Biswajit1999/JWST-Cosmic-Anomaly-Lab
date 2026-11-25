"""Check archival NIRCam slitless products for a spectrum at the brightening position."""

import json

import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from audit_jwst_module_b_spectra import RADIUS, order_positions
from check_jwst_transient_precursor import contains
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import fetch
from measure_jwst_epoch_candidate import verified_path

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

COVERAGE = OUT/"module_b_001_release_20261124.json"
STORAGE_MIB = 24_000


def main():
    release = json.loads(COVERAGE.read_text())
    sky = SkyCoord(*release["position_degrees"], unit="deg")
    listing = release["dawn"]["product"]
    exposures = pd.read_csv(verified_path(listing))
    modes = exposures[(exposures.extension == "rate") & exposures.pupil.isin(("GRISMR", "GRISMC"))]
    modes = modes[modes.footprint.map(lambda corners: contains(corners, sky))].sort_values("expstart")
    if len(modes) != 18:
        raise ValueError("Unexpected number of covering grism exposures in the pinned snapshot")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    products = []
    for row in modes.itertuples():
        name = row.file + "_x1d.fits"
        url = f"https://stpubdata.s3.us-east-1.amazonaws.com/jwst/public/{name[:7]}/{name[:13]}/{name}"
        product = fetch(cache, url, name, "observed_module_b_001_grism_x1d",
                        version="Public calibrated NIRCam slitless extraction")
        with fits.open(verified_path(product)) as h:
            hdr = h[0].header
            if (hdr["DETECTOR"], hdr["FILTER"]+"-"+hdr["PUPIL"]) != (row.detector, row.filter):
                raise ValueError(f"Unexpected grism metadata in {name}")
            order_groups = [order_positions(x, sky) for x in h if x.name == "EXTRACT1D"]
            if not order_groups:
                raise ValueError(f"Missing extraction extension in {name}")
            products.append(dict(product=product, visit=row.file[:13], filter=hdr["FILTER"], detector=hdr["DETECTOR"],
                                 disperser=hdr["PUPIL"], exposure_start_mjd=hdr["EXPSTART"], orders=order_groups))
        print(name, [x["nearest"] for x in order_groups], flush=True)
    write_json(OUT/"module_b_001_archived_grism.json", dict(
        executed_utc=utcnow(), candidate=release["candidate"], icrs_deg=release["position_degrees"],
        covering_exposure_listing=listing, release_sha256=text_digest(COVERAGE), match_radius_arcsec=RADIUS,
        products=products, order_position_helper_sha256=text_digest(ROOT/"scripts/audit_jwst_module_b_spectra.py"),
        script_sha256=text_digest(__file__),
        scope="Tests geometric direct-image footprints and published extraction positions. Dispersed wavelength"
              " coverage, blending, detection completeness and targeted spectrum extraction remain unmeasured.",
    ))


if __name__ == "__main__":
    main()

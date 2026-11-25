"""Check the newest public JWST follow-up holdings and TNS snapshot for JADES 103607."""

import json
from datetime import datetime, timezone
from urllib.parse import urlencode
from xml.etree import ElementTree

from astropy.coordinates import SkyCoord
from audit_jwst_brightening_followup import mast_followup
from audit_jwst_goods_s_prior_records import proposal
from check_jwst_transient_precursor import contains
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import fetch
from measure_jwst_epoch_candidate import verified_path
from vet_jwst_module_b import archive_query, tns_catalogue, tns_match

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

STORAGE_MIB = 25_500
OBSERVATION = 4
BRIGHT_EPOCH_END_MJD = 61323.10


def cloud_products(cache, tag):
    prefixes = {
        "visit_products": "jwst/public/jw12577/jw12577004001/",
        "stage3_products": "jwst/public/jw12577/L3/t/o004/",
    }
    listing = {}
    for label, prefix in prefixes.items():
        args = {"list-type": 2, "prefix": prefix, "max-keys": 1000}
        keys, pages = [], []
        while True:
            url = "https://stpubdata.s3.us-east-1.amazonaws.com/?"+urlencode(args)
            receipt = fetch(cache, url, f"module_b_001_o004_{label}_{tag}_{len(pages)+1}.xml",
                            "public_s3_observation_inventory", survey="MAST public S3",
                            version=f"GO 12577 observation {OBSERVATION}")
            root = ElementTree.parse(verified_path(receipt))
            page = [key.text for key in root.findall(".//{*}Key")]
            truncated = root.findtext("./{*}IsTruncated")
            if int(root.findtext("./{*}KeyCount")) != len(page) or truncated not in {"true", "false"}:
                raise ValueError(f"Invalid public S3 inventory at {prefix}")
            keys.extend(page)
            pages.append(receipt)
            if truncated == "false":
                break
            token = root.findtext("./{*}NextContinuationToken")
            if not token or token == args.get("continuation-token"):
                raise ValueError(f"Missing or repeated S3 continuation token at {prefix}")
            args["continuation-token"] = token
        listing[label] = dict(prefix=prefix, keys=keys, products=pages, truncated=False)
    return listing


def dawn_exposures(cache, sky, tag):
    data, receipt = archive_query(cache, "exposures", dict(
        ra=sky.ra.deg, dec=sky.dec.deg, arcmin=.013, output="csv"), f"module_b_001_release_dawn_{tag}.csv")
    covering = data.loc[data.extension == "rate"].copy()
    covering = covering[covering.footprint.map(lambda footprint: contains(footprint, sky))]
    if not covering.file.str.startswith("jw12577001001").any():
        raise ValueError("The DAWN inventory did not recover the October 2026 control visit")
    later = covering[covering.expstart > BRIGHT_EPOCH_END_MJD].sort_values("expstart")
    return dict(product=receipt, rows=len(data), covering_rows=len(covering),
                latest_covering_start_mjd=float(covering.expstart.max()),
                later_covering=json.loads(later[["file", "expstart", "detector", "filter"]].to_json(orient="records")))


def main():
    stamp = datetime.now(timezone.utc)
    tag = stamp.strftime("%Y%m%dT%H%M%SZ")
    prior = OUT/"module_b_001_native_epoch_flux.json"
    target = json.loads(prior.read_text())
    sky = SkyCoord(*target["icrs_deg"], unit="deg")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    plans = [proposal(cache, program, {target["candidate"]: sky}, stamp.strftime("%Y%m%d"))
             for program in (8060, 12577)]
    write_json(OUT/f"module_b_001_release_{stamp:%Y%m%d}.json", dict(
        executed_utc=utcnow(), candidate=target["candidate"], position_degrees=target["icrs_deg"],
        mast=mast_followup(sky), s3=cloud_products(cache, tag), dawn=dawn_exposures(cache, sky, tag),
        public_observing_plans=plans, tns=tns_match(tns_catalogue(cache), sky), native_receipt_sha256=text_digest(prior),
        script_sha256=text_digest(__file__),
    ))


if __name__ == "__main__":
    main()

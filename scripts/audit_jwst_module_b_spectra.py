"""Check the Module B positions against JADES and the F277W grism extractions."""

import json
from io import StringIO
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from check_jwst_transient_precursor import contains
from download_jwst_cross_detector import OUT, ROOT, download
from download_jwst_transient_grism_geometry import fetch
from measure_jwst_epoch_candidate import verified_path
from vet_jwst_module_b import archive_query

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

STORAGE_MIB = 16300
RADIUS = .20
JADES = {"001": 103607, "002": 117875}


def positional_counterpart(cache, position, source_id):
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"+urlencode({
        "-source": "V/159/gssample", "-out.max": 30,
        "-out": "NIRCamID,RA_ICRS,DE_ICRS,KronRad,zphot,zphotl68,zphotu68", "NIRCamID": source_id,
    })
    rec = download(cache, url, f"module_b_jades_{source_id}.tsv", "catalogue_source_lookup")
    table = pd.read_csv(StringIO(verified_path(rec).read_text(encoding="utf-8")), comment="#", sep="\t")
    candidates = table.apply(pd.to_numeric, errors="coerce")
    (row,) = list(candidates.loc[candidates.NIRCamID == source_id].itertuples())
    distance = position.separation(SkyCoord(row.RA_ICRS, row.DE_ICRS, unit="deg")).arcsec
    return dict(receipt=rec, nircam_id=source_id, icrs_deg=[float(row.RA_ICRS), float(row.DE_ICRS)],
                separation_arcsec=float(distance), photo_z=float(row.zphot), kron_radius_arcsec=float(row.KronRad))


def order_positions(hdu, sky):
    table = hdu.data
    located = np.isfinite(table["SOURCE_RA"]+table["SOURCE_DEC"]) & (table["SOURCE_ID"] > 0)
    pos = SkyCoord(table["SOURCE_RA"][located], table["SOURCE_DEC"][located], unit="deg")
    sep = sky.separation(pos).arcsec
    nearest = None
    if len(sep):
        index = int(np.argmin(sep))
        nearest = dict(source_id=int(table["SOURCE_ID"][located][index]), distance_arcsec=float(sep[index]),
                       icrs_deg=[float(pos[index].ra.deg), float(pos[index].dec.deg)])
    return dict(order=int(hdu.header["SPORDER"]), valid_positions=int(located.sum()),
                rows_within_match_radius=int((sep <= RADIUS).sum()), nearest=nearest)


def investigate_fading(cache, center):
    tar = SkyCoord(*center, unit="deg")
    data, listing = archive_query(cache, "exposures", dict(
        ra=center[0], dec=center[1], arcmin=.012, output="csv", filters="F277W-GRISMR,F277W-GRISMC"),
        "module_b_002_f277w_wfss.csv")
    matches = data[(data.extension == "rate") & data.file.str.startswith("jw07336")
                   & (data.detector == "NRCBLONG")].sort_values("expstart")
    products = []
    for row in matches.itertuples():
        name = row.file+"_x1d.fits"
        url = f"https://stpubdata.s3.us-east-1.amazonaws.com/jwst/public/{name[:7]}/{name[:13]}/{name}"
        receipt = fetch(cache, url, name, "observed_module_b_grism_x1d", version="GO 7336 public calibrated extraction")
        with fits.open(verified_path(receipt)) as h:
            hdr = h[0].header
            if (hdr["DETECTOR"], hdr["FILTER"], hdr["PUPIL"]) != ("NRCBLONG", "F277W", row.pupil):
                raise ValueError(f"Unexpected metadata on {name}")
            products.append(dict(product=receipt, filter=hdr["FILTER"], disperser=hdr["PUPIL"],
                                 geometric_direct_footprint=contains(row.footprint, tar),
                                 date_utc=f'{hdr["DATE-OBS"]}T{hdr["TIME-OBS"]}',
                                 extraction_orders=[order_positions(ext, tar) for ext in h if ext.name == "EXTRACT1D"]))
        print(name, products[-1]["extraction_orders"], flush=True)
    return dict(candidate_position_degrees=center, archive_listing=listing, products=products,
                positions_match_radius_arcsec=RADIUS,
                scope="Cross-match of published extraction coordinates; spectral coverage and sensitivity are untested")


def main():
    path = OUT/"module_b_vetting.json"
    vetting = json.loads(path.read_text())
    if len(vetting["candidates"]) != 2:
        raise ValueError("Require both completed imaging vetting records")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    counterparts, wfss = {}, None
    for row in vetting["candidates"]:
        pointing = row["candidate"].split("-")[1]
        sky = SkyCoord(*row["icrs_deg"], unit="deg")
        counterparts[pointing] = positional_counterpart(cache, sky, JADES[pointing])
        if row["sign"] < 0:
            wfss = investigate_fading(cache, row["icrs_deg"])
    write_json(OUT/"module_b_spectra.json", {
        "executed_utc": utcnow(), "direct_catalogue_matches": counterparts, "f277w_slitless": wfss,
        "vetting_sha256": text_digest(path), "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

"""Retrieve independent epochs and native detector exposures for the field B lead."""

import json
from urllib.parse import urlencode

import pandas as pd
from astropy.io import fits
from download_jwst_cross_detector import API, OUT, ROOT, download

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

STORAGE_MIB = 9400
EARLIER = {
    "F200W": [
        ("jw01180015001", "NRCA2"), ("jw01180013001", "NRCA2"),
        ("jw01180220001", "NRCB1"), ("jw08060003001", "NRCB3"),
        ("jw08060005001", "NRCB1"), ("jw05324027001", "NRCB1"),
    ],
    "F444W": [("jw01180015001", "NRCALONG"), ("jw01180013001", "NRCALONG")],
    "F356W": [("jw08060003001", "NRCBLONG"), ("jw08060005001", "NRCBLONG")],
}


def choose_candidate():
    matches = pd.read_csv(OUT / "field_b_matches.csv")
    eligible = matches[(matches.sign == -1) & matches.aperture_gate & (matches.f356w_aperture_gate == True)]  # noqa: E712
    if eligible.empty:
        raise ValueError("No fading source passes all three filter checks")
    return eligible.sort_values("minimum_score", ascending=False).iloc[0]


def native_exposures(cache, ra, dec, band, record):
    assoc = record["assoc_name"]
    filter_name = band+"-CLEAR"
    url = API+"/exposures?" + urlencode(dict(ra=ra, dec=dec, arcmin=.025, filters=filter_name, output="csv"))
    listing = download(cache, url, f"field_b_candidate_{band}_exposures.csv", "associated_exposures")
    exposures = pd.read_csv(ROOT / listing["path"]).query("assoc == @assoc and extension == 'rate'")
    if len(exposures) != record["exposure_count"]:
        raise ValueError(f"Unexpected number of detector exposures for {assoc}")
    records = []
    for e in exposures.sort_values("expstart").itertuples():
        name = e.file+"_cal.fits"
        url = "https://mast.stsci.edu/api/v0.1/Download/file?" + urlencode(dict(uri=f"mast:JWST/product/{name}"))
        path = cache.fetch(url, "cross_detector/"+name, survey="MAST JWST", version="calibrated native frame",
                           kind="observed_cal_exposure")
        with fits.open(path) as h:
            if h[0].header["FILTER"] != band or h[0].header["DETECTOR"] != record["detector"]:
                raise ValueError(f"Unexpected filter or detector for {name}")
        records.append(dict(file=name, path=path.relative_to(ROOT).as_posix(), url=url,
                            sha256=digest(path), expstart_mjd=float(e.expstart)))
        print(name, flush=True)
    return dict(listing=listing, products=records)


def main():
    best = choose_candidate()
    ra, dec = float(best.ra), float(best.dec)
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    visits_path = OUT / "mosaics.json"
    search = json.loads(visits_path.read_text())["images"]
    primary = [r for r in search if r["group"] == "b" and r["detector"] in (best.sw_detector, "NRCALONG")]
    receipt = dict(candidate=best.candidate, ra=ra, dec=dec, survey_visits=primary, older_epochs=[],
                   fitted_exposures={}, association_listings=[], cap_mib=STORAGE_MIB,
                   visits_sha256=text_digest(visits_path), shortlist_sha256=text_digest(OUT / "field_b_matches.csv"),
                   script_sha256=text_digest(__file__))

    def save():
        write_json(OUT / "field_b_candidate_products.json", {**receipt, "executed_utc": utcnow()})

    for band in ("F200W", "F444W"):
        older = next(r for r in primary if r["band"] == band+"-CLEAR" and r["epoch"] == "older")
        receipt["fitted_exposures"][band] = native_exposures(cache, ra, dec, band, older)
        save()
    for band, matches in EARLIER.items():
        filt = band+"-CLEAR"
        url = API+"/assoc_mosaic?" + urlencode(dict(coords=f"{ra:.8f},{dec:.8f},.02", filters=filt, output="csv"))
        listing = download(cache, url, f"field_b_candidate_{band}_assoc.csv", "associated_mosaics")
        receipt["association_listings"].append(listing)
        entries = pd.read_csv(ROOT / listing["path"])
        for visit, detector in matches:
            rows = entries.query("visit == @visit and detector == @detector and tile == 672 and status == 2")
            if len(rows) != 1:
                raise ValueError(f"Expected one match for {band} {visit} {detector}")
            row = rows.iloc[0]
            record = dict(band=filt, detector=detector, visit=visit, dateobs=row.dateobs, mjd=float(row.expstart),
                          assoc_name=row.assoc_name, exposure_count=int(row.nexp))
            for kind in ("sci", "wht"):
                url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
                record[kind] = download(cache, url, url.rsplit("/", 1)[1], "observed_candidate_epoch_"+kind)
            receipt["older_epochs"].append(record)
            save()
            print(band, visit, detector, flush=True)


if __name__ == "__main__":
    main()

"""Pin GOODS-S products for a search across NIRCam's two wavelength channels."""

from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
from astropy.coordinates import SkyCoord

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jwst_cross_detector"
STORAGE_MIB = 9000
CONTROL = (53.19541375266448, -27.77932804913642)
SEARCH_VISITS = {"older": "jw08060002001", "newer": "jw12577002001"}
CONTROL_VISITS = {"older": "jw08060001001", "newer": "jw12577001001"}
DETECTORS = {"F444W": ["NRCALONG"], "F356W": ["NRCALONG"], "F200W": ["NRCA1", "NRCA2", "NRCA3", "NRCA4"]}
API = "https://grizli-cutout.herokuapp.com"


def download(cache, url, name, kind):
    path = cache.fetch(url, "cross_detector/"+name, survey="DAWN JWST Archive",
                      version="public v7 mosaics or v4 spectra", kind=kind)
    return dict(url=url, path=path.relative_to(ROOT).as_posix(), sha256=digest(path))


def spectra(cache):
    url = API + "/nirspec_extractions?" + urlencode(
        dict(coords=f"{CONTROL[0]:.8f},{CONTROL[1]:.8f}", size=.167, output="csv"))
    query = download(cache, url, "control_nirspec_query.csv", "public_nirspec_metadata")
    rows = pd.read_csv(ROOT / query["path"])
    distance = SkyCoord(*CONTROL, unit="deg").separation(SkyCoord(rows.ra, rows.dec, unit="deg")).arcsec
    rows = rows.assign(separation_arcsec=distance)
    match = rows[(rows.separation_arcsec < .25) & (rows.grating == "PRISM_CLEAR") & (rows.grade >= 3)]
    assert len(match) == 4
    products = []
    for row in match.itertuples():
        url = f"https://s3.amazonaws.com/msaexp-nirspec/extractions/{row.root}/{row.file}"
        record = download(cache, url, row.file, "observed_extracted_nirspec")
        products.append({**record, "dataset": row.file, "separation_arcsec": row.separation_arcsec,
                         "msamet": row.msamet, "dja_z": row.z, "dja_grade": row.grade,
                         "dja_comment": row.comment if pd.notna(row.comment) else None})
    write_json(OUT / "control_spectra_products.json", {
        "executed_utc": utcnow(), "position_degrees": CONTROL, "query": query, "spectra": products,
        "script_sha256": text_digest(__file__),
    })


def mosaics(cache):
    records, searches = [], []
    for group, specs, query_point, visits in [
        ("control", {"F200W": ["NRCA3"]}, f"{CONTROL[0]:.8f},{CONTROL[1]:.8f},0.019", CONTROL_VISITS),
        ("b", DETECTORS, "53.187,-27.772,5.5", SEARCH_VISITS),
    ]:
        for band, detectors in specs.items():
            filt = band+"-CLEAR"
            url = API + "/assoc_mosaic?" + urlencode(dict(coords=query_point, filters=filt, output="csv"))
            listing = download(cache, url, f"{group}_{band}_associations.csv", "archive_association_metadata")
            searches.append(listing)
            rows = pd.read_csv(ROOT / listing["path"])
            for detector in detectors:
                for epoch, visit in visits.items():
                    match = rows.query("visit == @visit and status == 2 and tile == 672 and detector == @detector")
                    if len(match) != 1:
                        raise ValueError(f"Expected one match for {band} {visit} {detector}, found {len(match)}")
                    row = match.iloc[0]
                    record = dict(group=group, detector=detector, band=filt, epoch=epoch, dateobs=row.dateobs,
                                  detector_center=[float(row.ra), float(row.dec)], mjd=float(row.expstart),
                                  visit=visit, assoc_name=row.assoc_name, exposure_count=int(row.nexp))
                    for kind in ("sci", "wht"):
                        url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
                        record[kind] = download(cache, url, url.rsplit("/", 1)[1], "observed_epoch_"+kind)
                    records.append(record)
                    print(group, band, detector, visit, flush=True)
                    write_json(OUT / "mosaics.json", {
                        "executed_utc": utcnow(), "images": records, "associations": searches,
                        "script_sha256": text_digest(__file__), "raw_storage_cap_mib": STORAGE_MIB,
                    })


if __name__ == "__main__":
    download_cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    spectra(download_cache)
    mosaics(download_cache)

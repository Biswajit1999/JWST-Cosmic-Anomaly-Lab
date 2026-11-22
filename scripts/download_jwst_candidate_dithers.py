"""Retrieve calibrated detector exposures for brightening and fading checks."""

import json
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
from astropy.io import fits

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
out = root / "results/jwst_epoch_search"
target = json.loads((out / "target_followup_products.json").read_text())
primary = json.loads((out / "products.json").read_text())["images"]
budget_mib = 6200
cache = DownloadCache(root, budget_mib*1024**2)


def get_exposures(receipt, label):
    band, visit, assoc = (receipt[k] for k in ("band", "visit", "assoc_name"))
    params = urlencode(dict(ra=target["ra"], dec=target["dec"], arcmin=.025, filters=band, output="csv"))
    url = "https://grizli-cutout.herokuapp.com/exposures?" + params
    listing = cache.fetch(url, f"epoch_search/target_{label}_exposures.csv", survey="DJA",
                          version="epoch-search", kind="public_exposure_metadata")
    rows = pd.read_csv(listing).query("assoc == @assoc and extension == 'rate'")
    assert len(rows) == 3 and rows["file"].nunique() == 3
    records = []
    for r in rows.sort_values("expstart").itertuples():
        assert str(r.file).startswith(visit) and r.detector == "NRCALONG"
        name = r.file + "_cal.fits"
        download = "https://mast.stsci.edu/api/v0.1/Download/file?" + urlencode(dict(uri="mast:JWST/product/"+name))
        path = cache.fetch(download, "epoch_search/"+name, survey="MAST JWST",
                           version="calibrated exposure", kind="observed_cal_exposure")
        with fits.open(path, memmap=True) as h:
            assert h[0].header["FILTER"] == band.split("-")[0] and h[0].header["DETECTOR"] == "NRCALONG"
            assert h["SCI"].header["BUNIT"] == "MJy/sr"
        records.append(dict(file=name, path=path.relative_to(root).as_posix(), url=download, sha256=digest(path),
                            expstart_mjd=float(r.expstart), exposure_seconds=float(r.exptime)))
        print(name, path.stat().st_size, flush=True)
    return dict(band=band, visit=visit, assoc=assoc, products=records,
                listing_url=url, listing_sha256=digest(listing))


current = primary + target["images"]
groups = {}
for band, year, label in [
    ("F444W-CLEAR", "2026", "f444w"), ("F277W-CLEAR", "2026", "f277w"),
    ("F356W-CLEAR", "2026", "f356w"), ("F444W-CLEAR", "2025", "f444w_2025"),
]:
    match = next(r for r in current if r["band"] == band and r["dateobs"].startswith(year))
    groups[label] = get_exposures(match, label)
    write_json(out / "dither_products.json", dict(executed_utc=utcnow(), candidate=target["candidate"],
        target_ra=target["ra"], target_dec=target["dec"], groups=groups,
        script_sha256=text_digest(__file__), raw_budget_mib=budget_mib))

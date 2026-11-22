"""Retrieve fixed independent NIRCam visits from the public DJA association archive."""

from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
import yaml

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
out = root / "results/jwst_epoch_search"
config_path = root / "configs/jwst_epoch_search.yaml"
config = yaml.safe_load(config_path.read_text())
cache = DownloadCache(root, config["max_raw_storage_mib"] * 1024**2)
coord = ",".join(str(config[key]) for key in ["query_ra_deg", "query_dec_deg", "query_width_arcmin"])


def fetch(url, name, kind):
    path = cache.fetch(url, "epoch_search/" + name, survey="DAWN JWST Archive",
                       version="association mosaic v7.0", kind=kind)
    return dict(url=url, path=path.relative_to(root).as_posix(), sha256=digest(path))


images, listings = [], []
for band in config["filters"]:
    base = "https://grizli-cutout.herokuapp.com/assoc_mosaic?"
    listing = fetch(base+urlencode(dict(coords=coord, filters=band, output="csv")),
                    f"{band}_associations.csv", "archive_association_metadata")
    listings.append(listing)
    data = pd.read_csv(root / listing["path"])
    for epoch in ["older", "newer"]:
        visit = config[epoch+"_visit"]
        rows = data[(data.visit == visit) & (data.tile == config["tile"])
                    & (data["filter"] == band) & (data.status == 2)]
        if len(rows) != 1:
            raise ValueError(f"Ambiguous association for {band} {visit}: {len(rows)} rows")
        row = rows.iloc[0]
        assert row["file"].endswith("_sci.fits.gz") and row["instrument"] == "NIRCAM"
        entry = dict(visit=visit, epoch=epoch, band=band, assoc_name=row.assoc_name,
                     dateobs=row.dateobs, mjd=float(row.expstart), exposure_count=int(row.nexp),
                     exposure_seconds=float(row.exptime), pixel_scale_mas=float(row.pixscale_mas))
        for kind in ["sci", "wht"]:
            url = row["file"].replace("_sci.fits.gz", f"_{kind}.fits.gz")
            entry[kind] = fetch(url, url.rsplit("/", 1)[1], "observed_epoch_" + kind)
        images.append(entry)
        write_json(out / "products.json", dict(executed_utc=utcnow(), images=images, listings=listings,
                   config_sha256=text_digest(config_path)))
        print(band, visit, row.dateobs, row.nexp, flush=True)

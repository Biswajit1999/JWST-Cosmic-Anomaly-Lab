"""Retrieve fixed-coordinate additional-band and earlier-epoch checks."""

from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
import requests
import yaml

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
out = root / "results/jwst_epoch_search"
config = yaml.safe_load((root / "configs/jwst_epoch_search.yaml").read_text())
data = pd.read_csv(out / "matched_peaks.csv")
selected = data[(data.sign == 1) & data.compactness_gate].copy()
for band in ["f356w", "f444w"]:
    # Direct-aperture flux, not the rank of a small difference-image extremum.
    selected = selected[(selected[band+"_new_njy"] > 0)
                        & (selected[band+"_new_njy"]-selected[band+"_old_njy"] >= .25*selected[band+"_new_njy"])]
target = selected.sort_values("minimum_absolute_score", ascending=False).iloc[0]
ra, dec = float(target.ra), float(target.dec)
cache = DownloadCache(root, config["max_raw_storage_mib"] * 1024**2)


def fetch(url, name, kind):
    path = cache.fetch(url, "epoch_search/" + name, survey="Public archive",
                       version="pinned epoch-search follow-up", kind=kind)
    return dict(url=url, path=path.relative_to(root).as_posix(), sha256=digest(path))


spec = {"F277W-CLEAR": ["jw08060001001", "jw12577001001"],
        "F444W-CLEAR": ["jw01180009001"]}
receipts, association_receipts = [], []
for band, visits in spec.items():
    params = urlencode(dict(coords=f"{ra:.8f},{dec:.8f},0.03", filters=band, output="csv"))
    url = "https://grizli-cutout.herokuapp.com/assoc_mosaic?" + params
    listing = fetch(url, "followup_"+band+"_associations.csv", "archive_association_metadata")
    association_receipts.append(listing)
    rows = pd.read_csv(root / listing["path"])
    for visit in visits:
        found = rows[(rows.visit == visit) & (rows.status == 2) & (rows["filter"] == band)
                     & (rows.tile == config["tile"])]
        if len(found) != 1:
            raise ValueError(f"Ambiguous follow-up association for {visit} {band}: {len(found)}")
        row = found.iloc[0]
        rec = dict(band=band, visit=visit, dateobs=row.dateobs, mjd=float(row.expstart),
                   assoc_name=row.assoc_name, exposure_count=int(row.nexp))
        for kind in ["sci", "wht"]:
            url = row.file.replace("_sci.fits.gz", f"_{kind}.fits.gz")
            rec[kind] = fetch(url, url.rsplit("/", 1)[1], "observed_epoch_"+kind)
        receipts.append(rec)
        print(band, visit, row.dateobs, flush=True)
        write_json(out / "target_followup_products.json", dict(executed_utc=utcnow(), candidate=target.candidate,
            ra=ra, dec=dec, selection="Highest matched-score compact brightening with >=25% direct aperture change in both bands",
            images=receipts, associations=association_receipts, peak_table_sha256=text_digest(out / "matched_peaks.csv"),
            script_sha256=text_digest(__file__)))

# Public object query, distinct from prior published transient literature.
adql = f"SELECT TOP 30 main_id,ra,dec,otype FROM basic WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.0013888889))=1"
catalogues = {}
for label, url in [
    ("simbad", "https://simbad.cds.unistra.fr/simbad/sim-tap/sync?" + urlencode(
        dict(request="doQuery", lang="adql", format="csv", query=adql))),
    ("tns", "https://www.wis-tns.org/search?" + urlencode(dict(ra=ra, decl=dec, radius=5, coords_unit="arcsec", format="html"))),
]:
    try:
        record = fetch(url, label+"_target_query."+("csv" if label == "simbad" else "html"), "public_novelty_check")
        catalogues[label] = dict(status="RESPONSE_RETRIEVED", **record)
    except requests.RequestException as exc:
        catalogues[label] = dict(status="UNAVAILABLE", url=url, error=str(exc))
    write_json(out / "novelty_queries.json", dict(executed_utc=utcnow(), ra=ra, dec=dec, queries=catalogues))

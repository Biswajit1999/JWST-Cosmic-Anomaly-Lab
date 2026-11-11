import argparse
import sys

sys.path.insert(0, "src")
import pandas as pd

from cosmic_anomaly_lab.provenance import DownloadCache, write_json

parser = argparse.ArgumentParser(description="Download observed NIRCam pilot cutouts")
parser.add_argument("--all", action="store_true", help="Request all 20 pilot sources")
args = parser.parse_args()
cache = DownloadCache(".")
p = pd.read_csv("data/catalogues/pilot_selection.csv")
rows = pd.concat(
    [
        p[p["sample"] == "A"],
        p[(p.survey == "CEERS") & (p["sample"] == "B")].head(1),
        p[(p.survey == "JADES") & (p["sample"] == "B")].head(1),
    ]
)
result = []
if args.all:
    rows = p
output = "results/atlas_image_downloads.json" if args.all else "results/image_downloads.json"
for r in rows.itertuples():
    url = f"https://grizli-cutout.herokuapp.com/thumb?ra={r.ra}&dec={r.dec}&size=3&filters=f150w-clear,f277w-clear,f444w-clear&output=fits_weight"
    try:
        path = cache.fetch(
            url,
            "images/" + r.source_id + ".fits",
            survey=r.survey,
            version="DJA-cutout-live",
            kind="observed_image_cutout",
        )
        result.append(dict(source_id=r.source_id, path=path.name, status="downloaded"))
    except Exception as e:
        result.append(dict(source_id=r.source_id, status="failed", error=str(e)))
    print(result[-1], flush=True)
    write_json(output, result)

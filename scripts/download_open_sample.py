"""Download all prespecified fainter sources, preserving the earlier sample."""

from pathlib import Path

import pandas as pd
import yaml

from cosmic_anomaly_lab.batch_download import download_batch
from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
config = yaml.safe_load((root / "configs/open_search.yaml").read_text())
selection = config["selection"]
frame = pd.read_csv(root / "data/catalogues/observations.csv.gz")
previous = pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz")
eligible = (
    frame.loc[
        (frame.grating == selection["grating"])
        & (frame.grade == selection["quality_grade"])
        & frame.sn50.between(
            selection["minimum_archive_median_pixel_SN"],
            selection["maximum_archive_median_pixel_SN"],
            inclusive="left",
        )
        & frame.z_best.between(*selection["redshift_range"], inclusive="left")
        & ~frame.source_id.isin(previous.source_id)
    ]
    .sort_values(["sn50", "file"], ascending=[False, True])
    .drop_duplicates("source_id")
    .sort_values(["survey", "source_id"])
)
eligible.to_csv(root / "data/catalogues/open_selection.csv.gz", index=False)
write_json(
    root / "results/open_selection.json",
    dict(
        created_utc=utcnow(),
        config_sha256=text_digest(root / "configs/open_search.yaml"),
        parent_source_records=int(frame.source_id.nunique()),
        new_eligible_sources=len(eligible),
        previous_sources=len(previous),
        combined_sources=len(previous) + len(eligible),
        selected_by_survey=eligible.groupby("survey").size().to_dict(),
    ),
)
download_batch(
    root,
    [
        dict(
            url=r.url,
            filename="spectra/" + r.file,
            survey=r.survey,
            version="DJA-v4.4",
            kind="observed_spectrum",
        )
        for r in eligible.itertuples()
    ],
    budget_bytes=config["raw_storage_budget_MiB"] * 1024**2,
    maximum_bytes=config["maximum_product_MiB"] * 1024**2,
    workers=config["download_workers"],
    checkpoint="results/open_downloads.json",
)

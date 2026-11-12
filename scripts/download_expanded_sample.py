"""Download the fixed eligible sample, without stopping on a promising result."""

from pathlib import Path

import pandas as pd
import yaml

from cosmic_anomaly_lab.batch_download import download_batch
from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
config = yaml.safe_load((root / "configs/expanded_search.yaml").read_text())
selection = config["selection"]
frame = pd.read_csv(root / "data/catalogues/observations.csv.gz")
eligible = (
    frame.loc[
        (frame.grating == selection["grating"])
        & (frame.grade == selection["quality_grade"])
        & (frame.sn50 >= selection["minimum_archive_median_pixel_SN"])
        & frame.z_best.between(*selection["redshift_range"], inclusive="left")
    ]
    .sort_values(["sn50", "file"], ascending=[False, True])
    .drop_duplicates("source_id")
)
eligible = eligible.sort_values(["survey", "source_id"])
eligible.to_csv(root / "data/catalogues/expanded_selection.csv.gz", index=False)
write_json(
    root / "results/expanded_selection.json",
    dict(
        created_utc=utcnow(),
        config_sha256=text_digest(root / "configs/expanded_search.yaml"),
        parent_source_records=int(frame.source_id.nunique()),
        parent_observation_records=len(frame),
        eligible_sources=len(eligible),
        selected_by_survey=eligible.groupby("survey").size().to_dict(),
        scope="One high-SN PRISM observation per upstream source, not all JWST sources",
    ),
)
jobs = [
    dict(
        url=r.url,
        filename="spectra/" + r.file,
        survey=r.survey,
        version="DJA-v4.4",
        kind="observed_spectrum",
    )
    for r in eligible.itertuples()
]
download_batch(
    root,
    jobs,
    budget_bytes=config["raw_storage_budget_MiB"] * 1024**2,
    maximum_bytes=config["maximum_product_MiB"] * 1024**2,
    workers=config["download_workers"],
)

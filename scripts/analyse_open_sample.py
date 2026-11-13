"""Measure the fainter extension and scan combined spectra without class cuts."""

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy.cosmology import Planck18
from astropy.io import fits

from cosmic_anomaly_lab.expanded_search import screen_source
from cosmic_anomaly_lab.open_search import (
    band_measurements,
    normalized_shape,
    peer_mask,
    shape_rarity,
)
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum

root = Path(__file__).resolve().parents[1]
destination = root / "results/open_search"
destination.mkdir(parents=True, exist_ok=True)
config = yaml.safe_load((root / "configs/open_search.yaml").read_text())
shape_config = config["shape_search"]
triage = yaml.safe_load((root / "configs/expanded_search.yaml").read_text())["triage"]
new = pd.read_csv(root / "data/catalogues/open_selection.csv.gz")
previous = pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz")
with fits.open(root / "data/raw/instrument/prism_disp.fits") as hdus:
    resolution_wave = np.array(hdus[1].data["WAVELENGTH"])
    resolution_r = np.array(hdus[1].data["R"])
records = []
with gzip.open(destination / "new_source_fits.jsonl.gz", "wt", encoding="utf-8") as stream:
    for i, row in enumerate(new.itertuples(), 1):
        try:
            record, detail = screen_source(
                root / "data/raw/spectra" / row.file, row, resolution_wave, resolution_r, triage
            )
            stream.write(json.dumps(dict(**record, **detail), allow_nan=False) + "\n")
        except (OSError, ValueError, KeyError, np.linalg.LinAlgError) as exc:
            record = dict(
                source_id=row.source_id,
                file=row.file,
                survey=row.survey,
                z_best=float(row.z_best),
                status="FAILED",
                failure_reason=str(exc),
            )
        records.append(record)
        if i % 100 == 0:
            print(f"Fainter diagnostics {i}/{len(new)}", flush=True)
new_measurements = pd.DataFrame(records)
new_measurements.to_csv(destination / "new_measurements.csv", index=False)
combined = pd.concat([pd.read_csv(root / "results/expanded/measurements.csv"), new_measurements])
combined.to_csv(destination / "combined_measurements.csv", index=False)
selection = pd.concat([previous, new]).sort_values("source_id").reset_index(drop=True)
assert selection.source_id.is_unique
band_records = []
for i, row in enumerate(selection.itertuples(), 1):
    path = root / "data/raw/spectra" / row.file
    try:
        spec = read_spectrum(path)
        for mode, (low, high) in shape_config["rest_windows_um"].items():
            edges = np.geomspace(low, high, shape_config["bins_per_window"] + 1)
            measured = band_measurements(
                spec,
                row.z_best,
                edges,
                shape_config["minimum_bin_coverage"],
                shape_config["minimum_pixels_per_bin"],
            )
            band_records.append(
                dict(
                    source_id=row.source_id,
                    file=row.file,
                    survey=row.survey,
                    mode=mode,
                    z_best=float(row.z_best),
                    archive_sn50=float(row.sn50),
                    cosmic_age_Gyr=float(Planck18.age(row.z_best).value),
                    raw_sha256=digest(path),
                    edges_um=edges.tolist(),
                    **measured,
                )
            )
    except (OSError, ValueError, KeyError) as exc:
        for mode in shape_config["rest_windows_um"]:
            band_records.append(
                dict(source_id=row.source_id, mode=mode, status="FAILED", failure_reason=str(exc))
            )
    if i % 250 == 0:
        print(f"Broad spectral shapes {i}/{len(selection)}", flush=True)

rng = np.random.default_rng(shape_config["random_seed"])
scored, details = [], []
for mode in shape_config["rest_windows_um"]:
    available = []
    for item in [r for r in band_records if r["mode"] == mode]:
        if item["status"] == "MEASURED":
            sn = np.asarray(item["flux_uJy"]) / np.sqrt(np.diag(item["covariance_uJy2"]))
            item["median_band_SN"] = float(np.median(sn))
            if item["median_band_SN"] >= shape_config["minimum_median_band_SN"]:
                available.append(item)
            else:
                item["status"] = "LOW_BAND_SN"
    if not available:
        continue
    surveys = np.array([r["survey"] for r in available])
    redshifts = np.array([r["z_best"] for r in available])
    sn = np.array([r["archive_sn50"] for r in available])
    shapes = normalized_shape(np.array([r["flux_uJy"] for r in available]))
    for i, item in enumerate(available):
        mask = peer_mask(
            surveys,
            redshifts,
            sn,
            i,
            shape_config["maximum_normalized_redshift_gap"],
            shape_config["maximum_peer_SN_ratio"],
        )
        item["peer_count"] = int(mask.sum())
        if mask.sum() < shape_config["minimum_peers"]:
            item["status"] = "INSUFFICIENT_MATCHED_PEERS"
            continue
        rarity = shape_rarity(
            item["flux_uJy"],
            item["covariance_uJy2"],
            shapes[mask],
            rng,
            floor=shape_config["empirical_scatter_floor"],
            draws=shape_config["perturbations"],
        )
        item["status"] = "SCORED"
        record = {
            key: item[key]
            for key in [
                "source_id",
                "file",
                "survey",
                "mode",
                "z_best",
                "archive_sn50",
                "cosmic_age_Gyr",
                "median_band_SN",
                "peer_count",
            ]
        }
        record.update(
            {
                key: rarity[key]
                for key in [
                    "score",
                    "perturbation_q10",
                    "perturbation_q90",
                    "dominant_bin",
                ]
            }
        )
        record["dominant_rest_um"] = float(
            np.sqrt(
                item["edges_um"][rarity["dominant_bin"]]
                * item["edges_um"][rarity["dominant_bin"] + 1]
            )
        )
        scored.append(record)
        details.append(
            dict(
                **item, **rarity, peer_ids=[available[j]["source_id"] for j in np.flatnonzero(mask)]
            )
        )
scores = pd.DataFrame(scored)
scores["rank_percentile"] = scores.groupby("mode").perturbation_q10.rank(pct=True)
scores["previous_joint_priority"] = scores.source_id.isin(
    combined[combined.joint_followup_priority.fillna(False).astype(bool)].source_id
)
scores.to_csv(destination / "shape_scores.csv", index=False)
queue = (
    scores.sort_values(
        ["rank_percentile", "perturbation_q10", "source_id"], ascending=[False, False, True]
    )
    .drop_duplicates("source_id")
    .head(shape_config["inspection_count"])
)
queue.to_csv(destination / "shape_followup.csv", index=False)
for name, rows in [("band_measurements", band_records), ("shape_details", details)]:
    with gzip.open(destination / f"{name}.jsonl.gz", "wt", encoding="utf-8") as stream:
        for record in rows:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
statuses = pd.DataFrame(band_records).groupby(["mode", "status"]).size()
summary = dict(
    executed_utc=utcnow(),
    new_attempted=len(new),
    new_analysed=int((new_measurements.status == "ANALYSED").sum()),
    new_failed=int((new_measurements.status == "FAILED").sum()),
    new_joint_priorities=int(new_measurements.joint_followup_priority.fillna(False).sum()),
    combined_sources=len(selection),
    combined_analysed=int((combined.status == "ANALYSED").sum()),
    parent_sources=20144,
    fraction_of_parent=len(selection) / 20144,
    combined_joint_priorities=int(combined.joint_followup_priority.fillna(False).sum()),
    scored_source_mode_pairs=len(scores),
    unique_shape_scored_sources=int(scores.source_id.nunique()),
    shape_followup_sources=len(queue),
    shape_queue_overlap_joint_screen=int(queue.previous_joint_priority.sum()),
    shape_status_counts=[dict(mode=m, status=s, count=int(n)) for (m, s), n in statuses.items()],
    physical_classification="NOT_ESTABLISHED",
    cosmology="Astropy Planck18: flat LambdaCDM, H0=67.66 km/s/Mpc, Om0=0.30966",
    input_sha256={
        str(p.relative_to(root)).replace("\\", "/"): digest(p)
        if p.suffix == ".gz"
        else text_digest(p)
        for p in [
            root / "configs/open_search.yaml",
            root / "configs/expanded_search.yaml",
            root / "data/catalogues/open_selection.csv.gz",
            root / "data/catalogues/expanded_selection.csv.gz",
            root / "results/expanded/measurements.csv",
            root / "src/cosmic_anomaly_lab/open_search.py",
            Path(__file__),
            root / "src/cosmic_anomaly_lab/spectroscopy.py",
            root / "src/cosmic_anomaly_lab/expanded_search.py",
        ]
    },
)
write_json(destination / "summary.json", summary)
print(json.dumps(summary, indent=2), flush=True)

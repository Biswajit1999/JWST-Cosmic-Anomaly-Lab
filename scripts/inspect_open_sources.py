"""Inspect the fixed shape queue, including falsification-oriented stress checks."""

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from matplotlib.figure import Figure

from cosmic_anomaly_lab.open_search import band_measurements, normalized_shape, shape_rarity
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import (
    evidence_sheet,
    image_panel,
    read_images,
    save_figure,
)

root = Path(__file__).resolve().parents[1]
out = root / "results/open_search"
queue = pd.read_csv(out / "shape_followup.csv")
selection = pd.concat(
    [
        pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz"),
        pd.read_csv(root / "data/catalogues/open_selection.csv.gz"),
    ]
)
observations = pd.read_csv(root / "data/catalogues/observations.csv.gz")
config = yaml.safe_load((root / "configs/open_search.yaml").read_text())["shape_search"]
with gzip.open(out / "shape_details.jsonl.gz", "rt", encoding="utf-8") as stream:
    details = {(r["source_id"], r["mode"]): r for r in map(json.loads, stream)}
with gzip.open(out / "band_measurements.jsonl.gz", "rt", encoding="utf-8") as stream:
    bands = list(map(json.loads, stream))
lookup = {(r["source_id"], r["mode"]): r for r in bands}
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
rng = np.random.default_rng(20261113)
audit, prepared = [], []
shape_fig = Figure(figsize=(15, 12), layout="constrained")
shape_axes = shape_fig.subplots(4, 3).flat
for candidate, ax in zip(queue.itertuples(), shape_axes, strict=True):
    row = next(selection[selection.source_id == candidate.source_id].itertuples())
    item = details[(candidate.source_id, candidate.mode)]
    peers = normalized_shape(
        np.array([lookup[(sid, candidate.mode)]["flux_uJy"] for sid in item["peer_ids"]])
    )
    record = dict(
        source_id=row.source_id,
        survey=row.survey,
        source_number=int(row.srcid),
        mode=candidate.mode,
        z_best=float(row.z_best),
        assessment="UNKNOWN_OR_UNRESOLVED",
        physical_models_compared=False,
        novelty_audit="NOT_ESTABLISHED",
        peer_count=item["peer_count"],
        score=item["score"],
        score_without_dominant_bin=shape_rarity(
            np.delete(item["flux_uJy"], item["dominant_bin"]),
            np.delete(
                np.delete(item["covariance_uJy2"], item["dominant_bin"], axis=0),
                item["dominant_bin"],
                axis=1,
            ),
            normalized_shape(np.delete(peers, item["dominant_bin"], axis=1)),
            rng,
        )["score"],
    )
    # Fixed peers isolate band-placement sensitivity; +/-1% in (1+z) is a
    # stress test, not an estimate of the archive redshift uncertainty.
    spec = read_spectrum(root / "data/raw/spectra" / row.file)
    shifts = []
    for delta in [-0.01, 0.01]:
        z = row.z_best + delta * (1 + row.z_best)
        measured = band_measurements(spec, z, item["edges_um"])
        result = dict(fractional_shift_one_plus_z=delta, status=measured["status"])
        if measured["status"] == "MEASURED":
            values = shape_rarity(measured["flux_uJy"], measured["covariance_uJy2"], peers, rng)
            result.update(score=values["score"], perturbation_q10=values["perturbation_q10"])
        shifts.append(result)
    record["redshift_stress_fixed_peers"] = shifts
    other = [
        r
        for r in bands
        if r["mode"] == candidate.mode
        and r["status"] in ["SCORED", "INSUFFICIENT_MATCHED_PEERS"]
        and r["survey"] != row.survey
        and abs(r["z_best"] - row.z_best) / (1 + row.z_best) < 0.15
        and row.sn50 / 2 <= r["archive_sn50"] <= row.sn50 * 2
    ]
    record["other_survey_peer_count"] = len(other)
    record["other_survey_score"] = None
    if len(other) >= 30:
        record["other_survey_score"] = shape_rarity(
            item["flux_uJy"],
            item["covariance_uJy2"],
            normalized_shape(np.array([r["flux_uJy"] for r in other])),
            rng,
        )["score"]
    try:
        image_path = cache.fetch(
            f"https://grizli-cutout.herokuapp.com/thumb?ra={row.ra}&dec={row.dec}&size=3&filters=f150w-clear,f277w-clear,f444w-clear&output=fits_weight",
            "images/" + row.source_id + ".fits",
            survey=row.survey,
            version="DJA-cutout-live",
            kind="observed_image_cutout",
        )
        images, rgb, image_detail = read_images(image_path, row.ra, row.dec)
        record.update(
            image_status="RGB_AVAILABLE" if rgb is not None else "PARTIAL_COVERAGE",
            image_sha256=digest(image_path),
            image_diagnostics=image_detail,
        )
    except Exception as exc:
        images, rgb = {}, None
        record.update(image_status="UNAVAILABLE", image_error=str(exc))
    try:
        record.update(evidence_sheet(root, row, images, rgb, "figures/open_sources"))
    except (OSError, ValueError, KeyError) as exc:
        record["evidence_sheet_error"] = str(exc)
    prepared.append((row, images, rgb))
    alternatives = observations[
        (observations.source_id == row.source_id)
        & (observations.file != row.file)
        & (observations.grating == "PRISM")
        & (observations.grade == 3)
    ].sort_values(["sn50", "file"], ascending=[False, True])
    record["alternate_PRISM_products"] = len(alternatives)
    record["alternate_product_check"] = dict(status="NO_ELIGIBLE_ALTERNATE_PRODUCT")
    alternate_shape = None
    if len(alternatives):
        alt = alternatives.iloc[0]
        try:
            alternate_path = cache.fetch(
                alt.url,
                "spectra/" + alt.file,
                survey=alt.survey,
                version="DJA-v4.4",
                kind="observed_spectrum",
            )
            measured = band_measurements(
                read_spectrum(alternate_path), row.z_best, item["edges_um"]
            )
            check = dict(
                file=alt.file,
                raw_sha256=digest(alternate_path),
                archive_z_best=float(alt.z_best),
                status=measured["status"],
                interpretation="Same upstream source association and DJA reduction; "
                "not proven independent exposures or extraction",
            )
            if measured["status"] == "MEASURED":
                alternative = shape_rarity(
                    measured["flux_uJy"], measured["covariance_uJy2"], peers, rng
                )
                alternate_shape = np.array(alternative["normalized_flux"])
                check.update(
                    score=alternative["score"],
                    normalized_shape_RMS_difference=float(
                        np.sqrt(np.mean((alternate_shape - item["normalized_flux"]) ** 2))
                    ),
                )
            record["alternate_product_check"] = check
        except Exception as exc:
            record["alternate_product_check"] = dict(status="FAILED", reason=str(exc))
    edges = np.array(item["edges_um"])
    centers = np.sqrt(edges[:-1] * edges[1:])
    median, scatter = np.array(item["peer_median"]), np.array(item["peer_scatter"])
    flux = np.array(item["flux_uJy"])
    err = np.sqrt(np.diag(item["covariance_uJy2"])) / np.sqrt(np.mean(flux**2))
    ax.fill_between(
        centers, median - scatter, median + scatter, color="#d5dfe2", label="Peer MAD spread"
    )
    ax.plot(centers, median, color="#6b7b83", lw=1, label="Peer median")
    ax.errorbar(
        centers,
        item["normalized_flux"],
        yerr=err,
        fmt="o-",
        ms=3,
        color="#9d4b2d",
        lw=1,
        label="Target ± stat. error / RMS",
    )
    if alternate_shape is not None:
        ax.plot(centers, alternate_shape, color="#276e89", ls="--", label="Alternate product")
    ax.set(
        xlabel="Rest wavelength (µm)",
        ylabel="fν / RMS(fν)",
        title=f"{row.survey} {row.srcid} · z={row.z_best:.3f}\n"
        f"{candidate.mode} · score {item['score']:.2f} · {item['peer_count']} peers",
    )
    ax.grid(alpha=0.15)
    if len(audit) == 0:
        ax.legend(fontsize=7)
    audit.append(record)
    write_json(
        out / "inspection_audit.json",
        dict(
            executed_utc=utcnow(),
            sources=audit,
            script_sha256=text_digest(Path(__file__)),
            selection_sha256=text_digest(out / "shape_followup.csv"),
        ),
    )
    print(f"Inspected {len(audit)}/{len(queue)}: {row.source_id}", flush=True)
shape_fig.suptitle(
    "OBSERVED spectral-shape outliers · empirical comparison, not physical-model rejection",
    fontsize=14,
)
save_figure(shape_fig, root / "figures/open_shape_comparison")
atlas = Figure(figsize=(13, 15), layout="constrained")
for ax, (row, images, rgb) in zip(atlas.subplots(4, 3).flat, prepared, strict=True):
    image_panel(ax, images, rgb, row, None if rgb is not None else "F444W-CLEAR")
    ax.set_title(f"{row.survey} {row.srcid} · z={row.z_best:.3f}\n{row.source_id}", fontsize=9)
atlas.suptitle(
    "OBSERVED NIRCam images of the shape-follow-up queue\n"
    "R/G/B = F444W/F277W/F150W; F444W alone if RGB unavailable.\n"
    "Catalogue positions circled; 1″ bars; independent stretches, unmatched PSFs. No new class established.",
    fontsize=12,
)
save_figure(atlas, root / "figures/open_rgb_atlas")

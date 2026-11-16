"""Screen every DR4 record in two apertures; select a bounded diagnostic queue."""

from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_census import (
    LINES,
    aperture_minimum,
    contrast_score,
    normalize_filter,
)
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
config_path = root / "configs/jades_balmer_census.yaml"
config = yaml.safe_load(config_path.read_text())
out = root / "results/jades_balmer_census"
out.mkdir(parents=True, exist_ok=True)
cache = DownloadCache(root, budget_bytes=config["raw_storage_budget_MiB"] * 1024**2)
readme = cache.fetch(
    "https://jades.herts.ac.uk/DR4/Readme_DR4_catalogues.md",
    "metadata/jades_dr4_catalogues.md", survey="JADES", version="DR4",
    kind="catalogue_documentation",
)
catalogue = root / config["catalogue"]
records, invalid_labels = [], []
with fits.open(catalogue) as hdus:
    obs = hdus["Obs_info"].data
    ids = np.array(obs["Unique_ID"]).astype(str)
    assert len(set(ids)) == len(ids)
    for aperture in config["apertures"]:
        table = hdus[f"R1000_{aperture}pix"]
        table_ids = np.array(table.data["Unique_ID"]).astype(str)
        assert len(set(table_ids)) == len(table_ids) and set(table_ids) == set(ids)
        index = {sid: i for i, sid in enumerate(table_ids)}
        for i, sid in enumerate(ids):
            row = table.data[index[sid]]
            record = dict(
                dr4_id=sid, aperture_pixels=aperture, pid=int(obs[i]["PID"]),
                ra=float(obs[i]["RA_TARG"]), dec=float(obs[i]["Dec_TARG"]),
                z=float(obs[i]["z_Spec"]), z_flag=str(obs[i]["z_Spec_flag"]),
                field=str(obs[i]["Field"]),
            )
            for short, name in LINES.items():
                unit = u.Unit(table.columns[name + "_flux"].unit)
                assert unit.is_equivalent(u.erg / u.s / u.cm**2)
                scale = (1 * unit).to_value(u.erg / u.s / u.cm**2) / 1e-20
                flux = float(row[name + "_flux"]) * scale
                error = float(row[name + "_flux_err"]) * scale
                archive_snr = float(np.asarray(row[name + "_SNR"]).ravel()[0])
                valid = np.isfinite(flux) and np.isfinite(error) and error > 0
                label = normalize_filter(row[name + "_filter"])
                if label is None and valid:
                    invalid_labels.append(dict(dr4_id=sid, aperture=aperture, line=short,
                                               raw_value_repr=repr(row[name + "_filter"])))
                sn = flux / error if valid else np.nan
                record.update({
                    short + "_flux": flux, short + "_error": error,
                    short + "_flux_over_error": sn, short + "_archive_SNR": archive_snr,
                    short + "_filter": label,
                    short + "_usable": bool(valid and label is not None
                        and sn >= config["minimum_flux_over_error"]
                        and archive_snr >= config["minimum_archive_SNR"]),
                })
            for line, reference in config["reference_ratios"].items():
                eligible = (record[line + "_usable"] and record["hb_usable"]
                            and record["z_flag"] in config["redshift_flags"]
                            and np.isfinite(record["z"]) and record["z"] > 0)
                record[line + "_eligible"] = bool(eligible)
                record[line + "_over_hb"] = (
                    record[line + "_flux"] / record["hb_flux"] if eligible else np.nan
                )
                record[line + "_same_grating"] = bool(
                    eligible and record[line + "_filter"].removesuffix("_ext")
                    == record["hb_filter"].removesuffix("_ext"))
                for fraction in config["fractional_error_sensitivity"]:
                    for rho in config["statistical_correlation_sensitivity"]:
                        key = f"{line}_score_f{fraction:.2f}_rho{rho:.1f}"
                        score = contrast_score(record[line + "_flux"], record["hb_flux"],
                            record[line + "_error"], record["hb_error"], reference,
                            fraction, rho)
                        record[key] = float((-1 if line == "ha" else 1) * score) if eligible else np.nan
            records.append(record)
frame = pd.DataFrame(records)
frame.to_csv(out / "aperture_measurements.csv.gz", index=False)
write_json(out / "invalid_filter_labels.json", invalid_labels)

sources = frame[frame.aperture_pixels == 5][["dr4_id", "pid", "ra", "dec", "z", "z_flag", "field"]].copy()
sources = sources.set_index("dr4_id")
sensitivity = []
for line in config["reference_ratios"]:
    for aperture in config["apertures"]:
        by_id = frame[frame.aperture_pixels == aperture].set_index("dr4_id").reindex(sources.index)
        for suffix in ["over_hb", "eligible", "same_grating", "filter"]:
            sources[f"{line}_{suffix}_{aperture}px"] = by_id[f"{line}_{suffix}"]
    for fraction in config["fractional_error_sensitivity"]:
        for rho in config["statistical_correlation_sensitivity"]:
            key = f"{line}_score_f{fraction:.2f}_rho{rho:.1f}"
            pivot = frame.pivot(index="dr4_id", columns="aperture_pixels", values=key).reindex(sources.index)
            conservative = aperture_minimum(pivot.to_numpy().T)
            sensitivity.append(dict(line=line, fractional_error=fraction, statistical_correlation=rho,
                eligible_both_apertures=int(np.isfinite(conservative).sum()),
                both_scores_above_3=int((conservative > 3).sum()),
                interpretation="Uncalibrated diagnostic threshold, not detections"))
            if fraction == config["ranking_fractional_error"] and rho == config["ranking_statistical_correlation"]:
                sources[line + "_ranking_score"] = conservative

previous = pd.concat([pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz"),
                      pd.read_csv(root / "data/catalogues/open_selection.csv.gz")])
source_coords = SkyCoord(sources.ra.to_numpy() * u.deg, sources.dec.to_numpy() * u.deg)
previous_coords = SkyCoord(previous.ra.to_numpy() * u.deg, previous.dec.to_numpy() * u.deg)
ii, jj, sep, _ = search_around_sky(source_coords, previous_coords,
    config["positional_neighbour_radius_arcsec"] * u.arcsec)
links = pd.DataFrame(dict(dr4_id=sources.index.to_numpy()[ii],
    previous_source_id=previous.source_id.to_numpy()[jj], separation_arcsec=sep.arcsec))
links.to_csv(out / "previous_positional_neighbours.csv", index=False)
sources["previous_neighbour_count"] = links.groupby("dr4_id").size().reindex(sources.index, fill_value=0)
sources.to_csv(out / "source_diagnostics.csv", index=True)
queues = []
for line in config["reference_ratios"]:
    selected = sources[sources[line + "_ranking_score"] > 0].sort_values(
        [line + "_ranking_score", "dr4_id"], ascending=[False, True]).head(config["top_per_direction"]).copy()
    selected["direction"] = "low_ha_over_hb" if line == "ha" else "high_hg_over_hb"
    selected["rank"] = np.arange(1, len(selected) + 1)
    queues.append(selected.reset_index())
queue = pd.concat(queues, ignore_index=True)
queue.to_csv(out / "followup_queue.csv", index=False)
pd.DataFrame(sensitivity).to_csv(out / "sensitivity.csv", index=False)

fig = Figure(figsize=(13, 6), layout="constrained")
axes = fig.subplots(1, 2)
for ax, line, label in zip(axes, ["ha", "hg"], ["Hα/Hβ", "Hγ/Hβ"], strict=True):
    valid = sources[f"{line}_eligible_3px"] & sources[f"{line}_eligible_5px"]
    data = sources[valid]
    scatter = ax.scatter(data[f"{line}_over_hb_3px"], data[f"{line}_over_hb_5px"],
                         c=data.z, s=9, alpha=0.5, cmap="viridis")
    picked = queue[queue.direction == ("low_ha_over_hb" if line == "ha" else "high_hg_over_hb")]
    ax.scatter(picked[f"{line}_over_hb_3px"], picked[f"{line}_over_hb_5px"],
               facecolors="none", edgecolors="#c74b32", s=80, label="Six review priorities")
    limits = [min(data[f"{line}_over_hb_3px"].min(), data[f"{line}_over_hb_5px"].min()) * .85,
              max(data[f"{line}_over_hb_3px"].max(), data[f"{line}_over_hb_5px"].max()) * 1.15]
    ax.plot(limits, limits, color="gray", lw=.8)
    ref = config["reference_ratios"][line]
    ax.axvline(ref, color="#976540", ls="--", lw=.8)
    ax.axhline(ref, color="#976540", ls="--", lw=.8)
    ax.set(xscale="log", yscale="log", xlim=limits, ylim=limits,
           xlabel=f"{label}: 3-pixel extraction", ylabel=f"{label}: 5-pixel extraction",
           title=f"{len(data)} eligible records; shared photons")
    ax.legend(fontsize=8)
    fig.colorbar(scatter, ax=ax, label="Catalogue redshift")
fig.suptitle("JADES DR4 Balmer census · catalogue fits, not broad-line measurements\nDashed lines: fixed reference ratios, not physical exclusion boundaries", fontsize=12)
save_figure(fig, root / "figures/jades_balmer_apertures")

fig = Figure(figsize=(10, 5), layout="constrained")
ax = fig.subplots()
for line, label, color in [("ha", "Low Hα/Hβ", "#277e8e"), ("hg", "High Hγ/Hβ", "#b95739")]:
    s = pd.DataFrame(sensitivity)
    s = s[s.line == line]
    center = s[s.statistical_correlation == 0].sort_values("fractional_error")
    grouped = s.groupby("fractional_error").both_scores_above_3
    ax.fill_between(grouped.min().index * 100, grouped.min(), grouped.max(), color=color, alpha=.15)
    ax.plot(center.fractional_error * 100, center.both_scores_above_3, "o-", label=label, color=color)
ax.set(xlabel="Assumed independent fractional uncertainty per line (%)",
       ylabel="Records above diagnostic score 3 in BOTH apertures",
       title="Sensitivity of review flags, not discovery significance\nBands vary statistical cross-line correlation from −0.5 to +0.5")
ax.legend()
save_figure(fig, root / "figures/jades_balmer_sensitivity")
summary = dict(executed_utc=utcnow(), research_snapshot=config["research_snapshot"],
    catalogue_records=len(sources), aperture_records=len(frame),
    redshift_flag_counts=sources.z_flag.value_counts().to_dict(),
    both_apertures_eligible={line: int(np.isfinite(sources[line + "_ranking_score"]).sum()) for line in LINES if line != "hb"},
    invalid_filter_labels_with_valid_flux=len(invalid_labels),
    positional_neighbours_previous_records=int((sources.previous_neighbour_count > 0).sum()),
    no_neighbour_previous_records=int((sources.previous_neighbour_count == 0).sum()),
    queue_entries=len(queue), queue_distinct_dr4_ids=int(queue.dr4_id.nunique()),
    no_unique_object_union_claim=True, physical_classification="NOT_ESTABLISHED",
    inputs_sha256={config_path.relative_to(root).as_posix(): text_digest(config_path),
        config["catalogue"]: digest(catalogue),
        readme.relative_to(root).as_posix(): digest(readme),
        "scripts/analyse_jades_balmer_census.py": text_digest(Path(__file__)),
        "src/cosmic_anomaly_lab/balmer_census.py": text_digest(root / "src/cosmic_anomaly_lab/balmer_census.py"),
        **{f"data/catalogues/{name}_selection.csv.gz": digest(root / f"data/catalogues/{name}_selection.csv.gz") for name in ["expanded", "open"]}},
)
write_json(out / "summary.json", summary)
print(summary)
print(queue[["dr4_id", "z", "direction", "ha_ranking_score", "hg_ranking_score"]].to_string(index=False))

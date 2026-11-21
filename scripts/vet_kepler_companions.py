"""Crossmatch residual peaks and report event sensitivity for the replicated ones."""

import json
from pathlib import Path

import pandas as pd

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths
from cosmic_anomaly_lab.transits import phase_distance

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_companions"
roster = {r["kic"]: r for r in json.loads((out / "selection.json").read_text())["targets"]}
inputs = {}
tables = {}
for name in ["recent_koi", "tce"]:
    receipt = json.loads((out / (name+"_receipt.json")).read_text())
    assert digest(root / receipt["path"]) == receipt["sha256"]
    inputs[name] = receipt["sha256"]
    tables[name] = json.loads((root / receipt["path"]).read_text())
    assert len(tables[name]) == receipt["rows"]
peaks = pd.read_csv(out / "candidates.csv").query("sign == 'observed'")
measurements = pd.read_csv(out / "events.csv").query("sign == 'observed' and status == 'MEASURED'")
records = []
for row in peaks.itertuples():
    source = roster[row.kic]
    observed = measurements.query("kic == @row.kic and rank == @row.rank")
    matches = []
    for label, catalogue in [("pinned_koi", source["cumulative_koi"]), *tables.items()]:
        prefix = "tce_" if label == "tce" else "koi_"
        for obj in (o for o in catalogue if o["kepid"] == row.kic):
            period, epoch = obj.get(prefix+"period"), obj.get(prefix+"time0bk")
            if not period or epoch is None:
                continue
            duration = obj.get(prefix+"duration")
            tolerance = max(row.duration_days/2, duration/48 if duration else 0) + .05
            for numerator, denominator in [(1, 4), (1, 3), (1, 2), (1, 1), (2, 1), (3, 1), (4, 1)]:
                harmonic = numerator/denominator
                fractional_error = row.period_days/(period*harmonic)-1
                if abs(fractional_error) > .005:
                    continue
                phases = abs(phase_distance(observed.midpoint.to_numpy(), period/denominator, epoch))
                fraction = float((phases <= tolerance).mean()) if len(phases) else None
                matches.append(dict(catalogue=label, harmonic=harmonic, relative_period_error=fractional_error,
                    timing_overlap_fraction=fraction, measured_events=len(phases), reference=obj))
    timing_matches = [m for m in matches if m["timing_overlap_fraction"] is not None
                      and m["timing_overlap_fraction"] >= .5]
    leave_event_contrast, parity = [], []
    if row.replicated:
        later = observed.query("split == 'later'")
        for cycle in later.cycle:
            z = combine_event_depths(later[later.cycle != cycle].to_dict("records"))["contrast"]
            if z is not None:
                leave_event_contrast.append(z)
        for label in [0, 1]:
            r = combine_event_depths(later[later.cycle % 2 == label].to_dict("records"))
            parity.append(r)
    if not row.replicated:
        status = "BELOW_REPLICATION_GATE"
    elif timing_matches:
        status = "CATALOGUE_MATCH_OR_ALIAS"
    elif matches:
        status = "POSSIBLE_CATALOGUE_ALIAS"
    elif not row.wide_mask_replicated:
        status = "SENSITIVE_TO_KNOWN_TRANSIT_MASK"
    else:
        status = "UNMATCHED_REPLICATION_REQUIRES_VETTING"
    records.append(dict(kic=int(row.kic), role=row.role, rank=int(row.rank), period_days=row.period_days,
        status=status, min_later_contrast_excluding_one_event=min(leave_event_contrast, default=None),
        later_parity=parity, matches=matches))
write_json(out / "vetting.json", records)
table = pd.DataFrame(records).drop(columns=["later_parity", "matches"])
table.to_csv(out / "vetting.csv", index=False)
write_json(out / "vetting_receipt.json", dict(executed_utc=utcnow(), catalogue_sha256=inputs,
    script_sha256=text_digest(__file__), fitting_precedes_crossmatch=True, new_planets_confirmed=0,
    candidate_table_sha256=text_digest(out / "candidates.csv")))
print(table[table.status.ne("BELOW_REPLICATION_GATE")].to_string(index=False))

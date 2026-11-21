"""Compare fitted periods and transit phases with catalogue entries after the search."""

import json
from pathlib import Path

import pandas as pd

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transits import phase_distance

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_long_period"
summary = json.loads((out / "search_summary.json").read_text())
assert not summary["partial_run"]
receipts = json.loads((out / "catalogues.json").read_text())
catalogues = {}
for label in ["koi", "tce"]:
    record = receipts[label]
    path = root / record["path"]
    assert digest(path) == record["sha256"]
    catalogue = json.loads(path.read_text())
    assert len(catalogue) == record["rows"]
    catalogues[label] = catalogue

records = []
table = pd.read_csv(out / "search.csv")
for row in table[table.sign == "observed"].itertuples():
    if row.status != "EXECUTED":
        continue
    matched, coordinates = [], []
    for label, catalogue in catalogues.items():
        prefix = label + "_"
        for obj in (o for o in catalogue if o["kepid"] == row.kic):
            period, epoch = obj.get(prefix+"period"), obj.get(prefix+"time0bk")
            if period is None or period <= 0:
                continue
            coordinates.append(obj)
            for harmonic in [.25, 1/3, .5, 1, 2, 3, 4]:
                delta = row.period_days / (period*harmonic) - 1
                if abs(delta) <= .005:
                    phase = abs(float(phase_distance(row.epoch_bkjd, period*min(1, harmonic), epoch))) if epoch is not None else None
                    matched.append(dict(table=label, reference_period_days=period, harmonic=harmonic,
                                        relative_period_difference=delta, phase_distance_days=phase,
                                        same_epoch=phase is not None and phase <= max(.3, row.duration_days),
                                        reference=obj))
    same_epoch = any(m["same_epoch"] for m in matched)
    if not row.primary_replicated:
        status = "DID_NOT_REPLICATE"
    elif same_epoch:
        status = "CATALOGUE_SIGNAL_RECOVERY"
    elif matched:
        status = "PERIOD_ALIAS_REQUIRES_REVIEW"
    else:
        status = "NO_PERIOD_MATCH_REQUIRES_VETTING"
    records.append(dict(kic=int(row.kic), role=row.role, period_days=row.period_days,
                        primary_replicated=bool(row.primary_replicated), assessment=status,
                        same_target_catalogue_entries=len(coordinates), matches=matched))
write_json(out / "crossmatch.json", records)
write_json(out / "crossmatch_receipt.json", dict(executed_utc=utcnow(),
    catalogue_sha256={label: receipts[label]["sha256"] for label in catalogues},
    search_sha256=text_digest(out / "search.csv"), script_sha256=text_digest(Path(__file__)),
    missing_match_implies_novelty=False, validated_new_planets=0))
report = pd.DataFrame(records).drop(columns=["matches"])
report.to_csv(out / "crossmatch.csv", index=False)
print(report.to_string(index=False))

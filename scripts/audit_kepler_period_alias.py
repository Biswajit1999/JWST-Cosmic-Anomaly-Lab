"""Test the leading non-control peak by removing catalogue-matched transits."""

import json

import pandas as pd
from matplotlib.figure import Figure
from search_kepler_long_period import events_at, gate, load_curve, out, root, search

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths
from cosmic_anomaly_lab.transits import phase_distance
from cosmic_anomaly_lab.visual_evidence import save_figure

initial = pd.read_csv(out / "search.csv")
rank = initial.query("sign == 'observed' and role != 'confirmed_control'").dropna(subset=["replication_contrast"])
lead = rank.loc[rank.replication_contrast.idxmax()]
kic = int(lead.kic)
assessment = next(r for r in json.loads((out / "crossmatch.json").read_text()) if r["kic"] == kic)
match = next(m for m in assessment["matches"] if m["table"] == "koi" and m["same_epoch"])
reference = match["reference"]
period, epoch, duration = (reference["koi_period"], reference["koi_time0bk"], reference["koi_duration"]/24)
known = dict(period_days=period, epoch_bkjd=epoch, duration_days=duration)
data = json.loads((out / "products.json").read_text())
target = next(row for row in data["targets"] if row["kic"] == kic)
quarters = sorted([load_curve(p, kic) for p in target["products"]], key=lambda c: c["quarter"])
halves = {
    "discovery": [c for c in quarters if 1 <= c["quarter"] <= 8 and len(c["time"]) > 100],
    "replication": [c for c in quarters if c["quarter"] >= 9 and len(c["time"]) > 100],
}
masked, removed = {}, {}
for label, curves in halves.items():
    masked[label], removed[label] = [], 0
    for c in curves:
        # Mask the known event plus at least a half-duration margin on either side.
        keep = abs(phase_distance(c["time"], period, epoch)) > duration
        removed[label] += int((~keep).sum())
        masked[label].append(c | {key: c[key][keep] for key in ["time", "flux", "error"]})
new_fit = search(masked["discovery"])
result, append = {}, []
for label, curves, fit in [("original_peak", halves, lead), ("catalogue_planet", halves, known),
                           ("original_timing_after_mask", masked, lead), ("new_peak_after_mask", masked, new_fit)]:
    result[label] = dict(period_days=float(fit["period_days"]), epoch_bkjd=float(fit["epoch_bkjd"]))
    for split, half in curves.items():
        events = events_at(half, fit)
        result[label][split] = combine_event_depths(events)
        append.extend(dict(model=label, split=split, **r) for r in events)
    result[label]["replicated"] = all(gate(result[label][split]) for split in halves)

original_events = pd.DataFrame(append).query("model == 'original_peak' and status == 'MEASURED'").copy()
offsets = abs(phase_distance(original_events.midpoint, period, epoch))
original_events["distance_to_catalogue_transit_days"] = offsets
original_events.to_csv(out / "alias_overlap.csv", index=False)
pd.DataFrame(append).to_csv(out / "alias_audit_events.csv", index=False)
fig = Figure(figsize=(11, 5), layout="constrained")
axes = fig.subplots(1, 2)
all_curves = halves["discovery"] + halves["replication"]
for ax, fit, heading in zip(axes, [known, lead], ["Published 10.93-day planet", "32.79-day search peak"], strict=True):
    for curve in all_curves:
        phase = phase_distance(curve["time"], fit["period_days"], fit["epoch_bkjd"])
        use = abs(phase) <= .6
        ax.plot(phase[use]*24, (curve["flux"][use]-1)*1e6, ".", ms=2, color="#284c68", alpha=.23)
    ax.axvspan(-fit["duration_days"]*12, fit["duration_days"]*12, color="#dc8b3a", alpha=.22)
    ax.set(title=heading, xlabel="Hours from predicted transit", ylabel="Observed relative flux (ppm)")
fig.suptitle("Kepler-1437 · the search peak selects every third known transit", fontsize=12)
save_figure(fig, root / "figures/kepler_long_period/kepler1437_alias")
write_json(out / "alias_audit.json", dict(executed_utc=utcnow(), kic=kic, reference=reference,
    original_sampled_events=len(offsets), within_catalogue_duration=int((offsets < duration).sum()),
    largest_offset_minutes=float(offsets.max()*1440), removed_cadences=removed, variants=result,
    dataset_preserved=True, mask_half_width_days=duration,
    experiment="Post-search test prompted by a locked period and phase match",
    script_sha256=text_digest(__file__), search_sha256=text_digest(out / "search.csv")))
print(json.dumps(result, indent=2))

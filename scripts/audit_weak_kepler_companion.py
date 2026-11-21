"""Check the highest-ranked weak residual at its original, fixed timing."""

import json

import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.figure import Figure
from search_kepler_companions import events_at, out, root, split_curves

from cosmic_anomaly_lab.companion_search import prepare_companion_quarter
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths
from cosmic_anomaly_lab.visual_evidence import save_figure

peaks = pd.read_csv(out / "candidates.csv").query("role == 'search_host' and sign == 'observed'")
weight = peaks[["early_contrast", "later_contrast"]].min(axis=1, skipna=False)
fit = peaks.loc[weight.idxmax()]
kic = int(fit.kic)
star = next(s for s in json.loads((out / "selection.json").read_text())["targets"] if s["kic"] == kic)
target = next(s for s in json.loads((out / "products.json").read_text())["targets"] if s["kic"] == kic)
rows, masks, phase_stats = [], [], []
records = []
for column in ["PDCSAP_FLUX", "SAP_FLUX"]:
    for width in [1, 1.5, 2]:
        curves = []
        for p in target["products"]:
            path = root / p["path"]
            assert digest(path) == p["sha256"]
            with fits.open(path) as h:
                assert h[0].header["KEPLERID"] == kic and h[1].header["BJDREFI"] == 2454833
                assert h[1].header["TIMESYS"] == "TDB"
                cadence = h[1].header["TIMEDEL"]
                d = h[1].data
                t, y, e, masked, valid = prepare_companion_quarter(
                    d["TIME"], d[column], d[column+"_ERR"], d["SAP_QUALITY"], star["masks"], cadence,
                    mask_multiplier=width)
                q = int(h[0].header["QUARTER"])
                curves.append(dict(time=t, flux=y, error=e, cadence=cadence, quarter=q))
                masks.append(dict(quarter=q, photometry=column, mask_multiplier=width,
                                  masked=masked, quality_cadences=valid, retained=len(t)))
        halves = split_curves(sorted(curves, key=lambda c: c["quarter"]))
        for label, group in halves.items():
            events = events_at(group, fit)
            rows.extend(dict(split=label, photometry=column, mask_multiplier=width, **r) for r in events)
            records.append(dict(split=label, photometry=column, mask_multiplier=width,
                                **combine_event_depths(events)))
        if width == 1:
            for phase in np.linspace(.02, .98, 97):
                events = events_at(halves["later"], fit, phase=phase)
                phase_stats.append(dict(photometry=column, phase=float(phase), **combine_event_depths(events)))

pd.DataFrame(rows).to_csv(out / "weak_lead_events.csv", index=False)
pd.DataFrame(masks).to_csv(out / "weak_lead_mask_counts.csv", index=False)
pd.DataFrame(phase_stats).to_csv(out / "weak_lead_shifted_phases.csv", index=False)
df = pd.DataFrame(records)
df.to_csv(out / "weak_lead_sensitivity.csv", index=False)
fig = Figure(figsize=(11, 4.8), layout="constrained")
ax, control = fig.subplots(1, 2)
points = pd.DataFrame(rows).query("mask_multiplier == 1 and photometry == 'PDCSAP_FLUX' and status == 'MEASURED'")
for label, color in [("early", "#bb6c36"), ("later", "#366488")]:
    part = points[points.split == label]
    ax.errorbar(part.midpoint, part.depth*1e6, yerr=part.error*1e6, fmt=".", ms=4, lw=.6,
                alpha=.64, color=color, label=label)
ax.axhline(0, color="gray", lw=.8)
ax.set(xlabel="Predicted transit midpoint (BKJD)", ylabel="Signed depth (ppm)",
       title="Observed event fits at the original timing")
ax.legend()
for column, color in [("PDCSAP_FLUX", "#366488"), ("SAP_FLUX", "#bb6c36")]:
    part = df[(df.photometry == column) & (df.split == "later")]
    control.errorbar(part.mask_multiplier, part.depth*1e6, yerr=part.error*1e6, fmt="o-", lw=1,
                     color=color, label=column.replace("_FLUX", ""))
control.axhline(0, color="gray", lw=.8)
control.set(xlabel="Known-transit mask width multiplier", ylabel="Later-quarter depth (ppm)",
            title="Mask and photometry sensitivity")
control.legend()
fig.suptitle(f"Kepler-1850 · weak residual at {fit.period_days:.5f} days · below replication gate", fontsize=11)
save_figure(fig, root / "figures/kepler_companions/weak_lead")
write_json(out / "weak_lead_receipt.json", dict(executed_utc=utcnow(), kic=kic, rank=int(fit["rank"]),
    period_days=float(fit.period_days), epoch_bkjd=float(fit.epoch_bkjd), duration_days=float(fit.duration_days),
    selection="Largest minimum of early and later contrasts among observed search-host peaks",
    fixed_ephemeris=True, followup_is_blind=False,
    caution="SAP and PDCSAP share photons. Phase offsets share data and are not a false-alarm calibration.",
    qualitative_status="BELOW_PRIMARY_REPLICATION_GATE",
    script_sha256=text_digest(__file__), input_sha256={p["path"]: p["sha256"] for p in target["products"]}))
print(df.to_string(index=False))

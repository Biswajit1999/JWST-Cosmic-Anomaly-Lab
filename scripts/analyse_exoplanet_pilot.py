"""Execute a small unclassified-target Kepler screen with a known control."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.transits import flatten_segments, local_depth, phase_distance, scan_transits
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/exoplanet_pilot"
products = json.loads((out / "products.json").read_text())["products"]
koi = pd.read_json(root / "data/raw/exoplanets/koi_cumulative.json")
records, all_events, quarter_records, hashes, metadata = [], [], [], {}, []
for kid in sorted({p["kic"] for p in products}):
    arrays = []
    for p in [p for p in products if p["kic"] == kid and p["status"] == "DOWNLOADED"]:
        assert digest(root / p["path"]) == p["sha256"]
        hashes[p["path"]] = p["sha256"]
        with fits.open(root / p["path"]) as h:
            assert int(h[0].header["KEPLERID"]) == kid and int(h[0].header["QUARTER"]) == p["quarter"]
            assert h[1].header["BJDREFI"] == 2454833 and h[1].header["TIMESYS"] == "TDB"
            d = h[1].data
            ix, flux, error = flatten_segments(d["TIME"], d["PDCSAP_FLUX"], d["PDCSAP_FLUX_ERR"], d["SAP_QUALITY"])
            arrays.append((np.array(d["TIME"][ix]), flux, error, np.full(len(ix), p["quarter"])))
            metadata.append(dict(kic=kid, quarter=p["quarter"], rows=len(d), retained=len(ix),
                ra=h[0].header.get("RA_OBJ"), dec=h[0].header.get("DEC_OBJ"),
                catalogue_teff_K=h[0].header.get("TEFF"), catalogue_radius_solar=h[0].header.get("RADIUS"),
                crowding=h[1].header.get("CROWDSAP"), role=p["role"]))
    if len(arrays) != 2 or any(len(a[0]) < 100 for a in arrays):
        records.append(dict(kic=kid, status="INSUFFICIENT_TWO_QUARTER_DATA"))
        continue
    t, y, error, quarters = (np.concatenate([a[i] for a in arrays]) for i in range(4))
    order = np.argsort(t)
    t, y, error, quarters = [a[order] for a in (t, y, error, quarters)]
    result, events, power, efferr = scan_transits(t, y, error)
    period, epoch, duration = (result[k] for k in ["period_days", "epoch_bkjd", "duration_days"])
    all_events.extend(dict(kic=kid, **e) for e in events)
    quarter_snr = []
    for q in [2, 3]:
        use = quarters == q
        measurement = local_depth(t[use], y[use], efferr[use], period, epoch, duration)
        sn = measurement["depth"]/measurement["error"] if measurement["status"] == "MEASURED" else None
        quarter_snr.append(sn)
        quarter_records.append(dict(kic=kid, quarter=q, snr=sn, **measurement))
    known = koi[koi.kepid == kid]
    matches = []
    for row in known.itertuples():
        for harmonic in [.5, 1., 2.]:
            if abs(period/(row.koi_period*harmonic)-1) < .01:
                matches.append(dict(koi=row.kepoi_name, disposition=row.koi_disposition,
                                    catalogue_period=float(row.koi_period), harmonic=harmonic))
    screen = (result["depth_snr"] >= 8 and result["observed_events"] >= 3
              and all(s is not None and s >= 5 for s in quarter_snr)
              and result["odd_even_contrast"] < 3 and result["secondary_snr"] < 3
              and result["inverted_max_snr"] < .8*result["depth_snr"])
    role = next(p["role"] for p in products if p["kic"] == kid)
    result.update(kic=kid, role=role, retained_points=len(t), quarters=2, screen_pass=bool(screen),
                  koi_rows=len(known), matched_koi=json.dumps(matches),
                  status="KNOWN_SIGNAL_RECOVERY" if matches else ("FOLLOWUP_REQUIRED" if screen else "SCREEN_NOT_PASSED"),
                  physical_classification="NOT_VALIDATED_BY_THIS_ANALYSIS")
    records.append(result)
    pd.DataFrame(dict(time_bkjd=t, normalized_flux=y, formal_error=error, effective_error=efferr, quarter=quarters)).to_csv(out / f"kic{kid}_processed.csv.gz", index=False)
    fig = Figure(figsize=(11, 7), layout="constrained")
    ax = fig.subplots(2, 2)
    ax[0, 0].plot(t, y, ".", ms=1, alpha=.5)
    ax[0, 0].set(xlabel="BKJD (BJD TDB − 2454833)", ylabel="Normalized PDCSAP", title="Two observed quarters; no downward clipping")
    ax[0, 1].plot(power.period, power.power, lw=.5)
    ax[0, 1].set(xlabel="Trial period (days)", ylabel="BLS diagonal-noise S/N", title="Search ranking; not a false-alarm probability")
    phase = phase_distance(t, period, epoch)/duration
    for q, color in [(2, "#187b9a"), (3, "#bc6a38")]:
        use = (abs(phase)<3) & (quarters == q)
        ax[1, 0].plot(phase[use], y[use], ".", ms=2, alpha=.3, color=color, label=f"Q{q}")
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].set(xlabel="Time from transit / fitted duration", ylabel="Normalized flux", title=f"P={period:.6f} d; depth={result['depth']*1e6:.0f} ppm")
    good = [e for e in events if e["depth"] is not None]
    ax[1, 1].errorbar([e["midpoint_bkjd"] for e in good], [e["depth"]*1e6 for e in good],
                     yerr=[e["error"]*1e6 for e in good], fmt=".", ms=3, lw=.7)
    ax[1, 1].axhline(0, color="gray", lw=.6)
    ax[1, 1].set(xlabel="Event midpoint BKJD", ylabel="Local depth (ppm)", title="Individual events; effective diagonal errors")
    fig.suptitle(f"KIC {kid} · {role} · {result['status']}\nNo new planet is validated by this screening test", fontsize=12)
    save_figure(fig, root / f"figures/exoplanet_pilot/kic{kid}")
    pd.DataFrame(records).to_csv(out / "screening.csv", index=False)
    print(kid, result["status"], result["period_days"], result["depth_snr"], flush=True)
pd.DataFrame(all_events).to_csv(out / "individual_events.csv", index=False)
pd.DataFrame(quarter_records).to_csv(out / "quarter_depths.csv", index=False)
pd.DataFrame(metadata).to_csv(out / "input_metadata.csv", index=False)
write_json(out / "summary.json", dict(executed_utc=utcnow(), targets=len(records), input_sha256=hashes,
    script_sha256=text_digest(Path(__file__)), module_sha256=text_digest(root / "src/cosmic_anomaly_lab/transits.py"),
    koi_sha256=digest(root / "data/raw/exoplanets/koi_cumulative.json"),
    search_scope="Fixed 12-target convenience sample plus TrES-2 control, two quarters, 0.8–20 d",
    global_false_alarm_calibrated=False, centroid_vetting_complete=False,
    new_planets_validated=0, atmospheres_inferred_from_kepler=False))

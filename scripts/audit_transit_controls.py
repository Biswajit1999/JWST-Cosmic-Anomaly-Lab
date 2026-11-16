"""Expose limitations of box-model vetting; targeted known-candidate check."""

from pathlib import Path

import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.transits import local_depth

root = Path(__file__).resolve().parents[1]
out = root / "results/exoplanet_pilot"
events = pd.read_csv(out / "individual_events.csv")
control = events[(events.kic == 11446443) & events.depth.notna()]
means, variances = [], []
for parity in [0, 1]:
    group = control[control.event_number % 2 == parity]
    means.append(float(group.depth.mean()))
    variances.append(float(max(group.depth.var(ddof=1), (group.error**2).mean())/len(group)))
empirical = abs(means[0]-means[1])/np.sqrt(sum(variances))
catalogue = pd.read_json(root / "data/raw/exoplanets/koi_cumulative.json")
known = catalogue[catalogue.kepid == 11446424].iloc[0]
x = pd.read_csv(out / "kic11446424_processed.csv.gz")
p, epoch, duration = known.koi_period, known.koi_time0bk, known.koi_duration/24
targeted = local_depth(x.time_bkjd.to_numpy(), x.normalized_flux.to_numpy(), x.effective_error.to_numpy(), p, epoch, duration)
model = BoxLeastSquares(x.time_bkjd, x.normalized_flux, x.effective_error)
periods = 1/np.linspace(1/90, 1/20, 6000)
scan = model.power(periods, [.1, .16, .24, .36], objective="snr")
best = int(np.argmax(scan.power))
write_json(out / "control_audit.json", dict(executed_utc=utcnow(), script_sha256=text_digest(Path(__file__)),
    control_event_scatter_odd_even_contrast=float(empirical),
    control_original_box_odd_even_contrast=float(pd.read_csv(out / "screening.csv").query("kic==11446443").iloc[0].odd_even_contrast),
    empirical_diagnostic_is_not_validated_replacement=True,
    known_candidate=dict(kic=11446424, koi=known.kepoi_name, disposition=known.koi_disposition,
        period_days=p, duration_days=duration, epoch_bkjd=epoch,
        catalogue_candidate_radius_earth=float(known.koi_prad),
        targeted_measurement=targeted, targeted_snr=targeted["depth"]/targeted["error"],
        extended_scan_max_period=float(scan.period[best]), extended_scan_max_snr=float(scan.power[best]),
        selection="Follow-up after catalogue inspection, not blind discovery",
        validated=False)))
print('Empirical control odd/even contrast',empirical,'known candidate targeted S/N',targeted['depth']/targeted['error'])

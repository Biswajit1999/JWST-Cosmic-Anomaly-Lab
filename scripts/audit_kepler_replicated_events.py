"""Check whether replicated signals hinge on one event or a repeated seasonal window."""

from pathlib import Path

import numpy as np
import pandas as pd

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.transit_followup import combine_event_depths

root = Path(__file__).resolve().parents[1]
out = root / "results/kepler_long_period"
search = pd.read_csv(out / "search.csv")
data = pd.read_csv(out / "events.csv")
targets = search[(search.sign == "observed") & search.primary_replicated]
rows = []
for target in targets.itertuples():
    for split in ["discovery", "replication"]:
        subset = data[(data.kic == target.kic) & (data.sign == "observed") & (data.phase == 0)
                      & (data.split == split) & (data.status == "MEASURED")]
        events = subset.to_dict("records")
        rows.append(dict(kic=target.kic, split=split, exclusion="none", index=None,
                         **combine_event_depths(events)))
        for key, values in [
            ("event", subset.cycle),
            ("quarter", subset.quarter.unique()),
            ("roll_season", np.sort((subset.quarter % 4).unique())),
        ]:
            for index in values:
                if key == "event":
                    keep = subset.cycle != index
                elif key == "quarter":
                    keep = subset.quarter != index
                else:
                    keep = subset.quarter % 4 != index
                rows.append(dict(kic=target.kic, split=split, exclusion=key, index=int(index),
                                 **combine_event_depths(subset[keep].to_dict("records"))))
        for parity in [0, 1]:
            rows.append(dict(kic=target.kic, split=split, exclusion="keep_parity", index=parity,
                             **combine_event_depths(subset[subset.cycle % 2 == parity].to_dict("records"))))
frame = pd.DataFrame(rows)
frame.to_csv(out / "replicated_event_audit.csv", index=False)
write_json(out / "replicated_event_audit.json", dict(executed_utc=utcnow(),
    script_sha256=text_digest(Path(__file__)), event_sha256=text_digest(out / "events.csv"),
    replicating_observed_targets=int(len(targets)), rows=len(rows),
    limitations="Overlapping omissions are sensitivity checks and share observations. Quarter modulo four is a roll-season proxy."))
if len(frame):
    print(frame.groupby(["kic", "split", "exclusion"]).contrast.agg(["min", "max"]).to_string())

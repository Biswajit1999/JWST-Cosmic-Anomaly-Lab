"""Draw and retrieve a detection-independent sample of released G235M spectra."""

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import requests
from analyse_jades_cosmic_evolution import BIN_EDGES, OUT, ROOT
from inspect_jades_cosmic_evolution import spectrum_for_target

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

SEED = 94278
PER_STRATUM = 16
LEDGER = OUT / "forced_sample.json"


def spectrum_receipt(cache, entry):
    try:
        _, receipt = spectrum_for_target(cache, SimpleNamespace(**entry))
    except ValueError as exc:
        # MAST documents both GOODS-S Medium/JWST programme designations under one DR4 tier.
        if (entry["pid"] != 1180 or not entry["dr4_id"].startswith("goods-s-mediumjwst_")
                or not str(exc).endswith("failed the coordinate, programme or ID check")):
            raise
        _, receipt = spectrum_for_target(cache, SimpleNamespace(**(entry | {"pid": 1286})))
        receipt["association"] = "DR4_GOODS_S_MEDIUMJWST_HEADER_1286_CATALOGUE_1180"
    else:
        receipt["association"] = "CATALOGUE_PID_MATCH"
    return receipt


def main():
    targets = OUT / "observed_targets.csv"
    summary = OUT / "summary.json"
    summary_data = json.loads(summary.read_text(encoding="utf8"))
    if summary_data["target_table_sha256"] != text_digest(targets):
        raise ValueError("The parent target table does not match its summary")
    universe = pd.read_csv(targets)
    universe = universe[universe.exposure_s > 0].sort_values("dr4_id").copy()
    universe["redshift_bin"] = pd.cut(universe.z, BIN_EDGES, right=False, labels=False)
    rng = np.random.default_rng(SEED)
    selected, counts = [], []
    cols = ["dr4_id", "nircam_id", "pid", "field", "z", "ra", "dec", "redshift_bin", "m_uv", "exposure_s"]
    for (field, index), group in universe.groupby(["field", "redshift_bin"], sort=True):
        draw = group.iloc[np.sort(rng.choice(len(group), PER_STRATUM, replace=False))]
        selected.extend(dict(row, role="random_sample") for row in json.loads(draw[cols].to_json(orient="records")))
        counts.append(dict(field=field, redshift_bin=int(index), population=len(group), selected=len(draw)))
    selected.extend(dict(row, redshift_bin=row["bin"], role="sideband_control") for row in summary_data["example_targets"])
    cache = DownloadCache(ROOT, 28_000*1024**2)
    ledger = dict(executed_utc=utcnow(), seed=SEED, per_stratum=PER_STRATUM, strata=counts,
                  planned=len(selected), records=[], target_table_sha256=text_digest(targets),
                  catalogue_summary_sha256=text_digest(summary), script_sha256=text_digest(__file__),
                  archive_method_sha256=text_digest(ROOT / "scripts/inspect_jades_cosmic_evolution.py"),
                  dr4_tier_naming_documentation="https://archive.stsci.edu/hlsp/jades")
    for entry in selected:
        try:
            receipt = spectrum_receipt(cache, entry)
        except (ValueError, KeyError, OSError, requests.RequestException) as exc:
            record = dict(**entry, status="NO_ACCEPTED_PRODUCT", reason=str(exc))
        else:
            record = dict(**entry, status="ACCEPTED", **receipt)
        ledger["records"].append(record)
        write_json(LEDGER, ledger)
        print(len(ledger["records"]), entry["dr4_id"], record["status"], flush=True)
    if any(r["status"] != "ACCEPTED" for r in ledger["records"] if r["role"] == "sideband_control"):
        raise ValueError("An original example spectrum could not be recovered")


if __name__ == "__main__":
    main()

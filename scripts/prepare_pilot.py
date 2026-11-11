import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cosmic_anomaly_lab.catalogues import ingest, select_pilot
from cosmic_anomaly_lab.provenance import DownloadCache, write_json

root = Path(__file__).resolve().parents[1]
frame, sources = ingest(root)
config = yaml.safe_load((root / "configs/pilot.yaml").read_text())
pilot = select_pilot(frame, config=config)
pilot.to_csv(root / "data/catalogues/pilot_selection.csv", index=False)
print(pilot[["file", "survey", "field", "z_best", "sample"]].to_string(index=False), flush=True)
cache = DownloadCache(root, budget_bytes=config["raw_storage_budget_MiB"] * 1024**2)
downloads = []
for row in pilot.itertuples():
    try:
        path = cache.fetch(
            row.url,
            "spectra/" + row.file,
            survey=row.survey,
            version="DJA-v4.4",
            kind="observed_spectrum",
        )
        outcome = dict(file=row.file, status="downloaded", bytes=path.stat().st_size)
    except Exception as exc:
        outcome = dict(file=row.file, status="failed", error=str(exc))
    downloads.append(outcome)
    write_json(root / "results/spectrum_downloads.json", downloads)
    print(outcome, flush=True)

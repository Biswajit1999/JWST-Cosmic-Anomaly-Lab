import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cosmic_anomaly_lab.pipeline import run

run(Path(__file__).resolve().parents[1])

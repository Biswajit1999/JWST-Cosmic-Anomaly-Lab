"""Render the observed pilot atlas without changing the original analysis run."""

import argparse

from cosmic_anomaly_lab.visual_evidence import build_atlas

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", default=".")
args = parser.parse_args()
run = build_atlas(args.root)
print({k: run[k] for k in ["observed_sources", "rgb_sources", "measured_apertures"]})

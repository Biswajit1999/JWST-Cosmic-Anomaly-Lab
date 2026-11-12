import argparse

from cosmic_anomaly_lab.expanded_search import run_search

parser = argparse.ArgumentParser()
parser.add_argument("--root", default=".")
parser.add_argument("--limit", type=int, help="Smoke test only; does not overwrite the full run")
args = parser.parse_args()
print(run_search(args.root, args.limit))

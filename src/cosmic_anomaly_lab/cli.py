"""Small reproducible command surface; no implicit bulk archive download."""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="JWST Cosmic Anomaly Lab")
    parser.add_argument("--root", type=Path, default=Path.cwd())
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("analyse", help="Analyse already downloaded pilot products")
    commands.add_parser("ingest", help="Normalize local pinned catalogues and rebuild registry")
    commands.add_parser("verify", help="Verify every downloaded file against its manifest")
    search = commands.add_parser("sources", help="Search local SQLite registry")
    search.add_argument("--survey")
    search.add_argument("--limit", type=int, default=10)
    args = parser.parse_args()
    if args.command == "analyse":
        from .pipeline import run

        run(args.root)
    elif args.command == "ingest":
        from .catalogues import ingest

        ingest(args.root)
    elif args.command == "verify":
        from .provenance import digest

        records = json.loads((args.root / "data/manifests/downloads.json").read_text())
        latest = {r["path"]: r for r in records}
        failed = [
            r["path"]
            for r in latest.values()
            if not (args.root / r["path"]).exists() or digest(args.root / r["path"]) != r["sha256"]
        ]
        print(json.dumps({"products": len(latest), "failed": failed}))
        if failed:
            raise SystemExit(1)
    else:
        from .registry import Registry

        if not (args.root / "data/processed/registry.sqlite").exists():
            parser.error("Run ingest before searching the source registry")
        registry = Registry(args.root / "data/processed/registry.sqlite")
        result = registry.search(args.survey)
        print(json.dumps({"total": len(result), "rows": result[: max(0, args.limit)]}, indent=2))
        registry.close()


if __name__ == "__main__":
    main()

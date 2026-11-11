import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cosmic_anomaly_lab.literature import deduplicate, search
from cosmic_anomaly_lab.provenance import utcnow, write_json

records, logs = [], []
for provider in ["crossref", "openalex"]:
    for query in [
        "JWST JADES CEERS UNCOVER spectroscopy",
        "JWST little red dots alternative physical interpretations",
    ]:
        try:
            found, log = search(provider, query, 10)
            records.extend(found)
            logs.append(log)
        except Exception as exc:
            logs.append(
                dict(
                    provider=provider,
                    query=query,
                    status="failed",
                    error=str(exc),
                    accessed_utc=utcnow(),
                )
            )
        print(logs[-1], flush=True)
        write_json("literature/search_log.json", logs)
        write_json("literature/retrieved_metadata.json", deduplicate(records))

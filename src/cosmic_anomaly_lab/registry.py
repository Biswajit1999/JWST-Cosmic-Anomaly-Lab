"""SQLite source evidence and separate observation records; no positional auto-merging."""

import json
import sqlite3
from pathlib import Path
from uuid import NAMESPACE_URL, uuid5

from .provenance import utcnow

STATUSES = {
    "KNOWN_MODEL_CONSISTENT",
    "MULTIPLE_MODELS_VIABLE",
    "POORLY_EXPLAINED",
    "INSUFFICIENT_EVIDENCE",
    "LIKELY_SYSTEMATIC",
    "REQUIRES_FOLLOWUP",
}


def source_id(namespace, native_id):
    return "CAL-" + uuid5(NAMESPACE_URL, f"{namespace}:{native_id}").hex[:16]


class Registry:
    def __init__(self, path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(path)
        self.db.execute("PRAGMA foreign_keys=ON")
        self.db.executescript("""
          CREATE TABLE IF NOT EXISTS sources (
            source_id TEXT PRIMARY KEY, namespace TEXT NOT NULL, native_id TEXT NOT NULL,
            ra_deg REAL NOT NULL CHECK(ra_deg>=0 AND ra_deg<360),
            dec_deg REAL NOT NULL CHECK(dec_deg>=-90 AND dec_deg<=90),
            metadata_json TEXT NOT NULL, UNIQUE(namespace,native_id));
          CREATE TABLE IF NOT EXISTS observations (
            observation_id TEXT PRIMARY KEY, source_id TEXT NOT NULL REFERENCES sources,
            survey TEXT NOT NULL, field TEXT NOT NULL, product_url TEXT,
            metadata_json TEXT NOT NULL);
          CREATE INDEX IF NOT EXISTS obs_source ON observations(source_id);
          CREATE TABLE IF NOT EXISTS evidence (
            source_id TEXT NOT NULL REFERENCES sources, run_id TEXT NOT NULL,
            status TEXT NOT NULL, evidence_json TEXT NOT NULL, created_utc TEXT NOT NULL,
            PRIMARY KEY(source_id,run_id));
          CREATE TABLE IF NOT EXISTS source_links (
            source_a TEXT REFERENCES sources, source_b TEXT REFERENCES sources,
            relation TEXT CHECK(relation IN ('possible_counterpart','lensed_image_of','confirmed_counterpart')),
            evidence_json TEXT NOT NULL, PRIMARY KEY(source_a,source_b,relation));
        """)

    def add_source(self, namespace, native_id, ra_deg, dec_deg, **metadata):
        sid = source_id(namespace, str(native_id))
        self.db.execute(
            """INSERT INTO sources VALUES (?,?,?,?,?,?)
            ON CONFLICT(source_id) DO UPDATE SET ra_deg=excluded.ra_deg,
            dec_deg=excluded.dec_deg, metadata_json=excluded.metadata_json""",
            (
                sid,
                namespace,
                str(native_id),
                float(ra_deg),
                float(dec_deg),
                json.dumps(metadata, allow_nan=False),
            ),
        )
        return sid

    def add_observation(self, observation_id, sid, survey, field, url, **metadata):
        self.db.execute(
            """INSERT INTO observations VALUES (?,?,?,?,?,?)
            ON CONFLICT(observation_id) DO UPDATE SET metadata_json=excluded.metadata_json,
            product_url=excluded.product_url""",
            (observation_id, sid, survey, field, url, json.dumps(metadata, allow_nan=False)),
        )

    def add_evidence(self, sid, run_id, status, evidence):
        if status not in STATUSES:
            raise ValueError("Invalid evidence status")
        self.db.execute(
            "INSERT OR REPLACE INTO evidence VALUES (?,?,?,?,?)",
            (sid, run_id, status, json.dumps(evidence, allow_nan=False), utcnow()),
        )

    def search(self, survey=None):
        query = "SELECT DISTINCT s.* FROM sources s JOIN observations o USING(source_id)"
        if survey is not None:
            return self.db.execute(query + " WHERE o.survey=?", (survey,)).fetchall()
        return self.db.execute(query).fetchall()

    def close(self):
        self.db.commit()
        self.db.close()

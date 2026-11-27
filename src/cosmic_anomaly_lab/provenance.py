"""Budgeted, atomic downloads and append-only version manifests."""

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    # Windows indexers can briefly hold the destination after it is read.
    # Retain atomic replacement and fail after a bounded retry, never truncate.
    for attempt in range(5):
        try:
            temporary.replace(path)
            break
        except PermissionError:
            if attempt == 4:
                raise
            time.sleep(.1 * 2**attempt)


def text_digest(path):
    """Canonical UTF-8 text hash, insensitive to Windows/Git newline translation."""
    return hashlib.sha256(Path(path).read_text(encoding="utf-8").encode("utf-8")).hexdigest()


class DownloadCache:
    """Single-writer cache. Resume by completed product, restart interrupted products.

    Published checksums establish upstream integrity; SHA256 without a publisher
    checksum establishes only a reproducibility fingerprint. Refresh appends versions.
    """

    def __init__(self, root, budget_bytes=512 * 1024**2):
        self.root = Path(root).resolve()
        self.raw = self.root / "data/raw"
        self.raw.mkdir(parents=True, exist_ok=True)
        self.manifest = self.root / "data/manifests/downloads.json"
        self.records = json.loads(self.manifest.read_text()) if self.manifest.exists() else []
        self.budget = budget_bytes
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "JWST-Cosmic-Anomaly-Lab/0.1 (research pilot)"
        self.session.mount(
            "https://",
            HTTPAdapter(
                max_retries=Retry(
                    total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504]
                )
            ),
        )

    def fetch(
        self,
        url,
        filename,
        *,
        survey,
        version,
        kind,
        checksum=None,
        refresh=False,
        json_payload=None,
        form_payload=None,
        response=None,
    ):
        """Store one response, using an optional already-open GET for the same URL."""
        if json_payload is not None and form_payload is not None:
            raise ValueError("Specify only one POST request body")
        path = (self.raw / filename).resolve()
        if not path.is_relative_to(self.raw) or path == self.raw:
            raise ValueError("Download path must remain inside data/raw")
        records = [
            r
            for r in self.records
            if r["url"] == url
            and r["path"] == path.relative_to(self.root).as_posix()
            and r.get("request_json") == json_payload
            and r.get("request_form") == form_payload
        ]
        if path.exists() and not refresh and records:
            if digest(path) != records[-1]["sha256"]:
                raise ValueError(f"Cache checksum mismatch: {filename}")
            if checksum and digest(path, checksum.split(":")[0]) != checksum.split(":")[1]:
                raise ValueError("Published checksum mismatch")
            return path
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_suffix(path.suffix + ".part")
        part.unlink(missing_ok=True)
        used = sum(p.stat().st_size for p in self.raw.rglob("*") if p.is_file())
        available = self.budget - used
        if response is None:
            if json_payload is not None:
                response = self.session.post(url, json=json_payload, stream=True, timeout=(20, 90))
            elif form_payload is not None:
                response = self.session.post(url, data=form_payload, stream=True, timeout=(20, 90))
            else:
                response = self.session.get(url, stream=True, timeout=(20, 90))
        else:
            original = response.history[0].request if response.history else response.request
            expected_url = requests.Request("GET", url).prepare().url
            if (
                json_payload is not None or form_payload is not None
                or original.method != "GET" or original.url != expected_url
            ):
                raise ValueError("Prepared response does not match the cache GET")
        with response:
            response.raise_for_status()
            length = int(response.headers.get("Content-Length", "0"))
            if length > available:
                raise ValueError(f"Storage budget exceeded before download: {filename}")
            total = 0
            try:
                with part.open("wb") as stream:
                    for block in response.iter_content(1024 * 1024):
                        total += len(block)
                        if total > available:
                            raise ValueError("Storage budget exceeded during download")
                        stream.write(block)
                if checksum:
                    algorithm, expected = checksum.split(":", 1)
                    if digest(part, algorithm) != expected:
                        raise ValueError("Published checksum mismatch")
                record = dict(
                    url=url,
                    final_url=response.url,
                    path=path.relative_to(self.root).as_posix(),
                    bytes=total,
                    sha256=digest(part),
                    published_checksum=checksum,
                    accessed_utc=utcnow(),
                    etag=response.headers.get("ETag"),
                    last_modified=response.headers.get("Last-Modified"),
                    survey=survey,
                    version=version,
                    kind=kind,
                    previous_sha256=records[-1]["sha256"] if records else None,
                )
                if json_payload is not None:
                    record["request_json"] = json_payload
                if form_payload is not None:
                    record["request_form"] = form_payload
                part.replace(path)
                self.records.append(record)
                write_json(self.manifest, self.records)
            except Exception:
                part.unlink(missing_ok=True)
                raise
        return path

"""Bounded concurrent network reads with one serial cache/manifest writer."""

import hashlib
from concurrent.futures import ThreadPoolExecutor

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .provenance import DownloadCache, digest, utcnow, write_json


def read_remote(url, maximum_bytes):
    """Bound bytes in memory per worker, retry transient upstream responses."""
    with requests.Session() as session:
        session.mount(
            "https://",
            HTTPAdapter(
                max_retries=Retry(
                    total=2, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504]
                )
            ),
        )
        session.headers["User-Agent"] = "JWST-Cosmic-Anomaly-Lab/0.2 (public research)"
        with session.get(url, stream=True, timeout=(20, 60)) as response:
            response.raise_for_status()
            if int(response.headers.get("Content-Length", "0")) > maximum_bytes:
                raise ValueError("Product exceeds configured per-product byte cap")
            data = bytearray()
            for block in response.iter_content(65536):
                if len(data) + len(block) > maximum_bytes:
                    raise ValueError("Product exceeds configured per-product byte cap")
                data.extend(block)
            return data, dict(
                final_url=response.url,
                etag=response.headers.get("ETag"),
                last_modified=response.headers.get("Last-Modified"),
            )


def download_batch(
    root,
    jobs,
    *,
    budget_bytes,
    maximum_bytes=8 * 1024**2,
    workers=6,
    checkpoint="results/expanded_downloads.json",
    reader=read_remote,
):
    """Only the calling thread writes files/manifests; six bounded futures at most.

    Do not run another cache writer concurrently. A failed job is recorded and
    does not silently shrink the scientific denominator. Reruns verify and reuse
    completed products. Existing untracked files are rejected, not overwritten.
    """
    if not 1 <= workers <= 16 or maximum_bytes <= 0:
        raise ValueError("Invalid worker count or product cap")
    cache = DownloadCache(root, budget_bytes=budget_bytes)
    latest = {r["path"]: r for r in cache.records}
    used = sum(p.stat().st_size for p in cache.raw.rglob("*") if p.is_file())
    outcomes, pending = [], []
    for job in jobs:
        path = (cache.raw / job["filename"]).resolve()
        if not path.is_relative_to(cache.raw) or path == cache.raw:
            raise ValueError("Download path must remain inside data/raw")
        key = path.relative_to(cache.root).as_posix()
        if path.exists():
            record = latest.get(key)
            if not record or record["url"] != job["url"] or digest(path) != record["sha256"]:
                raise ValueError(f"Untracked or changed cache file: {key}")
            outcomes.append(dict(file=job["filename"], status="CACHED", sha256=record["sha256"]))
        else:
            pending.append((job, path, key))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for start in range(0, len(pending), workers):
            batch = pending[start : start + workers]
            futures = [pool.submit(reader, j["url"], maximum_bytes) for j, _, _ in batch]
            for (job, path, key), future in zip(batch, futures, strict=True):
                try:
                    data, response_meta = future.result()
                    if len(data) > maximum_bytes or used + len(data) > budget_bytes:
                        raise ValueError("Raw storage or per-product budget exceeded")
                    sha = hashlib.sha256(data).hexdigest()
                    path.parent.mkdir(parents=True, exist_ok=True)
                    part = path.with_suffix(path.suffix + ".part")
                    part.write_bytes(data)
                    part.replace(path)
                    cache.records.append(
                        dict(
                            url=job["url"],
                            path=key,
                            bytes=len(data),
                            sha256=sha,
                            published_checksum=None,
                            accessed_utc=utcnow(),
                            survey=job["survey"],
                            version=job["version"],
                            kind=job["kind"],
                            previous_sha256=None,
                            **response_meta,
                        )
                    )
                    used += len(data)
                    outcomes.append(dict(file=job["filename"], status="DOWNLOADED", sha256=sha))
                except Exception as exc:
                    outcomes.append(dict(file=job["filename"], status="FAILED", reason=str(exc)))
            write_json(cache.manifest, cache.records)
            write_json(cache.root / checkpoint, outcomes)
            print(
                f"Products {len(outcomes)}/{len(jobs)}; failures "
                f"{sum(r['status'] == 'FAILED' for r in outcomes)}; raw MiB {used / 1024**2:.1f}",
                flush=True,
            )
    write_json(cache.root / checkpoint, outcomes)
    return outcomes

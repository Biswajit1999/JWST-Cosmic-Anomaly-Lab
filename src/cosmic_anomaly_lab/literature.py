"""Auditable bibliographic searches; metadata retrieval is not a systematic review."""

import re

import requests

from .provenance import utcnow


def canonical_doi(value):
    if not value:
        return None
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", value.strip(), flags=re.I).lower()


def deduplicate(records):
    """Only merge explicit shared identifiers, not approximate title matches."""
    out = []
    for record in records:
        record = record.copy()
        record["doi"] = canonical_doi(record.get("doi"))
        record["arxiv"] = re.sub(r"v\d+$", "", record.get("arxiv") or "") or None
        match = next(
            (
                x
                for x in out
                if (record["doi"] and x.get("doi") == record["doi"])
                or (record["arxiv"] and x.get("arxiv") == record["arxiv"])
            ),
            None,
        )
        if match:
            match.setdefault("metadata_sources", []).extend(record.get("metadata_sources", []))
            for key, value in record.items():
                if match.get(key) is None and value is not None:
                    match[key] = value
        else:
            out.append(record)
    return out


def search(provider, query, limit=10):
    if provider == "crossref":
        url = "https://api.crossref.org/works"
        response = requests.get(
            url, params={"query.bibliographic": query, "rows": limit}, timeout=30
        )
        response.raise_for_status()
        records = [
            dict(
                title=" ".join(x.get("title", [])),
                doi=x.get("DOI"),
                publication_status=x.get("type"),
                year=(x.get("published", {}).get("date-parts", [[None]])[0][0]),
                metadata_sources=["crossref"],
            )
            for x in response.json()["message"]["items"]
        ]
    elif provider == "openalex":
        url = "https://api.openalex.org/works"
        response = requests.get(url, params={"search": query, "per-page": limit}, timeout=30)
        response.raise_for_status()
        records = [
            dict(
                title=x.get("title"),
                doi=x.get("doi"),
                year=x.get("publication_year"),
                publication_status=x.get("type"),
                metadata_sources=["openalex"],
            )
            for x in response.json()["results"]
        ]
    else:
        raise ValueError("Supported metadata providers: crossref and openalex")
    return records, dict(
        provider=provider,
        query=query,
        limit=limit,
        retrieved=len(records),
        accessed_utc=utcnow(),
        request_url=response.url,
        scope="Title/metadata retrieval only; returned papers not all relevance-reviewed",
    )

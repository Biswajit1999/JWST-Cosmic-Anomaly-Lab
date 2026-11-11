import sys
from concurrent.futures import ThreadPoolExecutor

import requests

sys.path.insert(0, "src")
from cosmic_anomaly_lab.provenance import utcnow, write_json

surveys = {
    "MAST": "https://archive.stsci.edu/missions-and-data/jwst",
    "JADES": "https://jades-survey.github.io/scientists/data.html",
    "CEERS": "https://ceers.github.io/releases.html",
    "UNCOVER": "https://jwst-uncover.github.io/DR4.html",
    "GLASS-JWST": "https://archive.stsci.edu/hlsp/glass-jwst",
    "COSMOS-Web": "https://cosmos.astro.caltech.edu/page/cosmosweb-dr",
    "DJA": "https://dawn-cph.github.io/dja/",
    "NEXUS": "https://nexus-jwst.github.io/",
    "RUBIES": "https://rubies-survey.github.io/",
    "PRIMER": "https://primer-jwst.github.io/",
    "FRESCO": "https://archive.stsci.edu/hlsp/fresco",
    "NGDEEP": "https://archive.stsci.edu/hlsp/ngdeep",
}


def probe(item):
    name, url = item
    try:
        r = requests.get(url, timeout=30)
        return dict(
            survey=name,
            url=url,
            http_status=r.status_code,
            bytes=len(r.content),
            accessed_utc=utcnow(),
            scope="landing_page_only; product availability is not inferred",
            final_url=r.url,
        )
    except Exception as e:
        return dict(survey=name, url=url, error=str(e), accessed_utc=utcnow())


results = list(ThreadPoolExecutor(max_workers=4).map(probe, surveys.items()))
write_json("data/manifests/survey_inventory.json", results)
for x in results:
    print(x)

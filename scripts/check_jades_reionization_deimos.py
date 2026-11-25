"""Check the GOODS-N Lyα leads against the published Keck/DEIMOS coordinates."""

import json
import re

from astropy.coordinates import SkyCoord
from check_jades_reionization_mmt import CELL, ROW, text
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

URL = "https://arxiv.org/html/2402.06070"
NUMBER = re.compile(r"^[+−-]?\s*\d+(?:\.\d+)?")


def number(cell):
    match = NUMBER.match(text(cell))
    if not match:
        raise ValueError("Missing numeric cell in the DEIMOS table")
    return float(match[0].replace("−", "-").replace(" ", ""))


def main():
    mmt = OUT / "mmt_literature_check.json"
    ids = {r["dr4_id"] for r in json.loads(mmt.read_text())["results"]}
    targets = [r for r in json.loads(LEDGER.read_text())["records"] if r["dr4_id"] in ids]
    cache = DownloadCache(ROOT, 28_000 * 1024**2)
    published = cache.fetch(
        URL,
        "metadata/tang_2402_06070.html",
        survey="Keck/DEIMOS",
        version="author arXiv HTML",
        kind="published_literature_table",
    )
    rows = []
    for html_row in ROW.findall(published.read_text(encoding="utf-8")):
        cells = CELL.findall(html_row)
        if len(cells) == 10 and (name := text(cells[0])).startswith("DEIMOS-"):
            rows.append(
                dict(
                    name=name,
                    ra=number(cells[1]),
                    dec=number(cells[2]),
                    line_redshift=number(cells[4]),
                )
            )
    if len(rows) != 19 or len({r["name"] for r in rows}) != len(rows):
        raise ValueError("Missing or duplicate Keck/DEIMOS emitter rows")
    sky = SkyCoord([r["ra"] for r in rows], [r["dec"] for r in rows], unit="deg")
    results = []
    for target in targets:
        deltas = SkyCoord(target["ra"], target["dec"], unit="deg").separation(sky).arcsec
        idx = int(deltas.argmin())
        results.append(
            dict(
                dr4_id=target["dr4_id"],
                matches_within_arcsec_1=int((deltas <= 1).sum()),
                nearest=dict(rows[idx], separation_arcsec=float(deltas[idx])),
            )
        )
    write_json(
        OUT / "deimos_literature_check.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            target_receipt_sha256=text_digest(mmt),
            published_table=receipt(published, URL),
            emitter_rows=len(rows),
            results=results,
            scope="Keck/DEIMOS emitter rows in Tang et al. 2024 Table 6 and its continuation. "
            "A coordinate miss does not establish coverage or literature priority.",
        ),
    )
    print(len(rows), results)


if __name__ == "__main__":
    main()

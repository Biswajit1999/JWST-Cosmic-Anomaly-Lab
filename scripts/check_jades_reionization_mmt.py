"""Check the leading GOODS-N Lyα targets against the observed MMT emitter table."""

import json
import re
from html import unescape

from astropy.coordinates import SkyCoord
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

URL = "https://arxiv.org/html/2509.18302"
ROW = re.compile(r"<tr\b[^>]*>(.*?)</tr>", re.IGNORECASE | re.DOTALL)
CELL = re.compile(r"<td\b[^>]*>(.*?)</td>", re.IGNORECASE | re.DOTALL)
SOURCE = re.compile(r"(?:z[5-8]-GN[DW]|Hu10|Stark11|Jung18)-\d+")


def text(value):
    return " ".join(unescape(re.sub(r"<[^>]+>"," ",value)).split())


def main():
    path = OUT/"crosschecks.json"
    old = json.loads(path.read_text())
    ids = {r["dr4_id"] for r in old["catalogue_unmatched"] if r["dr4_id"].startswith("goods-n")}
    targets = [r for r in json.loads(LEDGER.read_text())["records"] if r["dr4_id"] in ids]
    cache = DownloadCache(ROOT,28_000*1024**2)
    published = cache.fetch(URL,"metadata/prieto_lyon_2509_18302.html",survey="MMT/Binospec",
                            version="author arXiv HTML",kind="published_literature_table")
    entries = []
    for html_row in ROW.findall(published.read_text(encoding="utf8")):
        cells = CELL.findall(html_row)
        if len(cells) == 11 and SOURCE.fullmatch(name := text(cells[0])):
            entries.append(dict(name=name,ra=float(text(cells[1])),dec=float(text(cells[2])),
                                line_redshift=text(cells[3]),systemic_redshift=text(cells[4])))
    if len(entries) != 62 or len({r["name"] for r in entries}) != len(entries):
        raise ValueError(f"Parsed {len(entries)} coordinate rows, expected 62 unique emitters")
    sky = SkyCoord([r["ra"] for r in entries],[r["dec"] for r in entries],unit="deg")
    result = []
    for row in targets:
        distances = SkyCoord(row["ra"],row["dec"],unit="deg").separation(sky)
        closest = int(distances.argmin())
        result.append(dict(dr4_id=row["dr4_id"],matches_within_arcsec_1=int((distances.arcsec<=1.).sum()),
                           nearest=dict(**entries[closest],separation_arcsec=float(distances[closest].arcsec))))
    write_json(OUT/"mmt_literature_check.json",dict(executed_utc=utcnow(),script_sha256=text_digest(__file__),
        crosscheck_sha256=text_digest(path),published_table=receipt(published,URL),emitter_rows=len(entries),
        results=result,scope="Published MMT emitter table only. No inference about target allocation or other literature."))
    print(result,flush=True)


if __name__=="__main__":
    main()

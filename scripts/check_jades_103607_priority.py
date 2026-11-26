"""Cross-match JADES 103607 against the updated published JADES transient roster."""

import json
import re

from astropy.coordinates import SkyCoord
from check_jades_reionization_mmt import CELL, ROW, text
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from vet_jwst_module_b import tns_catalogue, tns_match

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

URL = "https://arxiv.org/html/2607.12018"
TNS_NAME = re.compile(r"(?:AT|SN)\s+(20\d{2}[a-z]{2,6})\b")
RADIUS = 5.0
WEB_QUERIES = (
    '"JADES 103607" "October 2026" transient GOODS-S JWST paper',
    'JADES "103607" galaxy redshift 53.18867 transient JWST NIRSpec',
    '"53.1886698" "27.8125394" transient JWST JADES public TNS',
    '"12577" "transient" JWST "October" "2026" JADES transient discovery DeCoursey Egami publication',
    '"October 2026" "GOODS-S" transient "JWST" JADES 2026 observing programme 12577 transient reports',
    '"03:32:45.28" "27:48:45" JADES transient GOODS-S ATel 2026',
    '"53.18867" "-27.81255" JWST GOODS-S transient supernova',
)


def main():
    target = OUT / "module_b_001_native_epoch_flux.json"
    position = json.loads(target.read_text())["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, 30_000 * 1024**2)
    listing = cache.fetch(
        URL,
        "metadata/jts2_2607_12018.html",
        survey="JADES Transient Survey",
        version="arXiv author HTML",
        kind="published_literature_table",
    )
    names = []
    for row in ROW.findall(listing.read_text(encoding="utf-8")):
        fields = CELL.findall(row)
        if fields and (match := TNS_NAME.match(text(fields[0]))):
            names.append(match[1])
    if len(names) != 83 or len(set(names)) != 83:
        raise ValueError(f"Expected 83 distinct TNS names in the JTS II tables, found {len(names)}")
    registry = tns_catalogue(cache, snapshot_name="tns_public_objects_verified_20261126.csv.zip")
    positions = registry[1][registry[1]["name"].isin(names)].copy()
    if len(positions) != len(names) or set(positions["name"]) != set(names):
        raise ValueError("The 83 published JTS II names did not resolve uniquely in the registry")
    positions["separation_arcsec"] = sky.separation(
        SkyCoord(positions["ra"], positions["declination"], unit="deg")
    ).arcsec
    nearest = positions.nsmallest(5, "separation_arcsec")
    prior = OUT / "goods_s_prior_records.json"
    past = json.loads(prior.read_text())
    matched_plans = []
    for program in past["public_apt"]:
        for cat in program["catalogues"]:
            for close in cat["positions"].values():
                for row in close["matches"]:
                    if row["ID"] == 103607 or str(row["ID"]) == "103607":
                        matched_plans.append(
                            dict(program=program["program"], catalogue=cat["name"], entry=row)
                        )
    write_json(
        OUT / "module_b_001_priority_20261126.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            target_input_sha256=text_digest(target),
            prior_registry_sha256=text_digest(prior),
            tns=tns_match(registry, sky),
            jts_ii=dict(
                table=receipt(listing, URL),
                source_count=len(positions),
                tns_names=names,
                radius_arcsec=RADIUS,
                matches=int((positions["separation_arcsec"] <= RADIUS).sum()),
                nearest=json.loads(nearest.to_json(orient="records")),
            ),
            prior_apt_catalogue_matches=matched_plans,
            public_web_search=dict(
                checked_utc_date="2026-11-26",
                queries=list(WEB_QUERIES),
                outcome="No earlier public report matching the position and October 2026 event was located in the results screened.",
            ),
            scope="JTS II reclassifies a pre-October 2026 sample. Its comparison checks for prior reports "
            "at the same sky position. The daily TNS catalogue has its own timestamp. "
            "Literature and catalogue coverage cannot exhaust unindexed or unpublished work.",
        ),
    )
    print(registry[0]["snapshot_label_utc"], len(registry[1]), len(positions))
    print(nearest[["name", "separation_arcsec"]].to_string(index=False))


if __name__ == "__main__":
    main()

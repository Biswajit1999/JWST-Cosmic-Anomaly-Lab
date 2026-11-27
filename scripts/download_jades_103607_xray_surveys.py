"""Pin eROSITA and ROSAT historical X-ray products for JADES 103607."""

import json
from urllib.parse import urlencode
from xml.etree import ElementTree

from astropy.coordinates import SkyCoord
from astropy.table import Table
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_xray_surveys_products"
POSITION = OUT / "module_b_001_native_november_fade.json"
ERODAT = "https://erosita.mpe.mpg.de/erodat"
HEASARC_TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
ROSAT = "https://heasarc.gsfc.nasa.gov/FTP/rosat/data/pspc/processed_data"
BANDS = {"024": [0.2, 2.3], "023": [2.3, 5.0]}
STORAGE_MIB = 34_200


def pin(cache, url, name, survey, version, kind):
    path = cache.fetch(url, f"xray_surveys/{name}", survey=survey, version=version, kind=kind)
    print(name, flush=True)
    return path, receipt(path, url)


def erosita(cache, ra, dec):
    tile_url = f"{ERODAT}/skyview/skytile_search_api/?" + urlencode(dict(RA=ra, DEC=dec, RAD=0))
    path, tile_product = pin(cache, tile_url, "erosita_tile.json", "eROSITA", "DR1", "tile_search")
    tile = json.loads(path.read_text(encoding="utf-8"))
    if not tile["success"]:
        raise ValueError("eROSITA tile lookup failed")
    field, = [str(row["srvmap"]).zfill(6) for row in tile["tiles"]
              if row["canonical"] and row["de_sky"]]
    base = f"{ERODAT}/data/download/{field[3:]}/{field[:3]}"

    def product(band, directory, stem):
        name = f"em01_{field}_{band}_{stem}_c010.fits.gz"
        return pin(cache, f"{base}/{directory}_010/{name}", name,
                   "eROSITA", "DR1 c010", stem)[1]

    limits, catalogues = {}, {}
    for band in BANDS:
        for release in ["DR1_eRASS1", "DR2_eRASSc3"]:
            name = f"{release}_{band}"
            args = urlencode(dict(ra=ra, dec=dec, band=band, dr_survey=release))
            limits[name] = pin(cache, f"{ERODAT}/upperlimit/service?{args}", name + ".json",
                               "eROSITA", release, "forced_aperture_limit")[1]
    for cat in ["DR1_Main", "DR1_Hard", "DR2_Main", "DR2_Hard"]:
        args = urlencode(dict(CAT=cat, RA=ra, DEC=dec, SR=0.05, VERB=3))
        catalogues[cat] = pin(cache, f"{ERODAT}/catalogue/SCS?{args}", cat + ".vot",
                              "eROSITA", cat, "cone_3_arcmin")[1]
    return dict(tile=tile_product, tile_id=field, limits=limits, catalogues=catalogues,
                events=product("020", "EXP", "EventList"), bands={
                    band: dict(counts=product(band, "EXP", "Image"),
                               exposure=product(band, "DET", "ExposureMap"),
                               background=product(band, "DET", "BackgrImage"))
                    for band in BANDS
                })


def rosat(cache, ra, dec):
    queries = {}
    sequences = []
    for catalog, radius, version in [
        ("rosmaster", 1.2, "pointed RDF"),
        ("rassmaster", 7.5, "RASS RDF"),
    ]:
        query = (f"SELECT TOP 100 * FROM {catalog} WHERE CONTAINS(POINT('ICRS',ra,dec),"
                 f"CIRCLE('ICRS',{ra},{dec},{radius}))=1")
        args = urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable", QUERY=query))
        path, queries[catalog] = pin(cache, f"{HEASARC_TAP}?{args}", catalog + ".vot",
                                     "ROSAT", version, "pointing_inventory")
        statuses = [
            node.attrib["value"] for node in ElementTree.parse(path).iter()
            if node.attrib.get("name", "").upper() == "QUERY_STATUS"
        ]
        if statuses != ["OK"]:
            raise ValueError(f"ROSAT {catalog} TAP status: {statuses}")
        table = Table.read(path)
        point = SkyCoord(ra, dec, unit="deg")
        separations = point.separation(SkyCoord(table["ra"], table["dec"], unit="deg")).deg
        chosen = list(range(len(table))) if catalog == "rosmaster" else [separations.argmin()]
        for index in chosen:
            row = table[index]
            sequence = row["DataLinkID"].lower()
            if catalog == "rosmaster" and not sequence.startswith("rp"):
                raise ValueError(f"Unhandled ROSAT observation {sequence}")
            group = "900000" if catalog == "rassmaster" else "100000"
            base = f"{ROSAT}/{group}/{sequence}"
            products = {}
            for suffix, kind in [
                ("bas.fits.Z", "standard_events"),
                ("mex.fits" if group == "900000" else "mex.fits.Z", "merged_exposure"),
            ]:
                filename = f"{sequence}_{suffix}"
                products[kind] = pin(cache, f"{base}/{filename}", filename,
                                     "ROSAT", version, kind)[1]
            sequences.append(dict(
                sequence=sequence.upper(), inventory=catalog,
                center_deg=[float(row["ra"]), float(row["dec"])],
                center_distance_deg=float(separations[index]), **products,
            ))
    return dict(queries=queries, sequences=sequences)


if __name__ == "__main__":
    position = json.loads(POSITION.read_text(encoding="utf-8"))["icrs_deg"]
    cache = DownloadCache(ROOT, STORAGE_MIB * 1024**2)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), icrs_deg=position,
        erosita=erosita(cache, *position), rosat=rosat(cache, *position),
        position_sha256=text_digest(POSITION), script_sha256=text_digest(__file__),
        raw_storage_cap_mib=STORAGE_MIB,
    ))

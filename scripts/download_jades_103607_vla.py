"""Fetch public VLA images and check archive coverage at JADES 103607."""

import json
from urllib.parse import urlencode, urljoin, urlsplit
from xml.etree import ElementTree

import pandas as pd
from astropy.table import Table
from astropy.time import Time
from bs4 import BeautifulSoup
from check_jades_103607_independent_2026 import END, SOURCE, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from requests import HTTPError

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

TAG = "module_b_001_vla_inventory"
BASE = "https://zenodo.org/records/18202043/files"
DOI = "10.5281/zenodo.18202043"
CADC = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca"
SKYVIEW = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
NRAO = "https://data-query.nrao.edu/tap/sync"
EVENT_CONE_DEG = 0.40
BUDGET = 34_200 * 1024**2
IMAGE = "single_pass_rp5_p3mJy.image.pbcor.tt0.fits"
CATALOGUE = "Adams_CDFS_StokesI.fits"
PRODUCTS = {
    IMAGE: ("ab7e694d712f797cbcd226c808dd4cbf", "observed_stokes_i_mosaic"),
    CATALOGUE: ("ad13e5e59d57ea43cbdb1eaf21a62d04", "published_stokes_i_components"),
}


def download_w_cdfs(cache):
    products = {}
    for filename, (md5, kind) in PRODUCTS.items():
        url = f"{BASE}/{filename}"
        path = cache.fetch(
            url, f"cross_detector/vla/{filename}", survey="NRAO VLA W-CDFS",
            version=DOI, kind=kind, checksum=f"md5:{md5}",
        )
        products[filename] = receipt(path, url)
        print(f"{filename}: {path.stat().st_size:,} bytes", flush=True)
    return dict(products=products, doi=DOI, observing_interval_utc=["2021-07-24", "2021-09-09"])


def query_status(path):
    xml = ElementTree.parse(path).getroot()
    statuses = [x.get("value") for x in xml.findall(".//{*}INFO") if x.get("name") == "QUERY_STATUS"]
    if statuses != ["OK"]:
        raise ValueError(f"Archive VOTable query did not succeed: {statuses}")


def optional(value):
    return None if str(value) == "--" else str(value)


def fetch_vlass(cache, position):
    ra, dec = position
    params = dict(POS=f"CIRCLE {ra} {dec} 0.01", COLLECTION="VLASS", MAXREC=100)
    url = f"{CADC}/sia/v2query?{urlencode(params)}"
    list_path = cache.fetch(
        url, f"metadata/{TAG}_cadc.vot", survey="CADC / VLASS",
        version="VLASS collection, SIA v2 coverage", kind="archive_geometry",
    )
    query_status(list_path)
    rows = Table.read(list_path, format="votable")
    if not len(rows) or len(rows) >= 100:
        raise ValueError("No usable public VLASS coverage inventory")
    products = []
    for row in sorted(rows, key=lambda v: v["obs_id"]):
        observation = str(row["obs_id"])
        link_url = str(row["access_url"])
        if urlsplit(link_url).hostname != urlsplit(CADC).hostname:
            raise ValueError("Unexpected VLASS DataLink provider")
        link_path = cache.fetch(
            link_url, f"metadata/vla/{observation}.vot", survey="CADC / VLASS",
            version="VLASS DataLink", kind="archive_product_links",
        )
        links = Table.read(link_path, format="votable")
        assets = {}
        for key, suffix in (("science", ".image.pbcor.tt0.subim.fits"), ("rms", ".image.pbcor.tt0.rms.subim.fits")):
            match = [link for link in links if str(link["access_url"]).endswith(suffix)]
            if not match and key == "rms":
                continue
            if len(match) != 1:
                raise ValueError(f"Unexpected number of VLASS {key} products for {observation}")
            product_url = str(match[0]["access_url"])
            filename = product_url.rsplit("/", 1)[-1]
            cut = urlencode(dict(ID=f"nrao:VLASS/{filename}", CIRCLE=f"{ra} {dec} 0.035"))
            cutout_url = f"{CADC}/caom2ops/sync?{cut}"
            try:
                item = cache.fetch(
                    cutout_url, f"cross_detector/vla/{observation}.{key}.fits",
                    survey="CADC / VLASS", version="Quick Look SODA cutout", kind=f"observed_{key}_image",
                )
            except HTTPError as exc:
                if key == "science":
                    raise
                assets["rms_request_http_status"] = exc.response.status_code
                continue
            assets[key] = receipt(item, cutout_url)
        products.append(dict(
            observation_id=observation, archive_mjd_range=[optional(row["t_min"]), optional(row["t_max"])],
            region=str(row["s_region"]), datalink=receipt(link_path, link_url), **assets,
        ))
        print(f"CADC {observation}: science cutout", flush=True)
    return dict(listing=receipt(list_path, url), products=products)


def fetch_nvss(cache, position):
    params = dict(survey="NVSS", position=",".join(map(str, position)), size=0.55, pixels=132, projection="Tan")
    url = f"{SKYVIEW}?{urlencode(params)}"
    rel = f"metadata/{TAG}_nvss.html"
    for refresh in (False, True):
        page = cache.fetch(
            url, rel, survey="NASA SkyView / NVSS", version="NVSS TAN reprojected image request",
            kind="cutout_service_response", refresh=refresh,
        )
        parsed = BeautifulSoup(page.read_bytes(), "html.parser")
        matches = [link["href"] for link in parsed.find_all("a", href=True) if link.get_text(strip=True) == "FITS"]
        if len(matches) != 1:
            raise ValueError("SkyView NVSS response did not identify one FITS image")
        href = urljoin(SKYVIEW, matches[0])
        if urlsplit(href).hostname != urlsplit(SKYVIEW).hostname:
            raise ValueError("Unexpected NVSS cutout provider")
        try:
            path = cache.fetch(
                href, "cross_detector/vla/nvss_jades103607_skyview.fits", survey="NASA SkyView / NVSS",
                version="NVSS original data, TAN cutout reprojected by SkyView", kind="observed_radio_image",
            )
        except HTTPError as exc:
            if refresh or exc.response.status_code not in (403, 404, 410):
                raise
            continue
        return dict(request=receipt(page, url), science=receipt(path, href), survey_observing_decade="1990s")
    raise RuntimeError("SkyView NVSS cutout did not resolve")


def event_search(cache, position):
    ra, dec = position
    begin, end = Time([START, END]).mjd
    fields = (
        "project_code,obs_publisher_did,target_name,t_min,t_max,s_ra,s_dec,s_fov,"
        "freq_min,freq_max,proprietary_status"
    )
    query = (
        f"SELECT TOP 400 {fields} FROM tap_schema.obscore WHERE "
        f"1=CONTAINS(POINT('ICRS',s_ra,s_dec),CIRCLE('ICRS',{ra},{dec},{EVENT_CONE_DEG:.2f})) "
        f"AND t_min >= {begin} AND t_min < {end} AND instrument_name='EVLA'"
    )
    url = f"{NRAO}?{urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=query))}"
    path = cache.fetch(
        url, f"metadata/{TAG}_event_nrao.csv", survey="NRAO VLA Archive",
        version="EVLA scan pointing cone and start time", kind="scan_pointing_metadata",
    )
    frame = pd.read_csv(path)
    if set(frame.columns) != set(fields.split(",")) or len(frame) >= 400:
        raise ValueError("Invalid or truncated VLA pointing result")
    return dict(
        query=receipt(path, url), radius_deg=EVENT_CONE_DEG, window_utc=[START, END],
        rows=json.loads(frame.to_json(orient="records")),
    )


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    cache = DownloadCache(ROOT, BUDGET)
    wide = download_w_cdfs(cache)
    vlass = fetch_vlass(cache, position)
    nvss = fetch_nvss(cache, position)
    event = event_search(cache, position)
    write_json(OUT / f"{TAG}.json", dict(
        executed_utc=utcnow(), script_sha256=text_digest(__file__),
        position_icrs_deg=position, position_source_sha256=text_digest(SOURCE),
        vlass=vlass, w_cdfs=wide, nvss=nvss, event_vla_pointing_search=event, raw_cache_budget_bytes=BUDGET,
        scope="The 2026 NRAO query selects scan pointing centers within 24 arcmin and scans "
        "starting in the date window; it is not an exact coverage test for wide beams beyond the "
        "cone. CADC SIA finds public VLASS footprint products, with FITS dates checked in the photometry.",
    ))


if __name__ == "__main__":
    main()

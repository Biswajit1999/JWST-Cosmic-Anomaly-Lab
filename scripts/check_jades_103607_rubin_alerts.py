"""Check public Rubin difference-image alerts around the October-November JWST event."""

import json
from urllib.parse import urlencode

from astropy.coordinates import SkyCoord
from astropy.time import Time
from check_jades_103607_independent_2026 import END, START
from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json

SOURCE = OUT / "module_b_001_native_november_fade.json"
TAG = "module_b_001_rubin_2026_alerts"
SERVER = "https://api.lsst.fink-portal.org/api/v1/conesearch"
COLUMNS = (
    "r:diaObjectId",
    "r:ra",
    "r:dec",
    "r:midpointMjdTai",
    "r:psfFlux",
    "r:psfFluxErr",
    "f:firstDiaSourceMjdTaiFink",
)


def cone(sky, radius, cache):
    query = {
        "ra": f"{sky.ra.deg:.10f}",
        "dec": f"{sky.dec.deg:.10f}",
        "radius": str(radius),
        "startdate": f"{START} 00:00:00",
        "stopdate": f"{END} 00:00:00",
        "kind": "across",
        "columns": ",".join(COLUMNS),
        "output-format": "json",
    }
    url = f"{SERVER}?{urlencode(query)}"
    path = cache.fetch(
        url,
        f"metadata/{TAG}_{radius}arcsec.json",
        survey="Rubin Observatory alerts via Fink",
        version="Public LSST Fink API v1",
        kind="observed_alert_registry",
    )
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError("Expected Fink cone-search result array")
    if len({row["r:diaObjectId"] for row in rows}) != len(rows):
        raise ValueError("Repeated Fink diaObject identifiers in one cone")
    decoded = []
    for row in rows:
        xy = [float(row["r:ra"]), float(row["r:dec"])]
        sep = float(sky.separation(SkyCoord(*xy, unit="deg")).arcsec)
        first, latest = (row["f:firstDiaSourceMjdTaiFink"], row["r:midpointMjdTai"])
        decoded.append(
            dict(
                dia_object_id=str(row["r:diaObjectId"]),
                position_icrs_deg=xy,
                separation_arcsec=sep,
                first_alert_utc=Time(first, format="mjd", scale="tai").utc.isot,
                latest_alert_utc=Time(latest, format="mjd", scale="tai").utc.isot,
                latest_psf_flux_njy=row.get("r:psfFlux"),
                latest_psf_flux_error_njy=row.get("r:psfFluxErr"),
            )
        )
    decoded.sort(key=lambda row: row["separation_arcsec"])
    return dict(radius_arcsec=radius, query=query, archive=receipt(path, url), objects=decoded)


def main():
    position = json.loads(SOURCE.read_text(encoding="utf-8"))["icrs_deg"]
    sky = SkyCoord(*position, unit="deg")
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    target = cone(sky, 5, cache)
    context = cone(sky, 120, cache)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=position,
            position_sha256=text_digest(SOURCE),
            search_window_utc=[START, END],
            target=target,
            surrounding_field=context,
            scope="Fink groups objects by diaObjectId and returns the latest alert record for "
            "objects active across the requested dates. Nearby alerts indicate operation in "
            "the field. The registry cannot establish target pixel coverage or a sensitivity limit.",
        ),
    )


if __name__ == "__main__":
    main()

"""Pin public G140M spectra for a stratified JADES Lyman-alpha pilot."""

import json
import re
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from analyse_jades_cosmic_evolution import CATALOGUE, MANIFEST, ROOT, select_sources
from astropy.coordinates import SkyCoord
from astropy.io import fits
from inspect_jades_cosmic_evolution import receipt

from cosmic_anomaly_lab.balmer_census import normalize_filter
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

OUT = ROOT / "results/jades_reionization"
LEDGER = OUT / "sample.json"
BIN_EDGES = (5.5, 6.5, 7.5, 9.)
MUV = (-22., -18.)
SEED = 518395
PER_FIELD_BIN = 18
MODE = "f070lp-g140m"


def targets(path):
    with fits.open(path) as hdus:
        obs, line = (hdus[name].data for name in ("Obs_info", "R1000_5pix"))
        fields = {
            "dr4_id": ("Unique_ID", str), "field": ("Field", str), "pid": ("PID", int),
            "ra": ("RA_TARG", float), "dec": ("Dec_TARG", float), "nircam_id": ("NIRCam_DR5_ID", int),
            "z": ("z_R1000", float), "z_spec": ("z_Spec", float), "z_lines": ("z_R1000n", float),
            "flag": ("z_Spec_flag", str), "m_uv": ("MUV", float), "exposure_s": ("tExp_G140M", float),
        }
        frame = pd.DataFrame({key: np.asarray(obs[name], dtype=kind) for key, (name, kind) in fields.items()})
        idx = pd.Index(np.asarray(line["Unique_ID"], dtype=str)).get_indexer(frame.dr4_id)
        if len(line) != len(frame) or np.any(idx < 0) or not frame.dr4_id.is_unique:
            raise ValueError("Target IDs fail the observation/line catalogue join")
        f, e, s = (np.asarray(line[f"O3_5007_{name}"][idx], dtype=float).reshape(-1) for name in ("flux", "flux_err", "SNR"))
        modes = [normalize_filter(mode) or "" for mode in line["O3_5007_filter"][idx]]
    frame = frame.assign(o3_snr=np.divide(f, e, out=np.full_like(f, np.nan), where=e > 0),
                         o3_archive_snr=s, o3_filter=modes)
    spectral = frame[(frame.flag == "A") & frame.z.ge(BIN_EDGES[0]) & frame.z.lt(BIN_EDGES[-1])
                     & frame.z_lines.ge(2) & ((frame.z-frame.z_spec).abs()/(1+frame.z)).le(.002)
                     & frame.o3_filter.str.startswith(("g235m_f170lp", "g395m_f290lp"))
                     & frame.o3_snr.ge(5) & frame.o3_archive_snr.ge(5)
                     & frame.m_uv.between(*MUV) & frame.exposure_s.gt(0)]
    sources, duplicate_audit = select_sources(spectral)
    sources = sources.sort_values("dr4_id").copy()
    sources["bin"] = pd.cut(sources.z, BIN_EDGES, right=False, labels=False)
    return sources, duplicate_audit


def listing_for(cache, tier):
    directory = f"https://jades.herts.ac.uk/DR4/{tier[:7]}/spectra/{MODE}/{tier}/"
    listing = cache.fetch(directory, f"metadata/dr4_{tier}_g140m_f070lp.html",
                          survey="JADES", version="DR4-listing", kind="archive_directory")
    links = re.findall(r'href="([^"]+)"', listing.read_text())
    index = {}
    for href in links:
        match = re.search(r"-(\d+)_" + MODE + r"_.*_x1d\.fits$", href)
        if match:
            index.setdefault(int(match[1]), []).append(href)
    return directory, listing, index


def acquire(cache, listings, row):
    tier, number = row["dr4_id"].rsplit("_", 1)
    number = int(number)
    if tier not in listings:
        listings[tier] = listing_for(cache, tier)
    directory, listing, index = listings[tier]
    links = index.get(number, [])
    if len(links) != 1:
        raise ValueError(f"{len(links)} listed G140M/F070LP 1D spectra")
    url = urljoin(directory, links[0])
    path = cache.fetch(url, f"independent/{Path(links[0]).name}", survey="JADES", version="DR4", kind="observed_spectrum")
    with fits.open(path) as h:
        header = h[0].header
        sep = SkyCoord(row["ra"], row["dec"], unit="deg").separation(
            SkyCoord(header["RA"], header["DEC"], unit="deg")
        ).arcsec
        pid, target = int(header["PROPOSID"]), int(header["TARGET"])
    product = dict(listing=receipt(listing, directory), spectrum=receipt(path, url),
                   header_pid=pid, header_target_number=target, offset_arcsec=float(sep))
    variant = tier == "goods-s-mediumjwst" and (row["pid"], pid) == (1180, 1286)
    if sep > .2 or target != number or not (pid == row["pid"] or variant):
        return dict(status="HEADER_MISMATCH", **product)
    association = "DR4_GOODS_S_MEDIUMJWST_HEADER_1286_CATALOGUE_1180" if variant else "CATALOGUE_PID_MATCH"
    return dict(status="ACCEPTED", association=association, **product)


def main():
    cache = DownloadCache(ROOT, 28_000*1024**2)
    path = cache.fetch(CATALOGUE, "jades_dr4.fits", survey="JADES DR4", version="v1.2.1", kind="observed_catalogue")
    catalogue = next(row for row in reversed(json.loads(MANIFEST.read_text())) if row["url"] == CATALOGUE)
    if digest(path) != catalogue["sha256"]:
        raise ValueError("Catalogue fingerprint mismatch")
    sources, duplicates = targets(path)
    OUT.mkdir(parents=True, exist_ok=True)
    roster = OUT / "eligible_targets.csv"
    sources.to_csv(roster, index=False)
    rng = np.random.default_rng(SEED)
    picked, counts = [], []
    for (field, redshift_bin), group in sources.groupby(["field", "bin"], sort=True):
        draws = group.iloc[np.sort(rng.choice(len(group), min(PER_FIELD_BIN, len(group)), replace=False))]
        picked.extend(draws.to_dict(orient="records"))
        counts.append(dict(field=field, bin=int(redshift_bin), available=len(group), selected=len(draws)))
    listings = {}
    ledger = dict(executed_utc=utcnow(), catalogue=catalogue, target_roster_sha256=text_digest(roster),
                  script_sha256=text_digest(__file__), strata=counts, duplicate_audit=duplicates,
                  planned=len(picked), redshift_bins=BIN_EDGES, uv_cut=MUV, seed=SEED,
                  naming_reference="https://archive.stsci.edu/hlsp/jades", records=[])
    for row in picked:
        try:
            record = acquire(cache, listings, row)
        except (ValueError, OSError, KeyError, requests.RequestException) as exc:
            record = dict(status="NO_ACCEPTED_PRODUCT", reason=str(exc))
        ledger["records"].append(dict(**row, **record))
        write_json(LEDGER, ledger)
        print(len(ledger["records"]), row["dr4_id"], record["status"], flush=True)


if __name__ == "__main__":
    main()

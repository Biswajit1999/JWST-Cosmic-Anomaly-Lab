"""Retrieve explicitly reviewed DR4 associations from actual archive listings."""

import re
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json

root = Path(__file__).resolve().parents[1]
config = yaml.safe_load((root / "configs/reduction_audit.yaml").read_text())
cache = DownloadCache(root, budget_bytes=config["raw_storage_budget_MiB"] * 1024**2)
reference = pd.read_csv(root / "data/catalogues/jades_dr4_reference.csv.gz")
selection = pd.concat(
    [
        pd.read_csv(root / "data/catalogues/expanded_selection.csv.gz"),
        pd.read_csv(root / "data/catalogues/open_selection.csv.gz"),
    ]
)
out = root / "results/reduction_audit"
records = []
for target in config["targets"]:
    r = reference[reference.Unique_ID == target["dr4_id"]]
    s = selection[selection.source_id == target["source_id"]]
    assert len(r) == len(s) == 1
    r, s = r.iloc[0], s.iloc[0]
    separation = (
        SkyCoord(r.ra_deg * u.deg, r.dec_deg * u.deg)
        .separation(SkyCoord(s.ra * u.deg, s.dec * u.deg))
        .arcsec
    )
    assert separation <= config["coordinate_tolerance_arcsec"]
    assert int(r.PID) == target["pid"]
    tier, number = target["dr4_id"].rsplit("_", 1)
    field = tier[:7]
    directory = f"https://jades.herts.ac.uk/DR4/{field}/spectra/clear-prism/{tier}/"
    record = dict(
        **target,
        association_separation_arcsec=float(separation),
        dja_file=s.file,
        archive_z=float(s.z_best),
        dr4_z=float(r.z_Spec),
        dr4_z_flag=r.z_Spec_flag,
        directory_url=directory,
    )
    try:
        listing = cache.fetch(
            directory,
            f"metadata/dr4_{tier}_prism.html",
            survey="JADES",
            version="DR4-listing-snapshot",
            kind="archive_directory",
        )
        links = re.findall(r'href="([^"]+)"', listing.read_text())
        matches = [
            link
            for link in links
            if (match := re.search(r"-(\d+)_clear-prism_", link))
            and int(match.group(1)) == int(number)
            and link.endswith("_x1d.fits")
        ]
        if len(matches) != 1:
            raise ValueError(f"Expected one listed x1d product, found {len(matches)}")
        name = matches[0]
        path = cache.fetch(
            urljoin(directory, name),
            f"independent/{name}",
            survey="JADES",
            version="DR4",
            kind="independently_reduced_observed_spectrum",
        )
        record.update(
            status="DOWNLOADED", path=path.relative_to(root).as_posix(), raw_sha256=digest(path)
        )
    except Exception as exc:
        record.update(status="UNAVAILABLE", reason=str(exc))
    records.append(record)
    write_json(
        out / "downloads.json",
        dict(
            executed_utc=utcnow(),
            config_sha256=text_digest(root / "configs/reduction_audit.yaml"),
            targets=records,
            not_associated=config["not_associated"],
        ),
    )
    print(target["dr4_id"], record["status"], flush=True)

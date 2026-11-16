"""Retrieve actual listed spectra and remeasure signed Balmer cores for the queue."""

import re
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import yaml
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_census import signed_core_integral
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/jades_balmer_census"
queue = pd.read_csv(out / "followup_queue.csv")
table = pd.read_csv(out / "aperture_measurements.csv.gz")
config = yaml.safe_load((root / "configs/jades_balmer_census.yaml").read_text())
cache = DownloadCache(root, budget_bytes=config["raw_storage_budget_MiB"] * 1024**2)
rest = dict(ha=.656452, hb=.486264, hg=.434165)
products, measurements = [], []
for target in queue.itertuples():
    tier, number = target.dr4_id.rsplit("_", 1)
    rows = table[table.dr4_id == target.dr4_id].set_index("aperture_pixels")
    labels = {str(rows.loc[a, line + "_filter"]).removesuffix("_ext")
              for a in [3, 5] for line in rest if pd.notna(rows.loc[a, line + "_filter"])}
    spectra = {}
    for label in sorted(labels):
        grating, filt = label.split("_")
        directory = f"https://jades.herts.ac.uk/DR4/{tier[:7]}/spectra/{filt}-{grating}/{tier}/"
        record = dict(dr4_id=target.dr4_id, grating_filter=label, directory_url=directory)
        try:
            listing = cache.fetch(directory, f"metadata/dr4_{tier}_{grating}_{filt}.html",
                                  survey="JADES", version="DR4-listing", kind="archive_directory")
            links = re.findall(r'href="([^\"]+)"', listing.read_text())
            matches = [link for link in links if link.endswith("_x1d.fits")
                       and (m := re.search(r"-(\d+)_" + re.escape(filt + "-" + grating) + "_", link))
                       and int(m.group(1)) == int(number)]
            if len(matches) != 1:
                raise ValueError(f"Expected one listed spectrum, found {len(matches)}")
            url = urljoin(directory, matches[0])
            path = cache.fetch(url, "independent/" + matches[0], survey="JADES", version="DR4",
                               kind="observed_spectrum")
            reference = read_dr4_spectrum(path, wavelength_range=(.6, 5.5))
            separation = SkyCoord(reference.metadata["RA"] * u.deg, reference.metadata["DEC"] * u.deg).separation(
                SkyCoord(target.ra * u.deg, target.dec * u.deg)).arcsec
            record.update(path=path.relative_to(root).as_posix(), sha256=digest(path), url=url,
                          separation_arcsec=float(separation), listing_sha256=digest(listing),
                          header_pid=int(reference.metadata["PROPOSID"]), catalogue_pid=target.pid,
                          header_target=str(reference.metadata["TARGET"]))
            if (separation > .2 or int(reference.metadata["TARGET"]) != int(number)
                    or int(reference.metadata["PROPOSID"]) != target.pid):
                raise ValueError("Target/PID/header-coordinate check failed")
            record.update(status="DOWNLOADED", path=path.relative_to(root).as_posix(),
                          sha256=digest(path), url=url, separation_arcsec=float(separation),
                          listing_sha256=digest(listing))
            for a in [3, 5]:
                for nods in [2, 3]:
                    spectra[(label, a, nods)] = read_dr4_spectrum(path, a, nods, wavelength_range=(.6, 5.5))
        except (OSError, ValueError, KeyError) as exc:
            record.update(status="ASSOCIATION_REVIEW_REQUIRED" if "header_pid" in record
                          else "UNAVAILABLE", reason=str(exc))
        products.append(record)
        write_json(out / "spectral_products.json", products)
        print(target.dr4_id, label, record["status"], flush=True)
    fig = Figure(figsize=(15, 8), layout="constrained")
    axes = fig.subplots(2, 3)
    for col, (line, wavelength) in enumerate(rest.items()):
        center = wavelength * (1 + target.z)
        for a, color in [(3, "#227f8e"), (5, "#b25835")]:
            label = str(rows.loc[a, line + "_filter"]).removesuffix("_ext")
            for nods in [3, 2]:
                spec = spectra.get((label, a, nods))
                for half in [400, 600]:
                    record = dict(dr4_id=target.dr4_id, line=line, aperture_pixels=a, nods=nods,
                                  grating_filter=label, fixed_z=float(target.z),
                                  center_um=float(center), beyond_previous_5p2um_limit=bool(center > 5.2),
                                  half_width_kms=half)
                    if spec is None:
                        record.update(status="NO_ACCEPTED_PRODUCT", flux=None, error=None)
                    else:
                        velocity = (spec.wavelength.value / center - 1) * 299792.458
                        flux = spec.flux.to_value(u.erg / u.s / u.cm**2 / u.AA,
                                                   equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
                        error = spec.error.to_value(u.erg / u.s / u.cm**2 / u.AA,
                                                     equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
                        record.update(signed_core_integral(velocity, flux, error, spec.good, center * 1e4, half))
                        if half == 600:
                            ax = axes[0 if nods == 3 else 1, col]
                            good = spec.good & (abs(velocity) <= 3000)
                            if good.any():
                                ax.errorbar(velocity[good], flux[good], yerr=error[good],
                                            fmt=".-", ms=2, lw=.5, color=color, alpha=.7,
                                            label=f"{a}px, {nods} nods")
                    measurements.append(record)
        for row in range(2):
            ax = axes[row, col]
            ax.axvspan(-600, 600, color="#f1dfc1", alpha=.35)
            ax.axhline(0, color="gray", lw=.6)
            ax.set(xlim=(-3000, 3000), xlabel="Velocity at catalogue z (km/s)",
                   ylabel="Fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)", title=f"{line.upper()} · {center:.3f} µm")
            handles, _ = ax.get_legend_handles_labels()
            if handles:
                ax.legend(fontsize=8)
            else:
                ax.text(.5, .6, "No accepted usable pixels", transform=ax.transAxes, ha="center")
    fig.suptitle(f"OBSERVED {target.dr4_id} · z={target.z:.5f}\nSigned line cores; apertures/nod combinations share photons; no physical classification", fontsize=12)
    save_figure(fig, root / "figures/jades_balmer_sources" / target.dr4_id)
pd.DataFrame(measurements).to_csv(out / "signed_core_measurements.csv", index=False)
write_json(out / "spectral_summary.json", dict(executed_utc=utcnow(), research_snapshot="2026-11-16",
    target_records=len(queue), listed_products=len(products),
    downloaded_products=sum(r["status"] == "DOWNLOADED" for r in products),
    association_review_products=sum(r["status"] == "ASSOCIATION_REVIEW_REQUIRED" for r in products),
    attempted_core_measurements=len(measurements),
    measured_cores=sum(r["status"] == "MEASURED_CORE" for r in measurements),
    analysis_wavelength_range_um=[.6, 5.5], flux_unit="1e-20 erg/s/cm2",
    aperture_and_nod_variants_independent=False, physical_classification="NOT_ESTABLISHED",
    input_sha256={"scripts/inspect_jades_balmer_queue.py": text_digest(Path(__file__)),
        "src/cosmic_anomaly_lab/balmer_census.py": text_digest(root / "src/cosmic_anomaly_lab/balmer_census.py"),
        "src/cosmic_anomaly_lab/reduction_audit.py": text_digest(root / "src/cosmic_anomaly_lab/reduction_audit.py"),
        "results/jades_balmer_census/followup_queue.csv": text_digest(out / "followup_queue.csv"),
        "results/jades_balmer_census/aperture_measurements.csv.gz": digest(out / "aperture_measurements.csv.gz")},
))

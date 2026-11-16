"""Check four DR4 priorities against released DJA medium-resolution reductions."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_census import signed_core_integral
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/multiwavelength"
out.mkdir(parents=True, exist_ok=True)
queue = pd.read_csv(root / "results/jades_balmer_census/followup_queue.csv")
apertures = pd.read_csv(root / "results/jades_balmer_census/aperture_measurements.csv.gz")
obs = pd.read_csv(root / "data/catalogues/observations.csv.gz")
dr4 = json.loads((root / "results/jades_balmer_census/spectral_products.json").read_text())
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
rest = dict(ha=.656452, hb=.486264, hg=.434165)
selection = {66336: ["G395M"], 89464: ["G395M"],
             10009106: ["G140M", "G235M"], 25734: ["G140M", "G235M"]}
products, cores, ratios = [], [], []
for number, gratings in selection.items():
    target = queue[queue.dr4_id.str.endswith("_" + str(number))].iloc[0]
    coordinate = SkyCoord(target.ra * u.deg, target.dec * u.deg)
    spectra = {}
    for grating in gratings:
        candidates = obs[(obs.srcid == number) & (obs.programme == target.pid)
                         & (obs.grating == grating)].copy()
        distances = coordinate.separation(SkyCoord(candidates.ra.to_numpy() * u.deg,
                                                    candidates.dec.to_numpy() * u.deg)).arcsec
        candidates = candidates[distances < .2]
        if len(candidates) != 1:
            raise ValueError(f"Ambiguous/missing product: {target.dr4_id} {grating}")
        row = candidates.iloc[0]
        path = cache.fetch(row.url, "followup/" + row.file, survey="DJA-JADES",
                           version=str(row.version), kind="observed_spectrum")
        spec = read_spectrum(path, wavelength_range=(.6, 5.5))
        separation = coordinate.separation(SkyCoord(spec.metadata["SRCRA"] * u.deg,
                                                     spec.metadata["SRCDEC"] * u.deg)).arcsec
        if separation > .2 or spec.metadata["GRATING"].upper() != grating:
            raise ValueError("Header identity mismatch")
        products.append(dict(dr4_id=target.dr4_id, grating=grating, url=row.url,
                             path=path.relative_to(root).as_posix(), sha256=digest(path),
                             separation_arcsec=float(separation), programme=int(target.pid),
                             source_number=number, status="VERIFIED"))
        spectra[grating] = spec
        print(row.file, "VERIFIED", flush=True)
    fig = Figure(figsize=(12, 4), layout="constrained")
    axes = fig.subplots(1, 3)
    for ax, (line, wavelength) in zip(axes, rest.items(), strict=True):
        aperture = apertures[(apertures.dr4_id == target.dr4_id) & (apertures.aperture_pixels == 5)].iloc[0]
        label = str(aperture[line + "_filter"]).removesuffix("_ext")
        grating = label.split("_")[0].upper()
        center = wavelength * (1 + target.z)
        alternatives = [("DJA", spectra.get(grating), "#126e82")]
        for p in dr4:
            if p["dr4_id"] == target.dr4_id and p["grating_filter"] == label and p["status"] == "DOWNLOADED":
                for nods, color in [(3, "#b94d36"), (2, "#9b8b4e")]:
                    alternatives.append((f"DR4 5px/{nods}nod", read_dr4_spectrum(
                        root / p["path"], 5, nods, wavelength_range=(.6, 5.5)), color))
        for name, spec, color in alternatives:
            if spec is None:
                continue
            v = (spec.wavelength.value / center - 1) * 299792.458
            unit = u.erg / u.s / u.cm**2 / u.AA
            f = spec.flux.to_value(unit, equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
            e = spec.error.to_value(unit, equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
            good = spec.good & (abs(v) <= 3000)
            if good.any():
                ax.errorbar(v[good], f[good], yerr=e[good], fmt=".-", ms=2, lw=.5,
                            color=color, alpha=.75, label=name)
            for half in [400, 600]:
                cores.append(dict(dr4_id=target.dr4_id, line=line, reduction=name,
                                  fixed_z=float(target.z), grating=grating,
                                  **signed_core_integral(v, f, e, spec.good, center * 1e4, half)))
        ax.axvspan(-600, 600, color="#d9c696", alpha=.2)
        ax.axhline(0, color="gray", lw=.6)
        ax.set(xlim=(-3000, 3000), xlabel="Velocity (km/s)", title=f"{line.upper()} | {center:.3f} µm",
               ylabel="Fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7)
    fig.suptitle(f"{target.dr4_id} | observed spectra, fixed z={target.z:.5f}\nNo flux rescaling; reductions share observations", fontsize=11)
    save_figure(fig, root / "figures/multiwavelength" / (target.dr4_id + "_reductions"))
data = pd.DataFrame(cores)
for (sid, reduction, half), block in data.groupby(["dr4_id", "reduction", "half_width_kms"]):
    b = block.set_index("line")
    for line in ["ha", "hg"]:
        if not {line, "hb"}.issubset(b.index):
            continue
        n, d = b.loc[line], b.loc["hb"]
        usable = n.status == d.status == "MEASURED_CORE"
        eligible = bool(usable and n.flux / n.error >= 5 and d.flux / d.error >= 5)
        ratios.append(dict(dr4_id=sid, reduction=reduction, half_width_kms=int(half), line=line,
                           eligible=eligible, numerator_snr=float(n.flux / n.error) if usable else None,
                           denominator_snr=float(d.flux / d.error) if usable else None,
                           ratio=float(n.flux / d.flux) if eligible else None,
                           error=float(np.hypot(n.error / d.flux, n.flux * d.error / d.flux**2)) if eligible else None))
data.to_csv(out / "dja_core_measurements.csv", index=False)
pd.DataFrame(ratios).to_csv(out / "dja_core_ratios.csv", index=False)
write_json(out / "dja_products.json", products)
write_json(out / "dja_summary.json", dict(executed_utc=utcnow(), target_records=4,
    products=len(products), core_attempts=len(data), measured_cores=int((data.status == "MEASURED_CORE").sum()),
    script_sha256=text_digest(Path(__file__)), raw_storage_budget_bytes=4096*1024**2,
    flux_unit="1e-20 erg/s/cm2", analysis_wavelength_range_um=[.6, 5.5],
    input_hash_convention="Raw bytes for gzip; canonical UTF-8/LF for repository text",
    input_sha256={p: (digest(root / p) if p.endswith(".gz") else text_digest(root / p)) for p in [
        "data/catalogues/observations.csv.gz", "results/jades_balmer_census/followup_queue.csv",
        "results/jades_balmer_census/aperture_measurements.csv.gz",
        "results/jades_balmer_census/spectral_products.json"]},
    independent_photons=False, physical_classification="NOT_ESTABLISHED"))

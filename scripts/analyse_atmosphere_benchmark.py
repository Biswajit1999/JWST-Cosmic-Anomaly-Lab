"""Published WASP-39 b spectrum: feature check and conditional bulk properties."""

import io
import json
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from astropy import constants as c
from matplotlib.figure import Figure

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/exoplanet_pilot"
product = json.loads((out / "atmosphere_product.json").read_text())
assert digest(root / product["path"]) == product["sha256"]
member = "G395H_paper_data/3_TRANSMISSION_SPECTRA/transit-spectrum-W39b-G395H-10pix_weighted-average.nc"
with zipfile.ZipFile(root / product["path"]) as z:
    with h5py.File(io.BytesIO(z.read(member)), "r") as h:
        wave, depth, error = [np.array(h[k]) for k in ["central_wavelength", "transit_depth", "transit_depth_error"]]
        declared_units = h["transit_depth"].attrs["units"].decode()
    # Read text models only. No pickle or downloaded notebook execution.
    models = {}
    for label, filename in [("PHOENIX", "fitted-model-spectrum-WASP-39b-PHOENIX-eq-Lothringer-Fig3.txt"),
                            ("PICASO", "fitted-model-spectrum-WASP-39b-PICASO-eq-Batalha-Fig3.txt")]:
        models[label] = np.loadtxt(io.BytesIO(z.read("G395H_paper_data/4_THEORY/"+filename)), delimiter=",")
valid = np.isfinite(wave+depth+error) & (error > 0)
assert valid.all() and np.all(np.diff(wave) > 0)
assert np.median(depth) > .01 and np.median(depth) < .04
# Stored values follow the release notebook's *100 -> percent convention.
# The HDF unit string is ambiguous: do NOT apply a Jupiter/Solar radius conversion.
pd.DataFrame(dict(wavelength_um=wave, fractional_depth=depth, fractional_error=error)).to_csv(out / "wasp39_transmission.csv", index=False)
band = (wave >= 4.20) & (wave <= 4.45)
design = np.column_stack([np.ones(len(wave)), wave-4.3])
results = []
fig = Figure(figsize=(12, 6), layout="constrained")
axes = fig.subplots(1, 2)
axes[0].errorbar(wave, depth*1e6, yerr=error*1e6, fmt=".", ms=2, lw=.5, color="#173c59", label="Published weighted spectrum")
for (label, model), color in zip(models.items(), ["#bc6a38", "#308c72"], strict=True):
    # Split the detector gap; do not draw invented data between detectors.
    for j, group in enumerate(np.split(np.arange(len(model)), np.where(np.diff(model[:, 0]) > .04)[0]+1)):
        axes[0].plot(model[group, 0], model[group, 1]*1e6, lw=.8, color=color, label=f"Published {label} model" if j == 0 else None)
axes[0].axvspan(4.20, 4.45, alpha=.12, color="gold")
axes[0].set(xlabel="Wavelength (µm)", ylabel="Transit depth (ppm)", xlim=(2.9, 5.15),
            title="Observed transmission spectrum; published models")
axes[0].legend(fontsize=7)
for name, ranges in [("reference_A", [(3.80, 4.00), (4.50, 4.60)]),
                     ("reference_B", [(3.90, 4.00), (4.50, 4.55)])]:
    side = np.logical_or.reduce([(wave >= a) & (wave <= b) for a, b in ranges])
    assert side.sum() >= 6 and band.sum() >= 6
    cov = np.linalg.inv((design[side].T/error[side]**2) @ design[side])
    beta = cov @ ((design[side].T/error[side]**2) @ depth[side])
    weights = (1/error[band]**2)
    weights /= weights.sum()
    cw = weights @ design[band]
    excess = weights @ depth[band]-cw @ beta
    unc = np.sqrt(np.sum((weights*error[band])**2)+cw @ cov @ cw)
    results.append(dict(reference=name, band_bins=int(band.sum()), continuum_bins=int(side.sum()),
                         excess_ppm=float(excess*1e6), diagonal_error_ppm=float(unc*1e6),
                         formal_contrast=float(excess/unc), molecule_detection_significance=False))
    axes[1].plot(wave, (design @ beta)*1e6, label=name, lw=1)
axes[1].errorbar(wave, depth*1e6, yerr=error*1e6, fmt=".", color="#173c59", ms=3, lw=.6)
axes[1].axvspan(4.20, 4.45, alpha=.14, color="gold")
axes[1].set(xlim=(3.75, 4.65), xlabel="Wavelength (µm)", ylabel="Transit depth (ppm)",
            title="4.2–4.45 µm excess above two smooth references")
axes[1].legend(fontsize=8)
fig.suptitle("WASP-39 b · known-planet atmospheric benchmark\nReanalysis of published JWST measurements; not a new planet or a new molecular detection", fontsize=12)
save_figure(fig, root / "figures/exoplanet_pilot/wasp39_atmosphere")
pd.DataFrame(results).to_csv(out / "wasp39_feature_contrasts.csv", index=False)
p = json.loads((root / "data/raw/exoplanets/wasp39_parameters.json").read_text())[0]
mass, radius = p["pl_bmassj"]*c.M_jup.value, p["pl_radj"]*c.R_jup.value
em, er = p["pl_bmassjerr1"]/p["pl_bmassj"], p["pl_radjerr1"]/p["pl_radj"]
rho = mass/(4*np.pi*radius**3/3)/1000
g = c.G.value*mass/radius**2
write_json(out / "wasp39_properties.json", dict(source="NASA Exoplanet Archive composite snapshot; heterogeneous literature measurements",
    retrieved_parameters=p, density_g_cm3=float(rho), density_diagonal_error=float(rho*np.hypot(em, 3*er)),
    gravity_m_s2=float(g), gravity_diagonal_error=float(g*np.hypot(em, 2*er)),
    independence_of_mass_radius_errors_assumed=True,
    new_abundances_retrieved=False, new_planet_discovered=False))
write_json(out / "atmosphere_summary.json", dict(executed_utc=utcnow(), released_spectral_bins=len(wave),
    source_doi=product["doi"], archive_member=member, raw_sha256=product["sha256"],
    script_sha256=text_digest(Path(__file__)), hdf_declared_depth_units=declared_units,
    depth_convention="Fractional depth verified against release notebook percent plotting; no Rjup/Rsun conversion",
    unit_metadata_ambiguity_documented=True, covariance_supplied=False,
    chemical_identification="Published attribution of 4.3 micron feature to CO2; this test measures only band contrast",
    known_benchmark=True))
print(pd.DataFrame(results).to_string(index=False))
print('Density and gravity',rho,g)

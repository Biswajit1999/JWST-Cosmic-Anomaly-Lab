"""Audit released DR4 intermediate extractions and leave-one-out combinations."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_census import signed_core_integral
from cosmic_anomaly_lab.exposure_audit import combine_extractions
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/exposure_audit"
out.mkdir(parents=True, exist_ok=True)
DownloadCache(root, budget_bytes=4096*1024**2).fetch(
    "https://jades.herts.ac.uk/DR4/Readme_DR4_spectra.md", "metadata/jades_dr4_spectra.md",
    survey="JADES", version="DR4", kind="spectrum_documentation")
dr4 = json.loads((root / "results/jades_balmer_census/spectral_products.json").read_text())
dja = json.loads((root / "results/multiwavelength/dja_products.json").read_text())
targets = pd.read_csv(root / "results/jades_balmer_census/followup_queue.csv")
unit = u.erg / u.s / u.cm**2 / u.AA
rest = dict(ha=.656452, hb=.486264, hg=.434165)
measurements, ratios, inventory, corrections, input_hashes, coverage_selection = [], [], [], [], {}, []
for number, grating, diagnostic in [(66336, "g395m_f290lp", "ha"),
                                   (10009106, "g140m_f070lp", "hg"),
                                   (89464, "g395m_f290lp", "hg")]:
    target = targets[targets.dr4_id.str.endswith("_" + str(number))].iloc[0]
    dp = next(p for p in dr4 if p["dr4_id"] == target.dr4_id and p["grating_filter"] == grating)
    jp = next(p for p in dja if p["dr4_id"] == target.dr4_id and p["grating"].lower() == grating.split("_")[0])
    for p in [dp, jp]:
        assert digest(root / p["path"]) == p["sha256"]
        input_hashes[p["path"]] = p["sha256"]
    with fits.open(root / dp["path"]) as hdus:
        hdu = hdus["INTERMEDIATE"]
        wave = (hdu.data["WAVELENGTH"] * u.Unit(hdu.columns["WAVELENGTH"].unit)).to_value(u.um)
        names = [n for n in hdu.columns.names if n.startswith("FLUX_")]
        flux, errors = [], []
        for name in names:
            error_name = name.replace("FLUX_", "ERR_", 1)
            flux.append((hdu.data[name] * u.Unit(hdu.columns[name].unit)).to_value(unit) / 1e-20)
            errors.append((hdu.data[error_name] * u.Unit(hdu.columns[error_name].unit)).to_value(unit) / 1e-20)
            inventory.append(dict(dr4_id=target.dr4_id, column=name, pointing_integration_region=name[5:],
                                  path=dp["path"], grating_filter=grating))
    flux, errors = np.array(flux), np.array(errors)
    variants = [(name, "individual", wave, f, e) for name, f, e in zip(names, flux, errors, strict=True)]
    full, err, contributors = combine_extractions(flux, errors)
    variants.append(("DIAGNOSTIC_IVAR_ALL", "combined", wave, full, err))
    for i, name in enumerate(names):
        f, e, _ = combine_extractions(flux, errors, omit=i)
        variants.append(("OMIT_" + name, "leave_one_out", wave, f, e))
    # Geometry/validity selection only: never select on ratio or positive flux.
    selected = []
    for i, name in enumerate(names):
        statuses = []
        for line in [diagnostic, "hb"]:
            center = rest[line]*(1+target.z)
            v = (wave / center - 1)*299792.458
            good = np.isfinite(flux[i]+errors[i]) & (errors[i] > 0) & (wave >= .6) & (wave <= 5.5)
            statuses.append(signed_core_integral(v, flux[i], errors[i], good, center*1e4, 600)["status"])
        accepted = all(s == "MEASURED_CORE" for s in statuses)
        coverage_selection.append(dict(dr4_id=target.dr4_id, column=name, accepted=accepted,
                                       numerator_status=statuses[0], denominator_status=statuses[1]))
        if accepted:
            selected.append(i)
    if selected:
        f, e, _ = combine_extractions(flux[selected], errors[selected])
        variants.append(("COMMON_COVERAGE", "common_coverage", wave, f, e))
    for label, reader, product in [("DR4_RELEASED", read_dr4_spectrum, dp), ("DJA_RELEASED", read_spectrum, jp)]:
        spec = reader(root / product["path"], wavelength_range=(.6, 5.5))
        f = spec.flux.to_value(unit, equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
        e = spec.error.to_value(unit, equivalencies=u.spectral_density(spec.wavelength)) / 1e-20
        variants.append((label, "released", spec.wavelength.value,
                         np.where(spec.good, f, np.nan), np.where(spec.good, e, np.nan)))
    for label, kind, w, f, e in variants:
        for line, wavelength in rest.items():
            center = wavelength * (1 + target.z)
            v = (w / center - 1) * 299792.458
            good = np.isfinite(f+e) & (e > 0) & (w >= .6) & (w <= 5.5)
            for half in [400, 600]:
                measurements.append(dict(dr4_id=target.dr4_id, variant=label, kind=kind, line=line,
                    fixed_z=float(target.z), **signed_core_integral(v, f, e, good, center*1e4, half)))
    with fits.open(root / jp["path"]) as hdus:
        table = hdus["SPEC1D"].data
        for line, wavelength in rest.items():
            center = wavelength * (1 + target.z)
            use = (abs((table["wave"] / center - 1) * 299792.458) <= 600)
            use &= np.isfinite(table["path_corr"]) & (table["path_corr"] > 0) & (table["npix"] > 0)
            corrections.append(dict(dr4_id=target.dr4_id, line=line,
                median_path_corr=float(np.median(table["path_corr"][use])) if use.any() else None,
                minimum_path_corr=float(np.min(table["path_corr"][use])) if use.any() else None,
                maximum_path_corr=float(np.max(table["path_corr"][use])) if use.any() else None,
                status="COLUMN_INSPECTION_ONLY_NOT_APPLIED"))
    fig = Figure(figsize=(12, 7), layout="constrained")
    axes = fig.subplots(2, 2)
    for col, line in enumerate([diagnostic, "hb"]):
        center = rest[line]*(1+target.z)
        for label, kind, w, f, e in variants:
            if kind == "leave_one_out":
                continue
            v = (w / center - 1) * 299792.458
            valid = np.isfinite(f+e) & (e > 0) & (abs(v) <= 1800)
            if not valid.any():
                continue
            color = {"DR4_RELEASED": "#b94d36", "DJA_RELEASED": "#126e82", "DIAGNOSTIC_IVAR_ALL": "black"}.get(label)
            ax = axes[0 if kind == "individual" else 1, col]
            window = abs(v) <= 1800
            ax.plot(v[window], np.where(valid[window], f[window], np.nan), lw=.8,
                    color=color, alpha=.7, label=label)
        for row in range(2):
            ax = axes[row, col]
            ax.axhline(0, color="gray", lw=.5)
            ax.axvspan(-600, 600, alpha=.12, color="gold")
            ax.set(xlim=(-1800, 1800), xlabel="Velocity at fixed catalogue z (km/s)",
                   ylabel="Fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)", title=f"{line.upper()} | {'individual extractions' if row == 0 else 'combinations'}")
            ax.legend(fontsize=6, ncol=2 if row == 0 else 1)
    fig.suptitle(f"{target.dr4_id}: released exposure-level audit\nIVAR combination has no clipping or covariance correction; not the survey stacking algorithm", fontsize=11)
    save_figure(fig, root / "figures/exposure_audit" / target.dr4_id)
data = pd.DataFrame(measurements)
for (sid, variant, kind, half), group in data.groupby(["dr4_id", "variant", "kind", "half_width_kms"]):
    group = group.set_index("line")
    line = "ha" if sid.endswith("_66336") else "hg"
    n, d = group.loc[line], group.loc["hb"]
    measured = n.status == d.status == "MEASURED_CORE"
    eligible = bool(measured and n.flux/n.error >= 5 and d.flux/d.error >= 5)
    ratios.append(dict(dr4_id=sid, variant=variant, kind=kind, half_width_kms=int(half), line=line,
        eligible=eligible, ratio=float(n.flux/d.flux) if eligible else None,
        numerator_snr=float(n.flux/n.error) if measured else None,
        denominator_snr=float(d.flux/d.error) if measured else None))
data.to_csv(out / "signed_cores.csv", index=False)
r = pd.DataFrame(ratios)
r.to_csv(out / "core_ratios.csv", index=False)
pd.DataFrame(inventory).to_csv(out / "intermediate_inventory.csv", index=False)
pd.DataFrame(corrections).to_csv(out / "dja_path_correction_columns.csv", index=False)
pd.DataFrame(coverage_selection).to_csv(out / "common_coverage_selection.csv", index=False)
summary = dict(executed_utc=utcnow(), targets=3, intermediate_extractions=len(inventory),
    spectrum_variants=len(data)//6, core_attempts=len(data), usable_cores=int((data.status == "MEASURED_CORE").sum()),
    input_sha256=input_hashes, script_sha256=text_digest(Path(__file__)),
    module_sha256=text_digest(root / "src/cosmic_anomaly_lab/exposure_audit.py"),
    documentation_sha256=digest(root / "data/raw/metadata/jades_dr4_spectra.md"),
    common_coverage_selection="Both lines satisfy coverage/pixel/sideband rules at 600 km/s; no S/N or ratio cut",
    independent_detection_trials=False, instrumental_cause_identified=False,
    physical_classification="NOT_ESTABLISHED")
write_json(out / "summary.json", summary)
print(r[r.kind != "individual"].to_string(index=False))

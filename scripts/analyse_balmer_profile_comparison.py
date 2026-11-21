"""Audit the Hγ/Hβ contrast with empirical, pixel-integrated line models."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_profiles import C_KMS, REST, fit_pair, pair_model
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/balmer_profiles"
out.mkdir(parents=True, exist_ok=True)
prev = root / "results/gn89464_followup"
legacy = json.loads((prev / "summary.json").read_text())
products = pd.read_csv(prev / "identity_checks.csv").query("status == 'ACCEPTED'").to_dict("records")
neighbor = json.loads((root / "results/spatial_context/neighbor_product.json").read_text())
products.append(dict(product="DJA_neighbor_G395M", source_id=38432, path=neighbor["path"]))
hashes = legacy["input_sha256"] | {neighbor["path"]: neighbor["sha256"]}
z = {89464: legacy["fixed_dr4_z"], 38432: 6.7170105817648125}
rows = []
fig = Figure(figsize=(11, 13), layout="constrained")
axes = fig.subplots(len(products), 2)
for i, p in enumerate(products):
    assert digest(root / p["path"]) == hashes[p["path"]]
    reader = read_dr4_spectrum if p["product"] == "DR4_G395M" else read_spectrum
    spectrum = reader(root / p["path"], wavelength_range=(.6, 5.5))
    wave = spectrum.wavelength.to_value(u.AA)
    unit = u.erg / u.s / u.cm**2 / u.AA
    flux, err = [v.to_value(unit, equivalencies=u.spectral_density(spectrum.wavelength)) / 1e-20
                 for v in [spectrum.flux, spectrum.error]]
    redshift = z[p["source_id"]]
    for width in [2400, 3000]:
        for order in [0, 1]:
            for motion in ["shared", "free"]:
                fit = fit_pair(wave, flux, err, spectrum.good, redshift,
                               window_kms=width, continuum_order=order, kinematics=motion)
                rows.append(dict(product=p["product"], source_id=p["source_id"], window_kms=width,
                                 continuum_order=order, kinematics=motion, **fit))
    # Plot the final fixed variant everywhere, including rejected fits and gaps.
    fit = rows[-1]
    velocity = [(wave / (rest * (1 + redshift)) - 1) * C_KMS for rest in REST]
    windows = [abs(v) <= 3000 for v in velocity]
    model = pair_model([wave[m] for m in windows], redshift, fit["parameters"]) if fit["status"] == "FIT" else None
    for j, (v, m, ax) in enumerate(zip(velocity, windows, axes[i], strict=True)):
        use = m & spectrum.good
        ax.step(v[m], np.where(spectrum.good[m], flux[m], np.nan), where="mid", color="#315676", lw=.9)
        ax.errorbar(v[use], flux[use], yerr=err[use], fmt=".", ms=2, lw=.5, color="#315676", alpha=.55)
        if model is not None:
            ax.plot(v[m], np.where(spectrum.good[m], model[j], np.nan), color="#c66c32", lw=1.2)
        ax.axvline(0, color="gray", lw=.5)
        ax.set(xlim=(-3000, 3000), ylabel="Flux density (10⁻²⁰ cgs)", title=
               f"{p['source_id']} · {p['product']} · {'Hγ' if j == 0 else 'Hβ'} · {fit['status']}")
        if i == len(products) - 1:
            ax.set_xlabel("Velocity relative to catalogue redshift (km/s)")
fig.suptitle("Observed spectra and empirical line models · fixed ±3000 km/s, linear continuum, free widths",
             fontsize=12)
save_figure(fig, root / "figures/balmer_profiles/profile_grid")
write_json(out / "fits.json", rows)
tab = pd.DataFrame(rows).drop(columns=["parameters", "parameter_errors", "covariance", "coverage"])
tab.to_csv(out / "fits.csv", index=False)
eligible = tab[tab.fit_quality_eligible == True]  # noqa: E712
summary = []
for (sid, product), group in tab.groupby(["source_id", "product"]):
    take = eligible[eligible["product"] == product]
    summary.append(dict(source_id=sid, product=product, attempts=len(group),
                        fits=int(group.status.eq("FIT").sum()), usable_ratios=len(take),
                        minimum=take.ratio.min() if len(take) else None,
                        maximum=take.ratio.max() if len(take) else None))
pd.DataFrame(summary).to_csv(out / "product_summary.csv", index=False)
write_json(out / "summary.json", dict(executed_utc=utcnow(), attempts=len(rows), products=summary,
    input_sha256=hashes, script_sha256=text_digest(Path(__file__)),
    module_sha256=text_digest(root / "src/cosmic_anomaly_lab/balmer_profiles.py"),
    model_scope="Observed Gaussian profiles including [O III]4363 nuisance; no LSF deconvolution",
    quality_selection="Both line S/N >= 5, fitted parameters interior, reduced chi2 <= 2",
    independent_trials=False, physical_classification="NOT_ESTABLISHED"))
print(pd.DataFrame(summary).to_string(index=False))

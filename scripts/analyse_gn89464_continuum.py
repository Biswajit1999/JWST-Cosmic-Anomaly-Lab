"""Fixed exploratory continuum/centering grid; every variant is retained."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from matplotlib import colormaps
from matplotlib.figure import Figure

from cosmic_anomaly_lab.balmer_census import signed_core_integral
from cosmic_anomaly_lab.continuum_audit import SIDEBANDS, bilateral_core
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import pixel_edges, read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/gn89464_followup"
out.mkdir(parents=True, exist_ok=True)
target = pd.read_csv(root / "results/jades_balmer_census/followup_queue.csv")
target = target[target.dr4_id.str.endswith("_89464")].iloc[0]
new = json.loads((out / "new_products.json").read_text())["products"]
dja = json.loads((root / "results/multiwavelength/dja_products.json").read_text())
dr4 = json.loads((root / "results/jades_balmer_census/spectral_products.json").read_text())
products = [("DJA_" + p["grating"], p, read_spectrum) for p in new]
products += [("DJA_G395M", next(p for p in dja if p["dr4_id"] == target.dr4_id), read_spectrum)]
products += [("DR4_G395M", next(p for p in dr4 if p["dr4_id"] == target.dr4_id
                              and p["grating_filter"] == "g395m_f290lp"), read_dr4_spectrum)]
cache = DownloadCache(root, budget_bytes=4096 * 1024**2)
records = {r["path"]: r for r in cache.records}
unit = u.erg / u.s / u.cm**2 / u.AA
rest = dict(hg=.434165, hb=.486264)
rows, legacy, identities, hashes, spectra = [], [], [], {}, {}
for label, p, reader in products:
    # Product records are also the reproducible retrieval recipe; single writer.
    record = records[p["path"]]
    path = cache.fetch(record["url"], p["path"].removeprefix("data/raw/"),
                       survey=record["survey"], version=record["version"], kind=record["kind"])
    assert digest(path) == p["sha256"]
    hashes[p["path"]] = p["sha256"]
    with fits.open(path) as h:
        if label.startswith("DJA"):
            header = h["SPEC1D"].header
            ra, dec, pid, sid = header["SRCRA"], header["SRCDEC"], int(header["PROGRAM"]), int(header["SOURCEID"])
        else:
            header = h[0].header
            ra, dec, pid, sid = header["RA"], header["DEC"], int(header["PROPOSID"]), int(header["TARGET"])
    sep = SkyCoord(ra*u.deg, dec*u.deg).separation(SkyCoord(target.ra*u.deg, target.dec*u.deg)).arcsec
    assert sep < .1 and pid == int(target.pid) and sid == 89464
    identities.append(dict(product=label, separation_arcsec=float(sep), pid=pid, source_id=sid,
                           status="ACCEPTED", path=p["path"]))
    s = reader(path, wavelength_range=(.6, 5.5))
    f = s.flux.to_value(unit, equivalencies=u.spectral_density(s.wavelength)) / 1e-20
    e = s.error.to_value(unit, equivalencies=u.spectral_density(s.wavelength)) / 1e-20
    spectra[label] = (s, f, e)
    for line, rw in rest.items():
        center0 = rw * (1 + target.z)
        v0 = (s.wavelength.value / center0 - 1) * 299792.458
        for half in [400, 600]:
            legacy.append(dict(product=label, line=line,
                **signed_core_integral(v0, f, e, s.good, center0*1e4, half)))
        for shift in [-200, -100, 0, 100, 200]:
            center = center0 * (1 + shift / 299792.458)
            v = (s.wavelength.value / center - 1) * 299792.458
            for half in [400, 600]:
                for name, windows in SIDEBANDS.items():
                    rows.append(dict(product=label, line=line, shift_kms=shift, half_width_kms=half,
                        sidebands=name, fixed_z=float(target.z),
                        **bilateral_core(v, f, e, s.good.copy(), center*1e4, half, windows)))
data = pd.DataFrame(rows)
ratios = []
for keys, group in data.groupby(["product", "shift_kms", "half_width_kms", "sidebands"]):
    indexed = group.set_index("line")
    n, d = indexed.loc["hg"], indexed.loc["hb"]
    measured = n.status == d.status == "MEASURED_CORE"
    eligible = bool(measured and n.flux/n.error >= 5 and d.flux/d.error >= 5)
    ratios.append(dict(zip(["product", "shift_kms", "half_width_kms", "sidebands"], keys, strict=True)) | dict(
        eligible=eligible, numerator_snr=float(n.flux/n.error) if measured else None,
        denominator_snr=float(d.flux/d.error) if measured else None,
        ratio=float(n.flux/d.flux) if eligible else None,
        formal_ratio_error=float(np.hypot(n.error/d.flux, n.flux*d.error/d.flux**2)) if eligible else None))
ratios = pd.DataFrame(ratios)
data.to_csv(out / "bilateral_cores.csv", index=False)
pd.DataFrame(legacy).to_csv(out / "legacy_cores.csv", index=False)
ratios.to_csv(out / "bilateral_ratios.csv", index=False)
pd.DataFrame(identities).to_csv(out / "identity_checks.csv", index=False)

fig = Figure(figsize=(12, 10), layout="constrained")
axes = fig.subplots(4, 2)
for row, label in enumerate(["DJA_G395M", "DR4_G395M", "DJA_G235M", "DJA_G395H"]):
    s, f, e = spectra[label]
    for col, (line, rw) in enumerate(rest.items()):
        v = (s.wavelength.value/(rw*(1+target.z))-1)*299792.458
        window = abs(v) <= 3400
        ax = axes[row, col]
        ax.step(v[window], np.where(s.good[window], f[window], np.nan), where="mid", lw=.8, color="#153b57")
        ax.fill_between(v[window], np.where(s.good[window], f[window]-e[window], np.nan),
                        np.where(s.good[window], f[window]+e[window], np.nan), step="mid", color="#153b57", alpha=.12)
        for a, b in SIDEBANDS["legacy"]:
            ax.axvspan(a, b, alpha=.13, color="#d78e24")
        for a, b in SIDEBANDS["near"]:
            ax.axvspan(a, b, alpha=.10, color="#3a9462")
        ax.axvspan(-600, 600, alpha=.07, color="gray")
        ax.axhline(0, lw=.5, color="gray")
        ax.set(xlim=(-3400, 3400), title=f"{label}: {line.upper()}",
               xlabel="Velocity relative to DR4 redshift (km/s)", ylabel="Fλ (10⁻²⁰ cgs / Å)")
        if label == "DJA_G395H" and line == "hb":
            ax.text(.02, .91, "Legacy blue sideband: 0 valid pixels", transform=ax.transAxes, fontsize=8, color="#a22c24")
fig.suptitle("GN-89464 observed spectra: continuum coverage matters\nOrange: legacy continuum; green: near continuum; gaps remain masked; shaded error: formal 1σ", fontsize=12)
save_figure(fig, root / "figures/gn89464_followup/continuum_coverage")

fig = Figure(figsize=(10, 5), layout="constrained")
ax = fig.subplots()
labels = ["DJA_G395M", "DR4_G395M", "DJA_G235M", "DJA_G395H", "DJA_PRISM"]
for i, label in enumerate(labels):
    rr = ratios[(ratios["product"] == label) & ratios.eligible]
    for j, (name, color) in enumerate(zip(SIDEBANDS, ["#176a8a", "#a85525", "#308849"], strict=True)):
        r = rr[rr.sidebands == name]
        if len(r):
            x = i + (j-1)*.20 + (r.shift_kms.to_numpy()/200)*.035
            ax.scatter(x, r.ratio, s=20, color=color, label=name if i == 0 else None, alpha=.7)
ax.axhline(.468, color="gray", ls="--", label="0.468 reference (not a universal boundary)")
counts = ratios.groupby("product").eligible.sum()
ax.set(xticks=range(5), xticklabels=[f"{s}\n{counts[s]}/30 eligible" for s in labels],
       ylabel="Hγ/Hβ signed core ratio", xlim=(-.5, 4.5),
       title="Fixed 30-variant grid per spectrum; dual S/N ≥5 required\nPoint estimates only; correlated variants, not confidence intervals", ylim=(.25, 1.2))
ax.legend(fontsize=8)
save_figure(fig, root / "figures/gn89464_followup/ratio_sensitivity")

# Actual released 2D SCI data, with signed intensities and explicit invalid pixels.
p = next(p for label, p, _ in products if label == "DJA_G395H")
with fits.open(root / p["path"]) as hdus:
    sci = hdus["SCI"]
    image = np.asarray(sci.data, float) * u.Unit(sci.header["BUNIT"]).to(u.uJy)
    weight = np.asarray(hdus["WHT"].data, float)
    wave = np.asarray(hdus["SPEC1D"].data["wave"], float)
    trace = float(sci.header["YTRACE"])
    assert image.shape == weight.shape and image.shape[1] == len(wave)
    assert np.allclose(wave, spectra["DJA_G395H"][0].wavelength.value)
    valid = np.isfinite(image) & np.isfinite(weight) & (weight > 0)
    valid &= spectra["DJA_G395H"][0].good[None, :]
fig = Figure(figsize=(12, 5), layout="constrained")
axes = fig.subplots(1, 2)
for ax, (line, rw) in zip(axes, rest.items(), strict=True):
    v = (wave/(rw*(1+target.z))-1)*299792.458
    window = abs(v) <= 3400
    data2d = np.where(valid[:, window], image[:, window], np.nan)
    limit = max(float(np.nanpercentile(abs(data2d), 98)), 1e-8)
    cmap = colormaps["RdBu_r"].copy()
    cmap.set_bad("#b8b8b8")
    mesh = ax.pcolormesh(pixel_edges(v[window]), np.arange(image.shape[0]+1)-.5-trace,
                         np.ma.masked_invalid(data2d), cmap=cmap, vmin=-limit, vmax=limit,
                         shading="flat", rasterized=True)
    ax.axhline(0, color="black", lw=.6, ls=":")
    for a, b in SIDEBANDS["legacy"]:
        ax.axvline(a, color="#946400", lw=.8, ls="--")
        ax.axvline(b, color="#946400", lw=.8, ls="--")
    ax.set(xlim=(-3400, 3400), ylim=(-8, 8), title=f"G395H: {line.upper()}",
           xlabel="Velocity relative to DR4 redshift (km/s)", ylabel="Spatial pixel relative to trace")
    fig.colorbar(mesh, ax=ax, label="Released SCI (µJy); panel-specific scale", shrink=.7)
fig.suptitle("GN-89464: actual JWST two-dimensional spectral images\nGray = invalid weight or masked 1D column; dashed boundaries = legacy continuum windows", fontsize=12)
save_figure(fig, root / "figures/gn89464_followup/observed_g395h")
summary = dict(executed_utc=utcnow(), target=target.dr4_id, fixed_dr4_z=float(target.z),
               catalogue_product_z_is_not_fitting_z=True, products=len(products),
               core_attempts=len(data), usable_cores=int((data.status == "MEASURED_CORE").sum()),
               ratio_attempts=len(ratios), eligible_ratios=int(ratios.eligible.sum()),
               input_sha256=hashes, script_sha256=text_digest(Path(__file__)),
               module_sha256=text_digest(root / "src/cosmic_anomaly_lab/continuum_audit.py"),
               selection="Exploratory grid defined after baseline coverage inspection, not blind confirmation",
               physical_classification="NOT_ESTABLISHED", independent_trials=False)
write_json(out / "summary.json", summary)
print(json.dumps(summary, indent=2))
print(ratios.groupby("product").agg(eligible=("eligible", "sum"), minimum=("ratio", "min"), maximum=("ratio", "max")))

"""Verify and display the prior GO 8060 spectroscopy of the fading lead."""

import json

import numpy as np
from astropy.io import fits
from audit_jwst_module_b_spectra import STORAGE_MIB
from download_jwst_cross_detector import OUT, ROOT, download
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import verified_path
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure


def main():
    path = OUT/"module_b_vetting.json"
    (lead,) = [r for r in json.loads(path.read_text())["candidates"] if r["sign"] < 0]
    entries = [r for r in lead["public_matches"]["nirspec"]["nearest"]
               if r["srcid"] == 2510132 and r["separation_arcsec"] < .10 and r["root"].startswith("gds-dts-egami")]
    if len(entries) != 2 or len({r["root"] for r in entries}) != 2:
        raise ValueError("Expected two matching GO 8060 extractions of the registered fading lead")
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    archive = "https://s3.amazonaws.com/msaexp-nirspec/extractions"
    records, spectra = [], []
    for entry in sorted(entries, key=lambda row: row["root"]):
        product = download(cache, f'{archive}/{entry["root"]}/{entry["file"]}', entry["file"],
                           "observed_prior_nirspec_transient")
        with fits.open(verified_path(product)) as h:
            data = h["SPEC1D"]
            hdr = data.header
            if hdr["SOURCEID"] != entry["srcid"] or (hdr["GRATING"], hdr["FILTER"]) != ("PRISM", "CLEAR"):
                raise ValueError("Unexpected spectrum metadata")
            if data.columns["wave"].unit != "um" or data.columns["flux"].unit != "uJy":
                raise ValueError("Unexpected units on the prior spectrum")
            wave, flux, err = (np.asarray(data.data[x], dtype=float) for x in ("wave", "flux", "err"))
            good = np.isfinite(wave+flux+err) & (err > 0)
            wave, flux, err = wave[good], flux[good], err[good]
            records.append(dict(entry=entry, product=product, date_utc=hdr["DATE-OBS"],
                                effective_exposure_seconds=hdr["EXPTIME"], msa=hdr["MSAMETFL"],
                                finite_samples=int(good.sum()), wavelength_um=[float(wave.min()), float(wave.max())]))
            spectra.append((wave, flux, err))

    figure = Figure(figsize=(9, 5.4), layout="constrained")
    axes = figure.subplots(2, 1, sharex=True, sharey=True)
    wavelength_limits = (0.70, 5.15)
    lower, upper = np.inf, -np.inf
    for ax, row, (wave, flux, err), colour in zip(axes, records, spectra, ("#265b80", "#92503b"), strict=True):
        visible = (wave >= wavelength_limits[0]) & (wave <= wavelength_limits[1])
        lower = min(lower, (flux[visible] - err[visible]).min())
        upper = max(upper, (flux[visible] + err[visible]).max())
        ax.axhline(0, color="#8a9295", lw=.6)
        ax.fill_between(wave, flux-err, flux+err, alpha=.16, color=colour, step="mid")
        ax.step(wave, flux, lw=.85, color=colour, where="mid")
        ax.set(ylabel="Flux density (µJy)", title=f'{row["date_utc"]} · DJA inspection grade {row["entry"]["grade"]}')
        ax.grid(alpha=.18)
    padding = 0.04 * (upper - lower)
    axes[-1].set(xlim=wavelength_limits, ylim=(lower - padding, upper + padding),
                 xlabel="Observed wavelength (µm)")
    figure.suptitle("Prior JWST/NIRSpec PRISM observations of the fading lead\n"
                    "GO 8060 source 2510132 · DAWN archival flux and per-sample uncertainty", fontsize=12)
    save_figure(figure, FIGURES/"module_b_002_prior_prism")
    write_json(OUT/"module_b_prior_prism.json", {
        "executed_utc": utcnow(), "candidate": lead["candidate"], "spectra": records,
        "scope": "Verification of prior spectroscopy; classification and redshift annotations are inherited from the archive",
        "vetting_sha256": text_digest(path), "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

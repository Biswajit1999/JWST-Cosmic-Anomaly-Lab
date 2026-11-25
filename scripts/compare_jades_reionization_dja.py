"""Measure the literature-unmatched Lyα responses in the DAWN reductions."""

import json

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from measure_jades_reionization import C_KMS, LYA, WIDTHS, best_model, line_models, spectral_arrays

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

INVENTORY = ROOT/"data/catalogues/observations.csv.gz"


def measure(data, zero, offset):
    models = line_models(data,zero,WIDTHS[0])
    out = {}
    for name, result in (("fixed", next((p for p in models if p["offset_kms"] == offset), None)),
                         ("best", best_model(models))):
        out[name] = (dict(offset_kms=result["offset_kms"], flux=result["flux"], error=result["error"],
                          snr=result["flux"]/result["error"], coverage=result["coverage"])) if result else None
    return out


def panel(ax,data,zero,label,color):
    wave,_,flux,err,good = data
    velocity = (wave/zero-1)*C_KMS
    show = (velocity > -1700)&(velocity < 4200)
    y = np.where(good,flux,np.nan)
    ax.step(velocity[show],y[show],where="mid",label=label,color=color,linewidth=1.)
    ax.fill_between(velocity[show],(y-err)[show],(y+err)[show],step="mid",color=color,alpha=.13,linewidth=0)


def main():
    crosschecks = OUT/"crosschecks.json"
    published = json.loads(crosschecks.read_text())
    rows = published["catalogue_unmatched"]
    sample = {r["dr4_id"]:r for r in json.loads(LEDGER.read_text())["records"]}
    obs = pd.read_csv(OUT/"measurements.csv").set_index("dr4_id")
    inventory = pd.read_csv(INVENTORY)
    inventory = inventory[(inventory.survey == "JADES") & (inventory.grating == "G140M")
                          & (inventory["filter"] == "F070LP")].copy()
    coords = SkyCoord(inventory.ra,inventory.dec,unit="deg")
    cache = DownloadCache(ROOT,28_000*1024**2)
    fig = Figure(figsize=(9.4,2.5*len(rows)),layout="constrained")
    axes = fig.subplots(len(rows),1,squeeze=False).ravel()
    comparisons = []
    for y,ax in zip(rows,axes,strict=True):
        sid = y["dr4_id"]
        r,c = sample[sid],obs.loc[sid]
        zero = LYA*(1+c.z)
        t = SkyCoord(r["ra"],r["dec"],unit="deg")
        separation = coords.separation(t).arcsec
        fits = inventory[separation <= .15].sort_values("root")
        if not len(fits):
            raise ValueError(f"No close DJA extraction for {sid}")
        jades = spectral_arrays(read_dr4_spectrum(ROOT/r["spectrum"]["path"],wavelength_range=(.65,5.4)))
        panel(ax,jades,zero,"JADES DR4", "#225c70")
        products = []
        for i,record in enumerate(fits.itertuples()):
            source = cache.fetch(record.url,f"spectra/{record.file}",survey="DAWN JWST Archive",version="v4",
                                 kind="observed_spectrum")
            spec = read_spectrum(source,wavelength_range=(.65,5.4))
            if (spec.metadata["GRATING"].upper(),spec.metadata["FILTER"].upper()) != ("G140M","F070LP"):
                raise ValueError(f"Wrong disperser in {record.file}")
            off = t.separation(SkyCoord(spec.metadata["SRCRA"], spec.metadata["SRCDEC"],unit="deg")).arcsec
            if off > .15:
                raise ValueError(f"DJA header position differs by {off} arcsec for {record.file}")
            if abs(record.z_best-c.z)/(1+c.z) > .002:
                raise ValueError(f"DJA redshift disagrees for {record.file}")
            data = spectral_arrays(spec)
            fit = measure(data,zero,c.best_velocity_kms)
            products.append(dict(product=receipt(source,record.url), dataset=record.dataset, root=record.root,
                                 header_separation_arcsec=float(off), programme=int(record.programme),
                                 archive_redshift=float(record.z_best), exposure_s=float(record.exptime), **fit))
            panel(ax,data,zero,record.root,("#b14e39","#785ba3")[i%2])
            print(sid,record.root,record.dataset,fit,flush=True)
        ax.axhline(0,color=".5",linewidth=.6)
        ax.axvline(0,color=".35",linestyle=":",linewidth=.7)
        ax.axvline(c.best_velocity_kms,color=".4",linestyle="--",linewidth=.7)
        ax.legend(frameon=False,fontsize=8,loc="upper right")
        ax.set(xlim=(-1700,4200),ylabel="Fλ (10⁻²⁰ cgs/Å)",title=f"{sid} · z={c.z:.3f}")
        comparisons.append(dict(dr4_id=sid, z=float(c.z), jades_measurement=measure(jades,zero,c.best_velocity_kms),
                                dja=products))
    axes[-1].set_xlabel("Velocity relative to systemic Lyα (km/s)")
    fig.suptitle("Observed spectra for four Lyα follow-up targets · JADES and DAWN reductions",fontsize=12)
    save_figure(fig,ROOT/"figures/jades_reionization_dja")
    write_json(OUT/"dja_comparison.json",dict(executed_utc=utcnow(), cohort_sha256=text_digest(crosschecks),
        script_sha256=text_digest(__file__), inventory_sha256=digest(INVENTORY), candidates=comparisons,
        scope="Distinct reduction products. Archive dataset names identify one reference exposure per extraction."
              " Inputs can overlap. Fixed offsets were set by the JADES spectra."))


if __name__=="__main__":
    main()

"""Audit JADES Lyman-alpha line responses across released integrations and slitlet groups."""

import json
import re

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from compare_jades_reionization_dja import measure as report_model
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import DENSITY_UNIT, receipt
from matplotlib.figure import Figure
from measure_jades_reionization import C_KMS, LYA, WIDTHS, profile, spectral_arrays
from scipy.stats import chi2

from cosmic_anomaly_lab.exposure_audit import combine_extractions
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import pixel_edges
from cosmic_anomaly_lab.visual_evidence import save_figure

DOC = "https://jades.herts.ac.uk/DR4/Readme_DR4_spectra.md"
STEM = "jades_reionization_exposure_audit"
CONTINUUM_WINDOWS = ((1800, 4800), (2200, 6800), (3600, 8500))
SIGMAS = (110, 150, 250)


def read_intermediate(path):
    with fits.open(path) as hdus:
        data = hdus["INTERMEDIATE"]
        unit = u.Unit(data.columns["WAVELENGTH"].unit)
        wave = (data.data["WAVELENGTH"]*unit).to_value(u.AA)
        names, flux, error = [], [], []
        for name in data.columns.names:
            match = re.fullmatch(r"FLUX_(\d+)_(\d+)_(\d+)", name)
            if match:
                names.append(tuple(map(int, match.groups())))
                for col, values in ((name, flux), ("ERR_"+name[5:], error)):
                    values.append(np.asarray(data.data[col], dtype=float)*u.Unit(data.columns[col].unit)
                                  .to(DENSITY_UNIT)/1e-20)
        if not names or len(set(names)) != len(names) or np.any(np.diff(wave) <= 0):
            raise ValueError(f"Unexpected intermediate data in {path.name}")
        dates = [hdus[0].header.get(key) for key in ("DATE-BEG", "DATE-END")]
        regions = sorted({int(i) for i in hdus["REGION"].data["COMPONENT"]}) if len(hdus["REGION"].data) else []
    return wave, names, np.array(flux), np.array(error), dates, regions


def arrays(wave, flux, error):
    good = np.isfinite(flux) & np.isfinite(error) & (error>0) & (wave>=6500) & (wave<=54000)
    return wave, pixel_edges(wave), flux, error, good


def diagnostic_fit(data, zero, offset, width, window, affine):
    wave, edges, flux, error, good = data
    v = (wave/zero-1)*C_KMS
    red = good & (edges[:-1] >= zero)
    core = red & (abs(v-offset) <= 4*width)
    side = red & (v>=window[0]) & (v<=window[1])
    shape = profile(edges,zero,offset,width)
    coverage = float(shape[core] @ np.diff(edges)[core])
    if coverage < .78 or core.sum() < 3 or side.sum() < 8:
        return None
    v_pivot = v/1000
    design = np.column_stack((np.ones(len(v)), v_pivot)) if affine else np.ones((len(v),1))
    w_design = design[side]/error[side,None]
    cov = np.linalg.inv(w_design.T @ w_design)
    beta = cov @ (w_design.T @ (flux[side]/error[side]))
    weights = shape[core]/error[core]**2
    weights /= weights @ shape[core]
    continuum_weights = weights @ design[core]
    value = float(weights @ flux[core]-continuum_weights @ beta)
    sigma = float(np.sqrt(np.sum((weights*error[core])**2)+continuum_weights @ cov @ continuum_weights))
    return dict(flux=value,error=sigma,snr=value/sigma,coverage=coverage,side_pixels=int(side.sum()))


def audit_target(target, observed, role):
    ref = target["spectrum"]
    path = ROOT/ref["path"]
    if digest(path) != ref["sha256"]:
        raise ValueError(f"Changed spectrum fingerprint for {path.name}")
    wave, names, f, e, dates, regions = read_intermediate(path)
    zero, offset = LYA*(1+observed.z), int(observed.best_velocity_kms)
    release = spectral_arrays(read_dr4_spectrum(path,wavelength_range=(.65,5.4)))
    combos = {}
    rows = []
    for i,(pointing,integration,region) in enumerate(names):
        fixed = report_model(arrays(wave,f[i],e[i]),zero,offset)["fixed"]
        rows.append(dict(pointing=pointing,integration=integration,region=region,fixed=fixed))
    groups = sorted({item[2] for item in names})
    subgroups, omitted = [], []
    name_region = np.array([r for _,_,r in names])
    for region in groups:
        for keep, kind, records in ((name_region==region,"group",subgroups),
                                     (name_region!=region,"omit_group",omitted)):
            if np.any(keep):
                flow, err, _ = combine_extractions(f[keep],e[keep])
                data = arrays(wave,flow,err)
                combos[kind,region] = data
                records.append(dict(region=region,integrations=int(keep.sum()),
                                    fixed=report_model(data,zero,offset)["fixed"]))
    flux, err, count = combine_extractions(f,e)
    combined = arrays(wave,flux,err)
    fixed_regions = [r["fixed"] for r in subgroups if r["fixed"] is not None]
    dispersion = None
    if len(fixed_regions) >= 2:
        y,s = (np.array([r[key] for r in fixed_regions]) for key in ("flux","error"))
        mean = float(np.average(y,weights=s**-2))
        stat = float(np.sum(((y-mean)/s)**2))
        dof = len(y)-1
        dispersion = dict(ivar_mean=mean,chi_squared=stat,dof=dof,
                          diagonal_noise_tail_probability=float(chi2.sf(stat,dof)))
    variants = []
    for sample,d in (("released",release),("intermediate_ivar",combined)):
        for lo,hi in CONTINUUM_WINDOWS:
            for affine in (False,True):
                for width in SIGMAS:
                    variants.append(dict(sample=sample,continuum_kms=[lo,hi],affine=affine,
                                         sigma_kms=width,fit=diagnostic_fit(d,zero,offset,width,(lo,hi),affine)))
    v = (wave/zero-1)*C_KMS
    core = (v >= 0) & (abs(v-offset) <= 4*WIDTHS[0])
    result = dict(dr4_id=target["dr4_id"],role=role,z=float(observed.z),offset_kms=offset,
                  dates_utc=dates,source=ref,slitlet_components_in_region_table=regions,
                  integrations=rows,subgroups=subgroups,leave_one_region_out=omitted,
                  released=report_model(release,zero,offset),intermediate_ivar=report_model(combined,zero,offset),
                  group_scatter_diagnostic=dispersion,continuum_and_width_sensitivity=variants,
                  contributors_per_core_bin=[int(x) for x in count[core]])
    return result, (zero,release,combined,combos)


def draw_spectrum(ax,data,zero,label,color,alpha):
    wave, _, flux, error, good = data
    v = (wave/zero-1)*C_KMS
    select = abs(v) <= 4000
    y = np.where(good,flux,np.nan)
    ax.step(v[select],y[select],where="mid",lw=.9,color=color,alpha=alpha,label=label)
    ax.fill_between(v[select],(y-error)[select],(y+error)[select],step="mid",color=color,alpha=.11,linewidth=0)


def exposure_figure(output, candidates):
    panels = [(r,d) for r,d in candidates if r["dr4_id"]=="goods-n-mediumjwst_58441"]
    (result,(zero,release,combined,combos)), = panels
    fig = Figure(figsize=(10,7.4),layout="constrained")
    ax,spectrum = fig.subplots(2,1,gridspec_kw={"height_ratios":(1,1.2)})
    y = 0
    for group,color,label in (("integrations","#6990a4","Individual integrations"),
                              ("subgroups","#b65844","Slitlet region combinations"),
                              ("leave_one_region_out","#7464a2","Omit one region")):
        start = y
        for r in result[group]:
            test = r["fixed"]
            key = (f"r{r['region']} / {r['pointing']}_{r['integration']}" if group=="integrations" else
                   f"{'omit ' if group=='leave_one_region_out' else ''}r{r['region']}")
            if test:
                ax.errorbar(y,test["flux"],yerr=test["error"],fmt="o",color=color,ms=4,capsize=2,
                            label=label if y==start else None)
            ax.text(y,-.055,key,rotation=35,ha="right",va="top",fontsize=7,transform=ax.get_xaxis_transform())
            y += 1
        y += .8
    ax.axhline(0,color=".4",lw=.7)
    ax.set(ylabel="Signed Lyα flux (10⁻²⁰ erg s⁻¹ cm⁻²)",xticks=[],xlim=(-.6,y-.3))
    ax.legend(frameon=False,loc="upper left",fontsize=8,ncols=3)
    for (kind,i),data in combos.items():
        if kind=="group":
            draw_spectrum(spectrum,data,zero,f"Region {i}",("#a54533","#7c5ca4","#527f8b")[(i-1)%3],.72)
    draw_spectrum(spectrum,release,zero,"Released JADES DR4","#20242a",1.)
    draw_spectrum(spectrum,combined,zero,"Diagnostic full combination","#308754",.85)
    spectrum.axhline(0,color=".4",lw=.6)
    spectrum.axvline(result["offset_kms"],color=".35",lw=.7,ls=":")
    spectrum.set(xlim=(-1500,4000),xlabel="Velocity relative to systemic Lyα (km/s)",
                 ylabel="Fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)")
    spectrum.legend(frameon=False,fontsize=8,ncols=2,loc="upper right")
    fig.suptitle("GOODS-N 58441 · released JWST G140M spectra and integration jackknife",fontsize=13)
    save_figure(fig,output)


def main():
    tracks = OUT/"crosschecks.json"
    verified = json.loads(tracks.read_text())
    choose = [(row["dr4_id"],"line_followup") for row in verified["catalogue_unmatched"]]
    choose.extend((row["dr4_id"],"published_control") for row in verified["examples"])
    rows = {r["dr4_id"]:r for r in json.loads(LEDGER.read_text())["records"]}
    obs = pd.read_csv(OUT/"measurements.csv").set_index("dr4_id")
    cache = DownloadCache(ROOT,28_000*1024**2)
    doc = cache.fetch(DOC,"metadata/jades_dr4_spectra.md",survey="JADES",version="DR4",kind="spectrum_documentation")
    audited = [audit_target(rows[sid],obs.loc[sid],role) for sid,role in choose]
    exposure_figure(ROOT/f"figures/{STEM}",audited)
    write_json(OUT/"exposure_audit.json",dict(
        executed_utc=utcnow(),script_sha256=text_digest(__file__),crosschecks_sha256=text_digest(tracks),
        measurement_sha256=text_digest(OUT/"measurements.csv"),documentation=receipt(doc,DOC),
        candidates=[item for item,_ in audited],
        scope="Released intermediate spectra, grouped by slitlet region. Formal errors assume diagonal noise."
              " Background nods can share observations. Diagnostic means use no clipping or survey path-loss recalibration.",
    ))
    for result,_ in audited:
        print(result["dr4_id"],result["released"]["fixed"]["snr"],
              [r["fixed"]["snr"] if r["fixed"] else None for r in result["subgroups"]],
              [r["fixed"]["snr"] if r["fixed"] else None for r in result["leave_one_region_out"]],flush=True)


if __name__=="__main__":
    main()

"""Audit GOODS-N 58441 calibration and exact exposure inventories across dispersers."""

import json
import re
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from audit_jades_reionization_exposures import arrays
from check_jades_reionization_prism import NOMINAL, PRIMARY, local_fit
from compare_jades_reionization_dja import measure as grating_fit
from download_jades_reionization_sample import OUT, ROOT
from inspect_jades_cosmic_evolution import DENSITY_UNIT, receipt
from matplotlib.figure import Figure
from measure_jades_reionization import C_KMS, LYA, spectral_arrays
from scipy.special import ndtr

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import pixel_edges, read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

UV_WINDOWS = ((1160,1200),(1250,1320),(1350,1500),(1550,1700))
STEM = "jades_reionization_58441_flux_scale"


def oxygen_doublet_check(product, z, offset):
    """Test the low-redshift [O III] identification with a simultaneous doublet fit."""
    spec = read_dr4_spectrum(ROOT/product["path"],wavelength_range=(.65,5.4))
    wave,edges,flux,error,good = spectral_arrays(spec)
    zero = LYA*(1+z)
    main = zero*(1+offset/C_KMS)
    weak = main*4960.295/5008.239
    use = good & (wave>weak-55) & (wave<main+55)
    checks = []
    for width in (110,150,220,300):
        line = [np.diff(ndtr((edges-c)/(c*width/C_KMS)))/np.diff(edges) for c in (main,weak)]
        for continuum_break in (False,True):
            terms = [np.ones(len(wave)),(wave-main)/100]
            if continuum_break:
                terms.append(ndtr((wave-zero)/(zero*width/C_KMS)))
            design = np.column_stack([*terms,*line])[use]/error[use,None]
            cov = np.linalg.inv(design.T @ design)
            fit = cov @ (design.T @ (flux[use]/error[use]))
            contrast = np.r_[np.zeros(len(terms)),1/2.98,-1.]
            diff = float(contrast@fit)
            sigma = float(np.sqrt(contrast@cov@contrast))
            checks.append(dict(sigma_kms=width,continuum_break=continuum_break,bins=int(use.sum()),
                               blended_main_flux=float(fit[-2]),weak_flux=float(fit[-1]),
                               weak_error=float(np.sqrt(cov[-1,-1])),
                               predicted_minus_measured_weak_flux=diff,contrast_error=sigma,formal_snr=diff/sigma))
    return dict(main_wavelength_a=main,weak_wavelength_a=weak,assumed_ratio=2.98,
                alternative_redshift=main/5008.239-1,checks=checks)


def uv_band(spec, window, z):
    edges = pixel_edges(spec.wavelength.to_value(u.AA))/(1+z)
    lo,hi = window
    span = np.maximum(0,np.minimum(edges[1:],hi)-np.maximum(edges[:-1],lo))
    use = spec.good & (span>0)
    coverage = float(span[use].sum()/(hi-lo))
    if coverage < .75 or use.sum()<3:
        return dict(coverage=coverage,mean_uJy=None,error_uJy=None)
    weights = span[use]/span[use].sum()
    mean = float(weights @ spec.flux.to_value(u.uJy)[use])
    err = float(np.linalg.norm(weights*spec.error.to_value(u.uJy)[use]))
    return dict(coverage=coverage,mean_uJy=mean,error_uJy=err)


def compare(jades, dja, items, z, offset, r, mode):
    record = dict(mode=mode)
    spectra = {}
    for survey, product, reader in (("jades",jades,read_dr4_spectrum),("dja",dja,read_spectrum)):
        path = ROOT/product["path"]
        if digest(path) != product["sha256"]:
            raise ValueError(f"Fingerprint changed for {path.name}")
        spectra[survey] = reader(path,wavelength_range=(.65,5.4))
        record[survey] = dict(product=product,uv=[uv_band(spectra[survey],w,z) for w in UV_WINDOWS])
    with fits.open(ROOT/dja["path"]) as h:
        t = h["SPEC1D"].data
        slit = h["SLITS"].data
        exposure_keys = []
        for f in slit["filename"]:
            match = re.match(r"(jw\d+_\d+_\d+)_nrs[12]_",str(f))
            if not match:
                raise ValueError(f"Unexpected exposure identifier {f}")
            exposure_keys.append(match[1])
        cut = (abs(t["wave"]*1e4/(1+z)-1217)<12) & (t["npix"]>0)
        if not cut.any() or len(t) != len(spectra["dja"].wavelength):
            raise ValueError("No usable DAWN pathloss samples near Lyα")
        corr = np.array(t["path_corr"],float)
        corr_lya = float(np.median(corr[cut]))
        record["dja"].update(disperser_input_files=sorted(set(map(str,slit["filename"]))),
                             detector_independent_exposure_keys=sorted(set(exposure_keys)),
                             lya_path_corr_median=corr_lya,
                             source_xpos_range=[float(np.min(slit["source_xpos"])),float(np.max(slit["source_xpos"]))])
    record["released_line_fluxes"] = items
    w,e,f,scr,good = spectral_arrays(spectra["dja"])
    available = np.isfinite(corr)&(corr>0)
    divisor = np.where(available,corr,np.nan)
    counterfactual = (w,e,f/divisor,scr/divisor,good&available)
    if mode=="PRISM":
        res,_ = local_fit(counterfactual,z,offset,r,*NOMINAL)
    else:
        res = grating_fit(counterfactual,LYA*(1+z),offset)["fixed"]
    record["dja"]["formal_division_by_pathcorr"] = res
    return record


def spatial_check(pr, sample):
    crosschecks = json.loads((OUT/"crosschecks.json").read_text())
    g = next(x for x in crosschecks["catalogue_unmatched"] if x["dr4_id"] == PRIMARY)["spectral_image"]
    listing = ROOT/pr["listing"]["path"]
    if digest(listing) != pr["listing"]["sha256"]:
        raise ValueError("PRISM archive listing changed")
    name = Path(pr["spectrum"]["path"]).name.replace("_x1d.fits","_s2d.fits")
    if name not in re.findall(r'href="([^"]+)"',listing.read_text()):
        raise ValueError("Expected PRISM s2d in the release directory")
    url = urljoin(pr["listing"]["url"],name)
    cache = DownloadCache(ROOT,28_000*1024**2)
    path = cache.fetch(url,f"independent/{name}",survey="JADES",version="DR4",kind="observed_spectral_image")
    products = [("PRISM",receipt(path,url)),("G140M",g)]
    fig = Figure(figsize=(9.8,6.2),layout="constrained")
    axes = fig.subplots(2,2,gridspec_kw=dict(width_ratios=(2.8,1)))
    output = []
    z,offset = pr["z"],pr["fixed_velocity_kms"]
    zero = LYA*(1+z)
    position = SkyCoord(sample["ra"],sample["dec"],unit="deg")
    for (mode,record),(img,profile) in zip(products,axes,strict=True):
        path = ROOT/record["path"]
        if digest(path) != record["sha256"]:
            raise ValueError("2D spectrum fingerprint differs")
        with fits.open(path) as f:
            meta = f[0].header
            sep = position.separation(SkyCoord(meta["RA"],meta["DEC"],unit="deg")).arcsec
            if sep>.2 or int(meta["TARGET"])!=sample["header_target_number"] or int(meta["PROPOSID"])!=sample["header_pid"]:
                raise ValueError("2D spectrum identity differs")
            w = f["WAVELENGTH"]
            wave = (w.data*u.Unit(w.header["BUNIT"])).to_value(u.AA)
            plane = f["FLUX"]
            flux = plane.data*u.Unit(plane.header["BUNIT"]).to(DENSITY_UNIT)/1e-20
            error = f["FLUX_ERR"].data*u.Unit(f["FLUX_ERR"].header["BUNIT"]).to(DENSITY_UNIT)/1e-20
            y = plane.header["CRVAL2"]+(np.arange(len(flux))+1-plane.header["CRPIX2"])*plane.header["CDELT2"]
        finite = np.isfinite(flux)&np.isfinite(error)&(error>0)
        signed = np.divide(flux,error,out=np.full_like(flux,np.nan),where=finite)
        v = (wave/zero-1)*C_KMS
        edge = 13000 if mode=="PRISM" else 2400
        use = abs(v-offset) <= edge
        img.pcolormesh(v[use],y,signed[:,use],shading="nearest",cmap="RdBu_r",vmin=-4,vmax=4)
        img.axhline(0,color=".35",linestyle=":",linewidth=.7)
        img.axvline(offset,color=".35",linestyle=":",linewidth=.7)
        img.set(xlabel="Velocity relative to systemic Lyα (km/s)",ylabel="Slit position (shutters)",title=f"Observed {mode} 2D")
        rows = []
        for i,row_y in enumerate(y):
            data = arrays(wave,flux[i],error[i])
            if mode=="PRISM":
                result,_ = local_fit(data,z,offset,pr["nominal_table_r"],*NOMINAL)
            else:
                result = grating_fit(data,zero,offset)["fixed"]
            rows.append(dict(row=i,slit_y=float(row_y),line=result))
        measured = [r for r in rows if r["line"] is not None]
        if measured:
            profile.errorbar([r["line"]["flux"] for r in measured],[r["slit_y"] for r in measured],
                             xerr=[r["line"]["error"] for r in measured],fmt="o",ms=3,capsize=1,color="#277887")
        profile.axhline(0,color=".35",ls=":",lw=.7)
        profile.axvline(0,color=".4",lw=.6)
        profile.set(xlabel="Signed Lyα core (10⁻²⁰ cgs)",title="Line fit by spatial row",ylim=img.get_ylim())
        output.append(dict(mode=mode,product=record,header_offset_arcsec=float(sep),rows=rows))
    fig.suptitle("GOODS-N 58441 · line location in two observed JWST dispersers\n"
                 "Image colour is signed flux / released pixel error (±4)",fontsize=12)
    save_figure(fig,ROOT/"figures/jades_reionization_58441_trace")
    return output


def main():
    prism_path, grating_path = (OUT/"prism_crosscheck.json",OUT/"dja_comparison.json")
    pr = next(r for r in json.loads(prism_path.read_text())["candidates"] if r["dr4_id"] == PRIMARY)
    gr = next(r for r in json.loads(grating_path.read_text())["candidates"] if r["dr4_id"] == PRIMARY)
    sample = next(r for r in json.loads((OUT/"sample.json").read_text())["records"] if r["dr4_id"] == PRIMARY)
    p = pr["source_diagnostics"]["dja"]
    g = gr["dja"][0]
    cases = [(pr["spectrum"],p["product"],[pr["nominal"],p["measurement"]],"PRISM"),
             (sample["spectrum"],g["product"],[gr["jades_measurement"]["fixed"],g["fixed"]],"G140M")]
    comparison = [compare(j,d,a,pr["z"],pr["fixed_velocity_kms"],pr["nominal_table_r"],mode)
                  for j,d,a,mode in cases]
    ids = [set(r["dja"]["detector_independent_exposure_keys"]) for r in comparison]
    common = sorted(ids[0]&ids[1])
    if common:
        raise ValueError("PRISM and G140M exposure lists overlap")
    spatial = spatial_check(pr,sample)
    doublet = oxygen_doublet_check(sample["spectrum"],pr["z"],pr["fixed_velocity_kms"])
    fig = Figure(figsize=(9.4,4.1),layout="constrained")
    left,right = fig.subplots(1,2,gridspec_kw=dict(width_ratios=(1.2,1)))
    limits = np.asarray(UV_WINDOWS)
    centers = limits.mean(axis=1)
    for row,c,mark in zip(comparison,("#287887","#914c99"),("o","s"),strict=True):
        for survey,shift in (("jades",-8),("dja",8)):
            v = np.array([np.nan if a["mean_uJy"] is None else a["mean_uJy"] for a in row[survey]["uv"]])
            e = np.array([np.nan if a["error_uJy"] is None else a["error_uJy"] for a in row[survey]["uv"]])
            left.errorbar(centers+shift,v,yerr=e,fmt=mark,color=c,mfc=c if survey=="jades" else "white",capsize=2,
                          label=f"{row['mode']} · {survey.upper()}")
    left.set(xlabel="Rest-frame wavelength (Å)",ylabel="Band-averaged flux density (μJy)",
             title="Observed continuum windows")
    left.legend(frameon=False,fontsize=8,ncols=2)
    left.axhline(0,color=".45",lw=.6)
    for i,row in enumerate(comparison):
        variants = [*row["released_line_fluxes"],row["dja"]["formal_division_by_pathcorr"]]
        for j,(r,c) in enumerate(zip(variants,("#267484","#b14b38","#cb765e"),strict=True)):
            if r:
                label = ["JADES release","DAWN release","DAWN ÷ path_corr"][j]
                right.errorbar(i+(j-1)*.16,r["flux"]/100,yerr=r["error"]/100,fmt="o",color=c,
                               mfc="white" if j==2 else c,capsize=3,label=label if i==0 else None)
    right.axhline(0,color=".45",lw=.6)
    right.set(xticks=[0,1],xticklabels=["PRISM","G140M"],ylabel="Lyα core flux (10⁻¹⁸ erg s⁻¹ cm⁻²)",
              title="Dependence on extraction and scaling")
    right.legend(frameon=False,fontsize=8)
    fig.suptitle("GOODS-N 58441 · flux-calibration diagnostics",fontsize=13)
    save_figure(fig,ROOT/f"figures/{STEM}")
    write_json(OUT/"58441_flux_scale.json",dict(executed_utc=utcnow(),script_sha256=text_digest(__file__),
        prism_sha256=text_digest(prism_path),grating_sha256=text_digest(grating_path),
        rest_continuum_windows_a=UV_WINDOWS,common_exposure_keys=common,comparisons=comparison,spatial_check=spatial,
        oxygen_doublet_alternative=doublet,
        scope="Exact DAWN slit filenames, collapsed over detector. Division by the tabulated correction"
              " is an illustrative scale test. It does not reproduce an uncorrected extraction or establish absolute flux."))
    for r in comparison:
        print(r["mode"],len(r["dja"]["detector_independent_exposure_keys"]),r["dja"]["lya_path_corr_median"],
              r["dja"]["formal_division_by_pathcorr"],r["jades"]["uv"],r["dja"]["uv"],flush=True)


if __name__=="__main__":
    main()

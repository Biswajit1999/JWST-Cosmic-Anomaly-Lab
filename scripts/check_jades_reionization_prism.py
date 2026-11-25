"""Check Lyα in separate released PRISM observations with a convolved UV break."""

import json
import re
from itertools import product
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from audit_jades_reionization_exposures import arrays, read_intermediate
from compare_jades_reionization_dja import INVENTORY
from download_jades_reionization_sample import LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from measure_jades_reionization import C_KMS, LYA, spectral_arrays
from scipy.ndimage import gaussian_filter1d
from scipy.special import ndtr

from cosmic_anomaly_lab.exposure_audit import combine_extractions
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import read_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

DISP = "https://jwst-docs.stsci.edu/files/216455634/216455648/1/1762452062676/jwst_nirspec_prism_disp.fits"
STEM = "jades_reionization_prism_check"
RESOLUTION_FACTORS = (1., 1.5, 2.)
SLOPES = (-2.5, -2., -1.5)
BREAK_OFFSETS_A = (-3., 0., 3.)
NOMINAL = (1.5, -2., 0.)
WINDOW_REST = (1100., 1400.)
PRIMARY = "goods-n-mediumjwst_58441"


def archive_spectrum(cache, entry):
    tier, number = entry["dr4_id"].rsplit("_",1)
    directory = f"https://jades.herts.ac.uk/DR4/{tier[:7]}/spectra/clear-prism/{tier}/"
    listing = cache.fetch(directory,f"metadata/dr4_{tier}_prism_clear.html",survey="JADES",
                          version="DR4-listing",kind="archive_directory")
    candidates = [name for name in re.findall(r'href="([^"]+)"',listing.read_text())
                  if re.search(rf"-0*{int(number)}_clear-prism_.*_x1d\.fits$",name)]
    if len(candidates) != 1:
        raise ValueError(f"Expected one listed PRISM spectrum for {entry['dr4_id']}, found {len(candidates)}")
    url = urljoin(directory,candidates[0])
    path = cache.fetch(url,f"independent/{Path(candidates[0]).name}",survey="JADES",version="DR4",
                       kind="observed_prism_spectrum")
    spec = read_dr4_spectrum(path,wavelength_range=(.65,5.4))
    sep = SkyCoord(entry["ra"],entry["dec"],unit="deg").separation(
        SkyCoord(spec.metadata["RA"],spec.metadata["DEC"],unit="deg")
    ).arcsec
    if sep > .2 or int(spec.metadata["TARGET"]) != int(number) or int(spec.metadata["PROPOSID"]) != entry["header_pid"]:
        raise ValueError(f"PRISM FITS identity check failed for {entry['dr4_id']}")
    return spec,dict(spectrum=receipt(path,url),listing=receipt(listing,directory),offset_arcsec=float(sep))


def design(edges, line_center, break_a, sigma, beta):
    """Bin-average both continuum bases after smoothing, and a unit-area line."""
    step = .4
    x = np.arange(edges[0]-8*sigma,edges[-1]+8*sigma+step,step)
    low = gaussian_filter1d((x < break_a).astype(float),sigma/step,mode="nearest")
    high = gaussian_filter1d(np.where(x >= break_a,(x/break_a)**beta,0.),sigma/step,mode="nearest")
    avg = []
    for basis in (low,high):
        integral = np.r_[0.,np.cumsum((basis[:-1]+basis[1:])*step/2)]
        avg.append(np.diff(np.interp(edges,x,integral))/np.diff(edges))
    avg.append(np.diff(ndtr((edges-line_center)/sigma))/np.diff(edges))
    return np.column_stack(avg)


def local_fit(data, redshift, offset, r, resolution_factor, beta, break_shift, minimum_area=.8):
    wave, edges, flux, error, valid = data
    unshifted = LYA*(1+redshift)
    center = unshifted*(1+offset/C_KMS)
    sigma = center/(2.354820045*r*resolution_factor)
    rescale = wave/(1+redshift)
    fit = (rescale >= WINDOW_REST[0]) & (rescale <= WINDOW_REST[1])
    indices = np.flatnonzero(fit)
    if not len(indices):
        return None,None
    region = slice(indices[0],indices[-1]+1)
    wave_data = wave[region]
    observed, uncertainty, good = flux[region],error[region],valid[region]
    if (good.sum() < 10 or np.sum(good & (wave_data < unshifted-2*sigma)) < 3
            or np.sum(good & (wave_data > unshifted+2*sigma)) < 3):
        return None,None
    edge_data = edges[indices[0]:indices[-1]+2]
    x = design(edge_data,center,(LYA+break_shift)*(1+redshift),sigma,beta)
    area = float(np.sum(np.diff(edge_data)[good]*x[good,2]))
    if area < minimum_area:
        return None,None
    weighted = x[good]/uncertainty[good,None]
    y = observed[good]/uncertainty[good]
    cov = np.linalg.inv(weighted.T @ weighted)
    theta = cov @ (weighted.T @ y)
    null,*_ = np.linalg.lstsq(weighted[:,:2],y,rcond=None)
    predicted,baseline = x@theta,x[:,:2]@null
    err = float(np.sqrt(cov[-1,-1]))
    value = float(theta[-1])
    total_stat = float(np.sum(((observed[good]-predicted[good])/uncertainty[good])**2))
    null_stat = float(np.sum(((observed[good]-baseline[good])/uncertainty[good])**2))
    result = dict(resolution_factor=resolution_factor,continuum_beta=beta,break_shift_rest_a=break_shift,
                  sigma_a=sigma,line_area_coverage=area,flux=value,error=err,snr=value/err,
                  fit_chi_squared=total_stat,delta_chi_squared=null_stat-total_stat,used_bins=int(good.sum()))
    return result,(wave_data,observed,uncertainty,good,predicted,baseline)


def source_diagnostics(cache, entry, products, data, z, offset, r):
    source = ROOT/products["spectrum"]["path"]
    two_nods = spectral_arrays(read_dr4_spectrum(source,nods=2,wavelength_range=(.65,5.4)))
    alternate_nods,_ = local_fit(two_nods,z,offset,r,*NOMINAL)
    wave,names,flux,err,_,_ = read_intermediate(source)
    group_id = np.array([row[2] for row in names])
    group_results = []
    for group in sorted(set(group_id)):
        keep = group_id == group
        mean,error,_ = combine_extractions(flux[keep],err[keep])
        measurement,_ = local_fit(arrays(wave,mean,error),z,offset,r,*NOMINAL)
        group_results.append(dict(region=int(group),integrations=int(keep.sum()),measurement=measurement))
    w,edges,f,error,good = data
    center = LYA*(1+z)*(1+offset/C_KMS)
    sigma = center/(2.354820045*r*NOMINAL[0])
    area = np.diff(ndtr((edges-center)/sigma))
    core = np.argsort(np.where(good,area,0.))[-3:][::-1]
    omitted = []
    for i in core:
        mask = good.copy()
        mask[i] = False
        measurement,_ = local_fit((w,edges,f,error,mask),z,offset,r,*NOMINAL,minimum_area=.40)
        omitted.append(dict(omitted_wavelength_a=float(w[i]),nominal_profile_area=float(area[i]),measurement=measurement))

    index = pd.read_csv(INVENTORY)
    index = index[(index.survey == "JADES") & (index.grating == "PRISM") & (index["filter"] == "CLEAR")]
    position = SkyCoord(entry["ra"],entry["dec"],unit="deg")
    close = position.separation(SkyCoord(index.ra,index.dec,unit="deg")).arcsec <= .15
    matches = index[close & (index.srcid == entry["header_target_number"])
                    & (index.programme == entry["header_pid"])]
    if len(matches) != 1:
        raise ValueError(f"Expected one target-matched DAWN PRISM product, got {len(matches)}")
    row = matches.iloc[0]
    path = cache.fetch(row.url,f"spectra/{row.file}",survey="DAWN JWST Archive",version="v4",kind="observed_spectrum")
    dja = read_spectrum(path,wavelength_range=(.65,5.4))
    off = position.separation(SkyCoord(dja.metadata["SRCRA"],dja.metadata["SRCDEC"],unit="deg")).arcsec
    if off>.15 or (dja.metadata["GRATING"].upper(),dja.metadata["FILTER"].upper()) != ("PRISM","CLEAR"):
        raise ValueError("DAWN PRISM header differs from archive inventory")
    dja_fit,_ = local_fit(spectral_arrays(dja),z,offset,r,*NOMINAL)
    return dict(two_nod_measurement=alternate_nods,slitlet_group_measurements=group_results,
                leave_one_core_bin_out=omitted,minimum_jackknife_profile_area=.40,
                dja=dict(product=receipt(path,row.url),inventory_sha256=digest(INVENTORY),archive_dataset=row.dataset,
                         header_offset_arcsec=float(off),measurement=dja_fit))


def main():
    input_path = OUT/"exposure_audit.json"
    parents = json.loads(input_path.read_text())
    catalogue = {r["dr4_id"]:r for r in json.loads(LEDGER.read_text())["records"]}
    cache = DownloadCache(ROOT,28_000*1024**2)
    dispersion = cache.fetch(DISP,"metadata/jwst_nirspec_prism_disp.fits",survey="STScI",
                             version="linked JDox curve",kind="instrument_reference")
    with fits.open(dispersion) as h:
        dl, resolution = (np.array(h[1].data[k],float) for k in ("WAVELENGTH","R"))
    if np.any(np.diff(dl) <= 0):
        raise ValueError("Instrument table wavelength must increase")
    records, plots = [], []
    for row in parents["candidates"]:
        sid = row["dr4_id"]
        entry = catalogue[sid]
        spec, product_receipts = archive_spectrum(cache,entry)
        z, offset = row["z"], row["offset_kms"]
        data = spectral_arrays(spec)
        r = float(np.interp(LYA*(1+z)/1e4,dl,resolution))
        variants = []
        nominal = None
        for factor,beta,shift in product(RESOLUTION_FACTORS,SLOPES,BREAK_OFFSETS_A):
            fit,pixels = local_fit(data,z,offset,r,factor,beta,shift)
            variants.append(dict(resolution_factor=factor,continuum_beta=beta,break_shift_rest_a=shift,measurement=fit))
            if (factor,beta,shift) == NOMINAL:
                nominal = fit
                if pixels is not None:
                    plots.append((row,pixels,fit))
        record = dict(dr4_id=sid,role=row["role"],z=z,fixed_velocity_kms=offset,**product_receipts,
                      nominal_table_r=r,nominal=nominal,variants=variants)
        if sid == PRIMARY:
            record["source_diagnostics"] = source_diagnostics(cache,entry,product_receipts,data,z,offset,r)
        records.append(record)
        print(sid,nominal,flush=True)
    show = [(row,x,fit) for row,x,fit in plots if row["dr4_id"] == PRIMARY
            or row["role"]=="published_control"]
    fig = Figure(figsize=(9.6,2.3*len(show)+.8),layout="constrained")
    axes = fig.subplots(len(show),1,squeeze=False).ravel()
    for (row,data,fit),ax in zip(show,axes,strict=True):
        wave,y,e,good,with_line,base = data
        rest = wave/(1+row["z"])
        ax.errorbar(rest[good],y[good],yerr=e[good],fmt="o",ms=3,color="#216e80",capsize=2,label="Observed PRISM")
        ax.step(rest,with_line,where="mid",color="#9c4430",label="Break + fixed line")
        ax.step(rest,base,where="mid",color="#7663a5",ls="--",label="Break only")
        ax.axvline(LYA,color=".4",ls=":",lw=.7)
        ax.set(ylabel="Fλ (10⁻²⁰ cgs/Å)",title=f"{row['dr4_id']} · z={row['z']:.3f} · formal line S/N {fit['snr']:.1f}")
        ax.legend(frameon=False,fontsize=8,ncols=3,loc="upper left")
    axes[-1].set_xlabel("Rest wavelength (Å)")
    fig.suptitle("Independent JWST disperser check · PRISM continuum and Lyα fits",fontsize=13)
    save_figure(fig,ROOT/f"figures/{STEM}")
    write_json(OUT/"prism_crosscheck.json",dict(
        executed_utc=utcnow(),exposure_audit_sha256=text_digest(input_path),script_sha256=text_digest(__file__),
        instrument_curve=receipt(dispersion,DISP),input_crosscheck_sha256=digest(OUT/"crosschecks.json"),
        candidates=records,nominal_scenario=dict(zip(["resolution_factor","continuum_beta","break_shift_a"],NOMINAL,strict=True)),
        scope="Separate disperser observations with line position fixed by G140M. Local continuum and Gaussian LSF models."
              " Quoted errors use released diagonal noise. No line-flux sign prior or fit-based redshift selection.",
    ))


if __name__=="__main__":
    main()

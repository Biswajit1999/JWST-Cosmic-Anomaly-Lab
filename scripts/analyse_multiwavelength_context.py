"""Enumerate positional candidates and retain signed HST/Spitzer measurements."""

import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from matplotlib.figure import Figure

from cosmic_anomaly_lab.multiwavelength import ab25_flux
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

root = Path(__file__).resolve().parents[1]
out = root / "results/multiwavelength"
targets = pd.read_csv(root / "results/jades_balmer_census/followup_queue.csv")
bands = {"f435w": .435, "f606w": .606, "f775w": .775, "f850lp": .9,
         "f125w": 1.25, "f140w": 1.4, "f160w": 1.6,
         "irac1": 3.6, "irac2": 4.5, "irac3": 5.8, "irac4": 8.0}
catalogues, xray, inputs, sizes = {}, {}, {}, {}
for field, short, xtable in [("GN", "goodsn", "cdfn2msnew"), ("GS", "goodss", "chandfs7ms")]:
    path = root / f"data/raw/ancillary/{short}_3dhst.v4.1.cats.tar.gz"
    with tarfile.open(path) as archive:
        member = next(m for m in archive.getmembers() if m.name.endswith(".cat.FITS"))
        with fits.open(io.BytesIO(archive.extractfile(member).read())) as hdus:
            catalogues[field] = Table(hdus[1].data).to_pandas()
    xpath = root / f"data/raw/ancillary/{xtable}.xml"
    xray[field] = Table.read(xpath, format="votable").to_pandas()
    sizes[field] = dict(photometric_records=len(catalogues[field]), xray_records=len(xray[field]))
    for p in [path, xpath]:
        inputs[p.relative_to(root).as_posix()] = digest(p)
neighbours, photometry, xneighbours, summary = [], [], [], []
for target in targets.itertuples():
    sky = SkyCoord(target.ra * u.deg, target.dec * u.deg)
    cat = catalogues[target.field]
    sep = sky.separation(SkyCoord(cat.ra.to_numpy() * u.deg, cat.dec.to_numpy() * u.deg)).arcsec
    indices = np.flatnonzero(sep <= .3)
    xcat = xray[target.field]
    xsep = sky.separation(SkyCoord(xcat.ra.to_numpy() * u.deg, xcat.dec.to_numpy() * u.deg)).arcsec
    xi = np.flatnonzero(xsep <= 2.)
    summary.append(dict(dr4_id=target.dr4_id, z=target.z, photometric_candidates=len(indices),
                        nearest_photometric_arcsec=float(sep.min()), xray_candidates=len(xi),
                        nearest_xray_arcsec=float(xsep.min())))
    for i in indices:
        row = cat.iloc[i]
        neighbours.append(dict(dr4_id=target.dr4_id, field=target.field, catalogue_id=int(row.id),
            ra=row.ra, dec=row.dec, separation_arcsec=float(sep[i]),
            use_phot=int(row.use_phot), star_flag=int(row.star_flag), flags=int(row["flags"]),
            contam_flag=int(row.contam_flag), association="POSITIONAL_CANDIDATE_ONLY"))
        for band, wavelength in bands.items():
            if "f_" + band not in row.index:
                continue
            f, e, w = (float(row[prefix + band]) for prefix in ["f_", "e_", "w_"])
            flux, error = ab25_flux(f, e, w)
            contamination = float(row.get(band + "_contam", np.nan))
            valid = flux is not None
            photometry.append(dict(dr4_id=target.dr4_id, catalogue_id=int(row.id),
                band=band, nominal_wavelength_um=wavelength,
                raw_flux_ab25=f, raw_error_ab25=e, relative_weight=w,
                flux_uJy=flux, error_uJy=error,
                use_phot=int(row.use_phot), contamination_ratio=contamination,
                usable_for_display=bool(valid and row.use_phot == 1),
                signal_to_noise=flux/error if valid else None,
                status="SIGNED_MEASUREMENT" if valid else "INVALID_OR_MISSING"))
    for i in xi:
        row = xcat.iloc[i]
        xneighbours.append(dict(dr4_id=target.dr4_id, catalogue_name=str(row["name"]),
            ra=float(row.ra), dec=float(row.dec), separation_arcsec=float(xsep[i]),
            positional_error_68pct_arcsec=float(row.error_radius),
            association="POSITIONAL_CANDIDATE_ONLY"))
pd.DataFrame(neighbours).to_csv(out / "photometric_neighbours.csv", index=False)
photo = pd.DataFrame(photometry)
photo.to_csv(out / "signed_photometry.csv", index=False)
pd.DataFrame(xneighbours, columns=["dr4_id", "catalogue_name", "ra", "dec", "separation_arcsec",
                                 "positional_error_68pct_arcsec", "association"]).to_csv(out / "xray_neighbours.csv", index=False)
pd.DataFrame(summary).to_csv(out / "association_summary.csv", index=False)
fig = Figure(figsize=(10, 7), layout="constrained")
for ax, number in zip(fig.subplots(2, 2).ravel(), [66336, 89464, 10009106, 25734], strict=True):
    target = targets[targets.dr4_id.str.endswith("_" + str(number))].iloc[0]
    block = photo[(photo.dr4_id == target.dr4_id) & photo.usable_for_display]
    for cid, group in block.groupby("catalogue_id"):
        for irac, color, label in [(False, "#126e82", "HST"), (True, "#b94d36", "Spitzer")]:
            sub = group[group.band.str.startswith("irac") == irac]
            ax.errorbar(sub.nominal_wavelength_um, sub.flux_uJy, yerr=sub.error_uJy,
                        fmt="o", ms=4, color=color, label=f"{label}; 3D-HST {cid}")
    if len(block):
        ax.legend(fontsize=7)
    else:
        ax.text(.5, .5, "No usable candidate photometry", transform=ax.transAxes, ha="center")
    ax.axhline(0, color="gray", lw=.7)
    ax.set(xscale="log", xlabel="Nominal observed wavelength (µm)", ylabel="Signed flux density (µJy)",
           title=f"{target.field} {number} | z={target.z:.3f}")
fig.suptitle("Hubble + Spitzer catalogue context for four JWST targets\nPositional candidates within 0.3 arcsec; total fluxes, not matched apertures", fontsize=12)
save_figure(fig, root / "figures/multiwavelength/hst_spitzer_context")
write_json(out / "context_summary.json", dict(executed_utc=utcnow(), target_records=len(targets),
    catalogue_rows=sizes, photometric_radius_arcsec=.3, xray_radius_arcsec=2.,
    targets_with_photometric_candidates=sum(s["photometric_candidates"] > 0 for s in summary),
    targets_with_xray_candidates=sum(s["xray_candidates"] > 0 for s in summary),
    photometric_neighbour_records=len(neighbours), photometric_band_rows=len(photo),
    valid_signed_band_rows=int((photo.status == "SIGNED_MEASUREMENT").sum()),
    display_usable_band_rows=int(photo.usable_for_display.sum()),
    negative_valid_flux_rows=int(((photo.flux_uJy < 0) & (photo.status == "SIGNED_MEASUREMENT")).sum()),
    association_is_confirmed_identity=False, xray_upper_limits_measured=False,
    physical_classification="NOT_ESTABLISHED", input_sha256=inputs,
    script_sha256=text_digest(Path(__file__))))
print(pd.DataFrame(summary).to_string(index=False))

"""Measure baseline-subtracted six-band NIRCam photometry of the GOODS-S transient."""

import json
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from download_jwst_cross_detector import API, OUT, ROOT, download
from download_jwst_transient_grism_geometry import STORAGE_MIB
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import pixel_scale, show_stamp
from vet_jwst_field_b_candidate import FIGURES, vet_mosaic

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

FADING_DATES = ("jw08060003001", "jw08060005001", "jw05324027001")
DETECTORS = {
    "F115W": [("jw08060002001", "NRCA2"), ("jw12577002001", "NRCA2")],
    "F150W": list(zip(FADING_DATES, ("NRCB3", "NRCB1", "NRCB1"), strict=True)),
    "F277W": [("jw08060002001", "NRCALONG"), ("jw12577002001", "NRCALONG")]
             + [(v, "NRCBLONG") for v in FADING_DATES],
}
FILTERS = ("F115W", "F150W", "F200W", "F277W", "F356W", "F444W")
RECENT = ("jw08060002001", *FADING_DATES)
FAINT = "jw12577002001"


def aperture_excess(bright, faint):
    values = []
    for measured, baseline in zip(bright["apertures"], faint["apertures"], strict=True):
        flux, noise = measured["flux_njy"], measured["local_mad_sigma_njy"]
        level, rms = baseline["flux_njy"], baseline["local_mad_sigma_njy"]
        values.append({
            "radius_arcsec": measured["radius_arcsec"],
            "excess_njy": float(flux-level) if flux is not None and level is not None else None,
            "local_noise_njy": float(np.hypot(noise, rms)) if noise is not None and rms is not None else None,
        })
    return values


def observed_colour(excesses, visit, first, second):
    keys = (visit, first), (visit, second)
    if not all(k in excesses for k in keys):
        return None
    rows = [excesses[k][-1] for k in keys]
    flux = np.array([r["excess_njy"] for r in rows], dtype=float)
    noise = np.array([r["local_noise_njy"] for r in rows], dtype=float)
    if not (np.isfinite(flux+noise).all() and np.all(flux > 0)):
        return None
    return dict(ab_colour=float(-2.5*np.log10(flux[0]/flux[1])),
                local_scatter_mag=float(2.5/np.log(10)*np.hypot(*(noise/flux))))


def plot(measurements, excesses, stamps):
    figure = Figure(figsize=(10.0, 7.0), layout="constrained")
    grid = figure.add_gridspec(2, 4, height_ratios=(1.5, 1))
    ax = figure.add_subplot(grid[0, :])
    colours = ("#204168", "#9c4f2d", "#66538d", "#35806a")
    for visit, colour in zip(RECENT, colours, strict=True):
        observed = [(i, excesses[visit, band][-1]) for i, band in enumerate(FILTERS)
                    if (visit, band) in excesses]
        time = next(x["dateobs"][:10] for x in measurements.values() if x["visit"] == visit)
        ax.errorbar([i for i, _ in observed],
                    np.array([m["excess_njy"] for _, m in observed], dtype=float),
                    yerr=np.array([m["local_noise_njy"] for _, m in observed], dtype=float),
                    fmt="o", color=colour, capsize=2, label=time)
    ax.set(xticks=range(len(FILTERS)), xticklabels=FILTERS, ylabel="Excess above October 2026 (nJy)")
    ax.grid(axis="y", alpha=.25)
    ax.legend(ncols=2)
    for i, band in enumerate(("F115W", "F277W")):
        images = [stamps[(visit, band)] for visit in (RECENT[0], FAINT)]
        limits = []
        for im, _, x, y, _, _ in images:
            limits.append(im[round(y)-7:round(y)+8, round(x)-7:round(x)+8].ravel())
        pair = np.concatenate(limits)
        upper = np.nanpercentile(pair, 99)
        for j, (im, wcs, x, y, _, _) in enumerate(images):
            p = figure.add_subplot(grid[1, i*2+j])
            show_stamp(p, im, x, y, pixel_scale(wcs), (-.08, upper), "gray", .24, half_arcsec=.58)
            p.set_title(f"{band} · {2025+j}", fontsize=9)
    figure.suptitle("GOODS-S transient · six-filter photometry\n"
                    "Observed NIRCam images · signed 0.24″ aperture excess", fontsize=13)
    save_figure(figure, FIGURES / "field_b_transient_colours")


def main():
    path = OUT / "field_b_transient_vetting.json"
    f150_path = OUT / "field_b_transient_f150w.json"
    initial = json.loads(path.read_text())
    f150 = json.loads(f150_path.read_text())
    sky = initial["position_degrees"]
    inputs = initial["observed_aperture_flux"]+f150["measurements"]
    measurements = {(x["visit"], x["band"]): x for x in inputs if x["visit"] in (*RECENT, FAINT)}
    cache = DownloadCache(ROOT, STORAGE_MIB*1024**2)
    receipts, associations, stamps = [], [], {}
    for band, visits in DETECTORS.items():
        f = band+"-CLEAR"
        url = API+"/assoc_mosaic?" + urlencode(dict(coords=f"{sky[0]:.8f},{sky[1]:.8f},.02", filters=f, output="csv"))
        query = download(cache, url, f"transient_colour_{band}_associations.csv", "archive_filter_coverage")
        associations.append(query)
        table = pd.read_csv(ROOT / query["path"])
        for visit, detector in visits:
            work = table.query("visit == @visit and detector == @detector and tile == 672 and status == 2")
            if len(work) != 1:
                raise ValueError(f"Expected one direct exposure association for {visit} {band}")
            v = work.iloc[0]
            record = dict(band=f, detector=detector, visit=visit, dateobs=v.dateobs, mjd=float(v.expstart),
                          assoc_name=v.assoc_name, exposure_count=int(v.nexp))
            for key in ("sci", "wht"):
                url = v.file.replace("_sci.fits.gz", f"_{key}.fits.gz")
                record[key] = download(cache, url, url.rsplit("/", 1)[-1], "observed_colour_epoch_"+key)
            receipts.append(record)
            stamps[visit, band], measurements[visit, band] = vet_mosaic(record, sky)
            print(band, visit, measurements[visit, band]["apertures"][-1]["flux_njy"], flush=True)

    excesses = {(v, b): aperture_excess(row, measurements[FAINT, b])
                for (v, b), row in measurements.items() if v in RECENT and (FAINT, b) in measurements}
    result = []
    for visit in RECENT:
        result.append(dict(
            visit=visit, date_utc=min(x["dateobs"] for x in measurements.values() if x["visit"] == visit),
            excesses={b: excesses[visit, b] for b in FILTERS if (visit, b) in excesses},
            f150w_minus_f356w=observed_colour(excesses, visit, "F150W", "F356W"),
            f200w_minus_f277w=observed_colour(excesses, visit, "F200W", "F277W"),
        ))
    plot(measurements, excesses, stamps)
    write_json(OUT / "field_b_transient_colours.json", {
        "executed_utc": utcnow(), "candidate": initial["candidate"], "position_degrees": sky,
        "reference_visit": FAINT, "raw_measurements": list(measurements.values()),
        "measurements_by_date": result, "archive_receipts": receipts, "association_queries": associations,
        "scope": "No aperture correction; local noise excludes shared baseline and systematic calibration effects",
        "original_photometry_sha256": text_digest(path), "early_f150w_sha256": text_digest(f150_path),
        "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

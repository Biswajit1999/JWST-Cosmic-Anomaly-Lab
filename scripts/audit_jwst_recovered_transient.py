"""Verify an earlier GOODS-S lead against NIRSpec and a separate NIRCam channel."""

import json

import numpy as np
from astropy.io import fits
from download_jwst_cross_detector import CONTROL, OUT, ROOT
from follow_jwst_fading_source import measure
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import (
    load_mosaic,
    pixel_scale,
    reference_sky,
    show_stamp,
    verified_path,
)

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

FIGURES = ROOT / "figures/jwst_cross_detector"


def main():
    spectra_path = OUT / "control_spectra_products.json"
    mosaics_path = OUT / "mosaics.json"
    source = json.loads(spectra_path.read_text())
    recs = {r["epoch"]: r for r in json.loads(mosaics_path.read_text())["images"] if r["group"] == "control"}
    images = [load_mosaic(recs[label], CONTROL, half=830) for label in ("older", "newer")]
    refs = reference_sky(images[0])
    older, newer = [measure(image, refs) for image in images]
    old_ref, new_ref = older[2], newer[2]
    valid = np.isfinite(old_ref+new_ref) & (old_ref > 40) & (new_ref > 0)
    calibration = (new_ref[valid]/old_ref[valid]).tolist()
    ratio = float(np.median(calibration)) if len(calibration) >= 6 else None
    recentered = [float(m[1]) if np.isfinite(m[1]) else None for m in (older, newer)]
    photometry = dict(older=older[0], newer=newer[0], usable_reference_count=len(calibration),
                      new_over_old_reference_ratio=ratio, reference_ratios=calibration,
                      raw_centroid_aperture_njy=recentered)
    if ratio and None not in recentered:
        photometry["relative_2025_excess_recentered"] = recentered[0]/(recentered[1]/ratio)-1
    figure = Figure(figsize=(10.2, 6.2), layout="constrained")
    grid = figure.add_gridspec(2, 4, height_ratios=[1, .75])
    ax = figure.add_subplot(grid[0, :])
    measures = []
    colours = ("#243851", "#4a82a7", "#9a4f2d", "#368772")
    for p, colour in zip(sorted(source["spectra"], key=lambda r: r["dataset"].split("_")[-2]), colours, strict=True):
        with fits.open(verified_path(p)) as h:
            table = h["SPEC1D"]
            if [col.unit for col in table.columns[:3]] != ["um", "uJy", "uJy"]:
                raise ValueError(f"Unexpected spectral units in {p['dataset']}")
            wave, flux, err = (np.asarray(table.data[c]) for c in ("wave", "flux", "err"))
            date = table.header["DATE-OBS"]
            measures.append(dict(dataset=p["dataset"], date_utc=date, start_mjd=float(table.header["MJD-BEG"]),
                                 aperture_pa_deg=float(h["SCI"].header["PA_APER"]),
                                 dja_z=p["dja_z"], dja_comment=p["dja_comment"],
                                 separation_arcsec=p["separation_arcsec"]))
        valid = np.isfinite(wave+flux+err) & (wave >= .8) & (wave <= 5.25) & (err > 0)
        ax.plot(wave[valid], flux[valid], color=colour, lw=1, alpha=.85,
                label=f"{date} · GO {int(p['dataset'].split('_')[-2])}")
    ax.set(xlabel="Observed wavelength (µm)", ylabel="Extracted Fν (µJy)", xlim=(.8, 5.25), ylim=(-.06, 1.45))
    ax.legend(fontsize=8, ncols=2)
    center_stamps = [g[0][round(g[3])-15:round(g[3])+16, round(g[2])-15:round(g[2])+16] for g in images]
    bounds = tuple(np.nanpercentile(np.concatenate([s.ravel() for s in center_stamps]), [4, 99]))
    for columns, (image, wcs, x, y, _, _), raw, date in zip(
            (slice(0, 2), slice(2, 4)), images, (older[0], newer[0]), ("2025-10-06", "2026-10-09"), strict=True):
        a = figure.add_subplot(grid[1, columns])
        show_stamp(a, image, x, y, pixel_scale(wcs), bounds, "gray", .24, half_arcsec=.47)
        a.set_title(f"F200W · {date} · {raw['flux_njy']:.1f} nJy", fontsize=10)
    figure.suptitle("ZFOURGE CDFS 18559 · archival spectroscopic recovery\n"
                    "Observed PRISM extractions and independent short-wave NIRCam imaging", fontsize=13)
    save_figure(figure, FIGURES / "recovered_control")
    write_json(OUT / "recovered_control.json", {
        "executed_utc": utcnow(), "position_degrees": CONTROL, "photometry_f200w": photometry,
        "spectroscopic_observations": sorted(measures, key=lambda m: m["start_mjd"]),
        "disposition": "PRIOR_ARCHIVE_TRANSIENT_RECOGNITION_AND_SPECTROSCOPIC_FOLLOWUP",
        "spectra_sha256": text_digest(spectra_path), "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

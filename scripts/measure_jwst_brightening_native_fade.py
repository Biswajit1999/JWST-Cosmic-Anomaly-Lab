"""Verify the November fade near JADES 103607 on the original detector pixels."""

import json
from warnings import catch_warnings, filterwarnings

import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import fetch
from matplotlib import colormaps
from matplotlib.figure import Figure
from measure_jwst_brightening_native import (
    apertures,
    compare,
    local_controls,
    native_metadata,
    positions_from_mosaic,
    stamps,
    zero_dq_core,
)
from measure_jwst_brightening_november import FILTERS
from measure_jwst_brightening_offset import BANDS
from measure_jwst_epoch_candidate import RADII, show_stamp
from measure_jwst_module_b_brightening import NEW, OLD
from search_jwst_cross_detector import window
from vet_jwst_field_b_candidate import FIGURES
from vet_jwst_module_b import native_alignment

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

NOVEMBER = "jw12577004001"
CLOUD = "https://stpubdata.s3.us-east-1.amazonaws.com/"
BAND_DETECTORS = {
    "F150W": ("04101", "nrca1"),
    "F200W": ("02101", "nrca1"),
    "F277W": ("02101", "nrcalong"),
    "F356W": ("04101", "nrcalong"),
}
LABELS = ("Oct 2025", "Oct 2026", "Nov 2026")
PAIRS = (("rise", 0, 1), ("november_minus_baseline", 0, 2), ("fade", 1, 2))
STEM = "module_b_001_native_november_fade"


def nova_products(cache, band, keys):
    activity, detector = BAND_DETECTORS[band]
    records = []
    for dither in range(1, 4):
        name = f"{NOVEMBER}_{activity}_{dither:05d}_{detector}_cal.fits"
        key = f"jwst/public/jw12577/{NOVEMBER}/{name}"
        if key not in keys:
            raise ValueError(f"Missing public November archive entry for {name}")
        item = fetch(
            cache,
            CLOUD + key,
            name,
            "observed_module_b_native",
            version="MAST public calibrated exposure",
        )
        header = fits.getheader(ROOT / item["path"])
        records.append(dict(item, file=name, expstart_mjd=float(header["EXPSTART"])))
    return dict(visit=NOVEMBER, products=records)


def measure_band(cache, sky, band, groups, base, reference_rows, keys, host_core_sigma):
    refs = positions_from_mosaic(base, sky, reference_rows)
    controls = local_controls(window(base, "sci", sky, 16.5)[1], sky)
    frames, audits, inventory = [], [], []
    for group in (*groups, nova_products(cache, band, keys)):
        visit = group["visit"]
        detector = (
            BAND_DETECTORS[band][1].upper()
            if visit == NOVEMBER
            else group["fits_headers"][0]["detector"]
        )
        metadata = native_metadata(group, visit, band, detector)
        with catch_warnings():
            filterwarnings(
                "ignore", message="'(datfix|obsfix)' made the change", category=FITSFixedWarning
            )
            batch, audit = native_alignment(base, group, sky)
        if len(batch) != 3:
            raise ValueError(f"Expected three covering dithers for {visit}/{band}")
        frames.extend(batch)
        audits.extend(audit)
        inventory.append(dict(visit=visit, products=group["products"], fits_headers=metadata))
    target, local, measurements, calibration = apertures(frames, refs, controls, sky)
    contrasts = {}
    for label, old, new in PAIRS:
        indices = np.r_[old * 3 : old * 3 + 3, new * 3 : new * 3 + 3]
        contrasts[label] = compare(
            target[indices], local[indices], [measurements[i] for i in indices]
        )
    result = dict(
        inventory=inventory,
        native_frames=[
            dict(**diagnostic, aperture_photometry=phot)
            for diagnostic, phot in zip(audits, measurements, strict=True)
        ],
        contrasts=contrasts,
        calibration=calibration,
    )
    if host_core_sigma is not None:
        trials = []
        for factor in (0.85, 1, 1.15):
            sigma = factor * host_core_sigma
            try:
                core = zero_dq_core(frames[3:], audits[3:], sigma, BANDS[band], measurements[3:])
            except ValueError as exc:
                trials.append(dict(fixed_sigma_arcsec=sigma, reason=str(exc)))
                continue
            # Compare every November fit to every October fit at the same width.
            for field in (
                "fitted_core_amplitude_njy_per_pixel",
                "without_brightest_pixel_njy_per_pixel",
            ):
                values = [c[field] for c in core["exposures"]]
                core[field + "_all_november_below_all_october"] = max(values[3:]) < min(values[:3])
            trials.append(core)
        result["fade_zero_dq_core_fits"] = trials
    return result, stamps(frames, audits) if band == "F200W" else None


def plot(bands, visits):
    fig = Figure(figsize=(11.5, 8.4), layout="constrained")
    grid = fig.add_gridspec(4, 3, height_ratios=(1.45, 1, 1, 1))
    ax = fig.add_subplot(grid[0, :])
    for i, (band, colour) in enumerate(
        zip(FILTERS, ("#1b6c8d", "#a04f29", "#695da7", "#417348"), strict=True)
    ):
        labelled = False
        for j, (label, _, _) in enumerate(PAIRS):
            result = bands[band]["contrasts"][label][-1]
            d, e = result["difference_njy"], result["local_difference_mad_njy"]
            if d is not None:
                ax.errorbar(
                    j + (i - 1.5) * 0.08,
                    d,
                    yerr=e,
                    fmt="o",
                    color=colour,
                    capsize=3,
                    label=band if not labelled else None,
                )
                labelled = True
    ax.set(
        ylabel="Signed 0.24″ flux change (nJy)",
        xticks=[0, 1, 2],
        xticklabels=[
            "October 2026 − October 2025",
            "November 2026 − October 2025",
            "November − October 2026",
        ],
    )
    ax.axhline(0, lw=0.7, color=".55")
    ax.legend(frameon=False, ncol=4)
    floor, top = np.nanpercentile(np.concatenate([d[0].ravel() for d in visits]), [8, 98.5])
    for index, (image, x, y, scale) in enumerate(visits):
        epoch, dither = divmod(index, 3)
        panel = fig.add_subplot(grid[epoch + 1, dither])
        show_stamp(
            panel,
            image,
            x,
            y,
            scale,
            (floor, top),
            colormaps["gray"].with_extremes(bad="#9b64a3"),
            0.24,
            half_arcsec=0.57,
        )
        panel.set_title(f"F200W · {LABELS[epoch]} · dither {dither + 1}", fontsize=9)
    fig.suptitle("JADES 103607 · rise and fade on native NIRCam detector pixels", fontsize=14)
    save_figure(fig, FIGURES / STEM)


def main():
    native = OUT / "module_b_001_native_epoch_flux.json"
    photometry = OUT / "module_b_001_six_filter_change.json"
    products = OUT / "module_b_products.json"
    release = OUT / "module_b_001_release_20261125.json"
    first, old, catalogue, followup = [
        json.loads(p.read_text()) for p in (native, photometry, products, release)
    ]
    mast_comparison = OUT / "module_b_001_mast_comparison.json"
    mast = json.loads(mast_comparison.read_text())["bands"]
    sky = old["icrs_deg"]
    if first["icrs_deg"] != sky or followup["position_degrees"] != sky:
        raise ValueError("Inconsistent forced photometry position")
    keys = set(followup["s3"]["visit_products"]["keys"])
    mosaics = [
        r
        for r in catalogue["images"]
        if r["pointing"] == "001"
        and r["detector"] == ("NRCB4" if r["band"] == "F200W" else "NRCBLONG")
    ]
    mosaics.extend(old["additional_products"])
    cache = DownloadCache(ROOT, 30_000 * 1024**2)
    bands, views = {}, None
    for band in FILTERS:
        group = first["native_band_measurements"][band]["sources"]
        if [r["visit"] for r in group] != [OLD, NEW]:
            raise ValueError("Unexpected order for the October visits")
        (reference,) = [r for r in mosaics if r["band"] == band and r["visit"] == OLD]
        refs = old["difference_photometry"][band]["local_flux_scale"]["references"]
        sigma = (
            mast[band]["astrometry"]["host_core"]["core_sigma_arcsec"] if band in BANDS else None
        )
        bands[band], view = measure_band(cache, sky, band, group, reference, refs, keys, sigma)
        if view:
            views = view
        print(
            band,
            {
                k: (r[-1]["difference_njy"], r[-1]["local_difference_mad_njy"])
                for k, r in bands[band]["contrasts"].items()
            },
            flush=True,
        )
    plot(bands, views)
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            icrs_deg=sky,
            script_sha256=text_digest(__file__),
            input_sha256={
                p.name: text_digest(p)
                for p in (native, photometry, products, release, mast_comparison)
            },
            native_aperture_script_sha256=text_digest(
                ROOT / "scripts/measure_jwst_brightening_native.py"
            ),
            epoch_labels=list(LABELS),
            aperture_radii_arcsec=list(RADII),
            bands=bands,
            scope="Calibrated native SCI pixels. Local reference positions and sky control apertures are fixed across visits. "
            "Aperture differences have no total flux or cross-detector PSF correction. "
            "Both the archived stage-3 and native analyses include the same November photons.",
        ),
    )


if __name__ == "__main__":
    main()

"""Identify catalogued sources that disperse onto the forced transient trace."""

import json
from contextlib import ExitStack

import asdf
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from download_jwst_cross_detector import OUT
from extract_jwst_transient_grism import BINS
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import show_stamp, verified_path
from measure_jwst_transient_grism_geometry import (
    direct_cutout,
    dispersed_trace,
    imaging_transform,
    pointing_header,
)
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

OVERLAP_RADIUS = 3
GRID = np.linspace(BINS[0], BINS[-1], 501)


def overlaps(record, calibrations, position, table):
    _, head = pointing_header(verified_path(record))
    catalogue_coords = table["sky_centroid"]
    pos = SkyCoord(*position, unit="deg")
    with ExitStack() as stack:
        refs = {k: stack.enter_context(asdf.open(verified_path(calibrations[v])))
                for k, v in record["references"].items()}
        (filter_offset,) = [row for row in refs["R_FILOFF"]["filters"]
                            if row["filter"] == "F356W" and row["pupil"] == "GRISMR"]
        wcs = imaging_transform(head, refs["R_DISTOR"]["model"], filter_offset)
        ox, oy = wcs.inverse(*position)
        x0, y0 = wcs.inverse(catalogue_coords.ra.deg, catalogue_coords.dec.deg)
        spec = refs["R_SPCWCS"]
        v = head["VELOSYS"]
        x, y, wave = dispersed_trace(spec, 1, ox, oy, v)
        cx, cy = np.interp(GRID, wave, x), np.interp(GRID, wave, y)
        rows, trace_offsets = [], {}
        aperture = (x0 >= -130) & (x0 <= 2164) & (y0 >= -130) & (y0 <= 2164)
        for order, filt, low, high in refs["R_WAVRAN"]["wavelengthrange"]:
            if filt != "F356W":
                continue
            look = aperture & (abs(y0-oy) <= 60) if order == 1 else aperture
            for i in np.flatnonzero(look):
                sx, sy, sw = dispersed_trace(spec, order, x0[i], y0[i], v)
                band = (sw >= low) & (sw <= high)
                separation = np.interp(cx, sx[band], sy[band], left=np.nan, right=np.nan)-cy
                mask = np.isfinite(separation) & (abs(separation) <= OVERLAP_RADIUS)
                if not mask.any():
                    continue
                star = table[i]
                coord = catalogue_coords[i]
                flux = star["aper_total_flux"]
                flux_ujy = float(flux)*1e6 if np.isfinite(flux) and not np.ma.is_masked(flux) else None
                entry = dict(catalogue_label=int(star["label"]), icrs=[float(coord.ra.deg), float(coord.dec.deg)],
                             direct_flux_microjy=flux_ujy, source_separation_arcsec=float(coord.separation(pos).arcsec),
                             order=order, diagnostic_grid_overlap=float(mask.mean()),
                             median_crossdisp_offset_pix=float(np.median(separation[mask])),
                             max_abs_crossdisp_offset_pix=float(np.max(abs(separation[mask]))))
                rows.append(entry)
                trace_offsets[order, entry["catalogue_label"]] = separation
    return rows, trace_offsets


def figure(source, pointing, projected_rows, traces):
    (image, x, y, scale), _ = direct_cutout(pointing, source)
    strongest = max(projected_rows[0], key=lambda row: row["direct_flux_microjy"] or 0)
    key = strongest["order"], strongest["catalogue_label"]
    from astropy.io import fits
    from astropy.wcs import WCS

    with fits.open(verified_path(pointing), memmap=True) as h:
        wcs = WCS(h["SCI"].header).celestial
        gx, gy = wcs.world_to_pixel_values(*strongest["icrs"])
        px, py = wcs.world_to_pixel_values(*source)
    dx, dy = (gx-px)*scale, (gy-py)*scale
    fig = Figure(figsize=(9.0, 4.1), layout="constrained")
    sky, detector = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1, 1.26]})
    show_stamp(sky, np.arcsinh(image/.4), x, y, scale, (-.1, 4.6), "gray", .24, half_arcsec=3.1)
    sky.scatter(dx, dy, marker="+", c="#ffc66e", s=120, linewidths=2)
    sky.annotate(f'Source {strongest["catalogue_label"]}', (dx, dy), xytext=(-5, -15),
                 textcoords="offset points", ha="right", color="#5e3b13", fontsize=8,
                 bbox={"facecolor": "#fff5df", "edgecolor": "none", "alpha": .9, "pad": 2})
    sky.annotate("Transient", (0, 0), xytext=(0, 13), textcoords="offset points", ha="center",
                 color="#174c61", fontsize=8,
                 bbox={"facecolor": "#dff6f9", "edgecolor": "none", "alpha": .9, "pad": 2})
    sky.set_title("Observed F356W direct image", fontsize=10)
    for i, row in enumerate(traces):
        if key in row:
            detector.plot(GRID, row[key], color="#a77833", alpha=.65, label="Catalogued trace" if i == 0 else None)
    detector.axhline(0, color="#21647d", lw=.85, label="Transient trace")
    detector.set(xlabel="Wavelength along transient trace (µm)", ylabel="Relative detector row (pixels)",
                 title="Six projected slitless traces")
    detector.grid(alpha=.17)
    detector.legend(fontsize=8, loc="upper right")
    fig.suptitle("Bright nearby source overlaps the GOODS-S transient spectrum", fontsize=12)
    save_figure(fig, FIGURES / "field_b_transient_grism_overlap")


def main():
    geom_path = OUT / "field_b_transient_grism_products.json"
    geom = json.loads(geom_path.read_text())
    audit_path = OUT / "field_b_transient_wfss.json"
    audit = json.loads(audit_path.read_text())
    catalogue = next(r for r in audit["direct_image_catalogues"] if "o006" in r["filename"])
    table = Table.read(verified_path(catalogue))
    if str(table["aper_total_flux"].unit) != "Jy":
        raise ValueError("Expected the published direct-image flux column in Jy")
    results, traces = [], []
    sky = geom["position_degrees"]
    for record in geom["rate_headers"]:
        if not record["rate_filename"].startswith("jw04549006001"):
            continue
        rows, offsets = overlaps(record, geom["calibrations"], sky, table)
        results.append({"filename": record["rate_filename"], "overlapping_sources": rows})
        traces.append(offsets)
        print(record["rate_filename"], rows, flush=True)
    direct = next(r for r in geom["direct_images"] if "o006" in r["path"])
    figure(sky, direct, [r["overlapping_sources"] for r in results], traces)
    write_json(OUT / "field_b_transient_grism_overlap.json", {
        "executed_utc": utcnow(), "candidate": geom["candidate"], "position_degrees": sky,
        "diagnostic_window_um": [float(GRID[0]), float(GRID[-1])], "crossdisp_tolerance_pix": OVERLAP_RADIUS,
        "archive_catalogue": catalogue, "overlap_results": results,
        "calculation_scope": "Centroid traces for first and second orders, not a contaminant flux model",
        "grism_products_sha256": text_digest(geom_path), "catalogue_audit_sha256": text_digest(audit_path),
        "script_sha256": text_digest(__file__),
    })


if __name__ == "__main__":
    main()

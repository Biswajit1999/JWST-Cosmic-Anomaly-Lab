"""Measure fixed spectral-line apertures for the possible JADES 103607 host."""

import json
from collections import defaultdict
from urllib.parse import urlsplit

import numpy as np
from astropy.io import fits
from audit_jwst_november_f200w_spike import NativeSections
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import CRDS, fetch
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import verified_path
from read_jwst_i2d_cutout import end_of_block
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

STEM = "module_b_001_grism_host_lines"
CONTROLS = (-34, -28, -21, -14, 14, 21, 28, 34)
HALF_ROWS = 50
JUMP_DET = 4


def rate_stripes(cache, row, start, stop):
    native = NativeSections(cache, row["file"])
    if native.etag != row["rate_etag"]:
        raise ValueError("Detector data no longer match the geometry audit")
    primary, cursor = native.header(0)
    if primary["FILENAME"] != row["file"] or primary["EXP_TYPE"] != "NRC_WFSS":
        raise ValueError("Unexpected grism detector metadata")
    arrays = {}
    for _ in range(8):
        hdr, base = native.header(cursor)
        width, height = hdr["NAXIS1"], hdr["NAXIS2"]
        kind = hdr["EXTNAME"]
        if kind in {"SCI", "DQ"}:
            expect = -32 if kind == "SCI" else 32
            if hdr["BITPIX"] != expect or hdr.get("BSCALE", 1) != 1:
                raise ValueError("Unsupported FITS scaling")
            if kind == "SCI" and (hdr["BUNIT"] != "DN/s" or hdr.get("BZERO", 0) != 0):
                raise ValueError("Unexpected grism rate units")
            path = native.read(base + start * width * 4, base + stop * width * 4)
            dtype = ">f4" if kind == "SCI" else ">i4"
            data = np.fromfile(path, dtype=dtype).reshape(stop - start, width)
            arrays[kind] = (
                data.astype(float) if kind == "SCI" else data.astype(np.int64) + hdr.get("BZERO", 0)
            )
        if set(arrays) == {"SCI", "DQ"}:
            return arrays, native.ranges
        cursor = base + end_of_block(width * height * abs(hdr["BITPIX"]) // 8)
    raise ValueError("Missing grism SCI or DQ data")


def line_cutout(pixels, x, y, offset, start):
    """Remove per-column cross-trace background and a local continuum."""
    cols = np.arange(round(x) - 13, round(x) + 14)
    rows = np.arange(round(y + offset) - 10, round(y + offset) + 11)
    dy = rows - y - offset
    if rows[0] < start or rows[-1] >= start + len(pixels):
        return None, None
    scl = pixels[np.ix_(rows - start, cols)]
    backgrounds = (abs(dy) >= 5) & (abs(dy) <= 10)
    usable_background = np.isfinite(scl[backgrounds]).sum(axis=0) >= 5
    corrected = np.full_like(scl, np.nan)
    corrected[:, usable_background] = scl[:, usable_background] - np.nanmedian(
        scl[backgrounds][:, usable_background], axis=0
    )
    core = abs(dy) <= 2.5
    n = np.isfinite(corrected[core]).sum(axis=0)
    trace = np.divide(
        np.nansum(corrected[core], axis=0) * core.sum(),
        n,
        out=np.full(len(cols), np.nan),
        where=n >= core.sum() - 1,
    )
    central_window = abs(cols - x) <= 2.5
    center = central_window & np.isfinite(trace)
    side = (abs(cols - x) >= 6) & np.isfinite(trace)
    if center.sum() < central_window.sum() - 1 or side.sum() < 10:
        return None, None
    linear = np.polynomial.polynomial.polyfit(cols[side] - x, trace[side], deg=1)
    line = trace[center] - np.polynomial.polynomial.polyval(cols[center] - x, linear)
    return float(line.mean() * central_window.sum() * 1000), corrected


def exposure(cache, row, flats):
    cases = row["positions"]["jades_103607"]["conditional_lines"]
    targets = [
        (z, label, line)
        for z, lines in cases.items()
        for label, line in lines.items()
        if line["intersects"]
    ]
    if not targets:
        return None, {}
    x1d = fits.getheader(verified_path(row["x1d"]))
    name = urlsplit(x1d["R_FLAT"]).netloc
    if not name.startswith("jwst_nircam_flat_") or not name.endswith(".fits"):
        raise ValueError("Invalid grism flat-field reference")
    if name not in flats:
        flats[name] = fetch(
            cache,
            CRDS + name,
            name,
            "calibration_reference",
            survey="JWST CRDS",
            version=x1d["CRDS_CTX"],
        )
    trace_y = [round(line["xy_pix"][1]) for _, _, line in targets]
    start, stop = max(4, min(trace_y) - HALF_ROWS), min(2044, max(trace_y) + HALF_ROWS + 1)
    arrays, ranges = rate_stripes(cache, row, start, stop)
    with fits.open(verified_path(flats[name]), memmap=False) as flat:
        divisor = flat["SCI"].data[start:stop].astype(float)
        good = (
            (flat["DQ"].data[start:stop] == 0) & ((arrays["DQ"] & ~JUMP_DET) == 0) & (divisor > 0)
        )
        science = np.divide(arrays["SCI"], divisor, out=np.full_like(divisor, np.nan), where=good)
    cuts, records = {}, []
    for z, label, line in targets:
        x, y = line["xy_pix"]
        cut = [line_cutout(science, x, y, offset, start) for offset in (0, *CONTROLS)]
        cuts[z, label] = cut[0][1]
        records.append(
            dict(
                z=float(z),
                line=label,
                wavelength_um=line["wavelength_um"],
                xy_pix=[x, y],
                aperture_mdn_s=[v for v, _ in cut],
            )
        )
    group = row["visit"]
    if row["filter"] == "F444W":
        group += "_" + row["file"].split("_")[1]
    return dict(
        file=row["file"],
        filter=row["filter"],
        group=group,
        flat=name,
        data_ranges=ranges,
        apertures=records,
    ), cuts


def robust_scatter(values):
    return float(1.4826 * np.median(abs(values - np.median(values))))


def summaries(rows):
    grouped = defaultdict(list)
    for row in rows:
        for a in row["apertures"]:
            if all(v is not None for v in a["aperture_mdn_s"]):
                grouped[a["z"], a["line"], a["wavelength_um"], row["filter"]].append(
                    (row["group"], a["aperture_mdn_s"])
                )
    output = []
    for (z, line, w, filt), entries in sorted(grouped.items()):
        by_group = {}
        for group in sorted({g for g, _ in entries}):
            signals = [a for g, a in entries if g == group]
            by_group[group] = dict(
                exposures=len(signals), combined_mdn_s=np.mean(signals, axis=0).tolist()
            )
        if len(by_group) != 2 or min(d["exposures"] for d in by_group.values()) < 2:
            continue
        means = np.mean([g["combined_mdn_s"] for g in by_group.values()], axis=0)
        center = float(np.median(means[1:]))
        output.append(
            dict(
                z=z,
                line=line,
                wavelength_um=w,
                filter=filt,
                groups=by_group,
                exposures=len(entries),
                aperture_mdn_s=float(means[0]),
                controls_mdn_s=means[1:].tolist(),
                control_median_mdn_s=center,
                empirical_scatter_mdn_s=robust_scatter(means[1:]),
            )
        )
    return output


def plot(records, example):
    fig = Figure(figsize=(10, 6.4), layout="constrained")
    ax, image = fig.subplots(2, 1, gridspec_kw={"height_ratios": [2.6, 1]})
    for i, r in enumerate(sorted(records, key=lambda d: (d["z"], d["wavelength_um"]))):
        value = r["aperture_mdn_s"] - r["control_median_mdn_s"]
        ax.scatter(
            i + np.linspace(-0.11, 0.11, 8),
            np.asarray(r["controls_mdn_s"]) - r["control_median_mdn_s"],
            s=12,
            color="#a5afb3",
            zorder=2,
        )
        colour = "#2d697d" if r["filter"] == "F356W" else "#ab5839"
        ax.errorbar(
            i,
            value,
            yerr=r["empirical_scatter_mdn_s"],
            fmt="o",
            ms=4,
            capsize=3,
            color=colour,
            zorder=3,
        )
    labels = [
        f"{r['line']}\n{r['z']}"
        for r in sorted(records, key=lambda d: (d["z"], d["wavelength_um"]))
    ]
    ax.axhline(0, color="#7e898b", lw=0.7)
    ax.set_xticks(range(len(labels)), labels, rotation=40, ha="right", fontsize=8)
    ax.set(
        ylabel="Aperture excess (mDN/s)",
        title="Fixed apertures with local control scatter · line label and trial redshift",
    )
    if example:
        stack = np.ma.median(np.ma.masked_invalid(example), axis=0).filled(np.nan) * 1000
        lim = np.nanpercentile(abs(stack), 97)
        detector = image.imshow(
            stack,
            origin="lower",
            cmap="RdBu_r",
            vmin=-lim,
            vmax=lim,
            extent=(-13.5, 13.5, -10.5, 10.5),
            aspect="auto",
        )
        fig.colorbar(detector, ax=image, pad=0.015, label="mDN/s per pixel")
        for axis in (image.axvline, image.axhline):
            axis(0, color="#182c2f", ls=":", lw=0.75)
    image.set(
        ylabel="Cross-trace pixels",
        xlabel="Dispersion pixels relative to predicted line",
        title="Observed FRESCO detector residuals · Paschen beta at z = 2.47",
    )
    fig.suptitle("JADES 103607 · pre-event NIRCam slitless host check", fontsize=13)
    save_figure(fig, FIGURES / STEM)


def main():
    path = OUT / "module_b_001_grism_host_coverage.json"
    geometry = json.loads(path.read_text())
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    flats, rows, example = {}, [], []
    for frame in geometry["exposures"]:
        if "positions" not in frame:
            raise ValueError("Coverage must be evaluated before measuring grism apertures")
        row, cuts = exposure(cache, frame, flats)
        if row is not None:
            rows.append(row)
            x = cuts.get(("2.47", "Pa beta"))
            if x is not None:
                example.append(x)
            print(
                row["file"],
                sum(a["aperture_mdn_s"][0] is not None for a in row["apertures"]),
                flush=True,
            )
    records = summaries(rows)
    if not records:
        raise ValueError("No line apertures passed the two-group and control-pixel tests")
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            geometry_sha256=text_digest(path),
            calibration_flats=flats,
            control_offsets_pix=list(CONTROLS),
            allowed_rate_dq_bits=JUMP_DET,
            exposures=rows,
            line_results=records,
            scope="Instrumental rates from flat-fielded native pixels with local background subtraction. "
            "Control scatter is descriptive. An absolute line flux, completeness, and global significance "
            "require contamination models and calibrated injection tests.",
        ),
    )
    plot(records, example)


if __name__ == "__main__":
    main()

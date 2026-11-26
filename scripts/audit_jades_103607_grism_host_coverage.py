"""Map archival slitless wavelength coverage at the JADES 103607 transient."""

import json
from contextlib import ExitStack
from urllib.parse import urlsplit

import asdf
import numpy as np
import requests
from astropy.io import fits
from audit_jwst_november_f200w_spike import NativeSections
from download_jwst_cross_detector import OUT, ROOT
from download_jwst_transient_grism_geometry import CRDS, fetch
from matplotlib.figure import Figure
from measure_jwst_epoch_candidate import verified_path
from measure_jwst_transient_grism_geometry import (
    coverage,
    dispersed_trace,
    imaging_transform,
    validate_pointing,
)
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

STEM = "module_b_001_grism_host_coverage"
FILTERS = {"F356W": (3.18, 3.94), "F444W": (3.90, 4.98)}
LINES = {
    "[S III] 9069": 0.90711,
    "[S III] 9531": 0.95332,
    "He I 10830": 1.08333,
    "Pa gamma": 1.09411,
    "Pa beta": 1.28216,
}
PHOTOMETRIC_SCENARIOS = (2.2459, 2.263, 2.47, 2.4858)
CRDS_KEYS = ("R_DISTOR", "R_FILOFF", "R_SPCWCS")


def frame_geometry(cache, row, coordinates, calibration_cache):
    x1d = row["product"]
    original = verified_path(x1d)
    header = fits.getheader(original)
    if tuple(header[k] for k in ("FILTER", "DETECTOR", "PUPIL")) != (
        row["filter"],
        row["detector"],
        "GRISMR",
    ):
        raise ValueError("Archived extraction header conflicts with the receipt")
    reference_names = {k: urlsplit(header[k]).netloc for k in CRDS_KEYS}
    refs = {}
    for key, name in reference_names.items():
        if not name.startswith("jwst_nircam_") or not name.endswith(".asdf"):
            raise ValueError(f"Invalid calibration reference for {key}")
        if name not in calibration_cache:
            calibration_cache[name] = fetch(
                cache,
                CRDS + name,
                name,
                "calibration_reference",
                survey="JWST CRDS",
                version=header["CRDS_CTX"],
            )
        refs[key] = calibration_cache[name]
    name = original.name.replace("_x1d.fits", "_rate.fits")
    native = NativeSections(cache, name)
    primary, first = native.header(0)
    sci, _ = native.header(first)
    if (
        primary["FILENAME"] != name
        or sci["EXTNAME"] != "SCI"
        or primary["DETECTOR"] != row["detector"]
    ):
        raise ValueError("Identity mismatch in rate image")
    positions = {}
    with ExitStack() as stack:
        m = {k: stack.enter_context(asdf.open(verified_path(v))) for k, v in refs.items()}
        (offset,) = [
            f
            for f in m["R_FILOFF"]["filters"]
            if (f["filter"], f["pupil"]) == (row["filter"], "GRISMR")
        ]
        wcs = imaging_transform(sci, m["R_DISTOR"]["model"], offset)
        order = m["R_SPCWCS"]
        checks = validate_pointing(
            dict(x1d_path=x1d["path"], x1d_sha256=x1d["sha256"]), wcs, order, sci["VELOSYS"]
        )
        for label, sky in coordinates.items():
            x0, y0 = wcs.inverse(*sky)
            x, y, wave = dispersed_trace(order, 1, x0, y0, sci["VELOSYS"])
            lo, hi = FILTERS[row["filter"]]
            physical = coverage(wave, x, y, sci["NAXIS1"], sci["NAXIS2"], lo, hi)
            cases = {}
            for z in PHOTOMETRIC_SCENARIOS:
                values = {k: rest * (1 + z) for k, rest in LINES.items()}
                cases[str(z)] = {}
                for key, w in values.items():
                    intersects = bool(
                        physical["first_um"] is not None
                        and physical["first_um"] <= w <= physical["last_um"]
                    )
                    cases[str(z)][key] = dict(
                        wavelength_um=w,
                        intersects=intersects,
                        xy_pix=[float(np.interp(w, wave, coord)) for coord in (x, y)]
                        if intersects
                        else None,
                    )
            positions[label] = dict(
                direct_xy=[float(x0), float(y0)],
                status="EVALUATED",
                nominal_band_um=[lo, hi],
                trace=physical,
                conditional_lines=cases,
            )
    return dict(
        x1d=x1d,
        rate_model_headers=native.ranges,
        rate_etag=native.etag,
        rate_url=native.url,
        science_start_date_utc=primary["DATE-BEG"],
        file=name,
        visit=row["visit"],
        detector=row["detector"],
        filter=row["filter"],
        calibration_names=reference_names,
        x1d_geometry_control=checks,
        positions=positions,
    )


def plot(records):
    fig = Figure(figsize=(10, 6.2), layout="constrained")
    a = fig.subplots()
    for i, row in enumerate(records):
        trace = row.get("positions", {}).get("transient", {}).get("trace")
        if trace and trace["first_um"] is not None:
            color = "#2d697d" if row["filter"] == "F356W" else "#ab5839"
            a.plot(
                [trace["first_um"], trace["last_um"]],
                [i, i],
                color=color,
                lw=3,
                solid_capstyle="butt",
            )
    for i, (line, rest) in enumerate(LINES.items()):
        w = rest * 3.47
        if 3.15 < w < 5.02:
            a.axvline(w, color="0.45", lw=0.7, ls="--")
            a.annotate(
                line,
                (w, 1.01),
                xycoords=("data", "axes fraction"),
                ha="left" if i == 3 else "right" if i == 2 else "center",
                fontsize=8,
            )
    a.set(
        xlim=(3.13, 5.03),
        xlabel="Observed wavelength (µm)",
        ylabel="Archive exposure",
        yticks=range(len(records)),
        yticklabels=[
            f"{r['filter']} · {r['detector']} · {r['file'].split('_')[1]}-{r['file'].split('_')[2][-2:]}"
            for r in records
        ],
    )
    a.invert_yaxis()
    a.tick_params(axis="y", labelsize=8)
    fig.suptitle(
        "JADES 103607 · archived JWST grism coverage and possible host lines at z = 2.47",
        fontsize=13,
    )
    save_figure(fig, FIGURES / STEM)


def main():
    path = OUT / "module_b_001_archived_grism.json"
    data = json.loads(path.read_text())
    c = {"transient": data["icrs_deg"], "jades_103607": [53.1886733, -27.8125526]}
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    calibrations, results = {}, []
    for row in data["products"]:
        try:
            outcome = frame_geometry(cache, row, c, calibrations)
        except (requests.RequestException, KeyError, ValueError) as exc:
            sample = row["product"]["path"].split("/")[-1].replace("_x1d.fits", "_rate.fits")
            outcome = dict(
                file=sample,
                filter=row["filter"],
                detector=row["detector"],
                visit=row["visit"],
                status="UNAVAILABLE",
                reason=str(exc),
            )
        results.append(outcome)
        trace = outcome.get("positions", {}).get("transient", {})
        print(
            outcome["file"],
            trace.get("trace", trace.get("reason", outcome.get("reason"))),
            flush=True,
        )
    write_json(
        OUT / f"{STEM}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            x1d_inventory_sha256=text_digest(path),
            coordinates_icrs_deg=c,
            calibrations=calibrations,
            conditional_redshifts=list(PHOTOMETRIC_SCENARIOS),
            rest_lines_um=LINES,
            exposures=results,
            scope="Geometric first-order coverage only. It does not measure line flux, contamination, or sensitivity. "
            "All listed observations predate the brightening.",
        ),
    )
    plot(results)


if __name__ == "__main__":
    main()

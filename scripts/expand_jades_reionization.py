"""Expand the JADES Lyα screen to every object in the eligible DR4 roster."""

import json

import numpy as np
import pandas as pd
import requests
from check_jades_reionization_prism import DISP, NOMINAL, archive_spectrum, local_fit
from download_jades_reionization_sample import BIN_EDGES, LEDGER, OUT, ROOT, acquire
from inspect_jades_reionization import literature_matches
from matplotlib.figure import Figure
from measure_jades_reionization import LYA, WIDTHS, best_model, line_models, spectral_arrays
from scipy.stats import binomtest

from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.visual_evidence import save_figure

RESULTS = ROOT / "results/jades_reionization_full"
SIGNAL_THRESHOLD = 5.0
PRISM_THRESHOLD = 3.5


def collect(cache):
    old = json.loads(LEDGER.read_text())
    roster_path = OUT / "eligible_targets.csv"
    if text_digest(roster_path) != old["target_roster_sha256"]:
        raise ValueError("Eligible roster changed since the pilot")
    roster = pd.read_csv(roster_path).sort_values("dr4_id")
    prior = {r["dr4_id"]: r for r in old["records"]}
    if len(roster) != 204 or not roster.dr4_id.is_unique or not set(prior) <= set(roster.dr4_id):
        raise ValueError("Unexpected set of JADES eligible targets")
    path = RESULTS / "archive_sample.json"
    archive = dict(
        executed_utc=utcnow(),
        roster_sha256=text_digest(roster_path),
        pilot_sha256=text_digest(LEDGER),
        script_sha256=text_digest(__file__),
        planned=len(roster),
        records=[],
    )
    resumed = {}
    if path.is_file():
        state = json.loads(path.read_text())
        if (
            state["roster_sha256"] == archive["roster_sha256"]
            and state["pilot_sha256"] == archive["pilot_sha256"]
        ):
            resumed = {r["dr4_id"]: r for r in state["records"] if r["status"] == "ACCEPTED"}
    listings = {}
    for row in roster.to_dict(orient="records"):
        sid = row["dr4_id"]
        if sid in prior and prior[sid]["status"] == "ACCEPTED":
            result = dict(prior[sid], cohort="pilot")
        elif sid in resumed:
            result = resumed[sid]
        else:
            try:
                record = acquire(cache, listings, row)
            except (KeyError, OSError, ValueError, requests.RequestException) as exc:
                record = dict(status="NO_ACCEPTED_PRODUCT", reason=str(exc))
            result = dict(**row, **record, cohort="additional")
            print("Archive", sid, result["status"], flush=True)
        archive["records"].append(result)
        write_json(path, archive)
    return path, archive["records"]


def fits_for(obj):
    row = {key: obj[key] for key in ("dr4_id", "field", "bin", "z", "m_uv", "cohort", "status")}
    row.update(above_screen=False, covered_templates=0)
    if obj["status"] != "ACCEPTED":
        return row
    path = ROOT / obj["spectrum"]["path"]
    if digest(path) != obj["spectrum"]["sha256"]:
        raise ValueError(f"Input checksum mismatch for {path.name}")
    zero = LYA * (1 + obj["z"])
    spec = read_dr4_spectrum(path, wavelength_range=(0.65, 5.4))
    data = spectral_arrays(spec)
    models = line_models(data, zero, WIDTHS[0])
    row["covered_templates"] = len(models)
    fixed = next((m for m in models if m["offset_kms"] == 350), None)
    if fixed:
        row.update(reference_flux=fixed["flux"], reference_error=fixed["error"])
    m = best_model(models)
    if m:
        score = m["flux"] / m["error"]
        row.update(
            best_offset_kms=m["offset_kms"],
            best_flux=m["flux"],
            best_error=m["error"],
            best_snr=score,
            above_screen=score >= SIGNAL_THRESHOLD,
        )
        if row["above_screen"]:
            for tag, ap, nod in (("three_pixel", 3, 3), ("two_nod", 5, 2)):
                other = read_dr4_spectrum(path, aperture=ap, nods=nod, wavelength_range=(0.65, 5.4))
                variants = line_models(spectral_arrays(other), zero, WIDTHS[0])
                value = next(
                    (
                        v["flux"] / v["error"]
                        for v in variants
                        if v["offset_kms"] == m["offset_kms"]
                    ),
                    None,
                )
                row[tag + "_snr"] = value
    return row


def prism_screen(cache, rows, inventory):
    from astropy.io import fits

    path = cache.fetch(
        DISP,
        "metadata/jwst_nirspec_prism_disp.fits",
        survey="STScI",
        version="linked JDox curve",
        kind="instrument_reference",
    )
    with fits.open(path) as hdus:
        wave, resolving = (np.array(hdus[1].data[name], float) for name in ("WAVELENGTH", "R"))
    receipts, views = [], []
    for row in rows:
        sid = row["dr4_id"]
        entry = inventory[sid]
        try:
            spec, products = archive_spectrum(cache, entry)
        except (KeyError, OSError, ValueError, requests.RequestException) as exc:
            receipts.append(dict(dr4_id=sid, status="NO_ACCEPTED_PRISM", reason=str(exc)))
            continue
        z, offset = row["z"], row["best_offset_kms"]
        r = float(np.interp(LYA * (1 + z) / 1e4, wave, resolving))
        best, data = local_fit(spectral_arrays(spec), z, offset, r, *NOMINAL)
        receipts.append(
            dict(
                dr4_id=sid,
                status="MEASURED" if best else "INSUFFICIENT_COVERAGE",
                g140m_snr=row["best_snr"],
                two_nod_snr=row["two_nod_snr"],
                result=best,
                **products,
            )
        )
        if best is not None:
            views.append((row, best, data))
        print("PRISM", sid, None if best is None else round(best["snr"], 2), flush=True)
    return receipts, views


def draw(summary, views):
    strong = sorted(views, key=lambda item: item[1]["snr"], reverse=True)[:4]
    fig = Figure(figsize=(9.4, 3.6 + 2.5 * len(strong)), layout="constrained")
    axes = fig.subplots(len(strong) + 1, 1, squeeze=False).ravel()
    census = axes[0]
    for cohort, marker, color, skew in (
        ("pilot", "o", "#297485", -0.08),
        ("additional", "s", "#a75038", 0.08),
    ):
        rows = [s for s in summary if s["cohort"] == cohort]
        for j, row in enumerate(rows):
            total, hits = row["usable"], row["formal_candidates"]
            if not total:
                continue
            val = hits / total
            interval = binomtest(hits, total).proportion_ci(method="wilson")
            x = row["bin"] + skew
            census.errorbar(
                x,
                val,
                yerr=[[val - interval.low], [interval.high - val]],
                fmt=marker,
                color=color,
                label=cohort.capitalize() if j == 0 else None,
                ms=6,
                capsize=3,
            )
            census.annotate(
                f"{hits}/{total}",
                (x, val),
                xytext=(8, 3),
                textcoords="offset points",
                ha="left",
                fontsize=8,
                color=color,
            )
    census.set(
        xticks=range(len(BIN_EDGES) - 1),
        xticklabels=[
            f"{lo:g} ≤ z < {hi:g}" for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:], strict=True)
        ],
        ylabel="Fraction above line screen",
        ylim=(0, 1.05),
        title="JADES DR4 · complete eligible Lyα target list",
    )
    census.legend(frameon=False)
    for row, ax in zip(strong, axes[1:], strict=True):
        lead, fit, data = row
        w, y, e, good, model, baseline = data
        rest = w / (1 + lead["z"])
        ax.errorbar(rest[good], y[good], yerr=e[good], fmt="o", ms=3, capsize=2, color="#246f81")
        ax.step(rest, model, where="mid", color="#a75038", label="Break + fixed line")
        ax.step(rest, baseline, where="mid", color="#7c60a3", ls="--", label="Break only")
        ax.set(
            title=f"{lead['dr4_id']} · z={lead['z']:.3f} · PRISM S/N {fit['snr']:.1f}",
            xlabel="Rest-frame wavelength (Å)",
            ylabel="Fλ (10⁻²⁰ cgs/Å)",
        )
        ax.legend(frameon=False, fontsize=8)
    save_figure(fig, ROOT / "figures/jades_reionization_full_screen")


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    cache = DownloadCache(ROOT, 28_000 * 1024**2)
    archive, entries = collect(cache)
    measured = pd.DataFrame([fits_for(row) for row in entries])
    art = RESULTS / "measurements.csv"
    measured.to_csv(art, index=False)
    measured = measured.set_index("dr4_id", drop=False)
    reference, matches = literature_matches(cache, entries, measured)
    published = {r["dr4_id"] for r in matches if r["matches"]}
    follow = measured[
        (measured.cohort == "additional")
        & measured.above_screen
        & measured.two_nod_snr.ge(SIGNAL_THRESHOLD)
        & measured.three_pixel_snr.ge(SIGNAL_THRESHOLD)
        & ~measured.index.isin(published)
    ].sort_values("best_snr", ascending=False)
    checks, views = prism_screen(
        cache, follow.to_dict(orient="records"), {t["dr4_id"]: t for t in entries}
    )
    aggregates = []
    for (cohort, index), group in measured.groupby(["cohort", "bin"], sort=True):
        aggregates.append(
            dict(
                cohort=cohort,
                bin=int(index),
                drawn=len(group),
                spectra=int(group.status.eq("ACCEPTED").sum()),
                usable=int(group.covered_templates.gt(0).sum()),
                formal_candidates=int(group.above_screen.sum()),
                published_table_recoveries=int(
                    (group.above_screen & group.index.isin(published)).sum()
                ),
            )
        )
    draw(aggregates, views)
    write_json(
        RESULTS / "screen_summary.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            archive_sha256=text_digest(archive),
            measurement_sha256=text_digest(art),
            counts=aggregates,
            literature_table=reference,
            catalogue_matches=matches,
            prism_checks=checks,
            prism_priorities=[
                r["dr4_id"]
                for r in checks
                if r.get("result") and r["result"]["snr"] >= PRISM_THRESHOLD
            ],
            scope="Exact DR4 eligible roster of 204. Fixed templates and release-noise significance."
            " PRISM priorities need exposure, 2D trace and literature vetting.",
        ),
    )
    print(aggregates, flush=True)


if __name__ == "__main__":
    main()

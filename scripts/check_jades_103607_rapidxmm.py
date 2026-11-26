"""Retrieve XMM-Newton archival aperture photometry at the JADES 103607 transient."""

import argparse
import json
from datetime import datetime
from math import isfinite
from statistics import median

from download_jwst_cross_detector import OUT, ROOT
from inspect_jades_cosmic_evolution import receipt
from matplotlib.figure import Figure
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import DownloadCache, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

SOURCE = OUT / "module_b_001_native_november_fade.json"
TAG = "module_b_001_rapidxmm"
SERVER = "https://nxsa.esac.esa.int/nxsa-sl/servlet/get-uls"
BANDS = {"band6": "0.2-2", "band7": "2-12", "band8": "0.2-12"}
EVENT_WINDOW = ("2026-10-01", "2026-11-27")
OBSERVATION_TYPES = {"pointed": 15, "slew": 30}


def valid_limits(row, band):
    return row[f"{band}_flags"] == 0 and all(
        value is not None and isfinite(value) and value > 0
        for value in (row[f"{band}_exposure"], row[f"{band}_ul_sigma3"])
    )


def summarize(observations):
    result = {}
    for kind in OBSERVATION_TYPES:
        rows = [row for row in observations if row["obstype"] == kind]
        bands = {}
        for band, energy in BANDS.items():
            usable = [row for row in rows if valid_limits(row, band)]
            rate = f"{band}_ul_sigma3"
            bands[energy] = dict(
                usable_rows=len(usable),
                by_detector={
                    instrument: sensitivity(
                        [row for row in usable if row["instrum"] == instrument], band, rate
                    )
                    for instrument in sorted({row["instrum"] for row in usable})
                },
            )
        result[kind] = dict(
            instrument_rows=len(rows),
            distinct_observation_ids=len({row["obsid"] for row in rows}),
            start_dates_utc=[rows[0]["start_date"], rows[-1]["start_date"]] if rows else None,
            bands_kev=bands,
        )
    return result


def sensitivity(rows, band, rate):
    values = [row[rate] for row in rows]
    exposures = [row[f"{band}_exposure"] for row in rows]
    best = min(rows, key=lambda row: row[rate])
    return dict(
        rows=len(rows),
        exposure_s_range=[min(exposures), max(exposures)],
        upper_count_rate_range_per_s=[min(values), max(values)],
        median_upper_count_rate_per_s=median(values),
        most_sensitive_row=dict(
            obsid=best["obsid"],
            start_date=best["start_date"],
            source_counts=best[f"{band}_src_counts"],
            expected_background=best[f"{band}_bck_counts"],
            exposure_s=best[f"{band}_exposure"],
        ),
    )


def plot(rows):
    fig = Figure(figsize=(8.5, 5.3), layout="constrained")
    axis = fig.subplots()
    markers = {"PN": "o", "M1": "s", "M2": "^"}
    for kind, radius in OBSERVATION_TYPES.items():
        for instrument, marker in markers.items():
            group = [
                row
                for row in rows
                if row["obstype"] == kind
                and row["instrum"] == instrument
                and valid_limits(row, "band6")
            ]
            if not group:
                continue
            axis.scatter(
                [datetime.fromisoformat(row["start_date"]) for row in group],
                [row["band6_ul_sigma3"] for row in group],
                marker=marker,
                s=23,
                alpha=0.7,
                label=f"{kind} {instrument} ({radius} arcsec)",
            )
    axis.axvspan(
        *[datetime.fromisoformat(day) for day in EVENT_WINDOW],
        color="0.85",
        label="2026 JWST search window",
    )
    axis.set(
        yscale="log",
        xlabel="XMM-Newton observation start (UTC)",
        ylabel="Soft-band 99.87% upper count rate (counts/s)",
        title="JADES 103607 position | XMM-Newton 0.2-2 keV archival coverage",
    )
    axis.grid(axis="y", alpha=0.2, which="both")
    axis.legend(frameon=False, fontsize=8, loc="center right")
    save_figure(fig, FIGURES / TAG)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="Query the live archive again.")
    args = parser.parse_args()
    ra, dec = json.loads(SOURCE.read_text())["icrs_deg"]
    url = f"{SERVER}?ra={ra:.10f}&dec={dec:.10f}"
    cache = DownloadCache(ROOT, 32_000 * 1024**2)
    path = cache.fetch(
        url,
        f"metadata/{TAG}.json",
        survey="XMM-Newton RapidXMM",
        version="ESA positional upper limit server",
        kind="observed_archive_aperture_photometry",
        refresh=args.refresh,
    )
    rows = json.loads(path.read_text())
    observation_keys = [(row["obstype"], row["obsid"], row["instrum"]) for row in rows]
    if not rows or len(set(observation_keys)) != len(rows):
        raise ValueError(
            "Expected nonempty RapidXMM results with one row per observation and detector"
        )
    rows.sort(key=lambda row: (row["start_date"], row["obsid"], row["instrum"]))
    start, stop = EVENT_WINDOW
    overlapping = sorted(
        {row["obsid"] for row in rows if row["start_date"] < stop and row["end_date"] >= start}
    )
    plot(rows)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=[ra, dec],
            position_sha256=text_digest(SOURCE),
            archive=receipt(path, url),
            aperture_radii_arcsec=OBSERVATION_TYPES,
            upper_posterior_probability=0.9987,
            event_window_utc_end_exclusive=EVENT_WINDOW,
            available_observations_overlapping_event_window=overlapping,
            summary=summarize(rows),
            observations=rows,
            scope="RapidXMM upper limits use nearby HEALPix cell centers and broad apertures. "
            "Limits are in detector counts per second. Aperture events could include nearby "
            "sources. Observation dates are scheduling bounds, not continuous target exposure. "
            "An empty time-window match describes this archive response only.",
        ),
    )


if __name__ == "__main__":
    main()

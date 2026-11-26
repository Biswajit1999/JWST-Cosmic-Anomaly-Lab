"""Summarize the fixed-position search across all 7 Ms CDF-S observing intervals."""

import json

import numpy as np
from astropy.time import Time
from measure_jades_103607_chandra_events import CAMPAIGNS, OUT, SOURCE, search_observations

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json

TAG = "module_b_001_chandra_full_history"


def main():
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    campaign_rows, best_trials, results = [], [], {}
    covered, short_trials = [], []
    for year in CAMPAIGNS:
        path = OUT / f"module_b_001_chandra_{year}_events.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        rows = data["observations"]
        expected, _ = search_observations(year)
        if (
            data["position_icrs_deg"] != source["icrs_deg"]
            or data["position_sha256"] != text_digest(SOURCE)
            or [row["journal_id"] for row in rows] != [row["journal_id"] for row in expected]
        ):
            raise ValueError(f"Incomplete or inconsistent CDF-S campaign {year}")
        usable = [row for row in rows if row["status"] == "measured"]
        intervals = [window for row in usable for window in row["windows"]]
        if data["number_of_tests"] != len(usable) + len(intervals):
            raise ValueError(f"Weakened trial accounting for campaign {year}")
        results[str(year)] = dict(path=path.name, sha256=text_digest(path))
        campaign_rows.append(
            dict(
                year=year,
                journal_intervals=len(rows),
                usable_intervals=len(usable),
                short_windows=len(intervals),
                full_acis_good_time_s=sum(row["observation"]["good_time_s"] for row in usable),
                short_window_good_time_s=sum(window["good_time_s"] for window in intervals),
                excluded_journal_ids=[
                    row["journal_id"] for row in rows if row["status"] != "measured"
                ],
                individual_photon_count=sum(row["observation"]["source_counts"] for row in usable),
            )
        )
        best_trials.append(dict(year=year, **data["best_trial"]))
        covered.extend(usable)
        short_trials.extend(intervals)

    if sum(row["journal_intervals"] for row in campaign_rows) != 102:
        raise ValueError("Expected 102 intervals in the published 7 Ms journal")
    number_of_tests = len(covered) + len(short_trials)
    best = min(best_trials, key=lambda row: row["binomial_background_tail"])
    dates = [min(row["start_utc"] for row in covered), max(row["end_utc"] for row in covered)]
    start, stop = Time(dates, format="isot", scale="utc")
    if any(a["end_utc"] > b["start_utc"] for a, b in zip(covered, covered[1:])):
        raise ValueError("Observing intervals overlap or are out of order")
    live = sum(row["short_window_good_time_s"] for row in campaign_rows)
    short_limits = [w["upper95_aperture_counts"] / w["good_time_s"] for w in short_trials]
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            position_icrs_deg=source["icrs_deg"],
            position_sha256=text_digest(SOURCE),
            inputs=results,
            campaigns=campaign_rows,
            journal_intervals_examined=102,
            usable_intervals=len(covered),
            combined_number_of_tests=number_of_tests,
            best_trial=best,
            best_bonferroni_p=min(1.0, best["binomial_background_tail"] * number_of_tests),
            measured_time_span_utc=dates,
            short_window_good_time_s=live,
            short_window_fraction_of_calendar=float(live / (stop - start).sec),
            short_window_95_aperture_count_rate_quantiles_per_s=dict(
                zip(
                    ("minimum", "median", "90_percentile", "maximum"),
                    map(float, np.quantile(short_limits, [0, 0.5, 0.9, 1])),
                    strict=True,
                )
            ),
            scope="Summed chip good time across nonoverlapping pointings for included windows. "
            "The calendar fraction describes sampled live time. It is not a flare recovery or "
            "recurrence probability. Count-rate bounds retain background, exclude response "
            "and aperture-fraction corrections and do not apply to short events diluted in a "
            "window. The p-value corrects for all full-interval and 10 ks tests at one position.",
        ),
    )


if __name__ == "__main__":
    main()

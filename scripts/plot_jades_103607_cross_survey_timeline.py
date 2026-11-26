"""Plot the actual date coverage and 2026 archival gaps at JADES 103607."""

import json
from datetime import datetime

from astropy.time import Time
from download_jwst_cross_detector import OUT
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.figure import Figure
from measure_jades_103607_hst_history import epoch_start
from vet_jwst_field_b_candidate import FIGURES

from cosmic_anomaly_lab.provenance import text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

PREFIX = "module_b_001_"
TAG = PREFIX + "cross_survey_timeline"
KEYS = (
    "hst_acs_goods_epochs",
    "hst_f160w_history",
    "spitzer_epochs",
    "vista_video_coverage",
    "chandra_full_history",
    "rapidxmm",
    "meerkat",
    "euclid_dr1_h_history",
    "native_november_fade",
    "independent_2026_coverage",
    "decam_2026_coverage",
    "rubin_2026_alerts",
)


def read_all():
    data, provenance = {}, {}
    for key in KEYS:
        path = OUT / f"{PREFIX}{key}.json"
        data[key] = json.loads(path.read_text(encoding="utf-8"))
        provenance[path.name] = text_digest(path)
    return data, provenance


def day_list(values):
    return sorted({str(value)[:10] for value in values})


def mjd_utc(value):
    return Time(value, format="mjd", scale="utc").isot


def collect(data, hashes):
    full_chandra = data["chandra_full_history"]
    chandra = []
    for meta in full_chandra["inputs"].values():
        path = OUT / meta["path"]
        if text_digest(path) != meta["sha256"]:
            raise ValueError("Chandra campaign changed after its combined audit")
        hashes[path.name] = meta["sha256"]
        chandra.extend(json.loads(path.read_text(encoding="utf-8"))["observations"])
    # The short Spitzer bands were observed in these independent epoch windows.
    # Exact input exposures are reported in the linked mosaic measurement.
    irac_windows = ("2004-08", "2005-08", "2006-02", "2011-09", "2012-03")
    radio = data["meerkat"]["contributing_date_range_utc"]
    history = {
        "Chandra ACIS": day_list(r["start_utc"] for r in chandra if r["status"] == "measured"),
        "XMM-Newton EPIC": day_list(
            r["start_date"] for r in data["rapidxmm"]["observations"] if r["band6_flags"] == 0
        ),
        "HST ACS F850LP": day_list(
            r["selected"]["first_input_date"] for r in data["hst_acs_goods_epochs"]["epochs"]
        ),
        "Spitzer IRAC": list(irac_windows),
        "HST CANDELS": day_list(
            epoch_start(r) for r in data["hst_f160w_history"]["coverage"] if r["center_weight"] > 0
        ),
        "VISTA VIDEO H-band pilot": day_list(
            r["start_utc"]
            for r in data["vista_video_coverage"]["h_band_pilot"]
            if r["status"] == "measured"
        ),
        "MeerKAT MIGHTEE coadd": radio,
        "Euclid DR1 H band": day_list(
            r["start_date_utc"]
            for r in data["euclid_dr1_h_history"]["visit_inventory"]
            if r.get("delta") is not None
        ),
    }
    time_window = data["independent_2026_coverage"]
    jwst_epochs = [
        mjd_utc(product["expstart_mjd"])
        for band in data["native_november_fade"]["bands"].values()
        for visit in band["inventory"]
        for product in visit["products"]
    ]
    dates = {
        "JWST NIRCam measured": day_list(
            day for day in jwst_epochs if "2026-10-01" <= day < "2026-11-27"
        ),
        "DECam covers target · unreleased": data["decam_2026_coverage"][
            "containing_observation_dates_utc"
        ],
        "Rubin alerts near target": day_list(
            r[label]
            for r in data["rubin_2026_alerts"]["surrounding_field"]["objects"]
            for label in ("first_alert_utc", "latest_alert_utc")
        ),
        "HST archives · footprints miss": day_list(
            mjd_utc(r["start_mjd"]) for r in time_window["mast_hst"] if r["contains"] is False
        ),
    }
    return history, dates


def write_timeline(history, recent):
    fig = Figure(figsize=(10.8, 6.6), layout="constrained")
    axes = fig.subplots(2, 1, gridspec_kw={"height_ratios": (1.65, 1)})
    for y, (name, values) in enumerate(history.items()):
        times = [
            datetime.fromisoformat(value if len(value) == 10 else f"{value}-15") for value in values
        ]
        if name == "MeerKAT MIGHTEE coadd":
            axes[0].plot(times, [y] * 2, "o--", color="#226e7f", ms=3, lw=1.15)
        else:
            axes[0].scatter(times, [y] * len(times), marker="|", s=100, color="#226e7f")
    axes[0].set_yticks(range(len(history)), list(history), fontsize=9)
    axes[0].set(
        xlim=(datetime(1999, 1, 1), datetime(2026, 7, 1)),
        title="Earlier observations at the position",
    )
    axes[0].xaxis.set_major_locator(YearLocator(4))
    axes[0].xaxis.set_major_formatter(DateFormatter("%Y"))
    axes[0].invert_yaxis()
    style = (("o", "#17687c"), ("D", "#9f522c"), ("^", "#6d597c"), ("x", "#585858"))
    for y, ((_, values), (mark, color)) in enumerate(zip(recent.items(), style, strict=True)):
        axes[1].scatter(
            [datetime.fromisoformat(v) for v in values],
            [y] * len(values),
            marker=mark,
            s=30,
            color=color,
        )
    axes[1].set_yticks(range(len(recent)), list(recent), fontsize=9)
    axes[1].set(
        xlim=(datetime(2026, 9, 30), datetime(2026, 11, 27)),
        title="2026 coverage and alert records",
        xlabel="UTC observation date",
    )
    axes[1].xaxis.set_major_locator(MonthLocator(bymonthday=(1, 9, 17, 25)))
    axes[1].xaxis.set_major_formatter(DateFormatter("%d %b"))
    axes[1].invert_yaxis()
    for axis in axes:
        axis.margins(y=0.18)
        axis.grid(axis="x", alpha=0.19)
    fig.suptitle("JADES 103607 transient | historical and contemporaneous records", fontsize=15)
    save_figure(fig, FIGURES / TAG)


def main():
    data, hashes = read_all()
    history, recent = collect(data, hashes)
    write_timeline(history, recent)
    write_json(
        OUT / f"{TAG}.json",
        dict(
            executed_utc=utcnow(),
            script_sha256=text_digest(__file__),
            input_result_sha256=hashes,
            historical_records=history,
            contemporaneous_records=recent,
            scope="Observation-start days or published Spitzer months. MeerKAT endpoints enclose "
            "the contributing observations in a single published coadd. Marks do not imply "
            "continuous exposure, equal depth or common flux limits. The 2026 lower panel includes "
            "public metadata for unreleased DECam images, nearby Rubin alerts, and HST products "
            "whose footprints miss the position.",
        ),
    )


if __name__ == "__main__":
    main()

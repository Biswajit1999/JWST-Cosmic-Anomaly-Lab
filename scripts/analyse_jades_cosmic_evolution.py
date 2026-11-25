"""Compare JWST oxygen-to-hydrogen line ratios at different cosmic ages."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
from astropy.io import fits
from matplotlib.figure import Figure
from scipy.stats import binomtest

from cosmic_anomaly_lab.balmer_census import LINES, normalize_filter
from cosmic_anomaly_lab.provenance import DownloadCache, digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.visual_evidence import save_figure

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jades_cosmic_evolution"
CATALOGUE = "https://jades.herts.ac.uk/DR4/Combined_DR4_external_v1.2.1.fits"
MANIFEST = ROOT / "data/manifests/downloads.json"
BIN_EDGES = (2.8, 3.5, 4.3, 5.2)
APERTURES = (5, 3)
FIELD_NAMES = {"GN": "GOODS-N", "GS": "GOODS-S", "both": "Both fields"}
GRATING = "g235m_f170lp"
FLUX_UNIT = 1e-20 * u.erg / u.s / u.cm**2
TRIALS = 4000
SEED = 1727
CASES = (
    ("primary", 5, 5, False),
    ("lower_snr", 5, 3, False),
    ("higher_snr", 5, 8, False),
    ("three_pixel", 3, 5, False),
    ("common_uv_range", 5, 5, True),
    ("grating_interior", 5, 5, False),
)


def load_catalogue(path):
    with fits.open(path) as hdus:
        data = hdus["Obs_info"].data
        columns = {
            "dr4_id": ("Unique_ID", str), "field": ("Field", str), "pid": ("PID", int),
            "tier": ("TIER", str), "nircam_id": ("NIRCam_DR5_ID", int),
            "ra": ("RA_TARG", float), "dec": ("Dec_TARG", float),
            "z": ("z_Spec", float), "flag": ("z_Spec_flag", str),
            "m_uv": ("MUV", float), "exposure_s": ("tExp_G235M", float),
        }
        obs = pd.DataFrame({key: np.asarray(data[col], dtype=typ) for key, (col, typ) in columns.items()})
        if not obs.dr4_id.is_unique or set(obs.field) != {"GN", "GS"}:
            raise ValueError("Unexpected DR4 target identifiers or field labels")
        for aperture in APERTURES:
            table = hdus[f"R1000_{aperture}pix"]
            row_index = pd.Index(np.asarray(table.data["Unique_ID"], dtype=str))
            order = row_index.get_indexer(obs.dr4_id)
            if len(row_index) != len(obs) or np.any(order < 0):
                raise ValueError("Line catalogue does not match observation IDs")
            for line, tag in (("O3_5007", "o3"), (LINES["hb"], "hb")):
                for field in ("flux", "flux_err", "SNR", "filter"):
                    col = f"{line}_{field}"
                    raw = table.data[col][order]
                    if field == "filter":
                        values = np.array([normalize_filter(item) or "" for item in raw])
                    else:
                        values = np.asarray(raw, dtype=float)
                        if field != "SNR":
                            values *= u.Unit(table.columns[col].unit).to(FLUX_UNIT.unit) / FLUX_UNIT.value
                    obs[f"{tag}_{field}_{aperture}"] = values
    return obs


def select_sources(obs):
    good = obs.flag.isin(("A", "B", "C")) & np.isfinite(obs.z) & (obs.z > 0)
    eligible = obs[good & (obs.nircam_id > 0)].copy()
    eligible["quality_order"] = eligible.flag.map({"A": 0, "B": 1, "C": 2})
    grouped = eligible.groupby(["field", "nircam_id"], sort=False)
    excluded = []
    duplicate_count = 0
    for key, group in grouped:
        if len(group) < 2:
            continue
        duplicate_count += 1
        sky = SkyCoord(group.ra, group.dec, unit="deg")
        spread = sky[:, None].separation(sky).arcsec.max()
        z_span = (group.z.max() - group.z.min()) / (1 + group.z.median())
        if spread > .25 or z_span > .01:
            excluded.append(dict(field=key[0], nircam_id=key[1], dr4_ids=group.dr4_id.tolist(),
                                 max_separation_arcsec=float(spread), fractional_z_span=float(z_span)))
    exclude_keys = {(x["field"], x["nircam_id"]) for x in excluded}
    if exclude_keys:
        eligible = eligible[~pd.MultiIndex.from_frame(eligible[["field", "nircam_id"]]).isin(exclude_keys)]
    exposure = eligible.exposure_s.fillna(0)
    ranked = eligible.assign(has_grating=exposure > 0, exposure_sort=exposure).sort_values(
        ["has_grating", "quality_order", "exposure_sort", "dr4_id"], ascending=[False, True, False, True]
    )
    sources = ranked.drop_duplicates(["field", "nircam_id"]).drop(
        columns=["has_grating", "quality_order", "exposure_sort"]
    )
    audit = dict(input_records=len(obs), quality_flags=obs.flag.value_counts().to_dict(),
                 secure_records=int(good.sum()), secure_records_without_dr5_id=int((good & (obs.nircam_id <= 0)).sum()),
                 repeat_dr5_groups=duplicate_count, excluded_discordant_groups=excluded,
                 retained_unique_sources=len(sources))
    return sources, audit


def measured(frame, aperture):
    covered = np.ones(len(frame), dtype=bool)
    snr = np.full(len(frame), np.inf)
    for line in ("o3", "hb"):
        flux = frame[f"{line}_flux_{aperture}"].to_numpy()
        err = frame[f"{line}_flux_err_{aperture}"].to_numpy()
        covered &= np.isfinite(flux) & np.isfinite(err) & (err > 0)
        covered &= frame[f"{line}_filter_{aperture}"].eq(GRATING).to_numpy()
        ratio = np.divide(flux, err, out=np.full(len(frame), np.nan), where=err > 0)
        snr = np.minimum(snr, np.minimum(ratio, frame[f"{line}_SNR_{aperture}"].to_numpy()))
    return covered, snr


def boot_median(values, rng):
    if len(values) < 5:
        return None, None
    sample = values[rng.integers(len(values), size=(TRIALS, len(values)))]
    draw = np.median(sample, axis=1)
    return np.percentile(draw, [2.5, 97.5]).tolist(), draw


def summarise(baseline):
    rng = np.random.default_rng(SEED)
    summaries, contrasts, examples = [], [], []
    for label, aperture, limit, uv_cut in CASES:
        case = baseline[baseline.m_uv.between(-21, -19)] if uv_cut else baseline
        if label == "grating_interior":
            case = case[case.z.between(2.90, 4.70)]
        for field in FIELD_NAMES:
            draws = {}
            for i, (lower, upper) in enumerate(zip(BIN_EDGES[:-1], BIN_EDGES[1:], strict=True)):
                rows = case[(case.z >= lower) & (case.z < upper)]
                if field != "both":
                    rows = rows[rows.field == field]
                observed = rows[rows.exposure_s > 0]
                fits_available, snr = measured(observed, aperture)
                detected = observed[fits_available & (snr >= limit)]
                values = np.log10(detected[f"o3_flux_{aperture}"] / detected[f"hb_flux_{aperture}"]).to_numpy()
                ci, draws[i] = boot_median(values, rng)
                median = float(np.median(values)) if len(values) else None
                if len(observed):
                    interval = binomtest(len(detected), len(observed)).proportion_ci(method="wilson")
                    detection = dict(fraction=len(detected) / len(observed), interval=[interval.low, interval.high])
                else:
                    detection = dict(fraction=None, interval=None)
                label_pairs = observed[[f"hb_filter_{aperture}", f"o3_filter_{aperture}"]].value_counts()
                summaries.append(dict(
                    case=label, field=field, bin=i, redshift=[lower, upper],
                    cosmic_age_gyr=Planck18.age([upper, lower]).to_value(u.Gyr).tolist(),
                    midpoint_age_gyr=float(Planck18.age((lower + upper) / 2).to_value(u.Gyr)),
                    secure_dr5_sources=len(rows), exposed_g235m=len(observed),
                    fits_available_same_grating=int(fits_available.sum()), two_line_detections=len(detected),
                    catalogue_grating_labels=[dict(hb=hb or None, o3=o3 or None, count=int(n))
                                              for (hb, o3), n in label_pairs.items()],
                    detection_fraction_in_exposed=detection, median_log10_ratio=median, bootstrap95=ci,
                    detected_median_m_uv=float(detected.m_uv.dropna().median()) if detected.m_uv.notna().any() else None,
                    detected_median_exposure_s=float(detected.exposure_s.median()) if len(detected) else None,
                ))
                if label == "primary" and field != "both" and i in (0, len(BIN_EDGES)-2) and len(values):
                    row = detected.iloc[np.argmin(abs(values - median))]
                    examples.append(dict(field=field, bin=i, dr4_id=row.dr4_id,
                                         nircam_id=int(row.nircam_id), pid=int(row.pid), z=float(row.z),
                                         ra=float(row.ra), dec=float(row.dec),
                                         log10_ratio=float(np.log10(row.o3_flux_5 / row.hb_flux_5))))
            lo, hi = draws[0], draws[len(BIN_EDGES)-2]
            result = dict(case=label, field=field, higher_minus_lower_redshift_dex=None, bootstrap95=None)
            if lo is not None and hi is not None:
                ends = [row for row in summaries if row["case"] == label and row["field"] == field]
                result.update(higher_minus_lower_redshift_dex=ends[-1]["median_log10_ratio"]-ends[0]["median_log10_ratio"],
                              bootstrap95=np.percentile(hi-lo, [2.5, 97.5]).tolist())
            contrasts.append(result)
    return summaries, contrasts, examples


def plot(summaries):
    fig = Figure(figsize=(9.4, 6.4), layout="constrained")
    ax, fraction = fig.subplots(2, 1, sharex=True)
    style = {"GN": ("#007089", "o", -.02), "GS": ("#ac4d35", "s", .02), "both": ("#242a34", "^", 0)}
    for field, (color, marker, offset) in style.items():
        rows = [s for s in summaries if s["field"] == field and s["case"] == "primary"]
        for j, row in enumerate(rows):
            t = row["midpoint_age_gyr"] + offset
            med, ci = row["median_log10_ratio"], row["bootstrap95"]
            d = row["detection_fraction_in_exposed"]
            if ci is not None:
                ax.errorbar(t, med, yerr=np.maximum(0, [[med-ci[0]], [ci[1]-med]]),
                            fmt=marker, color=color, capsize=3, label=FIELD_NAMES[field] if j == 0 else None)
            if d["interval"] is not None:
                p = d["fraction"]
                lo, hi = d["interval"]
                fraction.errorbar(t, p, yerr=[[p-lo], [hi-p]], fmt=marker, color=color, capsize=3)
    ax.set(ylabel="Median log₁₀([O III] 5007 / Hβ)", title="JADES DR4 · same-grating line ratios and detection fractions")
    ax.legend(frameon=False, ncols=3)
    fraction.set(xlabel="Age of the universe at observation (Gyr, Planck18)",
                 ylabel="Fraction with both lines ≥5 S/N", ylim=(0, 1))
    for axis in (ax, fraction):
        axis.grid(axis="y", linewidth=.6, alpha=.35)
    rows = [s for s in summaries if s["field"] == "both" and s["case"] == "primary"]
    redshift_axis = ax.secondary_xaxis("top")
    redshift_axis.set_xticks([s["midpoint_age_gyr"] for s in rows],
                             [f"z={sum(s['redshift'])/2:.2f}" for s in rows])
    save_figure(fig, ROOT / "figures/jades_cosmic_evolution")


def main():
    cache = DownloadCache(ROOT, 28_000*1024**2)
    path = cache.fetch(CATALOGUE, "jades_dr4.fits", survey="JADES DR4", version="v1.2.1", kind="observed_catalogue")
    inputs = json.loads(MANIFEST.read_text())
    receipt = next(row for row in reversed(inputs) if row["url"] == CATALOGUE)
    if receipt["sha256"] != digest(path):
        raise ValueError("Catalogue fingerprint mismatch")
    sources, audit = select_sources(load_catalogue(path))
    baseline = sources[sources.z.between(BIN_EDGES[0], BIN_EDGES[-1], inclusive="left")].copy()
    baseline["cosmic_age_gyr"] = Planck18.age(baseline.z.to_numpy()).to_value(u.Gyr)
    summaries, contrasts, examples = summarise(baseline)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "observed_targets.csv"
    baseline.to_csv(path, index=False)
    plot(summaries)
    times = (2.8, 3.5, 4.3, 5.2, 6, 8, 10, 14)
    write_json(OUT / "summary.json", dict(
        executed_utc=utcnow(), catalogue=receipt, audit=audit, target_table_sha256=text_digest(path),
        script_sha256=text_digest(__file__), bootstrap_draws=TRIALS, seed=SEED, grating_filter=GRATING,
        assumed_cosmology=dict(name=Planck18.name, H0_km_s_Mpc=Planck18.H0.value, Omega_m=Planck18.Om0),
        time_baseline=[dict(z=z, age_gyr=float(Planck18.age(z).to_value(u.Gyr)),
                            lookback_gyr=float(Planck18.lookback_time(z).to_value(u.Gyr))) for z in times],
        coverage_definition="Both lines have finite flux/positive errors and nominal G235M labels. Denominator is G235M exposure.",
        sample_definition="Unique secure DR5 counterparts. Repeat groups checked for coordinate and redshift agreement.",
        sensitivity_definition=dict(common_uv=(-21, -19), nominal_nonoverlap_z=(2.90, 4.70)),
        bins=summaries, endpoint_contrasts=contrasts, example_targets=examples,
    ))


if __name__ == "__main__":
    main()

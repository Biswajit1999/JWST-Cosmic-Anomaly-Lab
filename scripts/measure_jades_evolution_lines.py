"""Forced line measurements and artificial-line recovery in released JADES spectra."""

import json

import numpy as np
import pandas as pd
from analyse_jades_cosmic_evolution import BIN_EDGES, FIELD_NAMES, OUT, ROOT, boot_median
from analyse_jades_cosmic_evolution import measured as catalogue_coverage
from astropy import units as u
from astropy.cosmology import Planck18
from download_jades_evolution_sample import LEDGER
from inspect_jades_cosmic_evolution import DENSITY_UNIT, LINES_VACUUM_A
from matplotlib.figure import Figure
from scipy.special import ndtr

from cosmic_anomaly_lab.balmer_census import signed_core_integral
from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import pixel_edges
from cosmic_anomaly_lab.visual_evidence import save_figure

SPEED = 299792.458
HALF_WIDTH = 600
MASK_WIDTH = 800
MIN_SIDE_PIXELS = 2
SIGMAS_KMS = (125, 250)
FLUXES = (0, 10, 30, 100, 300)
LUMINOSITIES = (0, 1, 3, 10, 30)
DETECTION_SN = 5
TARGET_LINES = {"hb": LINES_VACUUM_A["Hβ"], "o3": LINES_VACUUM_A["[O III] 5007"]}
NEARBY_LINES_A = np.array([4687.02, 4712.69, 4741.49, 4862.64, 4923.31, 4960.29, 5008.24, 5017.08, 5199.35, 5201.70])
OFFSETS = np.r_[np.arange(-9600, -3999, 800), np.arange(4000, 9601, 800)]
RNG_SEED = 185191
STEM = "jades_evolution_forced_lines"


def line_window(wave, good, center, neighbours):
    v = (wave/center-1)*SPEED
    valid = good.copy()
    for near in neighbours:
        valid &= abs((wave/near-1)*SPEED) > MASK_WIDTH
    return v, valid


def measure(v, flux, err, valid, center):
    blue = valid & (v >= -3000) & (v <= -1500)
    red = valid & (v >= 2000) & (v <= 3000)
    if min(blue.sum(), red.sum()) < MIN_SIDE_PIXELS:
        return dict(status="MISSING_CONTINUUM_SIDE", flux=None, error=None)
    return signed_core_integral(v, flux, err, valid, center, HALF_WIDTH)


def injections(wave, flux, err, good, rest, line, redshift, dr4_id, field, index):
    edges = pixel_edges(wave)
    actual = NEARBY_LINES_A*(1+redshift)
    for offset in OFFSETS:
        center = rest*(1+redshift)*(1+offset/SPEED)
        if np.min(abs((actual/center-1)*SPEED)) < 1600:
            continue
        neighbours = actual[abs(actual/center-1)*SPEED < 4200]
        velocity, valid = line_window(wave, good, center, neighbours)
        base = measure(velocity, flux, err, valid, center)
        if base["flux"] is None:
            continue
        for sigma in SIGMAS_KMS:
            sigma_a = center*sigma/SPEED
            unit_line = np.diff(ndtr((edges-center)/sigma_a))/np.diff(edges)
            gain = measure(velocity, unit_line, err, valid, center)["flux"]
            yield dict(dr4_id=dr4_id, field=field, bin=index, redshift=redshift, line=line, offset_kms=int(offset),
                       wavelength_a=center, sigma_kms=sigma, base_flux=base["flux"], flux_error=base["error"],
                       response=gain)


def measure_target(entry):
    receipt = entry["spectrum"]
    path = ROOT / receipt["path"]
    if digest(path) != receipt["sha256"]:
        raise ValueError(f"Archive fingerprint mismatch for {path.name}")
    spec = read_dr4_spectrum(path, aperture=5, nods=3)
    wave = spec.wavelength.to_value(u.AA)
    density = u.spectral_density(spec.wavelength)
    flux = spec.flux.to_value(DENSITY_UNIT, equivalencies=density)/1e-20
    err = spec.error.to_value(DENSITY_UNIT, equivalencies=density)/1e-20
    z, sid, field, index = (entry[key] for key in ("z", "dr4_id", "field", "redshift_bin"))
    results, trials = {}, []
    for line, rest in TARGET_LINES.items():
        center = rest*(1+z)
        near = [LINES_VACUUM_A["[O III] 4959"]*(1+z)] if line == "o3" else []
        velocity, valid = line_window(wave, spec.good, center, near)
        results[line] = measure(velocity, flux, err, valid, center)
        if entry["role"] == "random_sample":
            trials.extend(injections(wave, flux, err, spec.good, rest, line, z, sid, field, index))
    if entry["role"] == "sideband_control":
        center = TARGET_LINES["o3"]*(1+z)
        velocity, _ = line_window(wave, spec.good, center, [])
        results["o3_uncorrected"] = signed_core_integral(velocity, flux, err, spec.good, center, HALF_WIDTH)
    return results, trials, (wave, flux, err, spec.good)


def signed_rows(entry, lines):
    row = {k: entry[k] for k in ("dr4_id", "field", "redshift_bin", "role", "z", "status")}
    for name, line in lines.items():
        flux, error = (line[k] for k in ("flux", "error"))
        row.update({f"{name}_status": line["status"], f"{name}_flux": flux, f"{name}_error": error,
                    f"{name}_snr": flux/error if flux is not None else None})
        if name == "hb":
            row["hb_formal_upper_3sigma"] = max(flux, 0)+3*error if flux is not None and flux/error < 3 else None
    return row


def statistics(measured, injection_data, catalogue):
    rng = np.random.default_rng(RNG_SEED)
    overview, contrast, recovery = [], [], []
    covered, snr = catalogue_coverage(catalogue, 5)
    detected = catalogue[covered & (snr >= DETECTION_SN)].set_index("dr4_id")
    catalogue_log_ratio = np.log10(detected.o3_flux_5/detected.hb_flux_5)
    measured = measured[measured.role == "random_sample"].copy()
    measured["catalogue_log_o3_hb"] = measured.dr4_id.map(catalogue_log_ratio)
    for field in FIELD_NAMES:
        subset = measured if field == "both" else measured[measured.field == field]
        metrics = {"dual_log_o3_hb": {}, "oxygen_selected_signed_hb_o3": {}}
        for index in range(len(BIN_EDGES)-1):
            rows = subset[subset.redshift_bin == index]
            available = rows[rows.hb_flux.notna() & rows.o3_flux.notna()]
            oxygen = available[available.o3_snr >= DETECTION_SN]
            joint = oxygen[oxygen.hb_snr >= DETECTION_SN]
            common = joint.dropna(subset=["catalogue_log_o3_hb"])
            paired_change = np.log10(common.o3_flux/common.hb_flux)-common.catalogue_log_o3_hb
            paired_ci, _ = boot_median(paired_change.to_numpy(), rng)
            row = dict(field=field, bin=index, selected=len(rows), spectral_products=int(rows.status.eq("ACCEPTED").sum()),
                       both_windows_measured=len(available), oxygen_5sigma=len(oxygen), both_5sigma=len(joint),
                       oxygen_only_hb_below_5sigma=len(oxygen)-len(joint),
                       oxygen_selected_hb_negative=int((oxygen.hb_flux < 0).sum()),
                       oxygen_selected_hb_below_3sigma=int((oxygen.hb_snr < 3).sum()),
                       catalogue_dual_detections_in_drawn_rows=int(rows.catalogue_log_o3_hb.notna().sum()),
                       paired_forced_minus_catalogue_dex=dict(n=len(common), median=float(paired_change.median())
                                                              if len(common) else None, bootstrap95=paired_ci))
            values = {"dual_log_o3_hb": np.log10(joint.o3_flux/joint.hb_flux),
                      "oxygen_selected_signed_hb_o3": oxygen.hb_flux/oxygen.o3_flux}
            for key, data in values.items():
                ci, draw = boot_median(data.to_numpy(), rng)
                estimate = float(data.median()) if len(data) else None
                row[key] = dict(n=len(data), median=estimate, interval95=ci)
                metrics[key][index] = (estimate, draw)
            overview.append(row)
        for key, by_bin in metrics.items():
            lo, a = by_bin[0]
            hi, b = by_bin[len(BIN_EDGES)-2]
            contrast.append(dict(field=field, metric=key, high_minus_low=(hi-lo if a is not None and b is not None else None),
                                 bootstrap95=np.percentile(b-a, [2.5, 97.5]).tolist()
                                 if a is not None and b is not None else None))
    distances = Planck18.luminosity_distance(injection_data.redshift.to_numpy()).to_value(u.cm)
    injection_data["flux_unit_per_luminosity_unit"] = 1e60/(4*np.pi*distances**2)
    for (index, line, sigma), group in injection_data.groupby(["bin", "line", "sigma_kms"], sort=True):
        blank_snr = group.base_flux/group.flux_error
        trials = [("flux", value, value) for value in FLUXES]
        trials.extend(("luminosity", value, value*group.flux_unit_per_luminosity_unit) for value in LUMINOSITIES[1:])
        for scale, level, flux in trials:
            found = (group.base_flux+group.response*flux)/group.flux_error >= DETECTION_SN
            per_source = group.assign(found=found).groupby("dr4_id").found.mean().to_numpy()
            sample = per_source[rng.integers(len(per_source), size=(4000, len(per_source)))]
            recovery.append(dict(bin=int(index), line=line, sigma_kms=int(sigma), scale=scale, level=level,
                                 galaxies=len(per_source),
                                 trials=len(found), detections=int(found.sum()),
                                 galaxy_mean_recovery=float(per_source.mean()),
                                 blank_snr_median=float(blank_snr.median()),
                                 blank_snr_robust_sigma=float(1.4826*(blank_snr-blank_snr.median()).abs().median()),
                                 bootstrap95=np.percentile(sample.mean(axis=1), [2.5, 97.5]).tolist()))
    return overview, contrast, recovery


def plot_recovery(recovery):
    fig = Figure(figsize=(9.5, 7.1), layout="constrained")
    axes = fig.subplots(2, 2, sharey=True)
    colors = ("#127089", "#8361a7", "#a64e35")
    units = (("flux", "Injected line flux (10⁻²⁰ erg s⁻¹ cm⁻²)"),
             ("luminosity", "Injected line luminosity (10⁴⁰ erg s⁻¹)"))
    for panels, (scale, xlabel) in zip(axes, units, strict=True):
        for ax, (line, label) in zip(panels, (("hb", "Hβ regions"), ("o3", "[O III] 5007 regions")), strict=True):
            for index, color in enumerate(colors):
                rows = [r for r in recovery if r["line"] == line and r["bin"] == index
                        and r["sigma_kms"] == SIGMAS_KMS[0] and r["scale"] == scale]
                x = [r["level"] for r in rows]
                y = [r["galaxy_mean_recovery"] for r in rows]
                ci = np.asarray([r["bootstrap95"] for r in rows])
                ax.errorbar(x, y, yerr=np.maximum(0., [y-ci[:, 0], ci[:, 1]-y]), fmt="o-", capsize=3, color=color,
                            label=f"{BIN_EDGES[index]} ≤ z < {BIN_EDGES[index+1]}")
            ax.set(xlabel=xlabel, title=label, ylim=(-.03, 1.03))
            ax.grid(axis="y", alpha=.3)
        panels[0].set_ylabel("Fraction recovered at formal S/N ≥5")
    axes[0, 1].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Artificial lines injected into observed JWST spectra · 125 km/s dispersion")
    save_figure(fig, ROOT / f"figures/{STEM}_recovery")


def plot_sideband(control, data):
    wave, flux, err, good = data
    rest = wave/(1+control["z"])
    keep = (rest >= 4945) & (rest <= 5065)
    y = np.where(good, flux, np.nan)
    fig = Figure(figsize=(9.2, 4.3), layout="constrained")
    ax = fig.subplots()
    ax.step(rest[keep], y[keep], where="mid", color="#216f82", label="Released G235M spectrum")
    ax.fill_between(rest[keep], (y-err)[keep], (y+err)[keep], step="mid", color="#216f82", alpha=.18)
    center = TARGET_LINES["o3"]
    companion = LINES_VACUUM_A["[O III] 4959"]
    ax.axvspan(center*(1-3000/SPEED), center*(1-1500/SPEED), color="#b47631", alpha=.18,
               label="Original blue continuum window")
    ax.axvspan(companion*(1-MASK_WIDTH/SPEED), companion*(1+MASK_WIDTH/SPEED), color="#ad4774", alpha=.14,
               label="[O III] 4959 mask")
    for wave, label in ((companion, "[O III] 4959"), (center, "[O III] 5007")):
        ax.axvline(wave, color="0.35", linestyle=":", linewidth=.8)
        ax.annotate(label, (wave, .96), xycoords=ax.get_xaxis_transform(), xytext=(4, 0), textcoords="offset points",
                    fontsize=9, va="top")
    ax.set(xlabel="Rest wavelength (vacuum Å)", ylabel="Fλ (10⁻²⁰ erg s⁻¹ cm⁻² Å⁻¹)",
           title=f"Observed JADES spectrum · {control['dr4_id']} · z={control['z']:.3f}", xlim=(4945, 5065))
    ax.margins(y=.15)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    save_figure(fig, ROOT / f"figures/{STEM}_continuum")


def main():
    ledger = json.loads(LEDGER.read_text(encoding="utf8"))
    if len(ledger["records"]) != ledger["planned"]:
        raise ValueError("Download ledger is incomplete")
    if ledger["target_table_sha256"] != text_digest(OUT / "observed_targets.csv"):
        raise ValueError("Download ledger does not match the target table")
    measurements, trials = [], []
    for entry in ledger["records"]:
        lines = {line: dict(status="NO_ACCEPTED_PRODUCT", flux=None, error=None) for line in TARGET_LINES}
        if entry["status"] == "ACCEPTED":
            lines, added, data = measure_target(entry)
            trials.extend(added)
            if entry["role"] == "sideband_control" and entry["field"] == "GN" and entry["redshift_bin"] == 0:
                plot_sideband(entry, data)
        measurements.append(signed_rows(entry, lines))
    measured, injection_data = pd.DataFrame(measurements), pd.DataFrame(trials)
    overview, contrast, recovery = statistics(measured, injection_data, pd.read_csv(OUT / "observed_targets.csv"))
    measures_path, injection_path = (OUT / f"{name}.csv" for name in ("forced_measurements", "injection_trials"))
    measured.to_csv(measures_path, index=False)
    injection_data.to_csv(injection_path, index=False)
    plot_recovery(recovery)
    write_json(OUT / "forced_summary.json", dict(
        executed_utc=utcnow(), input_sha256=text_digest(LEDGER), script_sha256=text_digest(__file__),
        measurements_sha256=text_digest(measures_path), trials_sha256=text_digest(injection_path),
        measurements=overview, high_minus_low=contrast, recovery=recovery, injected_fluxes=FLUXES,
        injected_luminosities=LUMINOSITIES[1:], luminosity_unit="1e40 erg s-1",
        luminosity_conversion="Planck18 luminosity distance, isotropic emission, with no attenuation or lens correction.",
        flux_unit="1e-20 erg s-1 cm-2", error_scope="Released diagonal variances including continuum-fit error.",
        injection_scope="Gaussian lines added to real off-line pixels, keeping noise and masks fixed. Errors bootstrap galaxies.",
        continuum_mask=dict(oiii_4959_half_width_kms=MASK_WIDTH, line_core_half_width_kms=HALF_WIDTH),
        target_summary_scope="Equal draw per field and redshift bin. Pooled values describe the pilot sample.",
    ))


if __name__ == "__main__":
    main()

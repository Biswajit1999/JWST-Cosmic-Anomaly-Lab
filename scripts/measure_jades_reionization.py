"""Measure signed Lyman-alpha fluxes and test recovery in JADES G140M spectra."""

import json

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.cosmology import Planck18
from download_jades_reionization_sample import BIN_EDGES, LEDGER, OUT, ROOT
from inspect_jades_cosmic_evolution import DENSITY_UNIT
from matplotlib.figure import Figure
from scipy.special import ndtr

from cosmic_anomaly_lab.provenance import digest, text_digest, utcnow, write_json
from cosmic_anomaly_lab.reduction_audit import read_dr4_spectrum
from cosmic_anomaly_lab.spectroscopy import pixel_edges
from cosmic_anomaly_lab.visual_evidence import save_figure

LYA = 1215.67
C_KMS = 299792.458
C_ANGSTROM_S = C_KMS*1e13
GRID = (150, 350, 550, 750)
WIDTHS = (150, 250)
EW_LEVELS = (10, 25, 50, 100)
CONTROLS = (-21000, -17000, -13000, -9000, 9000, 13000, 17000, 21000)
KNOWN_REST_LINES = (1215.67, 1238.82, 1242.80, 1260.42, 1302.17, 1334.53)
MIN_SN = 5
SEED = 453738
DRAWS = 3000


def spectral_arrays(spec):
    density = u.spectral_density(spec.wavelength)
    wave = spec.wavelength.to_value(u.AA)
    scale = 1e-20
    flux = spec.flux.to_value(DENSITY_UNIT, equivalencies=density)/scale
    err = spec.error.to_value(DENSITY_UNIT, equivalencies=density)/scale
    good = spec.good & np.isfinite(flux) & np.isfinite(err) & (err > 0)
    return wave, pixel_edges(wave), flux, err, good


def profile(edges, zero, shift, sigma_v):
    """Unit integrated flux in a red-side Gaussian, allowing for finite bins."""
    sigma = zero*sigma_v/C_KMS
    center = zero*(1+shift/C_KMS)
    limits = np.maximum(edges, zero)
    return np.diff(ndtr((limits-center)/sigma))/(ndtr(shift/sigma_v)*np.diff(edges))


def line_models(data, zero, width, omit=()):
    wave, edges, flux, err, good = data
    velocity = (wave/zero-1)*C_KMS
    valid = good & (edges[:-1] >= zero)
    for feature in omit:
        valid &= abs((wave/feature-1)*C_KMS) > 700
    side = valid & (velocity >= 2200) & (velocity <= 6800)
    if side.sum() < 8:
        return []
    norm = 1/err[side]**2
    continuum_var = 1/norm.sum()
    continuum_weights = norm*continuum_var
    continuum = continuum_weights @ flux[side]
    rows = []
    for shift in GRID:
        shape = profile(edges, zero, shift, width)
        core = valid & (abs(velocity-shift) <= 4*width)
        area = float(shape[core] @ np.diff(edges)[core])
        if core.sum() < 3 or area < .78:
            continue
        weights = shape[core]/err[core]**2
        weights /= weights @ shape[core]
        variance = np.sum((weights*err[core])**2)+weights.sum()**2*continuum_var
        rows.append(dict(offset_kms=shift, coverage=area, core=core, weights=weights,
                         side=side, continuum_weights=continuum_weights,
                         flux=float(weights @ (flux[core]-continuum)), error=float(np.sqrt(variance)),
                         continuum=float(continuum), continuum_error=float(np.sqrt(continuum_var))))
    return rows


def best_model(rows):
    return max(rows, key=lambda p: p["flux"]/p["error"]) if rows else None


def ew_flux_scale(z, m_uv):
    """Assume a flat F_nu continuum, normalised to the catalogue M_UV at 1500 Å."""
    d_l = Planck18.luminosity_distance(z).to_value(u.pc)
    fnu = 10**(-.4*(m_uv+48.6))*(10/d_l)**2*(1+z)
    flambda = fnu*C_ANGSTROM_S/(LYA*(1+z))**2/1e-20
    return flambda*(1+z)


def control_trials(data, target, zero, conversion):
    wave, edges, *_ = data
    avoid = np.array(KNOWN_REST_LINES)*(1+target["z"])
    for shift in CONTROLS:
        center = zero*(1+shift/C_KMS)
        if np.min(abs((avoid/center-1)*C_KMS)) < 1800:
            continue
        fits = line_models(data, center, WIDTHS[0], omit=avoid)
        if not fits:
            continue
        errors = np.array([row["error"] for row in fits])
        baseline = np.array([row["flux"] for row in fits])
        record = dict(dr4_id=target["dr4_id"], field=target["field"], bin=target["bin"],
                      z=target["z"], m_uv=target["m_uv"], pseudo_systemic_shift_kms=shift)
        yield dict(**record, ew_a=0, injected_offset_kms=None, dispersion_kms=None,
                   peak_snr=float(np.max(baseline/errors)), fixed_scale_flux=0.)
        good = data[4] & (edges[:-1] >= center)
        for feature in avoid:
            good &= abs((wave/feature-1)*C_KMS) > 700
        for width in WIDTHS:
            for offset in GRID:
                shape = profile(edges, center, offset, width)
                # Require usable pixels to contain at least 78% of the injected profile.
                area = shape[good] @ np.diff(edges)[good]
                if area < .78:
                    continue
                response = np.array([row["weights"] @ (shape[row["core"]]
                                     - row["continuum_weights"] @ shape[row["side"]]) for row in fits])
                for ew in EW_LEVELS:
                    injected = ew*conversion
                    peak = float(np.max((baseline+response*injected)/errors))
                    yield dict(**record, ew_a=ew, injected_offset_kms=offset, dispersion_kms=width,
                               peak_snr=peak, fixed_scale_flux=injected)


def measure_target(target):
    path = ROOT / target["spectrum"]["path"]
    if digest(path) != target["spectrum"]["sha256"]:
        raise ValueError(f"Spectrum fingerprint differs for {path.name}")
    data = spectral_arrays(read_dr4_spectrum(path, aperture=5, nods=3, wavelength_range=(.65, 5.4)))
    zero = LYA*(1+target["z"])
    conversion = ew_flux_scale(target["z"], target["m_uv"])
    models = line_models(data, zero, WIDTHS[0])
    best = best_model(models)
    confirmed_variant = None
    if best is not None and best["flux"]/best["error"] >= MIN_SN:
        spec_3px = read_dr4_spectrum(path, aperture=3, nods=3, wavelength_range=(.65, 5.4))
        variants = line_models(spectral_arrays(spec_3px), zero, WIDTHS[0])
        same = next((row for row in variants if row["offset_kms"] == best["offset_kms"]), None)
        confirmed_variant = same["flux"]/same["error"] if same else None
    base = dict(dr4_id=target["dr4_id"], field=target["field"], bin=target["bin"], z=target["z"],
                age_gyr=float(Planck18.age(target["z"]).to_value(u.Gyr)), m_uv=target["m_uv"],
                spectrum_status=target["status"], covered_templates=len(models), line_fit_snr_3px=confirmed_variant,
                continuum_proxy_flux_per_ew=conversion)
    reference = next((row for row in models if row["offset_kms"] == 350), None)
    if reference:
        base.update(reference_flux=reference["flux"], reference_error=reference["error"],
                    reference_signed_ew=reference["flux"]/conversion)
    if best:
        peak_snr = best["flux"]/best["error"]
        upper = max(0., max(row["flux"] for row in models))+3*max(row["error"] for row in models)
        base.update(best_velocity_kms=best["offset_kms"], best_snr=peak_snr, best_flux=best["flux"],
                    best_flux_error=best["error"], best_ew_a=best["flux"]/conversion,
                    nominal_ew_bound_a=upper/conversion if peak_snr < MIN_SN else None,
                    above_screen=peak_snr >= MIN_SN)
    else:
        base.update(above_screen=False)
    return base, list(control_trials(data, target, zero, conversion)), data


def bootstrap_means(values, rng):
    ix = rng.integers(len(values), size=(DRAWS, len(values)))
    return np.percentile(values[ix].mean(axis=1), [2.5, 97.5]).tolist()


def summarize(observed, trials):
    rng = np.random.default_rng(SEED)
    counts, recovery = [], []
    for i, (lo, hi) in enumerate(zip(BIN_EDGES[:-1], BIN_EDGES[1:], strict=True)):
        rows = observed[observed.bin == i]
        usable = rows[rows.covered_templates > 0]
        limits = usable.nominal_ew_bound_a.dropna()
        counts.append(dict(bin=i, redshift=[lo, hi], age_gyr=Planck18.age([hi,lo]).to_value(u.Gyr).tolist(),
                           drawn=len(rows), spectra=int(rows.spectrum_status.eq("ACCEPTED").sum()),
                           usable=len(usable), covered_all_offsets=int(usable.covered_templates.eq(len(GRID)).sum()),
                           formal_line_candidates=int(usable.above_screen.sum()),
                           negative_reference_flux=int((usable.reference_flux < 0).sum()),
                           median_nominal_ew_bound_a=float(limits.median()) if len(limits) else None))
    trials = trials.copy()
    trials["recovered"] = trials.peak_snr >= MIN_SN
    trials["offset_group"] = pd.cut(trials.injected_offset_kms, [0, 400, 800], labels=["150–350", "550–750"])
    trials["uv_group"] = pd.cut(trials.m_uv, [-22.01, -20, -17.99], labels=["bright", "faint"])
    for mode, filtered, axes in (
        ("pooled_offsets", trials[trials.dispersion_kms.eq(WIDTHS[0]) | trials.ew_a.eq(0)], ["bin", "ew_a"]),
        ("offset_width", trials[trials.ew_a > 0], ["bin", "ew_a", "dispersion_kms", "offset_group"]),
        ("uv_brightness", trials[trials.dispersion_kms.eq(WIDTHS[0]) | trials.ew_a.eq(0)], ["bin", "ew_a", "uv_group"]),
    ):
        for key, group in filtered.groupby(axes, observed=True, sort=True):
            frac = group.groupby("dr4_id").recovered.mean().to_numpy()
            record = dict(zip(axes, key, strict=True))
            recovery.append(dict(**record, mode=mode, targets=len(frac), trials=len(group),
                                 detections=int(group.recovered.sum()), recovery=float(frac.mean()),
                                 bootstrap95=bootstrap_means(frac, rng)))
    return counts, recovery


def plot(counts, recovery, observed, examples):
    fig = Figure(figsize=(9.3, 5.4), layout="constrained")
    ax = fig.subplots()
    for i, color in enumerate(("#196b80", "#8257a3", "#a44830")):
        rows = [r for r in recovery if r["mode"] == "pooled_offsets" and r["bin"] == i]
        xs = [r["ew_a"] for r in rows]
        y = np.array([r["recovery"] for r in rows])
        ci = np.array([r["bootstrap95"] for r in rows])
        ax.errorbar(xs, y, yerr=np.maximum(0., [y-ci[:,0], ci[:,1]-y]), color=color, marker="o", capsize=3,
                    label=f"{BIN_EDGES[i]:g} ≤ z < {BIN_EDGES[i+1]:g}")
    ax.set(xlabel="Injected rest-frame Lyα equivalent width (Å)", ylabel="Fraction above 5 formal S/N",
           ylim=(-.03,1.05), title="Lyα recovery in observed JADES G140M spectra")
    ax.grid(axis="y", alpha=.32)
    ax.legend(frameon=False)
    save_figure(fig, ROOT/"figures/jades_reionization_recovery")
    picks = []
    for count in counts:
        rows = observed[(observed.bin == count["bin"]) & observed.best_snr.notna()]
        if len(rows):
            picks.append(rows.loc[rows.best_snr.idxmax()])
    if not picks:
        return []
    fig = Figure(figsize=(9.3, 7.6), layout="constrained")
    for row, ax in zip(picks, fig.subplots(len(picks), 1, squeeze=False).ravel(), strict=True):
        wave, _, flux, err, good = examples[row.dr4_id]
        v = (wave/(LYA*(1+row.z))-1)*C_KMS
        take = (v >= -1700) & (v <= 4300)
        y = np.where(good, flux, np.nan)
        ax.step(v[take], y[take], where="mid", color="#156d82")
        ax.fill_between(v[take], (y-err)[take], (y+err)[take], step="mid", alpha=.22, color="#156d82")
        ax.axhline(0., color=".4", linewidth=.6)
        ax.axvline(0., color=".35", linestyle=":", linewidth=.8)
        ax.axvline(row.best_velocity_kms, color="#a44936", linestyle="--", linewidth=.8)
        ax.set(xlim=(-1700,4300), ylabel="Fλ (10⁻²⁰ cgs/Å)",
               title=f"{row.dr4_id} · z={row.z:.3f} · formal S/N {row.best_snr:.1f}")
    ax.set_xlabel("Velocity relative to systemic Lyα (km/s)")
    fig.suptitle("Observed JADES spectra · strongest Lyα template response in each bin")
    save_figure(fig, ROOT/"figures/jades_reionization_observed")
    return [row.dr4_id for row in picks]


def main():
    ledger = json.loads(LEDGER.read_text())
    if len(ledger["records"]) != ledger["planned"]:
        raise ValueError("Archive sample download is unfinished")
    observed, trials, examples = [], [], {}
    for target in ledger["records"]:
        if target["status"] != "ACCEPTED":
            observed.append(dict(dr4_id=target["dr4_id"], field=target["field"], bin=target["bin"],
                                 z=target["z"], spectrum_status=target["status"], covered_templates=0,
                                 above_screen=False))
            continue
        result, injected, data = measure_target(target)
        observed.append(result)
        trials.extend(injected)
        examples[target["dr4_id"]] = data
    obs, trial = pd.DataFrame(observed), pd.DataFrame(trials)
    counts, recovery = summarize(obs, trial)
    picks = plot(counts, recovery, obs, examples)
    observed_path, trial_path = OUT/"measurements.csv", OUT/"injections.csv.gz"
    obs.to_csv(observed_path, index=False)
    trial.to_csv(trial_path, index=False, compression={"method":"gzip", "mtime":0})
    lead = obs[obs.above_screen].sort_values("best_snr", ascending=False)
    write_json(OUT/"summary.json", dict(
        executed_utc=utcnow(), sample_sha256=text_digest(LEDGER), script_sha256=text_digest(__file__),
        observations_sha256=text_digest(observed_path), injections_sha256=digest(trial_path),
        counts=counts, recovery=recovery, candidate_ids=lead.dr4_id.to_list(), example_ids=picks,
        candidate_cut="maximum formal S/N >=5 in a fixed four-offset bank, primary dispersion 150 km/s",
        equivalent_width_scale="Planck18 luminosity distance, catalogue M_UV, flat F_nu",
        velocity_offsets_kms=GRID, injection_dispersions_kms=WIDTHS, injection_trials=len(trial),
        scope="Pilot in optically anchored, UV-limited targets, including nondetections. No neutral-fraction inference.",
    ))
    print(counts, flush=True)
    print(lead[["dr4_id", "z", "best_snr", "line_fit_snr_3px"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()

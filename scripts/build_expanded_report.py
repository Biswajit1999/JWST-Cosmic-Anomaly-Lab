"""Write a results-linked report and experiment ledger for the expanded search."""

import json
from pathlib import Path

import pandas as pd
import yaml

root = Path(__file__).resolve().parents[1]


def read(name):
    return json.loads((root / name).read_text(encoding="utf-8"))


summary = read("results/expanded/summary.json")
benchmark = read("results/glimpse_profile_comparison.json")
simulation = read("results/profile_simulation_validation.json")
inspection = read("results/expanded/inspection_audit.json")
cliff = read("results/cliff_sensitivity.json")
priorities = pd.read_csv(root / "results/expanded/followup_priority.csv")
survey_table = "\n".join(
    f"| {r['survey']} | {r['analysed']} | {r['joint_priorities']} |" for r in summary["by_survey"]
)
priority_table = "\n".join(
    f"| {r.source_id} | {r.survey} | {r.z_best:.3f} | {r.continuum_ratio:.2f} | "
    f"{r.minimum_broad_delta_bic:.1f} | {r.minimum_broad_sn:.1f} |"
    for r in priorities.head(12).itertuples()
)
inspection_lines = []
for r in inspection["sources"]:
    sid = r["source_id"]
    medium = r["medium_resolution_status"]
    if medium == "CONDITIONAL_PROFILE_FITS":
        medium = (
            f"{r['medium_grating']}: ΔBIC(narrow−Gaussian)={r['medium_delta_bic_narrow_minus_gaussian']:.1f}; "
            f"ΔBIC(Gaussian−exponential)={r['medium_delta_bic_gaussian_minus_exponential']:.1f}; "
            f"reduced χ²(G/E)={r['medium_fits']['gaussian']['reduced_chi2']:.2f}/"
            f"{r['medium_fits']['exponential']['reduced_chi2']:.2f}. "
            f"[Fit and residuals](../figures/expanded_sources/{sid}_medium.png)"
        )
    else:
        medium += ": " + r.get(
            "medium_reason", "No matching released medium-resolution product in this catalogue."
        )
    inspection_lines += [
        f"### {r['survey']} {r['source_number']} — {r['role']}",
        "",
        f"{sid}; z={r['redshift']:.4f}. Imaging: {r['image_status']}. {medium}",
        "",
        f"![Observed source evidence](../figures/expanded_sources/{sid}.png)",
        "",
    ]
bench_rows = []
for r in benchmark["comparisons"]:
    if r["lsf_scale"] == 1.7:
        bench_rows.append(
            f"| {r['window']} | {r['delta_bic_gaussian_minus_exponential']:.2f} | "
            f"{r['delta_bic_with_doubled_errors']:.2f} |"
        )
report = f"""# Expanded search for evidence of gas-enshrouded accretion

Research owner: Biswajit Jana. This is a measured research progress report, not a
black-hole-star discovery announcement. The conference posters summarize this
stage and can be regenerated when the evidence changes.

## Scope and coverage

The search completed diagnostics for **{summary["analysed"]:,} of {summary["eligible_denominator"]:,}
eligible upstream source records** ({summary["fraction_of_eligible_analysed"]:.1%}), with
{summary["failed"]} analysis failures. Eligibility was fixed before the bulk search:
PRISM, archive grade 3, median pixel S/N≥5, and 0.5≤z<12. The highest-S/N product
per upstream source ID was selected, with filename tie-breaking and no colour or
physical-class filter. We did not stop after an appealing source appeared.

The parent catalogue contains 20,144 upstream source records and 44,619
observations across four surveys. This run covers 13.1% of those source records,
not the majority of all public JWST objects. It covers the entire stated eligible
high-S/N subset. Targeting, release dates, redshift success and S/N selection
prevent interpreting these fractions as cosmic abundance or completeness.
Upstream IDs are not an independently verified count of physical galaxies.

| Assigned survey of selected product | Analysed | Joint follow-up priorities |
| --- | ---: | ---: |
{survey_table}

## Executed screening

- {summary["continuum_measured"]:,} usable continuum-ratio measurements.
- {summary["line_blends_measured"]:,} usable empirical Hα+[NII] blend comparisons.
- {summary["continuum_screen_pass"]} pass the continuum-only screen.
- {summary["broad_blend_screen_pass"]} pass the broad-blend-only screen.
- **{summary["joint_followup_priorities"]} pass both screens**, forming a follow-up queue.
- {summary["lsf_preference_changes"]} change the broad/narrow preference sign across LSF assumptions.
- {summary["joint_retained_with_doubled_line_errors"]} joint priorities retain the line criteria
  if line errors are doubled; this is sensitivity analysis, not a revised discovery threshold.

The continuum criterion requires red/blue fν>2 and its conditional 3σ lower
bound>1.5, with the lower-bound requirement repeated in shifted continuum windows.
The blend criterion requires ΔBIC(narrow−broad)>10 and broad-component S/N>5 under
each of three resolving-power scales (1, 1.3, 1.6). PRISM blends cannot separate
Hα, [NII], absorption, velocity structure and scattering. These cuts have no
calibrated physical-class false-positive rate or posterior probability.

![Expanded observed screening](../figures/expanded_search.png)

### Follow-up queue, not confirmed objects

| Upstream source ID | Survey | z | Continuum ratio | Minimum ΔBIC | Minimum blend S/N |
| --- | --- | ---: | ---: | ---: | ---: |
{priority_table}

The full queue and signed/missing measurements are preserved in
[follow-up priorities](../results/expanded/followup_priority.csv) and
[all measurements](../results/expanded/measurements.csv). Full per-source fit
details are in [compressed JSON Lines](../results/expanded/source_fits.jsonl.gz).
Selection does not establish novelty: some sources are already published, and
the complete queue has not undergone a source-by-source literature audit.

All eleven selected sources have poor absolute fits under the fixed-width PRISM
blend model: nominal reduced χ² ranges from approximately 4.9 to 4167. Relative
preference over an even poorer narrow model is not model adequacy. The queue is
therefore a list of spectra requiring better modelling, not eleven successful
physical-model fits.

### Literature and missed-benchmark checks

The queue includes already published sources. RUBIES-EGS-49140 and
UNCOVER-A2744-45924 are discussed in
[Wang et al. (2026)](https://arxiv.org/abs/2609.22623); RUBIES-EGS-9809 is described
as a dusty AGN host with stellar-population contributions in
[D'Eugenio et al. (2026)](https://arxiv.org/abs/2609.07825). These checks establish
that the queue is not a list of new discoveries. Later, deeper products in those
papers are not automatically present in our pinned DJA release.

The Cliff (RUBIES-UDS-154183; CAL-ebc2a6a997b3535d) is an especially important
benchmark from [de Graaff et al. (2025)](https://arxiv.org/abs/2503.16600). Our
continuum ratio is 8.46±2.12, but the fixed-width PRISM broad-line test rejects
it. A single 1500 km/s intrinsic sigma cannot represent all physically broad
profiles; it is also not the same quantity as FWHM. We retain this missed
published candidate as direct evidence of the screening method's incompleteness.
It is separately inspected and is not retroactively added to the eleven passes.

In its G395M spectrum, the full-window Gaussian-versus-exponential comparison
favours the exponential family: ΔBIC=56.8–69.6 over the LSF tests, remaining
10.7–13.9 when errors double. However, removing four central pixels within
±400 km/s reduces ΔBIC to 3.5–9.4 (−2.6 to −1.1 with doubled errors).
The signal is sensitive to the line core, where absorption is unmodelled.
This is conditional shape evidence, not a robust diagnosis from the wings alone.
See [all sensitivity cases](../results/cliff_sensitivity.json).

![The Cliff, higher-resolution profile comparison](../figures/expanded_sources/CAL-ebc2a6a997b3535d_medium.png)

## Published GLIMPSE-17775 benchmark

We downloaded the authors' stacked, slit-loss-corrected G395M spectrum and
verified its publisher MD5. This is reanalysis of their released spectrum, not
an independent detector-level reduction or a new source discovery.
Data: [Kokorev's public release](https://doi.org/10.5281/zenodo.17575391).
Context: [Kokorev et al.](https://arxiv.org/abs/2511.07515).

The comparison uses a linear local continuum, tied narrow nuisance components,
and either a Gaussian or a Gaussian convolved with symmetric exponential tails
for the broad component. Both are pixel integrated and convolved with a Gaussian
instrument response. We vary resolving power and use multi-start bounded fits.
These are empirical line shapes, not complete gas-envelope, stellar or AGN
forward models. Positive ΔBIC below favours the exponential shape.

| Window, resolving power ×1.7 | ΔBIC(G−E) | ΔBIC with doubled errors |
| --- | ---: | ---: |
{chr(10).join(bench_rows)}

**Hα is provisional:** the released flux and errors are finite, but NPIX is zero
throughout that window. The strict mask admits zero Hα pixels. The supplemental
finite-flux/error fit therefore requires resolution of the metadata discrepancy.
Even that exponential fit has large structured residuals (reduced χ²≈12.7), so
a large relative BIC improvement does not make it an adequate physical model.

**Paschen-δ uses the strict mask.** This implementation finds no clear evidence
for the extra exponential-tail parameter (ΔBIC≈−2.7); residuals remain above the
nominal error expectation. That differs from the published fit and is not a
refutation of it: the nuisance components, absorption treatment and fitting
assumptions have not been replicated in full. This disagreement is retained as
an unresolved methodological test, not hidden or counted as confirmation.

![Benchmark profile comparison and residuals](../figures/glimpse_profile_comparison.png)

## Validation and interpretation

Under a correctly specified isolated-line simulation, {simulation["gaussian"]["exponential_selected"]}/50
Gaussian controls and {simulation["exponential"]["exponential_selected"]}/50 exponential injections exceed
ΔBIC>10. The exact 95% intervals for these selection fractions are
{simulation["gaussian"]["binomial_interval_95"]} and
{simulation["exponential"]["binomial_interval_95"]}. This checks recovery under
known assumptions; it does not calibrate real-source discovery probabilities.
The simulations omit blending, absorption, covariance, calibration errors and
physical-model degeneracy. See [simulation records](../results/profile_simulation_validation.json).

**No new source is established as a black-hole star by this run.** The next
scientific requirements are independent extraction/metadata checks, realistic
line and continuum models, PSF-matched host/nucleus photometry, density and
fluorescence diagnostics, and comparisons against stellar, ordinary AGN and
composite alternatives. A negative or ambiguous result remains valid evidence.

## Bounded detailed inspection

We inspected the first twelve joint priorities (all eleven in this run) ordered by continuum lower bound
and blend ΔBIC, with up to four approximate same-survey comparison sources matched
in redshift and archive S/N. These controls are not certified normal galaxies.
The exact selection and higher-resolution outcomes are recorded in
[inspection audit](../results/expanded/inspection_audit.json).
Only one eligible medium-resolution counterpart per inspected source is used;
it shares the DJA processing and is not an independent reduction. Higher-resolution
fits use a nominal LSF ×1.3 and remain conditional on the adopted components.
The completed detailed set has eleven priorities, three comparison sources and
The Cliff as an added published benchmark. Seven sources yielded conditional
medium-resolution fits. In particular, JADES 10072503 prefers the narrow-only
model over a broad Gaussian (ΔBIC≈−11.6), warning that a PRISM screen can select
a source whose higher-resolution profile does not support the same conclusion.

{chr(10).join(inspection_lines)}
## Reproduction

```sh
python scripts/replay_manifest.py --budget-mib 4096
python scripts/download_expanded_sample.py
python scripts/analyse_expanded_sample.py
python scripts/analyse_glimpse_benchmark.py
python scripts/validate_line_profiles.py
python scripts/inspect_followup_sources.py
python scripts/audit_cliff_sensitivity.py
python scripts/build_expanded_report.py
cosmic-lab verify
```

Raw downloads use a finite 4 GiB cap and six bounded network workers with a single
manifest writer. All 2,647 requested PRISM products downloaded successfully.
The download cap does not authorize running multiple cache writers concurrently.
The earlier 20-source pilot and its manifest summaries remain historical records.
Use `scripts/replay_manifest.py --budget-mib 4096` to restore recorded bytes;
upstream checksum changes are rejected.
"""
(root / "reports/expanded_search.md").write_text(report, encoding="utf-8")
ledger = dict(
    stage="expanded-search-v1",
    executed_utc=summary["executed_utc"],
    scope="Full fixed high-SN eligible subset, not a complete JWST census",
    selection="configs/expanded_search.yaml",
    results="results/expanded/summary.json",
    benchmark="results/glimpse_profile_comparison.json",
    validation="results/profile_simulation_validation.json",
    inspection="results/expanded/inspection_audit.json",
    conclusion="Follow-up priorities; no new physical class established",
    remaining=[
        "Resolve GLIMPSE Halpha NPIX discrepancy",
        "Independent extraction and calibration",
        "Full competing physical forward models",
        "Realistic selection and false-positive validation",
        "Expand to fainter and uncertain-redshift sources without assuming the same completeness",
    ],
)
(root / "experiments/expanded_search.yaml").write_text(
    yaml.safe_dump(ledger, sort_keys=False), encoding="utf-8"
)
print("Wrote expanded research report and execution ledger")

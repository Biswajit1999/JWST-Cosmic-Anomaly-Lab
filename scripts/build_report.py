"""Build the pilot synthesis and experiment ledger from executed result files."""

import json
from pathlib import Path

import pandas as pd
import yaml

root = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((root / path).read_text(encoding="utf-8"))


run = read("results/run_summary.json")
catalogue = read("results/catalogue_summary.json")
spectra = pd.read_csv(root / "results/spectral_measurements.csv")
selection = pd.read_csv(root / "data/catalogues/pilot_selection.csv")
references = read("literature/curated_references.json")
sim = run["simulations"]
anomaly = run["anomaly"]
comparison = run["redshift_comparison"]
benchmark = spectra[spectra["sample"] == "A"].iloc[0]
lsf_switches = []
for path in (root / "reports/sources").glob("*.json"):
    evidence = json.loads(path.read_text())
    nominal = evidence["empirical_blend"]["fits"]
    if not nominal:
        continue
    preferred = min(nominal, key=lambda fit: fit["bic"])["model"]
    for scale, variant in evidence["lsf_sensitivity"].items():
        if variant["fits"]:
            alternate = min(variant["fits"], key=lambda fit: fit["bic"])["model"]
            if alternate != preferred:
                lsf_switches.append(
                    dict(
                        source_id=path.stem,
                        resolving_power_scale=scale,
                        nominal=preferred,
                        alternative=alternate,
                    )
                )
(root / "results/lsf_sensitivity_summary.json").write_text(
    json.dumps(lsf_switches, indent=2) + "\n"
)
lsf_switch_sources = len({r["source_id"] for r in lsf_switches})
manifest = read("data/manifests/downloads.json")
raw_bytes = sum(r["bytes"] for r in {r["path"]: r for r in manifest}.values())
survey_table = "\n".join(
    f"| {r['survey']} | {r['source_records']:,} | {r['observations']:,} |"
    for r in catalogue["by_survey"]
)
samples = selection["sample"].value_counts().to_dict()
reference_list = "\n".join(f"- [{r['title']}]({r['url']}) — {r['relevance']}" for r in references)
titles = [
    "Multi-survey public data discovery and catalogue construction",
    "Comparison of independent reductions for common sources",
    "Spectral diagnostics of diverse LRD candidates",
    "Compact galaxies versus obscured AGN",
    "Dense gas envelopes versus stellar atmospheres",
    "Balmer-break formation and interpretation",
    "Emission-line broadening mechanisms",
    "Gravitational lensing and luminosity uncertainty",
    "Photometric versus spectroscopic redshift disagreement",
    "Multi-band infrared morphology",
    "Multi-epoch variability",
    "Cross-survey source-population comparisons",
    "Unsupervised anomaly discovery",
    "Synthetic injection-recovery and false-positive testing",
    "Sources unexplained by all implemented physical families",
    "Forecasts for observations distinguishing physical interpretations",
]
questions = [
    "Can reproducible public source records be recovered in at least three surveys?",
    "Do independently reduced observations agree within measured uncertainties?",
    "Do multiple published LRD candidates share spectral diagnostics under uniform analysis?",
    "Can validated galaxy and AGN predictions be distinguished by the available data?",
    "Do gas-envelope or stellar-atmosphere models reproduce the line and continuum features?",
    "Can measured continuum discontinuities distinguish competing physical origins?",
    "Can instrumental, kinematic and scattering broadening be distinguished?",
    "Does lens-model uncertainty materially change the inferred intrinsic luminosity?",
    "Are photometric and spectroscopic redshift disagreements reproducible across reductions?",
    "Do wavelength-dependent sizes persist after PSF matching and deblending?",
    "Is intrinsic variability detectable after calibration and cadence effects?",
    "Do source properties remain consistent after controlling field and selection differences?",
    "Do outliers persist when their field is excluded from model training?",
    "Are known injected signals recovered and null false positives calibrated?",
    "Do all scientifically applicable validated physical models fail on the same sources?",
    "Which feasible measurement maximally separates competing predictive distributions?",
]
partial = {
    1: (
        "EXECUTED_PILOT",
        "results/catalogue_summary.json",
        None,
        "Three direct catalogue products; four survey subsets and five fields; complete survey selection functions remain unavailable.",
    ),
    2: (
        "PARTIAL",
        "results/redshift_comparison.json",
        "figures/catalogue_redshift_comparison.pdf",
        "Catalogue-level redshift comparison only; no independent spectral pixel-level reanalysis.",
    ),
    6: (
        "PARTIAL",
        "results/spectral_measurements.csv",
        "figures/observed_spectra.pdf",
        "Continuum proxy measurements only; no physical break-formation model comparison.",
    ),
    7: (
        "PARTIAL",
        "results/empirical_blend_fits.csv",
        "figures/observed_spectra.pdf",
        "Empirical narrow/broad profiles and LSF sensitivity; no separation of physical broadening mechanisms.",
    ),
    10: (
        "PARTIAL",
        "results/aperture_diagnostics.csv",
        "figures/observed_cutouts.pdf",
        "Registered cutouts and diagnostic apertures only; PSF matching and intrinsic morphology remain undone.",
    ),
    12: (
        "PARTIAL",
        "results/anomaly_feature_sample.csv.gz",
        "figures/archive_line_ratio_population.pdf",
        "Descriptive ratios and field holdouts; no completeness-corrected population inference.",
    ),
    13: (
        "EXECUTED_PILOT",
        "results/archive_ranking.csv",
        "figures/archive_line_ratio_population.pdf",
        "Statistical triage with field holdout; candidate vetting and broad physical feature coverage remain undone.",
    ),
    14: (
        "EXECUTED_BOUNDED_SIMULATION",
        "results/simulation_validation.json",
        "figures/simulation_validation.pdf",
        "Known-truth Gaussian population and idealized line recovery; performance does not transfer automatically to real galaxies.",
    ),
}
ledger = []
for number, (title, question) in enumerate(zip(titles, questions), 1):
    status, artifact, figure, limitation = partial.get(
        number, ("PLANNED", None, None, "Not executed in this pilot")
    )
    ledger.append(
        dict(
            id=f"EXP-{number:03d}",
            version=1,
            title=title,
            question=question,
            status=status,
            run_id=run["run_id"] if artifact else None,
            dataset="See product manifest and linked measurement table" if artifact else None,
            method="Documented in reports/multi_survey_investigation.md" if artifact else None,
            assumptions="See report, source evidence and data dictionary" if artifact else None,
            executed_code="src/cosmic_anomaly_lab/pipeline.py"
            if artifact and number != 1
            else "scripts/prepare_pilot.py"
            if number == 1
            else None,
            figures=[figure] if figure else [],
            measurements=artifact,
            uncertainties="Conditional statistics and explicit systematics in linked artifacts"
            if artifact
            else None,
            validation="tests/ and results/simulation_validation.json" if artifact else None,
            limitations=limitation,
            conclusion="See report; no new astrophysical class established" if artifact else None,
        )
    )
(root / "experiments/registry.yaml").write_text(
    yaml.safe_dump(ledger, sort_keys=False), encoding="utf-8"
)

report = f"""# A Multi-Survey Investigation of Unusual JWST Infrared Sources: Physical Modelling, Competing Interpretations, and the Search for Unexplained Astrophysical Populations

**Research owner: Biswajit Jana. Version 0.1 executed pilot.**

Run `{run["run_id"]}`; machine-recorded execution: `{run["executed_utc"]}`.
This report documents an initial investigation. It does not claim the full research
programme or all proposed experiments have been completed.

## Abstract

We constructed a persistent registry of {catalogue["ingested_source_records"]:,} upstream
source records and {catalogue["ingested_observations"]:,} spectral products from four JWST
surveys in five fields, retaining separate JADES and UNCOVER reference catalogues.
We downloaded and analysed {run["spectra_analysed"]} actual PRISM spectra, measured
{run["break_measurements"]} continuum-ratio proxies and {run["blend_fit_sources"]} local
emission-blend regions, and measured {run["aperture_measurements"]} diagnostic NIRCam
apertures. An error-aware two-feature search ranked {anomaly["ranked_sources"]} archive
sources using {anomaly["executed_folds"]} field/redshift holdouts. Eight nominal p<=0.01
outliers remain unvetted; none pass the reported 5% Benjamini-Hochberg adjustment.
Known-truth simulations give a {sim["anomaly"]["false_positive_rate"]:.1%} null false-positive
rate and {sim["anomaly"]["recovery_rate"]:.1%} recovery for the specified artificial injection.
Neither these simulations nor the observational ranking establishes a new physical class.

## Research motivation

Unusual infrared colour, a continuum break or a broad blend can arise from several
physical and observational causes. The platform therefore separates measurements,
published labels, model applicability and final interpretation. The pilot includes
one published unusual benchmark and independently selected comparison/archive objects.
GLIMPSE-17775 does not set the architecture or selection. The immediate falsifiable
tests concern archive accessibility, measurement recovery, catalogue agreement and
whether statistical outliers remain significant after accounting for many tests.

## Survey coverage and datasets

| Survey | Upstream source records | Spectral products |
| --- | ---: | ---: |
{survey_table}

Survey memberships overlap. The unique count of {catalogue["ingested_source_records"]:,}
uses DJA objid, not a new determination of independent physical objects. Lensed
images and blends have not been reconstructed into source-plane objects.
The full DJA snapshot contains {catalogue["dja_snapshot_rows"]:,} spectra; the four-survey
subset is retained. JADES DR4 contributes {catalogue["jades_reference_rows"]:,} reference
target rows and UNCOVER DR4.1 {catalogue["uncover_reference_rows"]:,}; they overlap the DJA
catalogue and are never added to its source count.

MAST public spectroscopy queries for programmes 1210, 1345 and 2561 returned
1,244, 5,253 and 1,187 observation rows respectively. Product inventories were
requested for one observation in each programme, not downloaded in bulk.
Release pages for additional collections are inventoried separately. A page
response is not evidence that every named product can be retrieved. Two guessed
landing endpoints failed (RUBIES and NGDEEP); RUBIES is actually released through
DJA, and NGDEEP has additional published archive paths, recorded in the inventory notes.

The raw cache contains {raw_bytes / 1024**2:.1f} MiB under a 512 MiB budget, including
catalogues, 20 spectra, three multi-band cutout files and a resolving-power table.
The DJA master catalogue matches its publisher MD5; all 28 products have local
SHA256 fingerprints. See [manifest](../data/manifests/downloads.json).

## Literature review

The curated references distinguish the survey/reduction publications from physical
interpretations. RUBIES, JADES and UNCOVER have different target-selection and
reduction histories. Published work on UNCOVER 45924 provides evidence beyond a
single low-resolution spectrum; this pilot does not reproduce all of that work.
Dense-gas interpretations and variability/X-ray non-detections motivate tests of
multiple explanations rather than adoption of one mechanism for every LRD.
Quasi-star scenarios are retained as hypotheses requiring quantitative forward models.

Crossref and OpenAlex were queried with two recorded queries each, retrieving ten
metadata records per query. Explicit identifiers deduplicate records; approximate
title matches do not automatically merge preprints and papers. The curated list
connects known arXiv/published versions where verified. ADS public access returned
an access error during research; no comprehensive ADS API search is claimed.
arXiv abstract pages and official release documentation were reviewed. Search logs
describe coverage; this is not a systematic or exhaustive literature review.

## Sample selection

Sample A has {samples.get("A", 0)} published-candidate benchmark (UNCOVER 45924).
Sample B has {samples.get("B", 0)} low-redshift comparison candidates; these are
operational comparison objects, **not independently certified normal-galaxy controls**.
Sample C has {samples.get("C", 0)} archive-selected objects. Grade 3 PRISM rows with
archive median pixel S/N>=5 are ordered by a fixed hash within survey and redshift
bins [0.5,2), [2,4), [4,7), [7,12). Repeated upstream objids are not counted twice.
There is no qualifying CEERS high-redshift entry for the requested final bin;
no replacement is invented. The selected redshift range is
{selection.z_best.min():.3f}-{selection.z_best.max():.3f}.

The parent spectroscopy is targeted and heterogeneous. Exposure time is recorded
per product; there is no common survey flux limit. Completeness is unestimated.
Redshift success, continuum S/N, wavelength coverage, upstream colour targeting,
public-release history and extraction quality all impose selection biases.
The line-ratio population additionally requires all three lines at S/N>=5 and a
line-fit redshift within 0.005(1+z) of the archive best estimate. It cannot measure
the incidence of weak-line, line-free, faint or otherwise excluded populations.

## Spectroscopy and imaging methods

FITS column units are validated and converted with Astropy. Nonfinite bins,
nonpositive errors, absent extraction pixels and flagged DQ are masked. Wavelengths
must increase strictly; 0.7-5.2 um is retained to limit edge/extended-order issues.
The published flux is used without an extra aperture/path-loss rescaling.
Released wavelength calibration is assumed; no independent lamp/sky calibration
is performed. Masked gaps remain gaps in figures.

The continuum ratio is inverse-variance mean fnu in rest 4000-4200 angstrom divided
by that in 3400-3600 angstrom. Each band needs at least three usable pixels and
mean S/N>3. This is a continuum proxy, not a unique physical Balmer-break measure.
Errors propagate from pixel variance with redshift held fixed.

The Halpha+[NII] region is fitted with a linear continuum and pixel-integrated
Gaussian blend templates at fixed redshift. Intrinsic sigma=100 and 1500 km/s are
template choices, not measured velocities. The LSF uses interpolated PRISM
resolving power with scale factors 1.0, 1.3 and 1.6; errors are also inflated by 1.5
for sensitivity checks. Linear coefficient covariance propagates to flux and EW.
The blend cannot separate Halpha, [NII], absorption, velocity shifts or scattering.

DJA cutouts supply WCS, SCI and inverse-variance WHT extensions for F150W, F277W
and F444W. Apertures have 0.3 arcsec radius; a clipped-median sky comes from
0.8-1.2 arcsec. PHOTFNU provides the unit conversion. The panels share sky positions
and angular scales but use independent display stretches. PSFs are unmatched;
drizzle covariance, aperture correction, segmentation and deblending are unfinished.
No intrinsic size or morphology difference is inferred.

![Observed spectra](../figures/observed_spectra.png)

![Observed cutouts](../figures/observed_cutouts.png)

## Physical models and competing interpretations

Implemented physical baselines are spherical blackbody Lnu, fully ionized-hydrogen
Thomson Eddington luminosity, efficiency-times-accretion-rest-energy luminosity,
and Planck18 luminosity-distance/redshift flux transformations. They use Astropy
constants and units. These are analytic approximations, not complete disc,
stellar-atmosphere, gas-envelope, SPS, nebular or lensing models. No blackbody is
fitted to a real source as proof of its physical class.

The hypothesis registry supports overlapping uncertain populations and checks
which observables a hypothesis requires. All full physical-family solvers are
honestly marked unimplemented. Empirical continuum/narrow/broad blend comparisons
are useful diagnostic alternatives, not a competition between AGN and stars.
The likelihood module supports covariance, a shared Gaussian calibration term,
and independent censored upper limits. Posterior sampling, physical prior studies
and reliable Bayesian evidences have not been performed.

## Simulation design and statistical inference

The seeded experiment (`20261110`) repeats a linear-continuum plus Gaussian line
recovery 500 times with known errors. The recovered amplitude pull mean is
{sim["line_recovery"]["mean_pull"]:.3f}, its standard deviation is
{sim["line_recovery"]["std_pull"]:.3f}, and nominal 95% coverage is
{sim["line_recovery"]["coverage_95"]:.1%}. This tests the conditional estimator under
its assumed model, not model misspecification or instrument systematics.

The anomaly experiment uses 1,500 synthetic training inliers, 1,000 calibration
inliers, 1,000 independent null controls and 500 shifted injections with varying
measurement errors. A 1.5 dex shift in one feature is intentionally large.
There are {sim["anomaly"]["false_positives"]}/1,000 false positives at p<=0.01
(95% binomial interval {sim["anomaly"]["false_positive_interval"][0]:.3%}-
{sim["anomaly"]["false_positive_interval"][1]:.3%}), and
{sim["anomaly"]["recovered"]}/500 injections recovered (95% interval
{sim["anomaly"]["recovery_interval"][0]:.1%}-{sim["anomaly"]["recovery_interval"][1]:.1%}).
There is no claim that real astrophysical contaminants have this distribution.

![Simulation validation](../figures/simulation_validation.png)

## Population results and anomaly detection

The measured continuum proxies span {spectra.break_ratio.min():.3f}-
{spectra.break_ratio.max():.3f} among the 16 usable cases. UNCOVER 45924 has
ratio {benchmark.break_ratio:.3f} +/- {benchmark.break_ratio_error:.3f} (conditional
statistical error), demonstrating a strong contrast between the specified windows.
It does not identify the emitting or absorbing physical medium.

The archive population includes {anomaly["feature_sources"]} sources. The two
features are log10([OIII] doublet/Hbeta) and log10([OIII] doublet/[OII]). Separate
[OIII] components are combined where necessary; their covariance is unavailable.
Robust intrinsic covariance is combined with each source's measurement covariance,
including the shared numerator. The calibration set is separate from training.
Each test field is held out in two redshift bins; a survey-independent reduction
is not implied. Full leave-one-survey-out validation remains additional work.

Eight nominal p<=0.01 objects are returned as an **unvetted triage list**. None
passes the reported 5% BH adjustment. In high-z folds, calibration samples of
64-91 objects have minimum p>0.01, so the absence of high-z flags is partly a
resolution limit of calibration, not evidence that the population is ordinary.
Population shift and shared fit errors further limit the p-value interpretation.
No cosmic abundance, physical discovery significance or real-data FPR is inferred.

![Observed line-ratio population](../figures/archive_line_ratio_population.png)

## Systematic uncertainty analysis

Within 0.3 arcsec, there are {comparison["candidate_pairs"]:,} positional reference
pairs; {comparison["unambiguous_robust_pairs"]:,} have unambiguous matches and qualifying
archive grades. Their median absolute normalized redshift difference is
{comparison["median_abs_normalized_difference"]:.6f}; {comparison["abs_normalized_difference_over_0p02"]}
exceed 0.02. These are candidates for checking, not 24 established incorrect
redshifts. Association, quality conventions, line identification and reduction
can all contribute. This comparison does not use published redshift errors to
test formal significance.

![Redshift comparison](../figures/catalogue_redshift_comparison.png)

LSF and error-inflation outputs are retained for each fitted source. In
{lsf_switch_sources} of the {run["blend_fit_sources"]} fitted sources, the preferred
empirical profile changes under the tested resolving-power assumptions. This
directly demonstrates why profile preference alone is not secure physical
broadening evidence. Absolute flux
calibration, correlated spectral noise, redshift uncertainty, aperture/path losses,
foreground contamination, alternate line identification, independent reductions,
blending and source-specific lensing remain material uncertainties. The 2D spectra
are downloaded but not systematically adjudicated. The shortlist is not classified
as physically unexplained until all those checks and a sufficient model library exist.

## Interpretation of unusual candidates

The published UNCOVER benchmark has diagnostic aperture fluxes of approximately
0.344, 2.487 and 4.016 microJy in F150W, F277W and F444W respectively, and a strong
measured continuum contrast. Independent stretches and unmatched
PSFs preclude a quantitative morphology conclusion. Published AGN interpretation
is recorded as literature evidence; it is not inferred solely from our broad
Gaussian template. All 20 source reports retain INSUFFICIENT_EVIDENCE for the
requested full physical comparison. Missing data do not imply a new class.

## Scientific limitations and future observational tests

The initial programme has delivered real metadata, spectra, cutouts and bounded
validation rather than complete physical modelling. It has not yet established
well-characterised controls across all requested populations, fitted complete
competing physical families, performed lens reconstructions or measured variability.
No entire sample is labelled as LRD, BH-star or anomalous.

The highest-priority next work is: inspect the 24 discrepant positional pairs
with their 2D and independent reductions; validate shortlist stability using
full line covariance and reduction alternatives; assemble independently characterised
control labels; build PSF-matched photometry; then add validated SPS, AGN and
gas/stellar-atmosphere grids before posterior comparisons. Hydrogen/helium ratios,
resolved line wings, long-wavelength continuum and multi-epoch observations are
qualitative discriminants. The generic separation function is tested, but no
sensitivity-based exposure forecast has been executed.

## Conclusions and experiment status

Public data support a broad, reproducible multi-survey investigation. Continuum
and line-ratio diversity is measurable in the pilot, while catalogue discrepancies
and selection/instrument limitations require further scrutiny. The implemented
search produces no multiplicity-adjusted physical discovery claim. No observation
has yet been established as unexplained by all meaningful physical models.

EXP-001 and EXP-013 ran as bounded pilots; EXP-014 completed its specified
synthetic experiment. EXP-002, 006, 007, 010 and 012 are partial. All other proposed
experiments are planned, including EXP-009 (photometric-vs-spectroscopic redshifts).
The current reference comparison is spec-z versus spec-z and does not complete
EXP-009. Every status and artifact is recorded in the experiment registry.

## Reproducibility and references

See README commands, the source-level JSON/Markdown reports, measurement CSVs,
input/code hashes in `results/run_summary.json`, the download manifest, dependency
lock, tests, and CI. Six figures are supplied as PNG and vector PDF. Source code
is MIT; dataset rights and acknowledgements remain with the original teams.

{reference_list}
"""
(root / "reports/multi_survey_investigation.md").write_text(report, encoding="utf-8")
print("Generated research report and versioned experiment registry.")

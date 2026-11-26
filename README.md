# JWST Cosmic Anomaly Lab

A provenance-first Python research platform for comparing ordinary, unusual, and
unclassified JWST sources across independent fields. **This is an executed
research programme in progress, not a completed physical classification.**

[Research record and report index](docs/research_record.md).

## Early galaxy evolution with JWST · 25 November 2026

The new JADES study compares oxygen and hydrogen emission in two GOODS
fields when the universe was 1.1–2.3 billion years old. Among 509 galaxies
with both lines detected in the same grating, the median ratio is about
10% higher in the earliest bin. The uncertainty includes zero change,
and the fraction of galaxies meeting the double-line cut varies by epoch.
The report records the selection checks and sets out follow-up studies
of reionization and early black-hole growth.
[Research goals and measurements](reports/jades_cosmic_evolution.md) ·
[Observed NIRCam cutouts and NIRSpec spectra](figures/jades_cosmic_evolution_observed.png).

We followed up with forced measurements in 96 public G235M spectra.
At the same artificial Hβ luminosity of 10⁴¹ erg s⁻¹, recovery fell
from 36% in the lowest-redshift bin to 5% in the highest. The follow-up
also documents and corrects contamination of the oxygen continuum window
in the original spectrum examples.
[Line recovery and continuum check](reports/jades_evolution_line_recovery.md).

## Lyman-alpha during reionization · 25 November 2026

We searched 93 public JADES G140M spectra from when the universe was
0.54–1.04 billion years old. The pilot recovered 18 published emitters
and found four line responses for follow-up. Tests with artificial lines
show that the sample's sensitivity varies across redshift bins.
[Reionization pilot](reports/jades_reionization.md) ·
[Observed JWST cutouts and spectra](figures/jades_reionization_unmatched.png).

The follow-up recovered a line consistent with Lyα in GOODS-N 58441
using separate G140M and PRISM exposure sets. The line position agrees
across the 2D traces, and the measured flux depends on the reduction.
The F090W excess at 58441 comes from a foreground star's diffraction
spike in the January 2025 SAPPHIRES imaging. Separate 2023 F090W
exposures show a compact source at 29–34 nJy. The literature audit
also recovered a published MMT detection at GOODS-N 38147.
[Exposure, calibration and identification checks](reports/jades_reionization_58441.md).

The expansion examined all 204 targets in the pilot's eligible JADES roster.
All 16 additional line hits match published emitters. Four CEERS objects
with both G140M and PRISM passed the positional and redshift checks.
The strongest response belongs to the previously published CEERS 1027.
[Full-sample results and observed CEERS spectra](reports/jades_reionization_full.md).

## Infrared rise and fade near JADES 103607 · 26 November 2026

A compact infrared source near **JADES 103607** brightened in five
JWST filters in October 2026 and largely returned to its earlier
brightness 46 days later. In F200W, the excess fell from **35.9 ± 4.3 nJy**
to **1.3 ± 4.0 nJy** above the October 2025 baseline. The October signal
appears in individual detector exposures and in two data reductions.
The November fade also appears on the original pixels of a different
NIRCam module. The fitted positions of the rise and fade agree within
9–21 mas in three filters. We found no earlier report of the event
in the public records checked through 26 November. Its physical cause
remains open.
[Analysis and archive receipts](reports/jwst_module_b.md) ·
[November photometry](results/jwst_cross_detector/module_b_001_november_photometry.json) ·
[Native detector images across three visits](figures/jwst_cross_detector/module_b_001_native_november_fade.png) ·
[Position of the extra light](figures/jwst_cross_detector/module_b_001_excess_position.png).

<img src="figures/jwst_cross_detector/module_b_001_november_photometry.png" width="810" alt="Observed JWST images of the October brightening and November fade near JADES 103607, with signed aperture photometry.">

The separate source-catalogue search of 40 public NIRCam products recovered
both the JADES 103607 fade and the known transient 2510042. The search
report includes the candidate measurements, observed cutouts, and a
native detector check that identified a saturation artefact.
[GOODS-S candidate search](reports/jwst_november_candidate_search.md).

The [additional Module A search](reports/jwst_module_a_pointing1.md)
covered 3.850 arcmin², recovered known transient 2510131 and traced a
three-filter false positive to a diffraction spike.

## Fading transient in JWST GOODS-S · 22 November 2026

A blind search recovered a fading infrared transient at **RA 53.162036°,
Dec −27.746384°**. It appears in six filters and was already bright on
30 September 2025 in a separate JWST programme. F200W fell from 415 nJy
in October 2025 to 49 nJy in October 2026, with intermediate measurements tracing the decline on
other detectors. The public Transient Name Server returned zero entries
within 5″. The public GO 8060 observing plan identifies the event as
transient **2510038**. A nearby galaxy overlaps the transient's predicted
slitless trace and accounts for a prominent hydrogen feature. The event's
physical classification remains open.
[Analysis and archive receipts](reports/jwst_cross_detector.md) ·
[Slitless spectrum overlap](reports/jwst_cross_detector.md#what-the-grism-frames-can-measure) ·
[Six-filter photometry and JWST images](figures/jwst_cross_detector/field_b_transient_colours.png) ·
[Hubble, Chandra and VLA comparison](reports/jwst_cross_detector.md#archival-x-ray-radio-and-optical-environment) ·
[MIRI and NIRSpec coverage](reports/jwst_cross_detector.md#coverage-before-brightening-and-in-the-transient-survey)

<img src="figures/jwst_cross_detector/field_b_transient_epochs.png" width="810" alt="JWST light curve and eight observed F200W cutouts showing the transient brighten and fade between 2022 and 2026.">

<img src="figures/jwst_cross_detector/field_b_transient_colours.png" width="810" alt="Six-filter JWST light curve and four observed F115W and F277W images of the transient in 2025 and 2026.">

## Earlier JWST GOODS-S discovery search · 22 November 2026

Six public JWST F444W visits show a candidate 2025 brightening at the known
galaxy ZFOURGE CDFS 18559. Its measured excess over 2026 persists after local
flux scaling and aperture recentering. Public NIRSpec records subsequently
identified a prior possible transient at that position. A separate apparent
brightening in three filters traced to the same detector pixel in all nine
exposures.
[Full analysis](reports/jwst_epoch_search.md) ·
[Native detector images](figures/jwst_epoch_search/target_dithers.png) ·
[NIRSpec follow-up](reports/jwst_cross_detector.md#recovery-of-an-earlier-transient-candidate)

<img src="figures/jwst_epoch_search/fading_lead_epochs.png" width="810" alt="Six real JWST F444W images of ZFOURGE CDFS 18559 from 2022 to 2026, with a photometry plot showing the measured 2025 peak.">

## Companion search in known Kepler systems · 21 November 2026

The search masked known transits across 236 Kepler light curves. Both
withheld known planets were recovered. None of 36 observed peaks from the
12 search systems passed replication. A weak residual on Kepler-1850 was
checked with two photometry products and three mask widths.
[Companion search and actual Kepler pixels](reports/kepler_companions.md)

## Transit discovery search · 21 November 2026

A 20–150 day search across 205 Kepler light curves recovered the known
planet Kepler-10 c. A 32.79-day peak on another star traced to every
third transit of the known Kepler-1437 b. A separate TESS residual search
found no repeating additional dip after fitting the dominant variability.
The report includes actual target pixels and signed event measurements.
[Transit discovery search](reports/transit_discovery_search.md)

## Kepler, TESS and Balmer profile follow-up · 21 November 2026

Across 13 Kepler quarters, candidate K08223.01 shows repeated dimming at its
catalogue timing, including 12 usable later events. Actual Kepler target-pixel
images accompany the aperture photometry. A separate TESS test uses 34 light
curves. Smooth variability models predict the second sector better than a
transit box for both search signals that repeat. The Kepler candidate remains
unvalidated. [Kepler and TESS results](reports/transit_followup.md)

Empirical Balmer profile fits retain the elevated Hγ/Hβ point estimates of
GN-89464 relative to GN-38432. The result depends on reduction and model
choices, and the physical explanation remains open. [Balmer profile results](reports/balmer_profiles.md)

## Actual images and a parallel exoplanet programme: 16 November 2026

The [spatial follow-up](reports/spatial_context.md) adds actual NIRCam images
for four targets. GN-89464 has a nearby distinct catalogue source, GN-38432,
whose hydrogen-line ratio is ordinary in the same measurement test. Their
centroids lie outside each other's released slit polygons; PSF contamination
and physical association remain unquantified.

The [first exoplanet pilot](reports/exoplanet_pilot.md) searches 24 Kepler
light curves from eleven search stars and one known control. TrES-2 is recovered;
no new planet is validated. Three synthetic injections are recovered, but a
control failure shows that automatic odd/even vetting needs calibration.
A published JWST spectrum of WASP-39 b supplies an atmospheric benchmark,
with explicit literature-derived physical properties and no new abundance claim.

![Actual source images and slit geometry](figures/spatial_context/gn89464_neighbor.png)

- [Actual Kepler host-star image](figures/exoplanet_pilot/tres2_observed_pixels.png)
- [JWST atmospheric spectrum](figures/exoplanet_pilot/wasp39_atmosphere.png)
- [Execution ledger](experiments/parallel_research.yaml)

## Additional gratings and a quality correction: 16 November 2026

The [GN-89464 follow-up](reports/gn89464_continuum.md) adds G235M/G395H spectra
and audits continuum coverage. Its elevated medium-resolution Hγ/Hβ ratios
survive a fixed 30-variant grid, but the additional gratings lack sufficient
dual-line signal to confirm them. The G395H Hβ measurement lacks its blue
continuum window. A two-sided coverage requirement also excludes eight earlier
GN-25734 Hα cores; 239 other previous cores remain usable with unchanged values.
This supersedes treating GN-25734's earlier two-nod ratios as accepted under
the stricter method. No physical classification is established.

- [Measurements and actual telescope images](reports/gn89464_continuum.md)
- [Execution ledger](experiments/gn89464_continuum.yaml)

![Actual high-resolution spectral images](figures/gn89464_followup/observed_g395h.png)

## Exposure-level research update: 16 November 2026

After preprint v1.0, the [new exposure audit](reports/exposure_audit.md) checks
21 individual-exposure/integration extractions for three Balmer priorities.
GS-10009106's combined ratio depends strongly on partial-coverage inputs; its
one extraction with usable coverage of both lines agrees closely with DJA.
GN-89464's elevated point estimates survive each single-extraction omission,
while GN-66336's inter-reduction discrepancy remains unresolved. These tests
do not identify a new black-hole star or establish the cause of the discrepancies.

- [Execution ledger](experiments/exposure_audit.yaml)
- [All signed core measurements](results/exposure_audit/signed_cores.csv)

![Exposure and coverage sensitivity](figures/exposure_audit/combination_sensitivity.png)

## Research paper and new validation: 16 November 2026

The [research preprint PDF](output/pdf/jwst_cosmic_anomaly_lab.pdf) consolidates
the executed project, with [editable LaTeX](paper/jwst_cosmic_anomaly_lab.tex),
[build instructions](paper/README.md) and a [LinkedIn draft](paper/linkedin_draft.md).
It is **not peer reviewed** and does not claim a new black-hole-star discovery.

Six additional DJA grating products change the interpretation of the earlier
Balmer queue: GN-66336 and GS-10009106's extreme DR4 ratios are not reproduced
at the same values; GN-89464 retains elevated Hγ/Hβ cores. Hubble/Spitzer
catalogues provide eight positional candidates for twelve targets. No Chandra
main-catalogue neighbour lies within 2 arcsec; no X-ray upper limit is inferred.

- [New measurements, telescope comparisons and limitations](reports/multiwavelength_validation.md)
- [Execution ledger](experiments/multiwavelength_validation.yaml)

![New cross-reduction check](figures/multiwavelength/reduction_ratio_comparison.png)

## Full JADES Balmer census: 16 November 2026

**5,190 DR4 target records** were audited in two extraction apertures. The
catalogue supports 1,216 dual-aperture Hα/Hβ comparisons and 448 Hγ/Hβ comparisons.
A twelve-record queue received 1D spectral checks: fifteen products passed
association checks, while one remains held for a programme-ID discrepancy.

There are 247 usable signed core measurements across 288 attempts. Three targets
retain measurable core pairs across all eight extraction/window variants, but
calibration, nod dependence, blending and stellar absorption remain competing
explanations. **No black-hole-star identification is established.**

- [New census, measurements and observed-image gallery](reports/jades_balmer_census.md)
- [Twelve-record follow-up queue](results/jades_balmer_census/followup_queue.csv)
- [Signed core measurements](results/jades_balmer_census/signed_core_measurements.csv)
- [Execution ledger](experiments/jades_balmer_census.yaml)

![Actual telescope spectral images](figures/jades_balmer_images/goods-n-mediumjwst_66336.png)

These 5,190 records overlap the earlier 4,731-record search; the counts are not
added together as unique objects. This is one catalogue census, not all archives.

## Wider high-resolution follow-up: 15 November 2026

**12 products for nine existing priorities** now have pixel-level coverage checks
and observed spectral sheets. Seven Hα products from six sources received
**42 conditional narrow/added-broad comparisons**. Four sources do not favour
the added broad Gaussian. Two already known sources improve substantially, but
their residuals still reject a satisfactory explanation by these simple models.

- [Results, individual spectra and model comparisons](reports/high_resolution_followup.md)
- [Full comparison table](results/high_resolution_coverage/halpha_comparisons.csv)
- [Execution ledger](experiments/high_resolution_coverage.yaml)

No new black-hole star is established; the total remains 4,731 analysed records.

## Independent-line control: 15 November 2026

The [OIII] doublet was checked in eleven extraction variants, independently of
the Hα fits. Twelve of 22 line/window checks are measurable; ten have no usable
G395M coverage. The observed widths constrain possible instrument blur but do
not uniquely separate it from gas motion. G235M widths also differ between
reductions, so the instrument-calibration ambiguity remains explicit.

- [New line-control report and observed spectra](reports/line_anchor.md)
- [Width measurements and conditional resolution bounds](results/line_anchor/width_constraints.csv)
- [Archive inventory for all 23 existing priorities](results/line_anchor/priority_followup_inventory.csv)

## Balmer-line benchmark: 15 November 2026

The known source JADES-GS-159717 now has explicit emission/absorption model
comparisons: **90 observed fits** across six extractions and five assumed
instrument resolutions, plus **60 synthetic controls**. The high-resolution
spectrum supports absorption structure within the tested model families.
Medium-resolution constraints depend strongly on instrument and extraction
assumptions; this is not a new black-hole-star discovery.
Synthetic controls also show that an incorrect instrument response can create
a false absorption preference even in emission-only spectra.

- [New benchmark report and observed model comparisons](reports/balmer_audit.md)
- [Resolution sensitivity and fit diagnostics](results/balmer_audit/model_preference.csv)
- [Execution ledger](experiments/balmer_audit.yaml)

## Reduction audit: 14 November 2026

Eight JADES DR4 spectra for seven shortlisted sources were compared across
32 extraction variants. Several unusual optical shapes persist, while the
GS-13905 ultraviolet anomaly weakens substantially in the separate reduction.
Its absolute UV flux also differs strongly. Actual 2D spectral images and
signed spatial measurements document the discrepancy; its origin remains unresolved.

**No new black-hole star is established.** This validates part of the existing
4,731-record search and does not add new unique objects.

- [New results, telescope spectral images and source comparisons](reports/reduction_audit.md)
- [Measurement table](results/reduction_audit/extraction_comparisons.csv)
- [Execution ledger](experiments/reduction_audit.yaml)

![Observed spectral images and reduction discrepancy](figures/uv_reduction_images.png)

## Open search: 13 November 2026

**4,731 source records analysed**, including 2,084 additional fainter spectra.
This covers 23.5% of the pinned parent catalogue and the entire stated grade-3,
S/N ≥ 2, 0.5 ≤ z < 12 PRISM subset. The added sample produces no additional
joint continuum/blend priorities; the earlier queue remains at 11.

The new optical/ultraviolet shape search does not require an LRD label, red
continuum, or broad hydrogen line. It preserves signed fluxes and missing
coverage, compares matched empirical peers, and checks noise, redshift,
dominant-band and survey dependence. **Rarity does not establish a new class.**

3,738 sources have enough coverage and matched peers for shape scoring. The
12-source follow-up queue includes seven from the fainter extension and has no
overlap with the previous 11. All 12 have observed evidence sheets; ten have
usable RGB cutouts. The method recovers published CEERS-3506 and JADES-GS-159717,
while alternate-product checks expose a fragile UV anomaly. A G395M/G395H
comparison also records the failure of a simple central-dip diagnostic on the
complex profile of JADES-GS-159717.

- [New research results, source images and spectra](reports/open_search.md)
- [Evolution, cosmology and tests for unfamiliar objects](reports/evolution_and_unknown_objects.md)
- [Run coverage and provenance](results/open_search/summary.json)
- [Execution ledger](experiments/open_search.yaml)

![Observed spectral-shape search](figures/open_search_overview.png)

The evolutionary calculation is an explicit set of conditional growth clocks,
not a mass measurement or a prediction of any source's future fate. Complete
competing physical forward models remain outstanding. The 14 November audit
above adds separate reductions for seven of the twelve priorities.

## Earlier high-S/N search: 12 November 2026

**2,647 / 2,647 eligible PRISM source records analysed**, across four surveys.
The fixed high-S/N subset is 13.1% of the 20,144-record parent catalogue, not a
complete JWST census. Eleven sources pass both empirical screens; all eleven,
three comparison sources and the published candidate The Cliff received detailed
image/spectrum inspection. Seven yielded conditional medium-resolution fits.

**No new black-hole star is established.** One PRISM priority favours a narrow
model at higher resolution, and The Cliff exposes a missed-candidate limitation. Its
exponential-profile preference weakens when the central pixels are excluded.
The GLIMPSE-17775 reanalysis also retains a quality-mask discrepancy and a
line-profile disagreement with the published analysis.

Start with the [expanded research report and source images](reports/expanded_search.md),
[measurements and coverage](results/expanded/summary.json), and
[execution ledger](experiments/expanded_search.yaml).

### Conference posters

Two A1 landscape posters of the **12 November high-S/N experiment**, with vector
plots and actual telescope images. They remain a dated research snapshot; the
new open-search results above extend that stage.

- [Project purpose and methods (PDF)](output/pdf/jwst_project_purpose_methods_A1.pdf)
- [Evidence and limitations (PDF)](output/pdf/jwst_evidence_limits_A1.pdf)
- [Poster previews and regeneration instructions](output/pdf/README.md)

![Project purpose and search methods](output/pdf/jwst_project_purpose_methods_A1.png)

The expanded download uses a finite **4 GiB** raw-data cap. Run:

```sh
python scripts/replay_manifest.py --budget-mib 4096  # restore all pinned inputs first
python scripts/download_expanded_sample.py
python scripts/analyse_expanded_sample.py
python scripts/analyse_glimpse_benchmark.py
python scripts/validate_line_profiles.py
python scripts/inspect_followup_sources.py
python scripts/audit_cliff_sensitivity.py
python scripts/build_expanded_report.py
```

## Original 20-source pilot

The pilot ingests 20,144 DJA source records and 44,619 observations from JADES,
CEERS, RUBIES, and UNCOVER, across GOODS-N, GOODS-S, EGS, UDS, and Abell 2744.
It also imports the independent JADES DR4 and UNCOVER DR4.1 reference catalogues.
These counts overlap: they are not counts of independent physical galaxies.

Twenty actual NIRSpec/PRISM spectra are analysed. Nine diagnostic apertures are
measured in NIRCam cutouts of three sources. The code produces error-aware
archive rankings, field-held-out calibration, synthetic injection/recovery tests,
source evidence reports, and figures labelled OBSERVED or SIMULATED.

**See the actual sources in the [image and spectrum gallery](reports/visual_gallery.md):
20 individual evidence sheets, 18 usable three-band colour composites, and
the observed 2D spectral images. Two sources have partial image-band coverage.**

![Observed NIRCam comparison](figures/observed_image_comparison.png)

**Have we identified a black-hole star? No.** The observational pilot is complete
within its stated scope; competing physical models have not been tested enough
to establish that interpretation. Read the [physical-test status](reports/black_hole_star_status.md).

The imaging extension is recorded separately in
[visual execution metrics](results/visual_atlas_run.json) and
[rendering methods](docs/visual_evidence.md). The original pilot's nine apertures
and run summary remain a historical record.

**For the original pilot, start with the [research report](reports/multi_survey_investigation.md),
[execution metrics](results/run_summary.json), and
[experiment registry](experiments/registry.yaml).**

![Observed pilot spectra](figures/observed_spectra.png)

## Installation and use

Python 3.12 or later. From the repository root:

```sh
python -m venv .venv
# Activate .venv using your shell's activation command.
python -m pip install -e '.[dev]'
pytest -q
ruff check src tests scripts
cosmic-lab --help
```

For the exact tested environment, `requirements-lock.txt` records installed
versions. Its Windows-only dependency is platform-marked for Linux CI.

```sh
python scripts/probe_archives.py       # ~227 MB of three catalogue products
python scripts/prepare_pilot.py        # ingest, select, download 20 spectra
python scripts/download_instrument.py  # resolving power reference
python scripts/download_cutouts.py     # three sources, three NIRCam filters
python scripts/download_cutouts.py --all  # request images for all 20 pilot sources
python scripts/build_visual_atlas.py   # images, 2D spectra, evidence sheets, gallery
cosmic-lab verify                     # SHA256 verification
cosmic-lab analyse                    # numerical analysis, figures, reports
python scripts/build_report.py        # human-readable synthesis of outputs
cosmic-lab sources --survey UNCOVER --limit 5
```

Run these commands from the root; `cosmic-lab --root PATH ...` can explicitly
select another checkout. Network-dependent scripts are separate from offline CI.
`python scripts/replay_manifest.py` restores the exact recorded download set,
refusing upstream changes that do not match recorded checksums.

## What exists

| Capability | Scope |
| --- | --- |
| Discovery | Public MAST metadata/product queries; pinned DJA, JADES, UNCOVER catalogue adapters |
| Provenance | Atomic bounded downloads, publisher checksum where available, SHA256 snapshots, version history |
| Registry | SQLite sources/observations/evidence/link tables; CSV exports; spherical candidate crossmatching |
| Spectroscopy | FITS units and quality masks, specutils interoperability, rest wavelengths, continuum ratio, pixel-integrated Gaussian blend fits, conditional EW/errors, LSF sensitivity |
| Imaging | WCS-aware aperture diagnostics and registered cutout panels; no PSF-corrected morphology |
| Inference | Gaussian covariance/calibration likelihood, independent censored upper-limit likelihood, empirical BIC comparisons |
| Physics | Explicit analytic blackbody, Eddington/accretion luminosity and cosmological flux baselines |
| Anomaly search | Complete-feature error-aware robust distance, shared-ratio covariance, holdout calibration and field validation |
| Validation | Seeded line recovery, synthetic anomaly FPR/recovery, unit/integration tests |
| Literature | Crossref/OpenAlex query logs, explicit-identifier deduplication, curated references and research gaps |

The complete mission is preserved in [the specification](docs/research_specification.md).
Only implemented functionality has Python modules. Configured hypotheses do not
imply implemented physical solvers. No source is declared a new class or a black-hole star.

## Data and reproducibility

Large data and SQLite are excluded from Git. Small catalogues, manifests,
measurement tables, source reports, and PNG/vector PDF figures are committed.
The default raw-data budget is 512 MiB; it includes cached products. Interrupted
downloads restart the incomplete file; completed files are verified and reused.
The cache is **single-writer**. A local hash is not a publisher authentication check.

`data/manifests/downloads.json` records exact URLs, versions, timestamps, sizes,
checksums, and ETags. `data/catalogues/sources.csv.gz` uses the upstream DJA objid
namespace, not newly inferred physical associations. JADES and UNCOVER reference
rows remain separate. Ambiguous positional matches are retained without merges.

The [data dictionary](docs/data_dictionary.md) explains units and missing values.
The [selection contract](configs/pilot.yaml) documents population and measurement
biases. Common DJA processing is **not independent reduction validation**.

## Scientific boundaries

Local continuum and blend fits cannot distinguish stellar populations, obscured
AGN, dense envelopes, or composites. Their statistical uncertainties omit several
systematics. PRISM broad-profile preference is not a black-hole mass measurement.
Image apertures lack PSF matching, aperture corrections, deblending and drizzle
covariance. A statistical ranking is a triage list, not a discovery catalogue.

Physical SPS/atmosphere/radiative-transfer grids, lens reconstructions, posterior
sampling/evidence, calibrated variability analyses, and quantitative observing
forecasts remain unimplemented. See the report and experiment registry for the
precise next work and for experiments that only ran in a limited pilot form.

## Credit and license

Research owner and code author: Biswajit Jana. Code: MIT. Upstream observations,
calibrations and catalogue rights remain with their creators; the code license
does not relicense their data. Cite JWST/STScI, survey teams, DJA/msaexp, Astropy,
astroquery, specutils, NumPy, SciPy, pandas, scikit-learn and Matplotlib as appropriate.
Survey and software references are in [the literature registry](literature/curated_references.json).


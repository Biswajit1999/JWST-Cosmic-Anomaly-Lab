# JWST Cosmic Anomaly Lab

A provenance-first Python research platform for comparing ordinary, unusual, and
unclassified JWST sources across independent fields. **This is an executed
research programme in progress, not a completed physical classification.**

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


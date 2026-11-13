# Have we tested the black-hole-star interpretation?

**Update, 13 November 2026:** the [open search](open_search.md) adds 2,084 fainter
sources for 4,731 analysed records in total. It also searches optical and UV
shapes without requiring the original LRD-like cuts, and adds an
[evolution/cosmology research note](evolution_and_unknown_objects.md).
No new physical class or black-hole-star identification is established.

**Update, 12 November 2026:** the [expanded search](expanded_search.md) now covers
2,647 eligible PRISM source records, 15 detailed source inspections, seven
conditional medium-resolution fits, GLIMPSE-17775 and The Cliff benchmark tests.
The table below documents the earlier pilot. The expanded work still does not
establish a new black-hole-star identification; it exposes both candidate signals
and important quality/selection limitations.

**Not yet to a standard that identifies a black-hole star.** We have executed an
observational pilot and visual evidence extension. We have not completed the
competing physical-model analysis requested in the research specification.
No source in this repository is a new black-hole-star discovery.

The central question remains whether unusual JWST sources are better explained
by gas-enshrouded accretion, stellar populations, ordinary or obscured AGN, or
composites. The broader multi-survey comparison and ordinary controls are needed
to avoid selecting only sources that already resemble the desired explanation.

## What has actually run

| Test | Executed here | What it establishes |
| --- | --- | --- |
| Observed spectra | 20 PRISM spectra with unit/quality checks | Released spectral data available for inspection |
| Continuum windows | 16 usable red/blue ratios | A measured continuum contrast, conditional on archive z and calibration |
| Local hydrogen-region profiles | 14 empirical Hα+[NII] fits with LSF sensitivity | Conditional broad/narrow blend preference; not an electron-scattering model |
| Images and spatial inspection | 20 evidence sheets; 18 RGB cutouts; two sources with partial band coverage; actual 2D spectra | Visual access to counterparts, nearby structure and rectified traces |
| Population triage | 857 archive line-ratio sources with field/redshift holdout calibration | Anomaly ranking; zero detections at the pilot's BH-adjusted p≤0.05 threshold |
| Known-truth validation | Idealized line and synthetic-population injections | Performance under those simulations, not completeness for real black-hole stars |

Here **BH-adjusted** means the Benjamini–Hochberg multiple-testing adjustment,
not a black-hole classification. The pilot's published LRD/AGN benchmark is
UNCOVER 45924 at archive z=4.4673. It is not GLIMPSE-17775. GLIMPSE-17775 was
absent from that original pilot; its released spectrum has now been analysed in
the expanded work, with an unresolved Hα quality-mask discrepancy explicitly retained.

## What a defensible physical test still needs

1. Validate counterpart, deblending, extraction and wavelength-dependent PSF/slit
   losses; use independent reductions and assess alternative redshifts. In lens
   fields, propagate magnification and possible differential magnification.
2. Assemble PSF-matched photometry, appropriate higher-resolution spectra and
   ancillary detections/limits. Test line profiles, absorption, line ratios and
   density/ionization diagnostics with their instrumental response and errors.
3. Fit validated stellar+nebular, AGN+dust, dense-envelope/reprocessing and
   composite forward models consistently. An analytic blackbody or a broad
   Gaussian is not a dense-gas radiative-transfer calculation.
4. Compare posterior predictions and held-out observables, including calibration
   nuisance terms, selection effects and model flexibility. Test whether each
   physically plausible competitor can explain the same observations.
5. Establish recovery and false-positive behaviour for realistic simulated
   alternatives before interpreting a mismatch as a new class.

These are outstanding research tasks, not results produced by the new gallery.

## Why images alone cannot settle it

NASA's account of GLIMPSE-17775 describes evidence from many spectral lines,
including line broadening and dense-gas diagnostics. Those are tests of physical
conditions that cannot be recovered from a red compact appearance alone.
See [NASA's GLIMPSE-17775 evidence summary](https://science.nasa.gov/missions/webb/nasa-webb-finds-strongest-evidence-yet-for-black-hole-stars/).

Competing evidence must also be retained: Ji et al. (2026) report far-UV line
features in Abell2744-QSO1 that challenge a fully covering envelope and suggest
open directions or clumpy gas. That is a source-specific constraint, not a
refutation of all enshrouded-accretion scenarios or a result for our benchmark.
See [the research preprint](https://arxiv.org/abs/2604.03370).

Read the [observed gallery](visual_gallery.md) alongside the
[original pilot report](multi_survey_investigation.md). The new images improve
access to evidence without changing the original scientific conclusions.

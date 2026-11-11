# CODEX MASTER PROMPT V2
## JWST Cosmic Anomaly Lab — Multi-Survey Astrophysical Discovery and Physical Modelling

**Repository:** `JWST-Cosmic-Anomaly-Lab`  
**GitHub account:** `Biswajit1999`  
**Research owner:** Biswajit Jana  
**Research domain:** Observational astrophysics, astronomical spectroscopy, computational physics, cosmology, statistical inference, and scientific machine learning.

### 1. Central research mission

Create, execute, validate, and maintain a substantial open-source Python research platform for investigating unusual astrophysical objects observed by the James Webb Space Telescope.

**Critical change from the previous specification: this project must NOT be restricted to GLIMPSE-17775, black hole stars, or previously identified Little Red Dots.**

GLIMPSE-17775 is one useful benchmark among many. It must not determine the architecture, source-selection criteria, scientific conclusions, or limits of the investigation.

Instead, search and analyse publicly available observations across numerous JWST programmes, fields, instruments, and source populations.

Our scientific objective is to answer four overarching questions:

1. What physical mechanisms produce the unusual infrared and spectroscopic properties of objects discovered by JWST?
2. Can established astrophysical models adequately explain these observations, including their uncertainties?
3. Are distinct observational populations being incorrectly grouped into a single physical class?
4. Are there reliable sources in public astronomical archives whose observable characteristics are not satisfactorily explained by the models currently available?

The outcome may support known physics, expose model degeneracies, challenge a published interpretation, or identify candidates for future observations.

Do not presuppose that the answer must involve a black hole, a new stellar class, or exotic physics.

**The goal is to investigate the Universe using observational evidence, not to generate evidence for a predetermined conclusion.**

---

### 2. Investigate multiple astronomical populations

Create a configurable source-class and hypothesis registry.

At minimum, investigate the following populations where relevant public data exist:

- Little Red Dots (LRDs)
- Black-hole-star candidates and gas-enshrouded accretors
- Supermassive-star candidates
- Direct-collapse black-hole candidates
- Ordinary and heavily obscured active galactic nuclei
- High-redshift quasars
- Compact star-forming galaxies
- Extreme emission-line galaxies
- Metal-poor galaxies and candidate Population III environments
- Strongly gravitationally lensed galaxies
- Sources with anomalous Balmer breaks
- Extreme red infrared sources
- Dust-obscured galaxies
- Variable accretion sources and transients
- Unresolved composite sources
- Foreground stellar or brown-dwarf contaminants
- Sources with uncertain or contradictory redshift determinations
- Objects rejected by standard photometric selection criteria
- Objects whose spectra are poorly explained by existing templates

Consider other physically motivated populations discovered during the literature review.

Do not assume these categories are mutually exclusive or observationally established as distinct physical classes.

Support multi-label and uncertain classifications.

Include an explicit `UNKNOWN_OR_UNRESOLVED` assessment, but do not treat it as proof of a new astrophysical class.

---

### 3. Discover and inventory observational datasets first

Before building a model around any particular object, perform a broad archive investigation.

Search:

**JWST and survey archives**

- MAST JWST: https://archive.stsci.edu/missions-and-data/jwst
- JADES: https://archive.stsci.edu/hlsp/jades
- JADES spectroscopic releases: https://jades-survey.github.io/scientists/data.html
- CEERS: https://ceers.github.io/releases.html
- UNCOVER: https://jwst-uncover.github.io/
- GLASS-JWST: https://archive.stsci.edu/hlsp/glass-jwst
- COSMOS-Web: https://cosmos.astro.caltech.edu/page/cosmosweb-dr
- DAWN JWST Archive: https://dawn-cph.github.io/dja/
- NEXUS: https://nexus-jwst.github.io/
- Public RUBIES spectroscopy
- PRIMER
- FRESCO
- NGDEEP
- Other relevant publicly released JWST observations

Use MAST programme searches, source coordinates, survey identifiers, and published catalogues to locate actual files.

Do not assume a named programme has a particular dataset available without checking.

**Ancillary astronomical resources**

Where justified, examine:

- HST imaging
- ALMA observations
- Chandra observations
- XMM-Newton observations
- VLA radio data
- Euclid public data
- Spitzer observations
- Gaia catalogues for foreground discrimination
- Ground-based deep imaging and spectroscopy
- Published gravitational lens mass models

Only include measurements relevant to the scientific questions and physically compatible with the target fields.

Record non-detections and upper limits instead of silently omitting them.

---

### 4. Build a multi-survey data discovery engine

Implement Python modules that automatically:

1. Search available public observations.
2. Retrieve catalogue metadata.
3. Identify spectroscopic and imaging products.
4. Distinguish raw, calibrated, and high-level science products.
5. Match objects across independent observations.
6. Download selected data products.
7. Verify file integrity and provenance.
8. Standardise units and metadata.
9. Construct searchable source-level records.
10. Track changes when archive catalogues or reductions are updated.

Use `astroquery`, `astropy`, `specutils`, `pandas`, `numpy`, and suitable archive APIs.

Create a persistent source registry with:

- Unique internal source ID
- Original survey source IDs
- Right ascension and declination
- Astrometric uncertainty
- Redshift and uncertainty
- Redshift method and quality flag
- Survey or field
- Observation programme IDs
- Instruments and spectral configurations
- Filter information
- Imaging and spectroscopy availability
- Crossmatched observations
- Published classifications
- Known possible lensing effects
- Data provenance
- Archive file references
- Source selection criteria
- Data-quality flags

Use established astronomical coordinate-matching conventions.

Avoid incorrect crossmatching in crowded or strongly lensed fields.

Maintain observation-level and source-level tables separately.

Support multiple observations of the same physical source and multiple images of the same lensed source.

Do not download the entire JWST archive indiscriminately.

Provide metadata-first discovery, selective downloads, caching, storage budgets, and resumable processing.

---

### 5. Construct a scientifically diverse research sample

The initial analysis sample must not consist exclusively of Little Red Dots.

Select sources across different surveys and physical categories.

Build three complementary samples:

**Sample A — Published unusual objects**

Include previously reported compact infrared sources, Little Red Dots, possible black-hole stars, and other unusual sources.

GLIMPSE-17775 may be included, but it must be treated as one object among other benchmark sources.

**Sample B — Astrophysical control populations**

Include conventional galaxies, AGN, emission-line galaxies, and other appropriate comparison populations with well-characterised observations.

**Sample C — Archive-discovered sources**

Select objects independently from broader public catalogues, without relying solely on previously assigned classifications or Little Red Dot colour cuts.

This sample is essential for discovering unexpected observational populations.

For every sample, document:

- Selection function
- Survey depth
- Wavelength coverage
- Redshift distribution
- Signal-to-noise requirements
- Completeness
- Selection biases
- Exclusion criteria

Use multiple independent fields to assess field-to-field variation.

Do not use the label "unknown" merely because an object is faint or has missing measurements.

---

### 6. Multi-object spectroscopy

Build a generalized spectroscopy pipeline applicable to different objects, surveys, and spectral resolutions.

Implement:

- FITS ingestion
- Wavelength calibration validation
- Flux-unit conversions
- Rest-frame conversion
- Spectral masking
- Continuum estimation
- Emission- and absorption-line detection
- Line profile fitting
- Equivalent-width measurements
- Balmer-break measurements
- Spectral slope estimates
- Velocity-width measurements
- Line ratios
- Uncertainty propagation
- Spectral residual analysis

Support different instruments, gratings, exposure depths, extraction methods, and spectral resolutions.

Convolve comparison models with the appropriate instrumental line-spread function.

Do not combine spectra from different instruments without assessing calibration and resolution differences.

Use line-identification hypotheses appropriate for the measured redshift.

Compare photometric and spectroscopic redshifts where both exist.

Investigate cases in which different identifications imply substantially different astrophysical interpretations.

For each object, generate a source-specific spectroscopic evidence report.

---

### 7. Multi-band imaging and photometry

Develop a parallel image-analysis pipeline using NIRCam, MIRI, and ancillary images.

Implement:

- FITS image reading
- WCS handling
- Astrometric alignment
- Background estimation
- Source detection
- Source segmentation
- Aperture photometry
- PSF fitting and PSF-matched photometry
- Photometric error calculations
- Source-size estimates
- Morphology measurements
- Multi-band image comparison
- Blending analysis
- Colour-colour diagrams
- Spectral energy distribution construction

Use published photometric catalogues as a reference when available.

Investigate whether an object changes appearance between infrared filters.

Test whether a supposedly compact source contains extended or blended components.

Check image calibration, PSF differences, deblending, detector artefacts, and aperture corrections before interpreting morphology differences as physical.

Generate multi-band source cutouts with consistent astrometry and scientifically meaningful scales.

---

### 8. Physical simulation programme

Create a substantial computational astrophysics package with modular forward models.

Simulate:

**A. Black-hole accretion**

- Black-hole mass
- Accretion rate
- Eddington ratio
- Radiative efficiency
- Disc temperatures
- Luminosity
- Emission spectra
- Characteristic variability timescales

**B. Dense gas envelopes**

- Gas density profiles
- Temperature structures
- Optical depths
- Ionisation
- Electron scattering
- Covering fractions
- Envelope geometry
- Spectral reprocessing

**C. Stellar and supermassive-star models**

- Stellar masses
- Effective temperatures
- Luminosities
- Atmospheric absorption
- Winds
- Emission and absorption-line behaviour

Use published stellar-atmosphere grids or established solvers where necessary.

**D. Galaxies and stellar populations**

- Stellar population synthesis
- Metallicity
- Star-formation history
- Nebular emission
- Dust attenuation
- Compact starburst populations
- Composite stellar and AGN contributions

**E. Gravitational lensing**

- Magnification
- Image multiplicity
- Distortion
- Source-plane reconstruction
- Differential magnification
- Lens-model uncertainty

Use real cluster mass models when available.

**F. Cosmological effects**

- Redshift
- Luminosity distance
- Cosmological dimming
- Intergalactic attenuation
- Observable flux and wavelength transformations

**G. Instrument-aware simulated observations**

Forward-model physical predictions into the actual filters, wavelength ranges, resolutions, and noise properties of the relevant observations.

**H. Population simulations**

Generate synthetic source populations with realistic parameter distributions, observational selection, and noise.

**I. Transient and variability models**

Where multi-epoch observations exist, test whether variability provides evidence favouring particular physical mechanisms.

All simulations must clearly distinguish approximate toy models from physically complete calculations.

Use tested physical constants and `astropy.units`.

Do not implement arbitrary physics without references.

---

### 9. Competing explanations for every object

Do not build a single BH* classifier.

Implement a general model-comparison system.

For each source, assess scientifically applicable alternatives, including:

- Accreting black hole
- Obscured AGN
- Accretion disc with unusual temperature structure
- Dense reprocessing envelope
- Supermassive star
- Compact galaxy
- Extreme starburst
- Composite source
- Gravitationally magnified source
- Incorrect redshift or foreground interloper
- Data-processing artefact

Only compare a hypothesis where its available forward model and the observations make the test meaningful.

Use Bayesian and frequentist diagnostics as appropriate.

Suggested libraries:

- SciPy
- emcee
- dynesty
- ArviZ
- PyMC or NumPyro where appropriate

Implement:

- Likelihood functions
- Physical parameter priors
- Posterior sampling
- Credible intervals
- Bayesian evidence when reliable
- Posterior predictive checks
- Prior-sensitivity studies
- Model residual diagnostics
- Simulated parameter-recovery tests

Account for upper limits, correlated noise, and calibration uncertainties where applicable.

Distinguish:

- Preferred model
- Statistically competitive alternatives
- Inconclusive evidence
- Missing observables
- Models rejected under documented assumptions

Never conclude that an object is a black-hole star simply because a gas-envelope model produces a plausible spectrum.

---

### 10. Search systematically for unexplained objects

This is one of the central scientific contributions.

Develop an open-set anomaly-discovery pipeline.

It must search beyond previously labelled exotic candidates.

Look for:

- Unusual infrared colour combinations
- Rare continuum shapes
- Unexpected spectral breaks
- Mixed absorption/emission signatures
- Extreme line ratios
- Unexpected emission-line widths
- Unexplained wavelength-dependent morphology
- Disagreement between spectroscopic and photometric measurements
- Sources with atypical variability
- Objects consistently misfitted by physical models
- Statistically significant outliers across independent observations

Begin with interpretable methods:

- PCA
- Robust distance metrics
- Isolation Forest
- Local Outlier Factor
- Gaussian mixture models
- Density-based clustering

Explore neural methods when a sufficiently representative dataset and validation framework exist:

- Autoencoders
- Variational autoencoders
- Contrastive learning
- Self-supervised spectral embeddings

Do not allow neural-network performance to substitute for physical interpretation.

Account for instrumental heterogeneity, survey selection functions, missing features, heteroscedastic errors, and redshift-dependent observables.

Develop simulated injection-recovery experiments and independent holdout testing.

Measure false-positive rates.

Before ranking a candidate as genuinely unusual, test:

1. Calibration and reduction errors.
2. Incorrect crossmatching.
3. Foreground contamination.
4. Alternative redshift solutions.
5. Strong emission-line contamination of broadband photometry.
6. Source blending.
7. Lensing effects.
8. Noise fluctuations.
9. Incomplete model libraries.
10. Known but rare astrophysical populations.

A high anomaly score is a reason for investigation, not evidence of a new physical class.

---

### 11. Search for physical mechanisms not adequately represented by current models

After testing established models, identify reproducible residual patterns that recur across independent sources.

Ask:

- Do multiple sources show the same unexplained spectral feature?
- Are apparent anomalies concentrated in a particular redshift range?
- Do anomalous objects share similar physical environments?
- Are discrepancies associated with particular observing configurations?
- Do competing models fail on the same observables?
- Can a new combination of known physical processes reproduce the evidence?
- Do any published models predict features not yet searched for in the data?

Investigate plausible mechanisms including:

- Early black-hole seed growth
- Gas-enveloped accretion
- Rapid stellar evolution
- Metal-poor stellar populations
- Nonstandard accretion-disc configurations
- Complex radiative transfer
- Extreme nebular emission
- Composite and unresolved systems
- Rare transients
- Strong or differential lensing

Review more speculative theoretical models only where they provide quantitative predictions that can be compared with observations.

Avoid proposing exotic explanations before conventional systematics and physical alternatives have been tested.

Any proposed new interpretation must produce at least one falsifiable observable prediction.

---

### 12. Cross-survey validation

Ensure that a claimed effect is not produced by one survey's processing or source-selection method.

Where possible:

- Repeat analyses using independent JWST fields.
- Compare independent reductions of the same observations.
- Compare photometric and spectroscopic evidence.
- Test different signal-to-noise thresholds.
- Repeat model fitting under alternate calibration assumptions.
- Check whether source ranking is stable across data-processing variants.
- Evaluate whether anomaly detection generalises to surveys excluded from training.

Use leave-one-survey-out validation where sample sizes permit.

Report whether physical conclusions remain stable after changing fields, templates, selection thresholds, and model priors.

---

### 13. Build an observational discriminant engine

For each unresolved physical classification, determine which measurements would best distinguish the leading explanations.

Examples include:

- Particular hydrogen or helium line ratios
- Line-wing morphology
- Balmer-break structure
- Continuum shape at longer infrared wavelengths
- Variability amplitude and timescale
- Relative compactness in different bands
- Expected X-ray or radio emission
- Lensing multiplicity or source reconstruction
- Constraints on source size and luminosity

Develop a software module that computes the expected separation between competing models in observable space.

Where reliable forward models and sensitivity information exist, estimate which additional observations would most reduce model ambiguity.

Distinguish proposals supported by quantitative forecasts from qualitative follow-up suggestions.

---

### 14. Build a source-level evidence database

Each analysed source must have a persistent scientific record.

Include:

- Source ID and aliases
- Survey and coordinates
- Redshift evidence
- Available data
- Data-quality assessment
- Spectral measurements
- Photometric measurements
- Morphological measurements
- Lensing information
- Competing model fits
- Posterior uncertainties
- Anomaly score and its calibration
- Possible systematic explanations
- Literature references
- Follow-up requirements
- Current interpretation

Generate machine-readable tables and source reports.

Classification statuses should include:

`KNOWN_MODEL_CONSISTENT`

`MULTIPLE_MODELS_VIABLE`

`POORLY_EXPLAINED`

`INSUFFICIENT_EVIDENCE`

`LIKELY_SYSTEMATIC`

`REQUIRES_FOLLOWUP`

Do not automatically translate `POORLY_EXPLAINED` into a new-object discovery.

---

### 15. Build a literature and research-gap database

Search ADS, arXiv, Crossref, OpenAlex, and relevant journal archives.

Cover theoretical predictions, observations, alternative interpretations, simulations, and follow-up studies.

Track:

- Paper metadata
- DOI or arXiv ID
- Publication status
- Target objects
- Relevant surveys
- Observational evidence
- Physical assumptions
- Conclusions
- Contradictory results
- Missing measurements
- Proposed future work

Deduplicate preprints and published versions.

Quantify literature coverage with documented search queries and dates.

Connect every scientific hypothesis and experimental result to the relevant literature.

Implement a research-gap registry that identifies:

- Unanswered observational questions
- Contradictory model predictions
- Unused public datasets
- Unreproduced published measurements
- Areas with insufficient validation
- Physical degeneracies requiring new observables

---

### 16. Proposed Python repository structure

Use a tested, installable `src`-layout package.

```text
JWST-Cosmic-Anomaly-Lab/
├── README.md
├── CITATION.cff
├── LICENSE
├── pyproject.toml
├── AGENTS.md
├── configs/
├── data/
│   ├── manifests/
│   ├── catalogues/
│   ├── raw/              # ignored by Git
│   ├── processed/        # ignored by Git
│   └── synthetic/        # ignored by Git
├── src/cosmic_anomaly_lab/
│   ├── archives/
│   ├── catalogues/
│   ├── crossmatching/
│   ├── imaging/
│   ├── photometry/
│   ├── spectroscopy/
│   ├── cosmology/
│   ├── physics/
│   ├── accretion/
│   ├── stellar_populations/
│   ├── radiative_transfer/
│   ├── lensing/
│   ├── simulations/
│   ├── inference/
│   ├── source_classification/
│   ├── anomaly_detection/
│   ├── population_analysis/
│   ├── model_discriminants/
│   ├── literature/
│   ├── validation/
│   ├── reporting/
│   └── cli.py
├── notebooks/
├── experiments/
├── figures/
├── results/
├── reports/
├── literature/
├── tests/
├── scripts/
└── .github/workflows/
```

Only create modules when implementing meaningful functionality.

Notebooks should call reusable package functions.

Keep downloaded large datasets outside ordinary Git history.

Use a manifest-based reproducibility system with URLs, checksums, archive identifiers, versions, and access dates.

---

### 17. Scientific experiments

Implement a versioned experiment registry.

Initial proposed experiments:

**EXP-001:** Multi-survey public data discovery and catalogue construction.

**EXP-002:** Comparison of independent reductions for common JWST sources.

**EXP-003:** Spectral diagnostics of diverse LRD candidates.

**EXP-004:** Compact galaxies versus obscured AGN.

**EXP-005:** Dense gas envelopes versus stellar atmospheres.

**EXP-006:** Balmer-break formation and interpretation.

**EXP-007:** Emission-line broadening mechanisms.

**EXP-008:** Gravitational lensing and luminosity uncertainty.

**EXP-009:** Photometric versus spectroscopic redshift disagreement.

**EXP-010:** Multi-band infrared morphology.

**EXP-011:** Multi-epoch variability where suitable data exist.

**EXP-012:** Cross-survey source-population comparisons.

**EXP-013:** Unsupervised anomaly discovery.

**EXP-014:** Synthetic injection-recovery and false-positive testing.

**EXP-015:** Sources unexplained by all implemented model families.

**EXP-016:** Forecasts for observations that distinguish competing interpretations.

Each experiment requires:

- Falsifiable research question
- Dataset
- Method
- Scientific assumptions
- Executed code
- Figures
- Measurements
- Uncertainties
- Validation
- Limitations
- Conclusions

Do not claim completion unless the experiment actually ran.

---

### 18. Deliver scientific results, not just source code

The project should generate:

1. A searchable, documented multi-survey source catalogue.
2. A curated observational dataset with verified provenance.
3. Reproducible individual-source analyses.
4. Population-level measurements.
5. Physical forward-model simulations.
6. Statistical comparisons of competing interpretations.
7. Validation and sensitivity experiments.
8. Quantified anomaly-search performance.
9. A ranked list of sources deserving further investigation.
10. Scientific reports explaining what the evidence supports and what remains uncertain.

All figures must identify whether data are observed or simulated.

All physical measurements must have units.

Scientific uncertainties must be propagated.

No results may be fabricated or labelled as executed if only planned.

A scientifically negative result is acceptable when it rules out a hypothesis under stated assumptions or establishes that the available evidence cannot distinguish alternatives.

---

### 19. Execution strategy

Do not spend the entire first development session on one source.

**Phase 0 — Repository and infrastructure**

Verify the authenticated GitHub account, local development environment, and requested author identity.

Set up packaging, testing, CI, documentation, and data provenance.

**Phase 1 — Multi-survey discovery**

Search at least three accessible public survey collections during the initial pilot, where available.

Create a real catalogue that spans multiple source types.

Retrieve a practical initial set of spectra and associated metadata, subject to storage, accessibility, and scientific quality.

Record the actual coverage; do not invent source counts.

**Phase 2 — Comparative observations**

Analyse several independently selected objects from at least two distinct survey fields, including normal control sources and unusual candidates.

GLIMPSE-17775 may be included but must not dominate this stage.

Generate reproducible spectral and photometric measurements.

**Phase 3 — Physical simulations**

Implement reusable physical baselines, instrument-aware simulations, and parameter-recovery tests.

**Phase 4 — Model comparison**

Fit multiple meaningful hypotheses and assess statistical degeneracies.

**Phase 5 — Population investigation**

Analyse variation in spectra, redshifts, morphology, and other available observable properties across source classes.

**Phase 6 — Blind anomaly search**

Evaluate sources that were not selected solely because they were already labelled unusual.

Calibrate false positives with controls and simulations.

**Phase 7 — Scientific synthesis**

Generate a research report identifying:
- Findings reproduced
- Plausible physical interpretations
- Hypotheses disfavoured
- Unresolved degeneracies
- Sources needing follow-up
- Possible unexplained populations
- Observational tests most likely to resolve the remaining questions

Keep each phase resumable.

Do not leave the repository full of incomplete scaffolding.

Prioritise validated implementations and actual scientific analysis.

---

### 20. GitHub authorship

Use the repository name:

`JWST-Cosmic-Anomaly-Lab`

Create it under:

`Biswajit1999`

Requested Git author:

`Biswajit Jana`

Use the owner's configured, verified Git email.

Do not invent an email address.

Do not attribute commits to another person.

Do not add AI co-author footers unless required by an applicable policy.

Accurately acknowledge external software, dataset creators, and scientific literature.

Use meaningful commits, run tests, and verify successful pushes.

Do not modify unrelated repositories.

If authentication is unavailable, complete the local implementation and report the exact limitation without claiming a successful remote operation.

---

### 21. Scientific integrity rules

Never:

- Treat a simulated source as an observed discovery.
- Claim every LRD contains a black hole.
- Treat an anomaly score as physical proof.
- Claim gravitational lensing is responsible without evidence.
- Infer precise physical parameters from insufficient data.
- Present photometric redshifts as spectroscopic confirmation.
- Ignore non-detections or selection effects.
- Fabricate measurements or literature citations.
- Treat approximate radiative-transfer models as physically complete.
- Claim a new astrophysical population without checking known alternatives.
- Use a single object to establish population-wide conclusions.
- Claim the entire research programme has been completed when only the software framework exists.

Always identify the assumptions that materially affect the interpretation.

---

### 22. Final research deliverables

Produce a research-grade repository and an accompanying report titled:

**"A Multi-Survey Investigation of Unusual JWST Infrared Sources: Physical Modelling, Competing Interpretations, and the Search for Unexplained Astrophysical Populations"**

Include:

- Abstract
- Research motivation
- Survey coverage
- Literature review
- Dataset descriptions
- Sample selection
- Spectroscopy and imaging methods
- Physical models
- Simulation design
- Statistical inference
- Population results
- Anomaly detection
- Systematic uncertainty analysis
- Interpretation of unusual candidates
- Scientific limitations
- Future observational tests
- Conclusions
- References

Provide publication-quality figures and machine-readable result tables.

At the end of each Codex execution, report:

- Actual GitHub URL and commit SHAs
- Data sources accessed
- Number of real sources successfully ingested
- Number of source spectra actually analysed
- Surveys represented
- Experiments completed
- Models implemented
- Simulations executed
- Validation results
- Figures generated
- Scientific findings supported by evidence
- Remaining unexplained observations
- Next highest-priority tasks

Do not invent completion metrics.

### 23. Begin execution

Start with the multi-survey data discovery engine and scientifically diverse source catalogue.

Verify datasets and archive access first.

Implement real-data processing before constructing elaborate speculative models.

Do not restrict initial development to GLIMPSE-17775.

Do not restrict the project to existing Little Red Dot catalogues.

Develop the ability to investigate known, unusual, and previously unclassified sources across independent observations.

**The ultimate objective is not merely to explain a newly reported black-hole-star candidate. It is to determine what JWST observations can establish about the physical nature of unusual cosmic sources, identify where current models fail, and discover which unanswered questions can be tested with data that already exist.**

Execute the work, validate it, and commit verified progress to GitHub.

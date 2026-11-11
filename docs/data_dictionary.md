# Data dictionary and measurement conventions

`data/catalogues/sources.csv.gz` has one row per upstream DJA v4.4 `objid` in the
four selected surveys. `source_id` is a deterministic UUID-derived internal key.
`ra_deg`, `dec_deg` are ICRS decimal degrees from the preferred upstream row.
`z_spec` is the upstream best redshift, not a new estimate; missing/sentinel
redshifts become empty CSV cells. `z_quality` uses the DJA grade convention.
`n_observations` counts released spectral products, not independent exposures.
Survey memberships overlap and their source counts must not be summed.

`observations.csv.gz` preserves one row per spectral file, original `srcid`,
`objid`, programme, instrument configuration, release root, extraction version,
exposure seconds, wavelength limits in micrometres, and the product URL.
`zline` is the redshift used for upstream line fitting; it can disagree with
`z_best`. `grade=3` indicates the archive's robust assessment. `valid` refers
to upstream redshift consistency and is not a complete detector-quality flag.

DJA `line_*` and `line_*_err` columns use **1e-20 erg/s/cm2**; the original column
documentation is downloaded and hashed. Archive `eqw_*` are observed-frame
angstroms. `phot_dr` is the nearest photometric counterpart distance in arcsec;
no such nearest neighbour is automatically accepted as a physical association.
`z_phot` is photometric evidence and is never spectroscopic confirmation.

JADES and UNCOVER reference files retain their own IDs and quality conventions.
JADES `z_Spec_flag` is not numerically converted to a DJA grade. The reference
comparison uses JADES A/B and UNCOVER quality >=2, and excludes ambiguous
positional candidates. The UNCOVER local-background advisory is retained.
UNCOVER redshift percentile columns and non-detections remain available.

## Derived results

- `spectral_measurements.csv`: continuum red/blue **fnu** ratio, dimensionless,
  using 4000-4200 / 3400-3600 rest-frame angstrom windows. Pixel errors propagate
  under independent Gaussian errors, conditional on the fixed archive redshift.
  Neither window is guaranteed free of stellar/nebular features. This is not a
  standard D4000 index or an isolated physical Balmer discontinuity.
- `empirical_blend_fits.csv`: local integrated Halpha+[NII] empirical Gaussian
  flux in erg/s/cm2, continuum in erg/s/cm2/angstrom, rest EW in angstrom,
  conditional uncertainties and BIC. Signed amplitudes are retained. The quoted
  `max(0,flux)+3*error` limit is an operational conservative Gaussian limit,
  not an exact bounded posterior quantile. Intrinsic profile width is an assumed
  template input in km/s, not a recovered gas velocity or virial mass.
- `aperture_diagnostics.csv`: flux density and diagonal error in microJy in a
  0.3 arcsec radius aperture. This is not total, PSF-matched photometry. Weights
  are taken as inverse variance in native image units, using `EXTVER=WHT`;
  their copied BUNIT is not interpreted as a flux image. `PHOTFNU` converts
  native science units to Jy. Aperture and clipped-sky uncertainty are propagated;
  drizzle covariance, crowded-field background bias and calibration remain open.
- `catalogue_counterpart_candidates.csv.gz`: positional pairs within 0.3 arcsec,
  including ambiguous ones. No merge occurs. The `ambiguous` flag is symmetric:
  either source with multiple neighbours is excluded from agreement statistics.
- `archive_ranking.csv`: two line-ratio features in dex, error-aware robust
  distance and held-out calibration p-value. Shared [OIII] numerator covariance
  is included; full spectral-fitting covariance is unavailable. Both [OIII]
  components are summed when the archive supplies separate components, assuming
  independent component errors. Reported adjusted p-values use Benjamini-Hochberg
  over the ranked sample and are descriptive under distribution shift.
- `anomaly_folds.json`: field-held-out train/calibration/test counts in each
  redshift bin. Minimum attainable p is 1/(Ncal+1), so small calibration sets
  cannot flag p<=0.01. This is limited power, not proof of no high-z anomalies.
- `simulation_validation.json`: explicitly simulated null controls, injections,
  line recovery and analytic physics; exact binomial confidence intervals.

All missing values are retained as null/empty rather than set to zero.
The catalogue imports observations rather than simulating them. Simulated
experiments cannot identify real discoveries or measure the archive's true FPR.

## Cache and version semantics

`downloads.json` is append-only for completed downloads. A repeat without refresh
verifies the cached SHA256 and makes no network request. Explicit refresh creates
a new manifest version with `previous_sha256`; old bytes are not retained. Replay
uses the newest record for each path and requires its SHA256. Cache writes are
single-process, and a download must fit the remaining storage budget including
existing files. Partial files are removed on failure, completed files persist.
Only the DJA master catalogue has a publisher MD5 in this pilot; other SHA256s
establish local reproducibility, not independent publisher authentication.


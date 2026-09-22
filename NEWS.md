# csdm 2.0.0.9000

This development version prepares the next major CRAN release. Several
corrections affect coefficient estimates, standard errors, and dependence-test
statistics. Analyses produced with `csdm` 1.0.1 should be re-estimated before
results are compared or reported.

## Estimation and inference

- Add `fullsample = TRUE` for CCE-based models. Cross-sectional averages are
  then calculated variable by variable from all finite observations in the
  selected sample before dynamic lag trimming.
- Construct model and cross-sectional-average lags by explicit time-grid
  matching, preventing lags from crossing gaps in a panel.
- Evaluate formulas, transformed cross-sectional-average variables, subsets,
  and missing-value policies consistently.
- Correct missing-value and leave-one-out cross-sectional averages.
- Use one eligible unit sample for mean-group coefficients and covariance, and
  require at least two eligible units.
- Correct fixed-weight mean-group covariance scaling and HC0--HC3 sandwich
  covariance calculations.
- Check structural identification after cross-sectional-average projection and
  report excluded units.
- Align CS-ARDL parameter components and covariance, and expose unit-level AR
  stability and long-run-ratio eligibility.
- Validate panel keys, retain original observation identity when sorting, and
  preserve named cross-sectional-average lag specifications.

## Dependence diagnostics

- Remove time periods containing no estimated residuals before assessing CD
  sample balance, while retaining the selected policy for partially observed
  periods.
- Implement the paper-defined pooled-variance CDw statistic and the
  correlation-scale screening term used by CDw+ for balanced samples.
- Correct CD-star unit-specific residual scales and validate PCA rank.
- Make randomized diagnostics opt-in, preserve seeded RNG state, and use
  pairwise samples for the classical CD statistic.
- Validate clustered covariance inputs and make residual transformations
  explicit.

## R interfaces

- Add standard extraction and update methods, original-row fitted and residual
  outputs, and stored-data model updates.
- Add `tidy()`, `glance()`, and `augment()` methods with inference and
  row-alignment checks.
- Document `cross_sectional_avg()` as a supported standalone data utility and
  distinguish it from averages configured by `csdm_csa()`.
- Retain evaluated subset, time-spacing, missing-value, and `pdata.frame` time
  information when models are updated.
- Reject unsupported model, trend, cross-sectional-average, long-run, and
  covariance specifications instead of silently storing or ignoring them.

## Deprecations

- Deprecate `csdm_pooled()` because pooled restrictions are not implemented.
- Deprecate `get_residuals()` in favor of the standard `residuals()` method.
- Deprecate `prepare_cd_input()` because its transformations are not used by
  `cd_test()` and can change the hypothesis represented by the residuals.
- Deprecate the exported `cluster_vcov()`, `sandwich_vcov()`, and
  `pooled_vcov()` utilities. Existing calls continue to work with a migration
  warning.

## Documentation and maintenance

- Rewrite the README, introductory vignette, and pkgdown navigation around the
  implemented API. Correct the description of `log_ngd`, update the CDw/CDw+
  explanations, and make the introductory examples runnable.
- Add reference fixtures, statistical validation scripts, platform checks, and
  Codecov reporting.
- Share unit-regression and sample-bookkeeping code across the MG, CCE, and DCCE
  engines, and remove obsolete internal helpers and imports.

# csdm 1.0.1

## Documentation and reference enhancements

- Add references for the implemented estimators and methods, including key
  papers and textbooks.
- Improve function documentation, estimator descriptions, assumptions, and
  documentation consistency.
- Add a link for reporting bugs.

# csdm 1.0.0

## Initial CRAN release

### Estimators

- Mean Group (MG)
- Common Correlated Effects (CCE)
- Dynamic CCE (DCCE)
- Cross-Sectionally Augmented ARDL (CS-ARDL)

### Inference and diagnostics

- Cross-sectional dependence (CD) tests
- Summary and printing methods

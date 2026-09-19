# Development

- Add a website-only CS-ARDL replication article with verified acquisition of
  the externally licensed CMPR data. The dataset is not distributed in the
  package or source repository.

- Implement `fullsample = TRUE` for CCE-based models. Cross-sectional averages
  are then calculated variable by variable from all finite observations in the
  selected sample before dynamic lag trimming.

- Rewrite the README, introductory vignette, and pkgdown navigation around the
  implemented API. Correct the bundled-data description of `log_ngd`, replace
  stale CDw/CDw+ explanations, and make all introductory examples runnable.

- Limit model, trend, CSA, long-run, and variance-covariance specifications to
  choices that are implemented. Unsupported names now fail at the relevant
  public entry point instead of being stored for later rejection.

- Deprecate `csdm_pooled()` because pooled restrictions are not implemented.
  Default fits now create their empty pooled metadata internally without warning.

- Remove an unused internal model-matrix helper and obsolete `stats` imports.

- Document `cross_sectional_avg()` as a supported standalone data utility and
  distinguish it from the model-term averages configured by `csdm_csa()`.

- Deprecate `get_residuals()` in favor of the standard `residuals()` method.
  `cd_test()` now uses an unexported residual-matrix adapter internally.

- Deprecate `prepare_cd_input()`. Its transformations are not used by
  `cd_test()` and can change the hypothesis represented by the residuals.

- Deprecate the exported `cluster_vcov()`, `sandwich_vcov()`, and
  `pooled_vcov()` matrix utilities. They are not used by `csdm()` estimators;
  existing calls continue to work with a migration warning.

- Remove time periods containing no estimated residuals before CD sample balance
  is assessed, while retaining the selected policy for partially observed periods.

- Model updates retain evaluated subset, time spacing, missing-value policy, and pdata.frame time indexes.

- Add tidy, glance, and augment methods with inference and row-alignment checks.

- Add explicit R model accessors, original-row outputs, and stored-data updates.

- Validate clustered covariance inputs and make residual transformations explicit.

- Implement paper-defined pooled-variance CDw and correlation-scale CDw+ screening for balanced samples.

- Correct CD-star unit-specific residual scales and validate PCA rank.

- Make randomized diagnostics opt-in, preserve seeded RNG state, and use pairwise classical CD samples.

- Align CS-ARDL parameter components and covariance; expose AR stability and ratio eligibility.

- Correct variance-of-mean scaling for fixed-weight MG covariance and document its assumptions.

- Correct HC0-HC3 sandwich meat and avoid dense hat matrices.

- Use one eligible unit sample for MG means and covariance; require at least two units.

- Check structural identification after CSA projection and report excluded units.

- Evaluate model formulas, transformed CSA variables, subset, and missing-value policies consistently.

- Share unit regression bookkeeping across MG, CCE, and DCCE engines.

- Correct missing-value and leave-one-out cross-sectional averages.

- Construct model and CSA lags by explicit time-grid matching.

- Validate panel keys and retain original observation identity when sorting.

- Preserve named CSA lag specifications and reject invalid lag orders.
- Reject unsupported estimation options instead of silently ignoring them.

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


# csdm 1.0.1

## Documentation and References enhancement

### References
- Added references for the implemented estimators and methods, including key papers and textbooks in the field

### Documentation
- Improved documentation for all functions, including detailed descriptions of the estimators, their assumptions, and
- Ensured consistency in the documentation style across different documents
- Added link for reporting bugs

# Development

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

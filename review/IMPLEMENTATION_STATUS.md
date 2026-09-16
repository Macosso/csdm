# Implementation status

Date: 2026-09-14. Branch: `codex/correctness-and-r-interface`.
Baseline: `78e661e18bc3ea0a3942d543b9cc71c4713c2942`.

The approved implementation milestones 1-19 are implemented. Milestone 20's local
Windows validation is complete; remote Windows/Linux CI and actual Stata execution
remain outstanding. The branch contains 23 focused commits, including three
additional fixes found during final validation. No version bump, push, or CRAN
submission was performed.

## Delivered

- Strict panel/specification validation, original observation identity, time-key
  lags, missing-value/leave-one-out CSA fixes, evaluated formulas, subset and NA policies.
- Shared unit fitting, structural identification after nuisance projection,
  common-unit MG estimation/covariance, and corrected HC/fixed-weight covariance.
- Aligned CS-ARDL components and covariance, ratio sample counts, AR-root metadata,
  and explicit levels-versus-ECM interpretation.
- Deterministic fitting, explicit diagnostic samples/RNG behavior, corrected
  CD-star scales, paper-defined balanced CDw/CDw+, and cluster/input validation.
- Standard R extraction and update methods, original-row outputs, and
  tidy/glance/augment integration.
- Migration and assumption documentation, executable vignette, and Windows/Linux
  CI configuration. Removed unused MASS dependency; generics and tibble support
  the tidy interface. numops was evaluated but not added: these fixes require
  identification checks and QR operations, not numerical clipping.

## Validation evidence

Environment: Windows, R 4.5.3, roxygen2 8.1.0. Dependency versions for the external
CCE comparison are retained in `plm-reference.txt`.

- Full testthat suite: passed. Includes independent lm references for MG/CCE and
  ARDL ratios, sandwich HC0-HC3 and clustered covariance comparisons, sample and
  lag invariants, formula/identification edge cases, and tidy/lmtest/modelsummary
  integration.
- `R CMD build --no-manual .`: passed, including vignette creation.
- `R CMD check --no-manual --output=.dev-check csdm_1.0.1.tar.gz`:
  **Status: OK**, zero errors, warnings, and notes. Includes installed tests,
  isolated namespace loading, examples, and vignette rebuilding.
- `review/validate-plm.R`: CCE slopes differ from plm::pcce by at most
  3.33e-16; covariance differs by at most 9.11e-18 (tolerance 1e-8).
- `review/validate-diagnostics.R`: 500 replications for each of four N/T
  settings under independent heterogeneous Gaussian errors and common-shock
  alternatives (seed 9142026). CDw null rejection 4.0-4.6%; CDw+ 4.4-5.0%.
  Common-shock rejection: CDw 47.6-74.4%, CDw+ 100%. Exact binomial intervals
  are retained in `diagnostic-calibration.csv`. These are selected designs,
  not a general calibration guarantee.
- `review/validate-cdstar.R`: 500 replications at N=50, T=100, seed 9142027.
  Heterogeneous nonproportional loadings: 5.0% null rejection (95% binomial
  interval 3.26-7.29%). Proportional loading/error-scale design: 100% rejection.
  Both results are retained in `cdstar-calibration.csv`. The latter has a
  zero limiting correction after standardization and diagnoses a substantive
  limitation; it is not counted as a passing calibration case.

The CD-star restriction follows the positive correction condition in equations
(6), (19), and (21) of [Pesaran and Xie](https://arxiv.org/pdf/2109.00408).
Finite-sample numerical rank checks cannot certify that asymptotic condition.
Documentation explicitly warns about near-zero corrections. A practical
screening rule would require further theoretical/statistical design; none was
invented to make the simulation pass.

## Check environment and reproducibility

Use an isolated development library in `.dev-library` plus the normal user R
library. Disable the repository's renv startup during package checks with
`R_PROFILE_USER=NUL` and `R_ENVIRON_USER=NUL` on Windows. Use `LC_ALL=C` and
`LANG=C` for this host's locale, and point `RSTUDIO_PANDOC` to the installed
Pandoc directory. Set `R_USER_CONFIG_DIR` to `.dev-check/config` so optional
reporting-package setup remains inside the workspace. Run `review/test.R`
with Rscript --vanilla for source tests; installed tests also run in R CMD check.
The final check log is copied to `review/windows-check.log`.

Baseline attempts exposed host startup/locale/library-path problems. Final checks
also caught missing generic imports, now fixed. These issues were resolved rather
than suppressed. The existing empty `.git/AUTO_MERGE.lock` was left untouched;
Git emitted a warning but created and verified every commit successfully.

## Remaining release evidence and scope

- No Stata executable was found. Obtain version-pinned Stata outputs/access before
  claiming parity. The diagnostic definitions intentionally follow primary
  equations where the reviewed Stata code differs.
- The GitHub Actions Windows/Linux matrix is configured, not executed locally
  or remotely. Remote execution requires a later push/PR.
- This was a full Windows source-package check with --no-manual, not --as-cran
  and not a PDF manual build. CRAN submission checks remain a release task.
- Near-degenerate CD-star correction remains a documented inference limitation.
  Broader error distributions, serial correlation, factor strengths and N/T ratios
  warrant further simulations before making broad diagnostic claims.
- New estimators, dynamic bias correction, pooled restrictions, forecasting, and
  unsupported covariance schemes remain deferred as approved.

Recommendation remains to maintain csdm as a compact, auditable implementation,
with these limitations and reference fixtures. The repairs justify continued
development; they do not establish full xtdcce2 coverage or remove the need for
independent release validation.

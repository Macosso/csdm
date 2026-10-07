# CD enhancement validation

Date: 2026-10-07. Branch: `feature/cd_enhancement`.
Baseline: `e562783` (`master`, csdm 2.0.0).
Environment: Windows, R 4.5.3, roxygen2 8.1.0.

## Confirmed interface

- Select numeric variables explicitly through `...`, using bare names, quoted
  names, or character vectors. Data methods require named controls after `...`.
- Plain data frames require distinct `id` and `time` column names. Indexed panel
  data frames use their stored indexes, including with `drop.index = TRUE`.
- Return one named `cd_test` result per variable in a `cd_test_list`, also for a
  single selection. Print a combined table of tests and sample dimensions.
- Default `reps = 1L`. For multiple draws, CDw uses the sum of its statistics
  divided by `sqrt(reps)`, following Juodis and Reese's equation (33). CDw+ adds
  the screening term once to that aggregate. Seeded calls restore RNG state.
- Leave the package version at 2.0.0; describe the changes in the development
  section of NEWS until the release version is chosen.

## Checks

- The full source testthat suite passed. Additional sample-exclusion checks
  passed after the final metadata change.
- Independent cross-product references verify equation (33), the single-draw
  calculation, normal p-values, and the CDw+ screening term.
- Data tests cover selection syntax, shuffled rows, structural gaps, independent
  missing samples, invalid keys, labeled/date time indexes, removed pdata.frame
  index columns, result printing, and RNG behavior.
- `R CMD build --no-manual .` passed, including vignette creation.
- `R CMD check --no-manual --output=.dev-check/cd-enhancement
  csdm_2.0.0.tar.gz` finished with **Status: OK**, with zero errors, warnings,
  and notes. All 465 installed-test assertions passed. Spelling, examples, S3 consistency, documentation,
  isolated namespace loading, and vignette rebuilding all passed.

The full package check used the normal execution context because the restricted
vignette subprocess encountered Windows Application Control DLL-loading blocks.
No Windows security settings or dependency installations were changed.

For this host, disable renv startup with `R_PROFILE_USER=NUL` and
`R_ENVIRON_USER=NUL`, set `LANG=C` and `LC_ALL=C`, include `.dev-library` in
`R_LIBS`, and set `R_USER_CONFIG_DIR` to `.dev-check/config`. Set
`RSTUDIO_PANDOC` to the existing directory
`C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools`.
Build/check logs are retained locally under `.dev-check`.

## Statistical smoke checks

Run `Rscript --vanilla review/validate-cd-reps.R` from the repository root.
The script uses seed 10072026 and 500 simulations for each combination of
panel dimensions, null/alternative, and weight draw count.

| N | T | Draws | CDw null rejection | CDw+ null rejection |
|---|---|---|---|---|
| 50 | 100 | 1 | 6.0% | 6.2% |
| 50 | 100 | 30 | 5.4% | 5.6% |
| 100 | 50 | 1 | 5.6% | 5.6% |
| 100 | 50 | 30 | 4.4% | 4.6% |

All null designs use independent Gaussian observations with heterogeneous unit
scales. Under an added common Gaussian shock, CDw+ rejected in all 500
simulations for each setting. Exact binomial confidence intervals and the CDw
alternative results are retained in `cd-reps-calibration.csv`.
The script passed its predeclared broad size and CDw+ power bounds. These checks
cover selected designs and do not establish calibration for every data process.

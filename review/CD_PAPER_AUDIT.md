# CD implementation and documentation audit

Date: 2026-10-07. Branch: `feature/cd_enhancement`.
Reviewed baseline: `aee1f29`. Package version: 2.0.0, development changes.

## Sources and equation mapping

| Test | Primary reference | Finding |
|---|---|---|
| Classical CD | [Pesaran's original diagnostic paper](https://docs.iza.org/dp1240.pdf), section 9, equation (31); [weak-dependence paper](https://docs.iza.org/dp6432.pdf) | The overlap-specific correlations, square-root overlap factors, and retained-unit normalization agree. |
| CDw | [Juodis-Reese, final published version](https://pure.uva.nl/ws/files/116332884/The_Incidental_Parameters_Problem_in_Testing_for_Remaining_Cross_Section_Correlation.pdf), p. 1197, equation (30) | The weighted covariance numerator and inverse pooled variance agree. |
| CDw+ | Same paper, equation (32) | The absolute-correlation screen and threshold agree; Fan-Liao-Yao provide the general enhancement principle. |
| Repeated CDw | Same paper, p. 1198, equation (33) | Sum independent draws and divide by the square root of the draw count. |
| CD* | [Pesaran-Xie working paper](https://www.ifo.de/DocDL/cesifo1_wp9234.pdf); [current author manuscript](https://arxiv.org/pdf/2109.00408), equations (12), (28)-(31), (38)-(40), Remark 6 | Correction algebra agrees. Standardization before PCA is a package variant. Fitted input needed correction. |

The published Juodis-Reese PDF was inspected visually. Its equation (29) is
the unscaled covariance expression; equation (30) includes the inverse pooled
variance. The previous vignette already contained that factor. Its compact pair
sum and different residual symbol were equivalent, rather than a different
statistic. The revised vignette uses explicit index bounds and identifies the
equation, residual demeaning, and variance separately. Since Rademacher weights
have square one, this variance is the mean squared demeaned input. Replacing it
with pair-specific standard deviations would implement a different statistic.

Repeated CDw+ adds one screen to the repeated CDw aggregate. This is the
explicitly agreed package extension combining equations (32) and (33).

## CD* decisions and corrected behavior

Two distinct choices were discussed and confirmed with the maintainer:

1. Preserve standardization before PCA and document it as an implementation
   variant. The reference PCA estimates loadings from unstandardized data;
   standardization can change the factor space for heterogeneous unit variances.
2. For fitted models, start CD* with the response minus the fitted economic and
   deterministic component, retaining the common-factor component. The previous
   code used full CCE residuals, which is a different input from the regression
   procedure described in the reference.

The implementation reconstructs that partial residual from the stored economic
design, response, unit coefficients, and original observation map. Economic
lags, intercepts, and trends are subtracted; CSA contributions are retained.
Only estimated observations are used. CD, CDw, and CDw+ continue to test full
fitted residuals. In `type = "all"`, CD* is computed separately, with its own
sample dimensions and exclusions recorded under `tests$CDstar`. Passing the
full residual matrix explicitly to CD* retains its literal matrix-input meaning.

The loading normalization is Gamma'Gamma/N = I. Factor-filtered residual
scales are unit-specific root mean squares. The denominator is the mean of
squared loading/scale adjustment terms, and the numerator adds the associated
square-root-time bias correction. The documentation now gives the complete
calculation. Zero requested PCs reduce to classical CD on the selected input.

The partial-residual construction aligns with the regression input in the
reference; it does not make the retained standardized-PCA variant an exact
reproduction of the paper, nor establish the paper's assumptions for every
MG, dynamic CCE, or CS-ARDL specification.

## Interpretation and references

All p-values use two-sided standard-normal approximations. There is no
implemented serial-correlation adjustment. Stationarity, factor strength,
residual properties, and panel asymptotics must be assessed for the application.
Ordinary CD can be biased after fitting common time parameters. Matrix/data
acceptance is an interface feature, not evidence that a particular series meets
the asymptotic conditions. For classical CD, the minimum overlap of two is
computational; skipped pairs or very short overlaps can affect calibration.

Corrected references:

- Pesaran-Xie's coauthor is Yimeng Xie. The 2021 source is a working paper,
  rather than the previously listed Econometric Reviews article. Retained that
  working-paper reference and added the [2026 Econometric Theory publication](https://doi.org/10.1017/S0266466625100212).
- Fan-Liao-Yao's third author is Jiawei Yao, as shown in the
  [published paper hosted by Fan](https://fan.princeton.edu/document/591).
- Added the verified Juodis-Reese and Fan-Liao-Yao DOIs.

## Verification

- Literal pair-overlap and balanced-panel references for classical CD.
- Literal triple-sum reference for CDw equation (30), including heterogeneous
  scales that distinguish it from a weighted-correlation statistic.
- Existing independent cross-product references for equation (33), normal
  p-values, the one-draw result, and adding the screen once.
- Independent eigenvector/loading normalization for two-factor CD*, applied to
  the agreed standardized input; verified invariance to positive unit scaling.
- Fitted partial residuals checked against direct response-minus-economic-term
  calculations and independent augmented regressions. Covered static CCE,
  shuffled data, missing values, sample subsets, dynamic terms, CS-ARDL, trends,
  MG, and numeric factor time labels in `pdata.frame` fits; other fitted test
  statistics retain their full-residual references.

These algebra checks verify computations; they are not universal size/power
guarantees for the retained PCA variant or arbitrary data processes.

Run `Rscript --vanilla review/validate-cdstar-fit.R` for the additional fitted
null simulation. It uses seed 10072027, N=50, T=100, 500 replications, one
Gaussian common factor, heterogeneous loadings/scales, a factor-dependent
regressor, and independent Gaussian idiosyncratic errors. The correctly
specified CCE fit is followed by one-PC CD* with standardized PCA.

| CD* input | Rejection at 5% | Exact 95% binomial interval |
|---|---|---|
| Corrected partial residuals | 5.8% | 3.92%-8.22% |
| Previous full CCE residuals | 100% | 99.26%-100% |

The partial-residual result passed the predeclared broad 12% smoke bound.
The full-residual result is a comparison, not an acceptance case. This single
design supports the input correction but does not establish general calibration
of the standardized-PCA variant. Results are in `cdstar-fit-calibration.csv`.

Final validation on Windows 11 with R 4.5.3:

- Source test suite passed, followed by a successful source-package build
  including the vignette.
- `R CMD check --no-manual` returned `Status: OK`: no check errors, warnings,
  or notes. Installed-package tests reported 502 passes, zero failures,
  zero warnings, and zero skips.
- Examples, spelling, S3 consistency, documentation, and rebuilding the
  vignette all passed. The generated HTML contains the corrected equation
  references, standardized-PCA qualification, partial-residual explanation,
  and publication DOIs.
- The dependency check printed a temporary CRAN index download failure but
  completed successfully using installed dependencies. This local check does
  not verify hosted GitHub Actions startup.

Local check logs are under `.dev-check/cd-paper-audit/` (ignored by Git).

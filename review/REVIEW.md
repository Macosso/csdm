# csdm implementation review

Reviewed 2026-09-14. Checkout: `78e661e18bc3ea0a3942d543b9cc71c4713c2942`, DESCRIPTION version 1.0.1.

**Recommendation: retain the package only with a focused correctness and interoperability programme. Do not retire it merely because another implementation exists; do not expand its estimator menu before repairing inference and sample handling.** Its clean-data MG/CCE/DCCE/CS-ARDL core is useful, but the current release has reproducible defects affecting estimates, uncertainty, and diagnostics.

Scope: all R source files, exported interface, documentation, existing tests and CI were inspected. `review/reproduce.R` sources this checkout and executes independent regression comparisons, edge cases, and a 200-replication diagnostic experiment under R 4.5.3. Results are in `review/reproduce-output.txt`. Package source was not changed. `testthat`, `Rdpack`, and `plm` are unavailable in this R library; I did not run the existing testthat suite, a full R CMD check, or Stata itself. Source-level comparisons with Stata are not a claim of end-to-end Stata replication. Competing packages were assessed from documentation, not audited for numerical correctness.

## What is working

- On the balanced, nonsingular reference panel, MG coefficients agree exactly with averages of separately fitted `lm()` coefficients. Its internal covariance agrees with `cov(B)/N` to floating-point precision.
- Static CCE agrees with independently constructed contemporaneous-CSA regressions.
- Dynamic augmented ARDL coefficients and unit-first long-run aggregation agree with independent regressions to approximately 1e-16 on the clean example.
- The adjustment sign and long-run algebra are correct: `adjustment_i = sum(phi_i)-1`; `theta_i = sum(beta_i)/(1-sum(phi_i))`. Averaging unit ratios is appropriate for the heterogeneous mean-group target; replacing this with a ratio of averages would change the estimand.
- The separated specification constructors and stable named coefficient/residual matrices provide a reasonable foundation. Documentation, examples, a test suite, and R CMD check CI already exist.

These are arithmetic/implementation checks, not demonstrations of consistency, coverage, or identification in arbitrary applications.

## Confirmed findings, ordered by impact

### 1. P1: CCE can report a structural slope that is not identified after augmentation

At [csdm_fit_engines.R:221](C:/Users/Joaoc/Documents/csdm/R/csdm_fit_engines.R:221), and the analogous DCCE regression, the structural regressors precede the CSA controls in an ordinary `lm()`. If a structural regressor lies in the CSA span, QR aliasing can discard the CSA column and retain the structural coefficient. The engine checks neither the rank of residualized structural regressors nor which nuisance controls were discarded.

Reproduction: set `x_it = sin(t)` for every unit. Then `x_it = mean_i(x_it)` exactly, yet CCE reports an `x` estimate of -0.03786 and SE 0.10923. The structural slope cannot be separated from the unrestricted unit loading on that identical CSA. This is distinct from harmless redundancy between two nuisance CSA columns.

Fix: establish a basis for the nuisance space first, project economic regressors and the response on its complement, and test structural identification there. Record aliased economic terms and enforce a documented exclusion policy. Also check usable observations, rank, and positive residual degrees of freedom before including each unit.

### 2. P1: Dynamic lags cross missing calendar periods

[csdm_fit_engines.R:402](C:/Users/Joaoc/Documents/csdm/R/csdm_fit_engines.R:402) and line 423 shift sorted observations, not time indexes. The CSA lag helper likewise shifts available time rows.

Reproduction: remove unit 1 at time 20. At time 21, the fitted DCCE regression uses time 19 as lag 1; it returns a finite residual when a calendar lag would be missing. Removing an observation and retaining an explicit missing row can therefore define different dynamics.

Fix: define a panel frequency/time-step contract and construct lags by `(id, time-k*delta)` matching. Validate duplicate, missing, and nonfinite indexes. Support standard panel/calendar representations or reject unsupported ones explicitly. Do not infer that every numeric sequence has step one.

### 3. P1: Coefficient means and covariance use different unit samples

[csdm_fit_engines.R:95](C:/Users/Joaoc/Documents/csdm/R/csdm_fit_engines.R:95) uses columnwise `na.rm=TRUE`; [csdm_internal_helpers.R:100](C:/Users/Joaoc/Documents/csdm/R/csdm_internal_helpers.R:100) removes any row with any missing coefficient before computing covariance.

Reproduction: make `z` constant for one of eight units. Its `x` estimate contributes to the reported mean, but not its variance. The reported `x` variance is 0.002058 versus 0.001559 for the variance of the eight available `x` estimates divided by eight. No unit is recorded as dropped. If `z` is unidentified everywhere, finite coefficients can remain while the whole covariance becomes NA and `z` becomes NaN.

Fix: choose a common-unit estimand by default or explicitly implement coefficient-specific inclusion and compatible cross-covariances. Report effective N by term. Pairwise covariance substitution alone is not a complete fix and can compromise positive semidefiniteness.

### 4. P1: CDw+ is severely miscalibrated on a basic null experiment

[utils_cd.R:202](C:/Users/Joaoc/Documents/csdm/R/utils_cd.R:202) compares `sqrt(T_ij)*abs(rho_ij)` with a threshold already scaled by `1/sqrt(T)`, then adds a large positive sum. In a balanced panel this effectively thresholds raw correlations at order `sqrt(log(N))/T`, not the documented correlation-scale threshold. Missing correlations can also contaminate the sum.

Reproduction: 200 independent Gaussian panels, N=20, T=100, seed 27182. At nominal 5%, CDw rejects 6.5%; CDw+ rejects **100%**. This is one finite-sample design, not a general Monte Carlo study, but it is sufficient to reject the current implementation as a trustworthy default diagnostic.

Fix: derive threshold and enhancement on consistent scales from the chosen paper/version, then test size and power across N/T regimes. Merely removing one square root is not adequate validation. The current [xtcd2 source](https://raw.githubusercontent.com/JanDitzen/xtdcce2/master/main/xtcd2.ado) also deserves scrutiny here: matching a reference implementation does not establish theoretical correctness.

### 5. P1: CD* collapses unit-specific residual scales to one scalar

[utils_cd.R:421](C:/Users/Joaoc/Documents/csdm/R/utils_cd.R:421) computes `sqrt(mean(res_defac^2))`. In the corresponding [xtcd2 implementation](https://raw.githubusercontent.com/JanDitzen/xtdcce2/master/main/xtcd2.ado), Mata `mean(res:^2)` produces a row vector of column means. This is established by the [Mata mean documentation](https://www.stata.com/manuals/m-5mean.pdf). R's `mean(matrix)` instead returns one scalar.

That scalar then cancels in the loading correction; the intended heterogeneity in residual scales is lost. Reproduction with one common factor: current statistic -2.83196, versus -1.36321 when the same calculation retains unit-specific residual scales. This changes the 5% conclusion. This reference calculation is an R translation of the relevant algebra, not a Stata execution.

Fix: preserve the unit scale vector with explicit columnwise operations, validate the full method against the paper and versioned Stata results, and constrain `n_pc` to feasible rank. Currently a 4-by-3 matrix with default `n_pc=4` fails with an obscure error. Unbalanced handling is also more limited than xtcd2.

### 6. P1: Exported HC covariance uses residuals instead of squared residuals

[utils_vcov.R:116](C:/Users/Joaoc/Documents/csdm/R/utils_vcov.R:116) calculates `crossprod(X, w * X)` although `w` contains signed, adjusted residuals. The sandwich meat requires squared adjusted residuals: for example `crossprod(X * w)`.

Reproduction: four intercept rows and residual vector (-1,-1,-1,-1) return variance -0.25 instead of +0.25. Even valid intercept-only OLS residuals (-1,1,-1,1) give zero instead of +0.25. HC0 through HC3 are affected.

This exported utility is not called by the main MG fitting path; it does **not** imply that every default `csdm()` standard error has this bug. Avoid allocating the full n-by-n hat matrix when computing leverage, and validate dimensions/rank.

### 7. P1: Exported `pooled_vcov()` returns dispersion rather than variance of the mean

[utils_vcov.R:184](C:/Users/Joaoc/Documents/csdm/R/utils_vcov.R:184) returns a weighted population covariance of unit coefficients without the required variance-of-mean scaling/correction. `pairwise` is accepted but unused.

Reproduction: coefficients (1,2,3,4) return 1.25; equal-weight MG variance is `var(1:4)/4 = 0.416667`. Unequal weights need an explicitly justified estimator; weighted coefficient dispersion is not automatically the uncertainty of their weighted average.

This is another exported utility defect separate from `.csdm_mg_vcov()`, whose complete-data formula is correct.

### 8. P1: CS-ARDL breaks the coefficient/covariance contract

[csdm_methods.R:371](C:/Users/Joaoc/Documents/csdm/R/csdm_methods.R:371) appends adjustment and long-run parameters to `coef()`, while `vcov()` at line 386 returns only the levels ARDL covariance. In the reproduction there are nine coefficients and a six-by-six covariance. `confint()` returns NA intervals for all three appended parameters despite their SEs appearing in the custom summary.

Fix: provide consistent parameter components in both methods. Either use explicit `component='levels'/'long_run'/'all'` with matching covariance or expose transformed effects through a separate accessor. Estimate the full covariance of the unit-level transformed vectors, including cross-covariances. The joint vector contains exact linear dependencies such as adjustment versus AR coefficients, so redundant joint Wald tests require explicit handling.

### 9. P1/P2: Accepted options silently do nothing

The fitting engines never dispatch on `vcov$type`: `mg`, `np`, `nw`, `wpn`, and `ols` give identical covariance. CSA `scope` and `cluster` do not alter aggregation. `pooled`, `fullsample`, and `mgmissing` are only stored. MG/CCE receive `lr` but ignore it. `subset`, `weights`, and `na.action` supplied through `...` are silently ignored.

Several arguments are honestly documented as reserved/stubs. The defect is accepting nondefault requests without rejecting them: a stored request is easily mistaken for a fitted restriction or covariance choice. Reproductions confirm ignored covariance, subset, weights, and MG lag requests.

Fix: reject unsupported nondefault settings and unused `...` until implemented. Do not invent estimator meanings for labels such as `nw` without a documented definition.

### 10. P2: Named CSA lag specifications are destroyed

[csdm_specs.R:48](C:/Users/Joaoc/Documents/csdm/R/csdm_specs.R:48) removes names with `as.integer()`, then subsets using the now-absent names. `csdm_csa(lags=c(y=2,x=1))$lags` becomes `integer(0)` and DCCE errors. Preserve names, reject duplicate/unknown names, define partial-name defaults, and reject fractional lags rather than silently truncating them. Clarify whether a length-one named lag applies to just that variable.

### 11. P2: Ordinary R formula behavior is broken

- `I(y^2) ~ x + z` fits internally but fails in R-squared construction: the engine requires a simple response column name there. Use the evaluated model response.
- `y ~ .` includes generated `.csdm_rowid__` as a regressor. Formula expansion must occur before bookkeeping columns are introduced.
- Default CSA selection uses raw `all.vars(formula)`. With transformed regressors it averages underlying columns rather than the transformed model variables. These are different factor proxies; neither should be substituted silently. Factors/interactions and non-syntactic names need a defined model-matrix policy.
- Dynamic formula restrictions are partly documented and explicitly rejected; they are limitations rather than all being silent bugs.

### 12. P2: Duplicate indexes produce inconsistent estimation and reporting

[csdm_internal_panel.R:3](C:/Users/Joaoc/Documents/csdm/R/csdm_internal_panel.R:3) does not reject duplicate id/time cells. Both rows enter regressions/averages, but residual and fitted matrices overwrite a cell. Reproduction: 321 input observations lead to 320 reported observations. Validate keys before fitting; distinguish observed, eligible, and estimated counts.

### 13. P2: CSA missing-value behavior contradicts its contract

[utils_avg.R:179](C:/Users/Joaoc/Documents/csdm/R/utils_avg.R:179) multiplies a missing own value by zero, which remains NA. For a time slice (NA,2,4), the leave-one-out mean for the first unit is NA rather than 3. Separately `na.rm=FALSE` still returns 3 for the ordinary mean although documented behavior is NA. Fix finite-value arithmetic before multiplication and honor the missing-value policy.

### 14. P2: Diagnostic defaults lose samples and alter the RNG

`cd_test(type='CD')` defaults to dropping every time with any missing residual before calculating pairwise correlations. This contradicts the documented claim that classical tests always use available pairs. One entirely dropped unit can make the default test impossible for all remaining units.

Every fit also invokes all randomized and PCA diagnostics in [csdm_internal_helpers.R:274](C:/Users/Joaoc/Documents/csdm/R/csdm_internal_helpers.R:274), even though only classic CD is printed. Fitting consumes `.Random.seed`; repeated fits can have different stored diagnostics and alter subsequent simulation/bootstrap draws. Errors/warnings are often suppressed.

Fix: diagnostics should be opt-in or lazily calculated; preserve RNG state for explicitly seeded randomized tests. Separate pairwise CD handling from balanced-PCA requirements and report actual diagnostic sample sizes and exclusions.

## Theoretical cautions and interpretation gaps

1. **CCE is not automatic endogeneity correction.** Its factor-proxy assumptions, structural identification, exogeneity conditions, and sufficient N/T must be stated clearly. It does not resolve arbitrary correlation between regressors and idiosyncratic shocks.
2. **DCCE needs a deliberate lag design.** Default `model='dcce'` adds neither a dependent-variable lag nor positive CSA lags. A name alone does not establish dynamic specification. Explain the role of growing CSA truncation orders and offer a documented heuristic/sensitivity workflow; no heuristic should be sold as a universal optimum. See the [xtdcce2 model documentation](https://janditzen.github.io/xtdcce2/).
3. **Short-run reporting is ambiguous.** The CS-ARDL summary labels untransformed levels coefficients “Short Run Est.” These are not the full ECM differenced coefficients for general ARDL(p,q). Label the table “levels ARDL coefficients”, or compute the actual ECM transformation and covariance.
4. **Long-run existence is not checked.** A denominator cutoff of 1e-8 prevents division by effectively zero but does not establish stability, cointegration, or well-behaved ratio inference. Report AR roots, denominator distributions, effective N, influential units, and sensitivity. Negative adjustment alone does not prove higher-order AR stability. Do not automatically exclude unstable units without explaining the resulting target population.
5. **Adjustment need not be dropped solely because its long-run denominator is zero.** Zero adjustment is meaningful even when a ratio is undefined. Separate eligibility rules for adjustment and long-run ratios.
6. **Absence of per-unit delta-method SEs does not by itself invalidate MG long-run SEs.** Cross-unit dispersion of correctly transformed unit estimates is a legitimate large-N approach under its assumptions. Missing joint covariance, unstable ratios, and absent finite-sample validation are the concrete problems.
7. **Classic CD after estimated common effects needs qualification.** Factor estimation can distort its null distribution. The [Pesaran–Xie paper](https://arxiv.org/abs/2109.00408) and weighted/bias-corrected diagnostics address this issue. Cross-sectionally demeaning residuals in `prepare_cd_input()` changes the tested object and can mechanically introduce dependence; it should not be recommended as generic harmless preprocessing.
8. **R-squared terminology and validation need tightening.** `R2_mg` uses a residual-variance adjustment, while tests contain older pooled-SSE comments and skip directly asserting the current R2_mg value. Distinguish ordinary within-unit, averaged, and adjusted CCE measures. Intercept-free formulas need their own definition. Do not claim exact Stata R2 equivalence without fixtures.

## Making this an R modelling package

The central change is an explicit fitted-object contract, not simply registering more S3 methods.

| Priority | Addition | Reason |
|---|---|---|
| Immediate | Matching `coef()` / `vcov()` components; working `confint()` | Reliable extraction and downstream inference |
| Immediate | `nobs()`, `model.frame()`, `model.matrix()`, explicit `terms()` and `formula()` behavior | Preserve the evaluated estimation problem |
| Immediate | Explicit `fitted()` / `residuals()` with vector, long-data, and matrix formats | Restore original row order and missingness while retaining panel diagnostics |
| Immediate | `na.action`, `subset`, contrasts, factor levels, formula environments and original row identifiers | Normal modelling semantics and reproducible `update()` |
| Next | `broom::tidy()`, `glance()`, `augment()`; integration tests with modelsummary and lmtest | Tables and inference workflows should use public methods |
| Next | `predict(newdata)` with explicit unit/factor information requirements | Conditional fitted values differ from forecasts of future common factors |
| Next | Per-unit coefficient/covariance accessors and sample/alias diagnostics | Heterogeneity is the point of MG estimation |
| Later | Diagnostic plots, lag/CSA sensitivity and long-run influence plots | Help assess model adequacy |

Observed details matter: `formula(fit)` already works through the default; `fitted(fit)` currently works by partial matching `fitted_xb`. Conversely `nobs(fit)` errors, `df.residual(fit)` is NULL, and `model.frame(fit)` returns the string `"mg"` because `$model` stores an estimator label instead of a model frame. Do not call all unregistered methods broken, and do not depend on accidental partial matching.

There is no unique OLS residual df or likelihood automatically appropriate for every heterogeneous estimator. Define meaningful aggregate and per-unit df. Add `logLik`, AIC/BIC, `anova`, sandwich scores/bread, emmeans, and marginaleffects support only where the estimand, likelihood, covariance, and prediction semantics justify them. A fabricated likelihood or ordinary-OLS sandwich interface would be worse than an explicit limitation. `plm::pmg()` means mean-group estimation, not the Pesaran–Shin–Smith pooled-mean-group estimator.

## Features worth adding after repairs

The [xtdcce2 documentation](https://janditzen.github.io/xtdcce2/) describes broader coverage than this package: pooled/partially pooled restrictions, IV, long-run alternatives, bias corrections, regularized CSA and specification diagnostics. Prioritize:

1. A versioned R-versus-Stata validation gallery and machine-readable fixtures. Include unit coefficients, covariance, long-run effects, sample membership, and diagnostics; disclose intentional differences.
2. Half-panel jackknife and recursive mean adjustment for dynamic finite-T bias, with simulation evidence.
3. Per-regressor lag orders, time-aware lag/difference syntax, genuine CSA scope/fullsample handling, and effective-sample summaries.
4. CCEP/partial pooling and CS-DL; a correctly defined direct ECM/PMG implementation if there is sufficient maintenance capacity. These are new estimators, not aliases or inverse-variance averaging shortcuts.
5. CSA lag/set sensitivity, static rank-condition classification, and regularized CCE. Respect methods' static/dynamic applicability.
6. Suitable bootstrap inference, dependence-strength and slope-heterogeneity diagnostics. Reuse established packages where possible rather than implementing every panel test.
7. IV-CCE only after the basic contract and identification diagnostics are reliable.

Implementation maintenance: unify the three largely duplicated fitting loops, construct/evaluate each model frame once, retain sample provenance, collect per-unit outputs before binding, use QR rather than normal-equation inverses, compute leverage without a dense hat matrix, and use thin SVD rather than a full T-by-T eigendecomposition when appropriate. These changes should follow correctness fixtures, not precede them.

## Validation required before the next release

Existing tests predominantly check object shape, names, signs, finite values, and printed output. Some compare residual R-squared against manual calculations, which is useful. They do not currently establish inference calibration or Stata equivalence; the CDw+ test only checks finiteness.

- Convert each reproduction into a regression test with an independent expected result.
- Add balanced/unbalanced, gaps, duplicates, missing units, insufficient T, CSA-spanned regressors, rank deficiency, factors, transformed formulas, and index-class cases.
- Compare HC covariance to an established implementation on true OLS residuals and test invariance to residual sign reversal.
- Verify MG and transformed MG covariance against independent unit-level calculations, including missing coefficient policies.
- Add property checks: row-order invariance; explicit missing row versus absent time consistency; exact coefficient/covariance name agreement; deterministic fitting; `y = fitted + residual` on the estimation sample.
- Test CD size/power across N/T, factor strength, heteroskedasticity, and missingness using fixed independent seeds and Monte Carlo tolerances. Keep substantial simulations separate from fast CRAN checks.
- Run package-installed namespace tests, the existing test suite, examples, R CMD check, and at least Windows/Linux CI. Sourcing files does not verify namespace/dependency packaging.

## Community value and retirement decision

Static MG/CCE are already available in [plm](https://cran.r-project.org/web/packages/plm/refman/plm.html). As of this review, [CRAN dcce 0.4.2](https://cran.r-project.org/package=dcce), published 2026-05-05, advertises MG/CCE/DCCE, CS-ARDL, CS-DL, pooling, bias-related tools/diagnostics, and broom-compatible methods. This is documented overlap, not evidence that those implementations are correct or superior.

Consequently, “CCE and DCCE are available in R” is no longer sufficient differentiation. A credible role for csdm is **a small, transparent, carefully validated implementation of dynamic heterogeneous panels, with excellent R methods and documented Stata comparisons**. Its limited dependencies and readable R code are assets for teaching, auditability, and independent verification. A new estimator theorem is not required for useful econometric software.

My recommendation is a bounded rehabilitation: publish a corrective release with explicit limitations, fix the confirmed statistical defects, establish reference fixtures, then complete the R interface. Decide on expansion only after comparing the corrected package with alternatives on the same designs. No usage/download or maintenance-capacity evidence was available, so I would not infer community impact from CRAN publication alone.

If maintaining those correctness guarantees is not feasible, a managed deprecation or collaboration with another project is preferable: publish the known limitations, preserve reproducibility, validate migration examples, and identify which replacement covers which use case. Abrupt retirement without migration guidance would not help existing users. Leaving demonstrably unreliable diagnostics unchanged would also not be an acceptable maintenance strategy.

# csdm correctness and R interface implementation plan

Status: proposed, awaiting user review. No package implementation changes are authorized until this plan is approved.

Branch: `codex/correctness-and-r-interface`.
Baseline: `78e661e18bc3ea0a3942d543b9cc71c4713c2942`, version 1.0.1.
Evidence: [review](C:/Users/Joaoc/Documents/csdm/review/REVIEW.md), [reproductions](C:/Users/Joaoc/Documents/csdm/review/reproduce.R).

## Scope proposed for approval

Repair confirmed correctness defects, establish theoretically justified inference, and complete the essential R modelling contract. Preserve correct clean-data estimator behavior. Each item below is a separate focused commit, including its relevant tests and documentation; split an item further if review would otherwise be difficult.

New estimators (CS-DL, pooled/partially pooled CCE, PMG, IV-CCE), bias corrections, and broad new diagnostics are a later milestone requiring a separate plan. Unsupported requests will be explicitly rejected in this milestone rather than silently approximated. New-data forecasting is also deferred pending its own statistical/API design.

## Development rules

- First establish an installed-package test baseline; distinguish pre-existing failures from regressions.
- For each fix: reproduce the defect, implement the smallest coherent change, run focused tests, update affected documentation and NEWS, then commit. Do not combine unrelated fixes or unrelated formatting changes.
- Keep functions small and names idiomatic to the existing R package. Prefer explicit matrix dimensions and names; preserve formula environments and panel indexes.
- Use short comments for statistical intent, non-obvious invariants, or references. Explain public behavior in roxygen documentation, not long implementation comments.
- Use QR/SVD-based calculations where appropriate; avoid arbitrary clipping or generalized inverses that conceal nonidentification.
- Investigate `numops` for repeated numerical operations once work is approved. It is optional: verify its actual API, numerical semantics, CRAN availability, and dependency impact before proposing its use. Do not add a dependency solely for an isolated clamp.
- Preserve user work. Stage explicit files or hunks; inspect staged changes before committing. Never use broad staging that captures unrelated files.
- No rewriting existing history, force pushes, remote publication, CRAN submission, or automatic version bump is included.
- If a statistical choice, compatibility decision, or instruction remains unclear, ask the user before implementing that part. Continue only independent, already-approved work.

## Decisions to approve with the plan

These are proposals, not assumptions. They determine implementation and tests.

1. **Missing/aliased economic coefficients:** default to a common set of units that identify every requested economic coefficient, with warnings and recorded reasons for exclusions. Reject a fit with fewer than two eligible units. Keep `mgmissing=TRUE` explicitly unsupported until a separately justified coefficient-specific covariance design is agreed. Distinguish redundant nuisance columns from unidentified economic coefficients.
2. **Time semantics:** introduce an explicit numeric time-step argument, proposed name `time_step`, with default 1 for annual/integer-step panels. Validate that indexes lie on that grid and build lags by time-key matching. Non-unit grids require an explicit value. This milestone supports numeric time and numeric-like `pdata.frame` indexes; Date/year-quarter/year-month support requires a later explicit frequency design. Never compress gaps or silently infer a frequency.
3. **CS-ARDL extraction:** retain the existing default `coef()` output (levels, adjustment, long-run), make `vcov()` match it, and add matching `component = c("all", "levels", "adjustment", "long_run")` arguments. Use a common eligible unit set within each returned component, with explicit counts; if eligibility differs for a combined component, recompute its mean and covariance on the combined common set and make that selection visible. Do not silently mix samples. An alternative is to restrict all components to one common set; resolve this choice during plan review if that is preferred.
4. **Diagnostics:** stop calculating randomized/PCA diagnostics inside ordinary fitting. Keep the existing classic-CD summary field, using explicit pairwise handling on eligible units and reporting its sample. Additional diagnostics run only when requested. Seeded randomized calls preserve the caller's RNG state; unseeded explicit randomized calls use the normal R RNG behavior.
5. **Formula contract:** support ordinary evaluated model responses and design terms for MG/CCE. For default CCE augmentation, use evaluated response/design variables rather than silently substituting underlying raw variables; document generated CSA names and factor/interaction handling. Expand `.` against user data before generated columns exist. Keep dynamic transformed-term restrictions explicit until lagging evaluated terms is separately designed.
6. **Unsupported options and weights:** implement `subset` and `na.action`; reject nondefault unimplemented pooling, CSA scope/fullsample, covariance settings, and estimation weights. Retain working weights in the CSA utility. Do not silently turn weighted CSA calculations into weighted unit regressions or weighted MG aggregation.
7. **Interface compatibility:** preserve matrix residual/prediction defaults for existing users, add explicit vector/long formats with original-row mapping, and make standard generics reliable. Do not fabricate likelihoods or an OLS interpretation of heterogeneous-model degrees of freedom.

## Commit sequence and acceptance criteria

### A. Baseline and safe input handling

1. `test: establish estimator reference fixtures`
   - Set up an isolated development library, install needed testing/documentation dependencies, and run existing tests and R CMD check.
   - Add independent, clean-data MG/CCE/ARDL references, checking coefficients, covariance, fitted values, residuals, and long-run ratios.
   - Keep review artifacts out of the package build through a narrowly scoped build-ignore entry if needed. Keep large simulation output separate from fast tests.
   - Acceptance: baseline outcomes recorded; reference tests pass without changing estimators.

2. `fix: validate specifications and reject ignored options`
   - Preserve named CSA lags; reject fractional, negative, nonfinite, duplicate, and unknown lag specifications.
   - Validate spec classes/structure at the fitting boundary, with concise actionable errors.
   - Reject unused `...` and nondefault unimplemented settings; do not reject supported settings.
   - Acceptance: formerly ignored requests now fail clearly; all supported examples remain executable.

3. `fix: validate panel keys and preserve observation identity`
   - Reject duplicate/missing/nonfinite panel keys; validate id/time names and internal-name collisions.
   - Preserve original row identifiers before sorting. Safely normalize supported pdata.frame indexes.
   - Acceptance: duplicates cannot enter estimation; shuffled inputs preserve estimates and recover original-order outputs.

4. `fix: construct panel and CSA lags on the time grid`
   - Implement one tested key-based lag helper with explicit time-step validation.
   - Apply it to response, regressor, and CSA lags; use actual grid distance for trends.
   - Acceptance: absent periods and explicit missing rows produce the same lag availability; unit boundaries and all-panel gaps are respected.

5. `fix: honor missing values in cross-sectional averages`
   - Repair leave-one-out arithmetic for missing own observations and `na.rm=FALSE` propagation.
   - Validate weight lengths/values and zero denominators; protect generated columns against collisions.
   - Acceptance: independently computed weighted/unweighted and leave-one-out examples agree, including all-missing and zero-weight slices.

### B. Estimation sample, identification, and formulas

6. `refactor: share unit regression and sample bookkeeping`
   - Extract common model-frame, unit-fit, output-collection, and exclusion-recording helpers from duplicated engines.
   - Avoid repeated model-frame evaluation and rbind growth; retain unit fit metadata needed for inference and extraction.
   - Acceptance: clean-data outputs unchanged within numerical tolerances; no new statistical behavior in this refactor commit.

7. `fix: respect R model frames and evaluated formula terms`
   - Implement subset/NA handling; preserve terms, contrasts, levels, formula environments, and evaluated responses.
   - Expand dot formulas before bookkeeping and generate default CSAs from the approved evaluated-variable policy.
   - Use evaluated responses in fit statistics; define intercept-free R-squared behavior.
   - Acceptance: transformed response, factors/interactions, non-syntactic names, dot formulas, and NA restoration have meaningful tests. Explicitly unsupported dynamic formulas error clearly.

8. `fix: enforce identification after CSA projection`
   - Build the nuisance-space basis first; test rank of residualized economic regressors with a documented numerical tolerance.
   - Record rank, residual degrees of freedom, terms not identified, and exclusion reasons per unit.
   - Acceptance: a regressor fully spanned by CSA is not reported as identified; duplicate nuisance columns do not unnecessarily invalidate economic estimates; insufficient-T units are handled explicitly.

9. `fix: align mean-group estimates and covariance samples`
   - Apply the approved common-unit policy consistently to estimates, covariance, and reported counts.
   - Expose observed/eligible/estimated N and per-unit effective T; eliminate silent NaN aggregate results.
   - Acceptance: estimates and covariance match independent calculations on the exact same unit set, including rank-deficient examples.

### C. Inference and diagnostics

10. `fix: correct HC covariance calculations`
    - Square adjusted residuals in sandwich meat; compute leverage without allocating an n-by-n matrix.
    - Validate dimensions, rank, and leverage-one/zero-df cases instead of masking invalid calculations.
    - Acceptance: HC0-HC3 match a trusted implementation on actual OLS residuals; residual-sign invariance holds.

11. `fix: define and correct weighted mean-group covariance`
    - Write down the equal/unequal-weight estimand, normalization, finite-sample correction, and missingness policy before coding.
    - Implement only the agreed valid cases of `pooled_vcov`; reject unsupported missing/weight combinations and give `pairwise` explicit behavior.
    - Acceptance: equal weights match sample covariance divided by N, unequal-weight cases match independently derived references, and invalid weights fail.
    - Gate: ask the user if the weighted estimand or compatibility policy cannot be resolved from the approved design; do not substitute weighted coefficient dispersion.

12. `fix: align CS-ARDL parameter covariance and reporting`
    - Construct unit-level transformed parameter vectors and matching covariance for approved components.
    - Keep zero adjustment reportable when long-run ratios are undefined; expose exclusions and denominator diagnostics.
    - Rename the levels table to avoid claiming a full ECM transformation; record stability diagnostics without silently selecting a different population.
    - Acceptance: coefficient/covariance names match; confint agrees with reported SEs; exact parameter dependencies are documented and tested.

13. `fix: make dependence diagnostics explicit and reproducible`
    - Remove automatic weighted/PCA work; separate pairwise CD handling from PCA sample requirements.
    - Record diagnostic N/T, overlaps, excluded units, and missing results with reasons.
    - Respect the approved RNG policy and validate diagnostic arguments.
    - Acceptance: ordinary fitting does not change RNG state; one missing unit cannot silently erase a valid diagnostic sample.

14. `fix: correct CD-star unit scale calculations`
    - Check equations against the primary paper and versioned xtcd2 code, retaining columnwise scales with explicit sweep operations.
    - Validate feasible PCA counts and degenerate residuals; clearly state supported balanced/unbalanced behavior.
    - Acceptance: deterministic algebra references agree; heterogeneous-scale simulations and fixtures support the correction.

15. `fix: implement a theoretically validated CDw-plus statistic`
    - Resolve the exact paper-defined threshold, enhancement scaling, and reference distribution before implementation.
    - Treat the Stata implementation as a comparison target, not proof of correctness; document any justified divergence.
    - Acceptance: independent statistic checks and predeclared size/power simulations across relevant N/T designs pass. Use Monte Carlo uncertainty, not arbitrary exact rejection counts.
    - Gate: if the named method or derivation remains ambiguous, present the alternatives to the user. Do not ship an unverified replacement or silently relabel it.

16. `fix: validate remaining diagnostic and covariance edge cases`
    - Clarify prepare_cd_input demeaning semantics; honor explicit residual-type requests.
    - Validate cluster vectors, missing cluster IDs, one-cluster and insufficient-df cases; compare one/two-way covariance on supported full-rank designs.
    - Acceptance: clear errors for unsupported inputs and reference agreement for supported inputs. No unrelated new covariance estimator is added.

### D. R modelling interface and release readiness

17. `feat: provide reliable model extraction methods`
    - Add explicit nobs, model.frame, model.matrix, terms, formula, fitted, and residual accessors using retained metadata.
    - Add original-row vector/long formats without changing existing matrix defaults; define df access honestly and test update behavior.
    - Acceptance: generics return expected types and samples; no reliance on accidental partial matching; installed-namespace tests pass.

18. `feat: add tidy model summaries and augmentation`
    - Add tidy/glance/augment methods following a deliberate dependency strategy.
    - Test modelsummary/lmtest compatibility where applicable; expose component and sample choices rather than hiding them.
    - Acceptance: displayed inference agrees with coef/vcov/confint; augmentation aligns to original rows and does not claim unsupported predictions.

19. `docs: document estimator assumptions and compatibility changes`
    - Consolidate user-facing assumptions, sample policies, diagnostic limitations, supported options, and migration examples.
    - Document ordinary versus adjusted fit statistics and current long-run interpretation.
    - Keep earlier commits' documentation intact; this commit is an overview/migration pass, not delayed documentation for fixes.

20. `test: complete cross-platform and statistical validation`
    - Run the full installed-package suite, examples, documentation generation, and R CMD check.
    - Add/verify Windows and Linux CI and separate slower statistical validation scripts.
    - If Stata is available, add actual version-pinned reference outputs; otherwise clearly mark Stata execution as outstanding and ask for access/fixtures before claiming parity.
    - Acceptance: no unexplained regressions or package-check issues; statistical validation results recorded with seeds, versions, and tolerances.

## Completion criteria

All confirmed defects have a tested fix or an explicit, approved unsupported behavior. Coefficient and covariance samples agree. Unsupported options never silently succeed. R methods expose the actual fitted model. Statistical diagnostics have derivations and numerical evidence. Each commit is independently reviewable and package code has concise documentation. The final handoff includes commit summaries, validation results, compatibility changes, and any remaining evidence gaps.

## Review requested

Please approve or amend the scope and the seven proposed design decisions above before implementation. In particular, confirm the common-unit exclusion policy, explicit time-step interface, CS-ARDL component/sample behavior, and deferral of new estimators and forecasting. Any unresolved statistical/API choice remains a question, not an implementation assumption.

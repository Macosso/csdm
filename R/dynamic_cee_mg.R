#' Dynamic Common Correlated Effects Mean Group (CCE-MG) Estimator
#'
#' @description
#' Levels ARDL(1) with Common Correlated Effects (CCE) augmentation and
#' Mean Group (MG) averaging. Includes p lags of cross-sectional averages (CSAs)
#' of y and regressors. Designed for heterogeneous panels with latent common factors
#' (Chudik & Pesaran; Pesaran & Smith).
#'
#' @param formula y ~ x1 + x2 + ...
#' @param data data.frame, tibble, or `plm::pdata.frame`.
#' @param id,time Character names of unit and time identifiers (required if not a `pdata.frame`).
#' @param p Integer, number of CSA lags (default 3).
#' @param na.action NA handler (default `na.omit`).
#' @param vcov_unit One of c("OLS","HC"); per-unit covariance for reporting only.
#' @param cd_method One of c("CDstar","CD","none"); residual cross-section test to run in summary.
#' @details
#' \strong{Model and purpose.}
#' This function implements a heterogeneous dynamic panel regression with latent common factors
#' using the \emph{Common Correlated Effects} (CCE) approach and \emph{Mean Group} (MG) averaging.
#' For each cross-section \eqn{i=1,\dots,N} and time \eqn{t=1,\dots,T}, we estimate an ARDL(1) in levels:
#' \deqn{ y_{it} = \alpha_i + \phi_i y_{i,t-1} + x_{it}'\beta_i \;+\; z_t'\gamma_i \;+\; u_{it}, }
#' where \eqn{x_{it}} are unit-specific regressors and \eqn{z_t} collects \emph{cross-sectional averages} (CSAs)
#' of \eqn{y} and \eqn{x} and their lags (lag-augmented CCE).
#' The CSAs proxy the unobserved common factors and their dynamics; coefficients \eqn{(\phi_i,\beta_i,\gamma_i)}
#' are allowed to be heterogeneous across units. The reported MG coefficients are averages of the
#' unit-level estimates for the subset of “economic” parameters (intercept, \eqn{\phi_i}, and the
#' \eqn{\beta_i} on user-supplied regressors), excluding the CSA controls.
#'
#' \strong{CCE augmentation and lagging.}
#' The CCE idea is to include at each \eqn{t} the contemporaneous CSAs \eqn{ \bar y_t, \bar x_t } to span the
#' factor space induced by pervasive shocks. With weakly exogenous regressors and dynamic factors, it is
#' recommended to \emph{augment} the controls with \eqn{p \ge 1} lags of those CSAs to better track factor
#' persistence. Practically, the function:
#' \enumerate{
#'   \item Computes time-wise means \eqn{ \bar y_t = N_t^{-1}\sum_{i\in \mathcal{I}_t} y_{it}} and
#'         \eqn{ \bar x_t = N_t^{-1}\sum_{i\in \mathcal{I}_t} x_{it}}, where \eqn{ \mathcal{I}_t} is the set of
#'         units observed at \eqn{t} (works with unbalanced panels).
#'   \item Forms \eqn{ \{\bar y_t, \bar x_t\}} on a \emph{time-level} table, then adds lags
#'         \eqn{ \{\bar y_{t-\ell}, \bar x_{t-\ell}\}_{\ell=1}^p } by time only.
#'   \item Merges those CSA controls back to the panel and estimates each unit equation by OLS.
#' }
#' This ordering is essential: lagging must be done \emph{by time}, not within the id–time rows, to avoid
#' mixing information across units.
#'
#' \strong{Dynamic specification (ARDL(1)).}
#' The default specification includes one lag of \eqn{y_{it}} per unit (an ARDL(1) in levels). Higher-order
#' own lags could be added in extensions, but ARDL(1) is a common baseline in DC-CCE applications. Deterministic
#' components (unit intercepts) are included by default; unit-specific trends can be added in future versions.
#'
#' \strong{Mean Group (MG) estimation and inference.}
#' Let \eqn{\hat\theta_i} be the vector of unit-level coefficients of interest (intercept, lag of \eqn{y}, and
#' user regressors; CSA controls are \emph{nuisance} and are not averaged). The MG estimator is
#' \eqn{\bar\theta = N^{-1}\sum_{i=1}^N \hat\theta_i}. Its variance is computed from the cross-sectional dispersion:
#' \deqn{ \widehat{\mathrm{Var}}(\bar\theta) = \frac{1}{N(N-1)} \sum_{i=1}^N
#'        \left(\hat\theta_i - \bar\theta\right)\left(\hat\theta_i - \bar\theta\right)'. }
#' Reported standard errors are the square-roots of the diagonal entries of this matrix. This is the usual MG
#' large-\eqn{N} inference, robust to slope heterogeneity across units. Per-unit OLS (or HC) SEs can be provided
#' for diagnostics, but MG uncertainty is driven by cross-section dispersion of the unit estimates.
#'
#' \strong{Residual cross-sectional dependence.}
#' Even after CCE filtering, residuals may exhibit \emph{weak} cross-sectional dependence.
#' The summary method can run a dependence test; for best small-sample behavior, a bias-corrected CD test
#' (often called CD\*) is recommended. If only the classic CD test is available, interpret results cautiously
#' in factor-driven panels.
#'
#' \strong{Unbalanced panels and missing data.}
#' CSAs at time \eqn{t} are computed using only units observed at \eqn{t}; lags of CSAs are constructed on the
#' time index and then merged. Per-unit regressions are estimated on the rows that are complete for that unit
#' (with respect to \eqn{y}, regressors, and all selected CSA controls), so some units may be dropped if
#' \eqn{T_i} becomes too small after lagging. MG averaging automatically adapts to the set of units actually
#' estimated (effective \eqn{N}).
#'
#' \strong{Choosing the number of CSA lags \code{p}.}
#' The parameter \code{p} controls how aggressively we capture factor dynamics. Rules of thumb include small
#' fixed values (e.g., 1–3) for moderate \eqn{T}, or data-driven selection via information criteria in larger samples.
#' Excessively large \code{p} with small \eqn{T} can over-fit and reduce usable observations.
#'
#' \strong{What is averaged (and what is not).}
#' The \code{coef} and \code{se} returned by the function pertain to the intercept, the lag of \eqn{y}, and the
#' user-specified regressors. The CSA controls (contemporaneous and lagged averages) are treated as high-dimensional
#' nuisance covariates used to span the factor space; their coefficients are \emph{not} summarized in the MG output.
#'
#' \strong{Relation to alternative estimators.}
#' \itemize{
#'   \item \emph{MG vs PMG/DFE:} MG allows full slope heterogeneity and averages post-estimation.
#'         Pooled Mean Group (PMG) imposes long-run parameter homogeneity but allows heterogeneous adjustment
#'         dynamics; Dynamic Fixed Effects (DFE) imposes common slopes. MG is the least restrictive but can be less
#'         efficient if homogeneity holds.
#'   \item \emph{CCE vs standard FE:} Without CSA controls, correlated latent factors bias FE/OLS estimates when
#'         regressors load on the factors. CCE is designed to mitigate this by projecting on the space spanned by CSAs.
#' }
#'
#' \strong{Numerical notes and diagnostics.}
#' Collinearity can arise when \code{p} is large relative to \eqn{T} or when regressors are highly persistent and
#' strongly correlated with CSAs. Inspect per-unit condition numbers, and consider reducing \code{p} or the set
#' of regressors. After estimation, inspect: (i) dispersion of unit-level coefficients, (ii) residual cross-section
#' tests (CD/CD\*), and (iii) sensitivity of MG estimates to \code{p}.
#'
#' \strong{Limitations.}
#' The current version uses ARDL(1) and time-invariant \code{p} for all units. It does not (yet) include unit-specific
#' trends, higher-order dynamics, endogenous regressors with external instruments, or small-\eqn{T} finite-sample
#' corrections. These can be added in future releases.
#'
#' @references
#' \itemize{
#'   \item Chudik, A. and M. H. Pesaran (2015). \emph{Common correlated effects estimation with weakly exogenous regressors}.
#'   \item Pesaran, M. H. and R. Smith (1995). \emph{Estimating long-run relationships from dynamic heterogeneous panels}.
#'   \item Pesaran, M. H., Y. Shin, and R. P. Smith (1999). \emph{Pooled mean group estimation of dynamic heterogeneous panels}.
#'   \item Ditzen, J. (2021). \emph{Cross-sectional dependence in panel data models: concepts and estimation} (overview).}
#' (For residual dependence testing, see bias-corrected CD tests.)
#' @return An object of class 'dynamic_cce_mg' with:
#'   - coef: MG coefficients (named vector; intercept, lag_y, own regressors only)
#'   - se: MG standard errors from cross-section dispersion
#'   - vcov_mg: MG variance-covariance matrix
#'   - coef_i, vcov_i: per-unit estimates and vcov
#'   - residuals: pdata.frame(id,time,residual)
#'   - r.squared: per-unit R^2
#'   - meta: list with N,T,p, regressors, CSA spec
#' @importFrom dplyr group_by ungroup arrange mutate select across all_of sym summarise left_join everything n
#' @importFrom purrr map map_dbl map2
#' @importFrom tidyr nest unnest unnest_wider pivot_longer pivot_wider
#' @importFrom stats lm model.frame as.formula na.omit complete.cases coef vcov
#' @export
dynamic_cce_mg <- function(formula, data, id = NULL, time = NULL,
                           p = 3, na.action = na.omit,
                           vcov_unit = c("OLS","HC"),
                           cd_method = c("CDstar","CD","none")) {
  vcov_unit <- match.arg(vcov_unit)
  cd_method <- match.arg(cd_method)

  # ---- indices ----
  if (inherits(data, "pdata.frame")) {
    idx <- names(plm::index(data)); id <- idx[1]; time <- idx[2]
    df <- as.data.frame(data)
  } else {
    if (is.null(id) || is.null(time)) stop("Provide 'id' and 'time' when data is not a pdata.frame.")
    df <- as.data.frame(data)
  }

  # ---- variables ----
  mf <- model.frame(formula, data = df)
  vars <- all.vars(formula)
  y <- vars[1]; x_vars <- vars[-1]
  if (length(x_vars) == 0L) stop("At least one regressor is required on the RHS.")

  # keep only needed cols
  df <- df |>
    dplyr::select(!!sym(id), !!sym(time), dplyr::all_of(c(y, x_vars)))

  # ---- build CSAs at time-level (correct) ----
  # time-level table of CSAs
  csa_time <- df |>
    dplyr::group_by(!!sym(time)) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(c(y, x_vars)),
                                   ~ mean(.x, na.rm = TRUE),
                                   .names = "cs_{.col}"),
                     .groups = "drop") |>
    dplyr::arrange(!!sym(time))

  # add p lags of CSAs by time only
  add_lags <- function(tab, vnames, L) {
    for (l in seq_len(L)) {
      for (v in vnames) {
        tab[[paste0(v, "_lag", l)]] <- dplyr::lag(tab[[v]], l)
      }
    }
    tab
  }
  csa_cols <- names(csa_time)[startsWith(names(csa_time), "cs_")]
  if (p > 0L) csa_time <- add_lags(csa_time, csa_cols, p)

  # ---- merge CSAs back to panel ----
  df <- dplyr::left_join(df, csa_time, by = time)

  # ---- add ARDL(1) lag of y by unit ----
  df <- df |>
    dplyr::arrange(!!sym(id), !!sym(time)) |>
    dplyr::group_by(!!sym(id)) |>
    dplyr::mutate(!!paste0("lag_", y) := dplyr::lag(!!sym(y), 1L)) |>
    dplyr::ungroup()

  # regressors = lag_y + own x + contemporaneous CSAs + CSA lags
  avg_names <- csa_cols
  lag_labs  <- if (p > 0L) {
    as.vector(outer(avg_names, seq_len(p), FUN = function(v, l) paste0(v, "_lag", l)))
  } else character(0)
  regressors <- unique(c(paste0("lag_", y), x_vars, avg_names, lag_labs))

  # filter rows with complete cases for y and regressors (use provided na.action)
  df2 <- na.action(df[, c(id, time, y, regressors)])

  # ---- nest by unit and fit ----
  nested <- df2 |>
    tidyr::nest(.by = !!sym(id)) |>
    dplyr::rename(.unit_id = !!sym(id)) # preserve id column

  fit_formula <- stats::as.formula(paste(y, "~", paste(regressors, collapse = " + ")))

  get_vcov <- function(fit) {
    if (vcov_unit == "HC") {
      if (!requireNamespace("sandwich", quietly = TRUE))
        stop("Install 'sandwich' for HC vcov.")
      sandwich::vcovHC(fit, type = "HC1")
    } else {
      stats::vcov(fit)
    }
  }

  nested <- df2 |>
    tidyr::nest(.by = !!sym(id)) |>
    dplyr::rename(.unit_id = !!sym(id)) |>
    dplyr::mutate(
      results = purrr::map2(data, .unit_id, ~ {
        di  <- .x
        uid <- .y

        cols   <- c(y, regressors)
        mask   <- stats::complete.cases(di[, cols, drop = FALSE])
        di_fit <- di[mask, , drop = FALSE]

        if (nrow(di_fit) == 0L) {
          return(list(
            model = NULL,
            coef  = setNames(numeric(0), character(0)),
            vc    = NULL,
            r2    = NA_real_,
            residuals = setNames(
              data.frame(numeric(0), numeric(0), numeric(0)),
              c(id, time, "residual")
            )
          ))
        }

        fit <- stats::lm(fit_formula, data = di_fit)
        r   <- stats::residuals(fit)

        # build residual frame using the group id (uid), not di_fit[[id]]
        df_res <- data.frame(uid, di_fit[[time]], r, check.names = FALSE)
        names(df_res) <- c(id, time, "residual")

        list(
          model = fit,
          coef  = stats::coef(fit),
          vc    = get_vcov(fit),
          r2    = summary(fit)$r.squared,
          residuals = df_res
        )
      }),
      model     = purrr::map(results, "model"),
      coef      = purrr::map(results, "coef"),
      vc        = purrr::map(results, "vc"),
      r2        = purrr::map_dbl(results, "r2"),
      residuals = purrr::map(results, "residuals")
    ) |>
    dplyr::select(-results)



   # note: some units may have NA coef for some terms due to collinearity or insufficient data

  # ---- align coefficients across units ----
  all_terms <- unique(unlist(purrr::map(nested$coef, names)))
  coef_mat <- do.call(rbind, lapply(nested$coef, function(b) {
    out <- setNames(rep(NA_real_, length(all_terms)), all_terms)
    out[names(b)] <- b
    out
  }))
  rownames(coef_mat) <- nested$.unit_id

  # MG average on the "MG terms": intercept, lag_y, and user regressors only (exclude CSA controls)
  mg_terms <- intersect(c("(Intercept)", paste0("lag_", y), x_vars), colnames(coef_mat))
  B <- coef_mat[, mg_terms, drop = FALSE]
  b_mg <- colMeans(B, na.rm = TRUE)

  # Cross-section dispersion VCOV (Pesaran–Smith MG)
  # Use only units with all terms present for each column to compute dispersion robustly
  complete_rows <- stats::complete.cases(B)
  Bc <- B[complete_rows, , drop = FALSE]
  N_eff <- nrow(Bc)
  if (N_eff < 2L) stop("Not enough complete unit estimates to compute MG dispersion.")
  dev <- sweep(Bc, 2L, b_mg, "-")
  V_mg <- crossprod(as.matrix(dev)) / (N_eff * (N_eff - 1))
  se_mg <- sqrt(diag(V_mg))

  # residuals and R2
  res_df <- tidyr::unnest(nested, residuals)
  r2 <- setNames(nested$r2, nested$.unit_id)

  out <- list(
    coef      = b_mg,
    se        = setNames(se_mg, names(b_mg)),
    vcov_mg   = V_mg,
    coef_i    = setNames(split(coef_mat, row(coef_mat)), NULL), # raw matrix is often sufficient
    vcov_i    = setNames(nested$vc, nested$.unit_id),
    coef_df   = data.frame(unit = nested$.unit_id, B, check.names = FALSE),
    r.squared = r2,
    residuals = plm::pdata.frame(res_df, index = c(id, time)),
    call      = formula,
    p         = p,
    id        = id,
    time      = time,
    meta      = list(
      N = length(unique(df2[[id]])),
      T = length(unique(df2[[time]])),
      regressors = regressors,
      csa_vars = c(y, x_vars),
      csa_lags = p,
      cd_method = cd_method,
      vcov_unit = vcov_unit
    )
  )
  class(out) <- "dynamic_cce_mg"
  out
}

#' @export
print.dynamic_cce_mg <- function(x, ...) {
  cat("Dynamic CCE–MG (ARDL(1) + CCE controls)\n")
  cat("Call:", deparse(x$call), "\n")
  cat("CSA lags (p):", x$p, "\n\n")
  cat("Mean-group coefficients:\n")
  print(round(x$coef, 6))
  cat("\nStd. errors (cross-section dispersion):\n")
  print(round(x$se, 6))
  invisible(x)
}

#' @export
summary.dynamic_cce_mg <- function(object, ...) {
  mg <- data.frame(
    estimate = object$coef,
    std.error = object$se,
    z = object$coef / object$se,
    row.names = names(object$coef)
  )
  mg$p.value <- 2 * (1 - pnorm(abs(mg$z)))
  mg$significance <- cut(mg$p.value,
                         breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                         labels = c("***","**","*","."," "),
                         right = FALSE)

  CD <- NULL
  if (object$meta$cd_method != "none") {
    if (object$meta$cd_method == "CDstar" && exists("cd_star")) {
      CD <- cd_star(object$residuals, var = "residual")  # TODO: implement CD* per Pesaran & Xie (2022)
    } else if (object$meta$cd_method == "CD" && exists("cd")) {
      CD <- cd(object$residuals, var = "residual")
    }
  }

  out <- list(
    call = object$call,
    coefficients = mg,
    r.squared = mean(object$r.squared, na.rm = TRUE),
    CD = CD,
    meta = object$meta
  )
  class(out) <- "summary.dynamic_cce_mg"
  out
}

#' @export
print.summary.dynamic_cce_mg <- function(x, digits = 6, ...) {
  cat("Dynamic CCE–MG (heterogeneous ARDL(1) with CCE)\n")
  cat("Call:", deparse(x$call), "\n\n")

#  print(within(round(x$coefficients, digits), { }), right = FALSE)

  x$coefficients <- x$coefficients |>
    dplyr::mutate(
      estimate = formatC(estimate, format = "f", digits = digits),
      std.error = formatC(std.error, format = "f", digits = digits),
      z = formatC(z, format = "f", digits = digits),
      p.value = formatC(p.value, format = "f", digits = 4)
    )
  print(x$coefficients)
  print("")
  cat("\nMean R^2 across units:", format(round(x$r.squared, 3), nsmall = 3), "\n")
  if (!is.null(x$CD)) {
    cat("Cross-section dependence test:", x$CD$method, "\n",
        "stat =", round(x$CD$statistic, 2), " p-value =", round(x$CD$p.value, 4), "\n")
  }
  invisible(x)
}

#' @export
vcov.dynamic_cce_mg <- function(object, ...) {
  object$vcov_mg
}

#' @export
coef.dynamic_cce_mg <- function(object, ...) {
  object$coef
}

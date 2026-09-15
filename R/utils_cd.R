# utils_cd.R - Cross-sectional dependence tests for panel models

#' Cross-sectional dependence (CD) tests for panel residuals
#'
#' Computes Pesaran CD, CDw, CDw+, and CD* tests for cross-sectional dependence
#' in panel residuals. The implementation supports residual matrices or fitted
#' \code{csdm_fit} objects and provides consistent handling of unbalanced panels.
#'
#' @param object A \code{csdm_fit} model object or a numeric matrix of residuals (N x T).
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class \code{cd_test} with fields \code{tests}, \code{type},
#'   \code{N}, \code{T}, \code{na.action}, \code{excluded_units},
#'   \code{excluded_times}, \code{kept_times}, and \code{call}. The \code{tests}
#'   list contains one or more test results, each with \code{statistic} and
#'   \code{p.value}.
#'
#' @details
#' ## Notation
#'
#' Let \eqn{E} be the residual matrix with \eqn{N} cross-sectional units and \eqn{T}
#' time periods. For each unit pair \eqn{(i,j)}, let \eqn{T_{ij}} be the number of
#' overlapping time periods and \eqn{\rho_{ij}} the pairwise correlation.
#'
#' ## Test statistics
#'
#' \describe{
#'   \item{CD (Pesaran, 2015)}{
#'     \deqn{CD = \sqrt{\frac{2}{N(N-1)}} \sum_{i<j} \sqrt{T_{ij}} \, \rho_{ij}}
#'   }
#'   \item{CDw (Juodis and Reese, 2022)}{
#'     Independent Rademacher weights \eqn{w_i \in \{-1,1\}} are applied by
#'     unit. For a balanced residual panel, the statistic is
#'     \deqn{CD_W = \left(\frac{1}{NT}\sum_{i,t}w_i^2 e_{it}^2\right)^{-1}
#'     \sqrt{\frac{2}{TN(N-1)}}
#'     \sum_t\sum_{i<j}w_i e_{it}w_j e_{jt}.}
#'     The first factor is the inverse pooled residual variance. One set of
#'     random weights is drawn per call; use \code{seed} for reproducibility.
#'   }
#'   \item{CDw+ (Juodis and Reese, 2022; Fan, Liao, and Yao, 2015)}{
#'     The power-enhanced statistic is
#'     \deqn{CD_{W+} = CD_W + \sum_{i<j}|\rho_{ij}|
#'     1\left\{|\rho_{ij}| > 2\sqrt{\log(N)/T}\right\}.}
#'     Here \eqn{\rho_{ij}} is the ordinary residual correlation, without
#'     multiplication by \eqn{\sqrt{T}}. The nonnegative screening term is
#'     asymptotically zero under the conditions of the null hypothesis.
#'   }
#'   \item{CD* (Pesaran and Xie, 2021)}{
#'     CD is computed on residuals after removing \code{n_pc} principal components
#'     from \eqn{E}. This provides a bias-corrected test under multifactor errors.
#'   }
#' }
#'
#' CD* requires a nondegenerate bias-correction denominator. Near-zero
#' denominators can produce severe size distortions, including proportional
#' loading/error-scale designs after standardization. Numerical rank checks
#' do not establish the validity of the asymptotic approximation.
#'
#' ## Missing data and balance
#'
#' Time periods containing no finite residual for any retained unit are outside
#' the effective residual sample and are always removed before balance is
#' assessed. Partially observed periods are handled according to \code{na.action}.
#'
#' \describe{
#'   \item{CD}{Uses pairwise-complete observations by default. Each pairwise
#'   correlation uses available overlaps.}
#'   \item{CDw, CDw+}{Require a balanced sample; explicitly select complete times if desired.}
#'   \item{CD*}{Requires a balanced panel. Explicitly setting \code{na.action = "drop.incomplete.times"}
#'   removes any time period with missing observations. With \code{na.action = "pairwise"},
#'   CD* returns \code{NA} and a warning when missing values are present.}
#' }
#'
#' @references
#' \insertRef{Pesaran2015}{csdm}
#'
#' \insertRef{Pesaran2021}{csdm}
#'
#' \insertRef{JuodisReese2021}{csdm}
#'
#' \insertRef{FanLiaoYao2015}{csdm}
#'
#' \insertRef{PesaranXie2021}{csdm}
#'
#' @examples
#' # Simulate independent and dependent panels
#' set.seed(1)
#' E_indep <- matrix(rnorm(100), nrow = 10)
#' E_dep <- matrix(rnorm(10), nrow = 10, ncol = 10, byrow = TRUE)
#'
#' # Compute all tests
#' cd_test(E_indep, type = "all")
#' cd_test(E_dep, type = "CD")
#'
#' # Specific test with parameters
#' cd_test(E_indep, type = "CDstar", n_pc = 2)
#'
#' # From a fitted csdm model
#' data(PWT_60_07, package = "csdm")
#' df <- PWT_60_07
#' ids <- unique(df$id)[1:10]
#' df_small <- df[df$id %in% ids & df$year >= 1970, ]
#' fit <- csdm(
#'   log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   data = df_small,
#'   id = "id",
#'   time = "year",
#'   model = "cce",
#'   csa = csdm_csa(vars = c("log_rgdpo", "log_hc", "log_ck", "log_ngd"))
#' )
#' cd_test(fit, type = "all")
#'
#' @export
cd_test <- function(object, ...) {
  UseMethod("cd_test")
}

#' @rdname cd_test
#' @param type Which test(s) to compute: one of \code{"CD"}, \code{"CDw"}, \code{"CDw+"},
#'   \code{"CDstar"}, or \code{"all"} (default: \code{"CD"}).
#' @param n_pc Number of principal components for CD* (default 4).
#' @param seed Integer seed for weight draws. Seeded calls restore the caller's RNG state; NULL uses the current RNG stream.
#' @param min_overlap Minimum number of overlapping time periods required for a unit
#'   pair to be included in CD/CDw/CDw+ (default 2).
#' @param na.action How to handle missing data: \code{"drop.incomplete.times"}
#'   removes time periods with any missing observations to create a balanced panel for CD*;
#'   \code{"pairwise"} (default) uses pairwise correlations for CD; unbalanced CDw/CDw+ requests error and CD* warns.
#' @export
#' @method cd_test default
cd_test.default <- function(object,
                            type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                            n_pc = 4L,
                            seed = NULL,
                            min_overlap = 2L,
                            na.action = c("pairwise", "drop.incomplete.times"),
                            ...) {
  type <- match.arg(type)
  na.action <- match.arg(na.action)

  # Validate input
  if (!is.matrix(object) || !is.numeric(object)) {
    stop("cd_test.default: object must be a numeric matrix (N x T).")
  }
  if (nrow(object) < 2 || ncol(object) < 2) {
    stop("cd_test: At least 2 units and 2 time periods required.")
  }

  min_overlap <- .csdm_integer(min_overlap, "min_overlap")
  if (min_overlap < 2L) stop("'min_overlap' must be at least two.")
  if (any(is.infinite(object))) stop("Residuals may contain NA, but not infinite values.")
  if (!is.null(seed)) {
    seed <- .csdm_integer(seed, "seed")
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
  }
  usable <- apply(object, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) >= min_overlap && stats::sd(x) > 0
  })
  if (sum(usable) < 2L) stop("At least two nonconstant units with sufficient observations are required.")
  excluded_units <- which(!usable)
  E <- object[usable, , drop = FALSE]

  # Periods with no estimated residuals are not part of the residual sample.
  empty_times <- colSums(is.finite(E)) == 0L
  excluded_times <- if (is.null(colnames(E))) which(empty_times) else colnames(E)[empty_times]
  if (any(empty_times)) E <- E[, !empty_times, drop = FALSE]
  if (ncol(E) < 2L) stop("At least two time periods with residual observations are required.")

  # Handle missing data
  if (na.action == "drop.incomplete.times" && anyNA(E)) {
    # Drop time periods with any missing observations
    complete_times <- colSums(is.na(E)) == 0
    n_dropped <- sum(!complete_times)
    if (sum(complete_times) < 2) {
      stop("cd_test: After dropping incomplete time periods, fewer than 2 periods remain.")
    }
    E <- E[, complete_times, drop = FALSE]
    if (n_dropped > 0) {
      message(sprintf("cd_test: Dropped %d incomplete time period%s (%.1f%%). Balanced panel: %d units x %d periods.",
                      n_dropped, if(n_dropped > 1) "s" else "",
                      100 * n_dropped / ncol(object), nrow(E), ncol(E)))
    }
  }

  N <- nrow(E)
  Tt <- ncol(E)

  # Convert to T x N for correlation computation (Stata convention)
  data_tn <- t(E)

  out <- list()

  # 1. CD (classic Pesaran)
  if (type %in% c("CD", "all")) {
    cd_res <- .cd_compute_classic(data_tn, N, Tt, min_overlap)
    out$CD <- list(
      statistic = cd_res$statistic,
      p.value = cd_res$p.value,
      N = N,
      T = Tt,
      pairs_used = cd_res$pairs_used
    )
  }

  if (type %in% c("CDw", "CDw+", "all")) {
    if (anyNA(E)) stop("CDw/CDw+ currently require a balanced sample; use explicit drop.incomplete.times or classical CD.")
    if (!is.null(seed)) set.seed(seed)
    w <- sample(c(-1, 1), N, replace = TRUE)
    centered <- sweep(data_tn, 2L, colMeans(data_tn))
    variance <- mean(centered^2)
    weighted <- sweep(centered, 2L, w, "*")
    pair_sum <- sum(rowSums(weighted)^2 - rowSums(weighted^2)) / 2
    statistic <- sqrt(2 / (Tt * N * (N - 1))) * pair_sum / variance
    out$CDw <- list(statistic = statistic,
      p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
      N = N, T = Tt, pairs_used = choose(N, 2))
    if (type %in% c("CDw+", "all")) {
      correlation <- stats::cor(centered)
      rho <- abs(correlation[upper.tri(correlation)])
      threshold <- 2 * sqrt(log(N) / Tt)
      enhancement <- sum(rho[rho > threshold])
      statistic <- statistic + enhancement
      out$CDw_plus <- list(statistic = statistic,
        p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE),
        N = N, T = Tt, pairs_used = choose(N, 2),
        threshold = threshold, enhancement = enhancement)
    }
  }

  # 4. CD* (bias-corrected with PCA factor removal)
  if (type %in% c("CDstar", "all")) {
    # CD* requires balanced panel
    if (anyNA(E)) {
      warning("cd_test: CD* requires a balanced panel; returning NA due to missing values.")
      out$CDstar <- list(
        statistic = NA_real_,
        p.value = NA_real_,
        N = N,
        T = Tt,
        n_pc = as.integer(n_pc)
      )
    } else {
      cdstar_res <- .cd_compute_star(data_tn, N, Tt, n_pc)
      out$CDstar <- list(
        statistic = cdstar_res$statistic,
        p.value = cdstar_res$p.value,
        N = N,
        T = Tt,
        n_pc = as.integer(n_pc)
      )
    }
  }

  res <- list(
    tests = out,
    type = if (type == "all") names(out) else type,
    N = N,
    T = Tt,
    na.action = na.action,
    excluded_units = excluded_units,
    excluded_times = excluded_times,
    kept_times = colnames(E),
    call = match.call()
  )
  class(res) <- "cd_test"
  res
}

#' @rdname cd_test
#' @export
#' @method cd_test csdm_fit
cd_test.csdm_fit <- function(object,
                              type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                              n_pc = 4L,
                              seed = NULL,
                              min_overlap = 2L,
                              na.action = c("pairwise", "drop.incomplete.times"),
                              ...) {
  type <- match.arg(type)
  na.action <- match.arg(na.action)
  E <- get_residuals(object, type = "auto", strict = TRUE)
  cd_test.default(E, type = type, n_pc = n_pc, seed = seed,
                  min_overlap = min_overlap, na.action = na.action, ...)
}

#' @rdname cd_test
#' @param x An object of class \code{cd_test}.
#' @param digits Number of digits to print (default 3).
#' @export
#' @method print cd_test
print.cd_test <- function(x, digits = 3, ...) {
  if (is.null(x$tests) || length(x$tests) == 0) {
    cat("cd_test: no results\n")
    return(invisible(x))
  }

  cat("Cross-sectional dependence tests\n")
  if (!is.null(x$N) && !is.null(x$T)) {
    cat(sprintf("N = %d, T = %d\n", x$N, x$T))
  }
  cat("\n")

  tests <- x$tests
  type_all <- is.character(x$type) && length(x$type) > 1
  if (!type_all && is.character(x$type) && length(x$type) == 1) {
    type_all <- identical(x$type, "all")
  }
  type_map <- c(
    "CD" = "CD",
    "CDw" = "CDw",
    "CDw+" = "CDw_plus",
    "CDstar" = "CDstar"
  )

  get_num <- function(test, key) {
    if (!is.null(test[[key]])) as.numeric(test[[key]]) else NA_real_
  }

  if (type_all) {
    test_names <- names(tests)
    display_names <- gsub("CDw_plus", "CDw+", test_names, fixed = TRUE)
    stat <- vapply(tests, get_num, key = "statistic", FUN.VALUE = numeric(1))
    pval <- vapply(tests, get_num, key = "p.value", FUN.VALUE = numeric(1))
    out <- data.frame(
      statistic = stat,
      p.value = pval,
      row.names = display_names,
      stringsAsFactors = FALSE
    )
  } else {
    key <- type_map[[as.character(x$type[1])]]
    if (is.null(key) || is.null(tests[[key]])) {
      key <- names(tests)[1]
    }
    cd_test <- tests[[key]]
    test_name <- as.character(x$type[1])
    out <- data.frame(
      statistic = get_num(cd_test, "statistic"),
      p.value = get_num(cd_test, "p.value"),
      row.names = test_name,
      stringsAsFactors = FALSE
    )
  }

  fmt_num <- function(x) {
    ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
  }
  out$statistic <- fmt_num(out$statistic)
  out$p.value <- fmt_num(out$p.value)
  print(out, row.names = TRUE, right = TRUE)
  invisible(x)
}

# Internal helper: compute classic CD statistic with pairwise overlap
.cd_compute_classic <- function(data_tn, N, Tt, min_overlap = 2L) {
  # data_tn: T x N matrix
  # Returns: list(statistic, p.value, pairs_used)

  corr_mat <- stats::cor(data_tn, use = "pairwise.complete.obs")

  # Count valid pairs (those with sufficient overlap)
  pairs_used <- 0
  cd_sum <- 0

  for (i in seq_len(N - 1)) {
    for (j in (i + 1):N) {
      ok <- is.finite(data_tn[, i]) & is.finite(data_tn[, j])
      T_ij <- sum(ok)
      if (T_ij >= min_overlap) {
        rho_ij <- corr_mat[i, j]
        if (is.finite(rho_ij)) {
          cd_sum <- cd_sum + sqrt(T_ij) * rho_ij
          pairs_used <- pairs_used + 1
        }
      }
    }
  }

  if (pairs_used == 0) {
    warning("cd_test: No unit pairs have sufficient overlap for CD computation.")
    return(list(statistic = NA_real_, p.value = NA_real_, pairs_used = 0L))
  }

  # The Pesaran normalization uses the retained number of units.
  cd_stat <- sqrt(2 / (N * (N - 1))) * cd_sum
  cd_p <- 2 * stats::pnorm(abs(cd_stat), lower.tail = FALSE)

  list(statistic = cd_stat, p.value = cd_p, pairs_used = as.integer(pairs_used))
}

# Internal helper: compute CD* with PCA factor removal
.cd_compute_star <- function(data_tn, N, Tt, n_pc) {
  n_pc <- .csdm_integer(n_pc, "n_pc")
  if (n_pc >= min(N, Tt - 1L)) stop("'n_pc' must be below min(N, T - 1).")
  scales <- apply(data_tn, 2L, stats::sd)
  if (any(!is.finite(scales) | scales <= 0)) stop("CD* requires nonconstant complete unit series.")
  data_std <- sweep(sweep(data_tn, 2L, colMeans(data_tn)), 2L, scales, "/")
  if (n_pc == 0L) return(.cd_compute_classic(data_std, N, Tt))
  decomposition <- svd(data_std, nu = n_pc, nv = 0L)
  if (decomposition$d[n_pc] <= decomposition$d[1L] * 1e-7) stop("Requested factors exceed numerical rank.")
  factors <- cbind(1, decomposition$u[, seq_len(n_pc), drop = FALSE])
  beta <- qr.coef(qr(factors), data_std)
  residual <- data_std - factors %*% beta
  sigma <- sqrt(colMeans(residual^2))
  if (any(!is.finite(sigma) | sigma <= sqrt(.Machine$double.eps))) stop("CD* residual scales are degenerate.")
  loadings <- beta[-1L, , drop = FALSE]
  gamma <- sweep(loadings, 1L, sqrt(rowMeans(loadings^2)), "/")
  phi <- rowMeans(sweep(gamma, 2L, sigma, "/"))
  a <- as.numeric(1 - t(sweep(gamma, 2L, sigma, "*")) %*% phi)
  correction <- mean(a^2)
  if (!is.finite(correction) || correction <= sqrt(.Machine$double.eps)) stop("CD* bias correction is degenerate.")
  cd <- .cd_compute_classic(residual, N, Tt)$statistic
  statistic <- (cd + sqrt(Tt / 2) * (1 - correction)) / correction
  list(statistic = statistic, p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE))
}

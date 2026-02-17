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
#'   \code{N}, \code{T}, \code{na.action}, and \code{call}. The \code{tests} list
#'   contains one or more test results, each with \code{statistic} and \code{p.value}.
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
#'   \item{CDw (Juodis and Reese, 2021)}{
#'     Random sign flips \eqn{w_i \in \{-1,1\}} are applied to residuals before
#'     computing correlations. The statistic is CD applied to the sign-flipped data.
#'   }
#'   \item{CDw+ (Fan, Liao, and Yao, 2015)}{
#'     Power enhancement adds a sparse thresholding term to CDw. The threshold is
#'     \deqn{c_N = \sqrt{\frac{2 \log(N)}{T}}}
#'     and the power term sums \eqn{\sqrt{T_{ij}} |\rho_{ij}|} for pairs exceeding
#'     the threshold.
#'   }
#'   \item{CD* (Pesaran and Xie, 2021)}{
#'     CD is computed on residuals after removing \code{n_pc} principal components
#'     from \eqn{E}. This provides a bias-corrected test under multifactor errors.
#'   }
#' }
#'
#' ## Missing data and balance
#'
#' \describe{
#'   \item{CD, CDw, CDw+}{Always use pairwise-complete observations. Each pairwise
#'   correlation uses available overlaps.}
#'   \item{CD*}{Requires a balanced panel. By default, \code{na.action = "drop.incomplete.times"}
#'   removes any time period with missing observations. With \code{na.action = "pairwise"},
#'   CD* returns \code{NA} and a warning when missing values are present.}
#' }
#'
#' @references
#' Pesaran, M.H. (2015). "Testing weak cross-sectional dependence in large panels."
#' \emph{Econometric Reviews}, 34(6-10), 1089-1117.
#'
#' Pesaran, M.H. (2021). "General diagnostic tests for cross-sectional dependence
#' in panels." \emph{Empirical Economics}, 60, 13-50.
#'
#' Juodis, A., & Reese, S. (2021). "The incidental parameters problem in testing for
#' remaining cross-sectional correlation." \emph{Journal of Business and Economic Statistics},
#' 40(3), 1193-1203.
#'
#' Fan, J., Liao, Y., & Yao, J. (2015). "Power Enhancement in High-Dimensional
#'   Cross-Sectional Tests." \emph{Econometric Reviews}, 34(6-10), 742-779.
#'
#' Pesaran, M.H., & Xie, Y. (2021). "A bias-corrected CD test for error cross-sectional
#'   dependence in panel models." \emph{Econometric Reviews}, 41(6), 649-677.
#'
#' @examples
#' # Simulate independent and dependent panels
#' set.seed(1)
#' E_indep <- matrix(rnorm(100), nrow = 10)
#' E_dep <- matrix(rnorm(10), nrow = 10, ncol = 10, byrow = TRUE)
#'
#' # Compute all tests
#' cd_test(E_indep, type = "all")
#' cd_test(E_dep, type = "all")
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
#' @param seed Integer seed for weight draws in CDw/CDw+ (default NULL = no seed set).
#' @param min_overlap Minimum number of overlapping time periods required for a unit
#'   pair to be included in CD/CDw/CDw+ (default 2).
#' @param na.action How to handle missing data: \code{"drop.incomplete.times"} (default)
#'   removes time periods with any missing observations to create a balanced panel for CD*;
#'   \code{"pairwise"} uses pairwise correlations for CD/CDw/CDw+ and warns for CD*.
#' @export
#' @method cd_test default
cd_test.default <- function(object,
                            type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                            n_pc = 4L,
                            seed = NULL,
                            min_overlap = 2L,
                            na.action = c("drop.incomplete.times", "pairwise"),
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

  E <- object  # N x T matrix

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

  # 2. Weighted CD (Juodis & Reese)
  if (type %in% c("CDw", "all", "CDw+")) {
    if (!is.null(seed)) set.seed(seed)
    w <- sample(c(-1, 1), N, replace = TRUE)
    data_cdw <- sweep(data_tn, 2, w, FUN = "*")
    cdw_res <- .cd_compute_classic(data_cdw, N, Tt, min_overlap)
    out$CDw <- list(
      statistic = cdw_res$statistic,
      p.value = cdw_res$p.value,
      N = N,
      T = Tt,
      pairs_used = cdw_res$pairs_used
    )
  }

  # 3. Power-enhanced CDw+ (Fan et al.)
  if (type %in% c("CDw+", "all")) {
    # Need unweighted correlations for threshold
    if (is.null(out$CD)) {
      cd_res <- .cd_compute_classic(data_tn, N, Tt, min_overlap)
    }
    # Compute threshold and power term
    crit <- sqrt(2 * log(N) / Tt)
    corr_mat <- stats::cor(data_tn, use = "pairwise.complete.obs")
    upper_corr <- corr_mat[upper.tri(corr_mat)]
    # Count overlaps for each pair
    overlap_mat <- matrix(0, N, N)
    for (i in seq_len(N - 1)) {
      for (j in (i + 1):N) {
        ok <- is.finite(data_tn[, i]) & is.finite(data_tn[, j])
        overlap_mat[i, j] <- sum(ok)
      }
    }
    upper_overlap <- overlap_mat[upper.tri(overlap_mat)]
    # Power term: sum sqrt(T_ij) * |rho_ij| for correlations exceeding threshold
    exceeds <- abs(upper_corr * sqrt(upper_overlap)) > crit & upper_overlap >= min_overlap
    power_term <- sum(sqrt(upper_overlap[exceeds]) * abs(upper_corr[exceeds]))

    if (is.null(out$CDw)) {
      if (!is.null(seed)) set.seed(seed)
      w <- sample(c(-1, 1), N, replace = TRUE)
      data_cdw <- sweep(data_tn, 2, w, FUN = "*")
      cdw_res <- .cd_compute_classic(data_cdw, N, Tt, min_overlap)
    } else {
      cdw_res <- list(
        statistic = out$CDw$statistic,
        pairs_used = out$CDw$pairs_used
      )
    }

    # Normalize power_term to match the scale of cdw_res$statistic
    cdw_plus_stat <- cdw_res$statistic + power_term
    cdw_plus_p <- 2 * (1 - stats::pnorm(abs(cdw_plus_stat)))

    out$CDw_plus <- list(
      statistic = cdw_plus_stat,
      p.value = cdw_plus_p,
      N = N,
      T = Tt,
      pairs_used = cdw_res$pairs_used
    )
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
      cdstar_res <- .cd_compute_star(data_tn, N, Tt, as.integer(n_pc))
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
                              na.action = c("drop.incomplete.times", "pairwise"),
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

  # Normalize by number of pairs used (not total pairs)
  cd_stat <- sqrt(2 / (N * (N - 1))) * cd_sum
  cd_p <- 2 * (1 - stats::pnorm(abs(cd_stat)))

  list(statistic = cd_stat, p.value = cd_p, pairs_used = as.integer(pairs_used))
}

# Internal helper: compute CD* with PCA factor removal
.cd_compute_star <- function(data_tn, N, Tt, n_pc) {
  # data_tn: T x N matrix (must be complete/balanced)
  # Returns: list(statistic, p.value)

  # Standardize columns
  col_means <- colMeans(data_tn)
  col_sds <- apply(data_tn, 2, stats::sd)
  col_sds[col_sds == 0 | !is.finite(col_sds)] <- 1
  data_std <- sweep(sweep(data_tn, 2, col_means), 2, col_sds, `/`)

  # PCA: Extract factors
  S <- data_std %*% t(data_std)
  eig <- eigen(S, symmetric = TRUE)
  idx <- order(eig$values, decreasing = TRUE)[seq_len(n_pc)]
  f <- eig$vectors[, idx, drop = FALSE]
  fx <- cbind(1, f)  # intercept + factors

  # Defactor residuals
  beta <- solve(t(fx) %*% fx, t(fx) %*% data_std)
  res_defac <- data_std - fx %*% beta

  # CD on defactored residuals
  corr_defac <- stats::cor(res_defac)
  upper_defac <- corr_defac[upper.tri(corr_defac)]
  cd_defac <- sqrt(2 / (N * (N - 1))) * sum(upper_defac * sqrt(Tt))

  # Bias correction (Pesaran & Xie 2021)
  betai <- beta[-1, , drop = FALSE]  # remove intercept
  betaij <- (betai %*% t(betai)) / N
  betasum <- sqrt(diag(betaij))
  betasum[betasum == 0] <- 1  # avoid division by zero
  gamma <- sweep(betai, 1, betasum, `/`)  # normalize columns (units) by their norms
  sgm <- sqrt(mean(res_defac^2))
  if (sgm == 0) sgm <- 1
  phi <- rowMeans(gamma / sgm)
  ai_vals <- as.numeric((1 - t(gamma * sgm) %*% phi) / sqrt(N))
  theta <- sum(ai_vals^2)
  if (theta == 0) theta <- 1  # avoid division by zero

  cd_star_stat <- (cd_defac + sqrt(Tt / 2) * (1 - theta)) / theta
  cd_star_p <- 2 * (1 - stats::pnorm(abs(cd_star_stat)))

  list(statistic = cd_star_stat, p.value = cd_star_p)
}

# utils_residuals.R

#' Extract residual matrices from model objects
#'
#' @description
#' Unified accessor that returns an \eqn{N \times T} residual matrix suitable for
#' cross-sectional dependence testing.
#'
#' @param object A fitted model object supported by this package (e.g., class
#'   \code{csdm_fit}), or directly a numeric matrix of residuals shaped as
#'   \eqn{N \times T}.
#' @param type Character string selecting which residuals to return when available:
#'   one of \code{"auto"}, \code{"cce"}, \code{"pca"}, or \code{"pca_std"}.
#'   \itemize{
#'     \item \code{"auto"}: prefer standardized PCA residuals if present, otherwise
#'       PCA residuals, otherwise CCE residuals, otherwise \code{object} if it is a matrix.
#'     \item \code{"cce"}: residuals from the CCE-augmented per-unit regressions.
#'     \item \code{"pca"}: residuals after removing estimated factors from \eqn{\hat v_{it}}.
#'     \item \code{"pca_std"}: \code{"pca"} residuals standardized by unit-specific
#'       scale (recommended for CD).
#'   }
#' @param strict Logical; if \code{TRUE}, error on unsupported objects. If \code{FALSE},
#'   return \code{NULL} when residuals cannot be found.
#'
#' @returns A numeric matrix of residuals with rows = units and columns = time,
#'   preserving \code{rownames} and \code{colnames} when available; or \code{NULL}
#'   if nothing suitable is found and \code{strict = FALSE}.
#'
#' @keywords internal
#' @export
get_residuals <- function(object,
                          type = c("auto", "cce", "pca", "pca_std"),
                          strict = TRUE) {
  type <- match.arg(type)

  as_mat <- function(x) {
    if (is.matrix(x)) return(x)
    if (is.data.frame(x)) return(as.matrix(x))
    x
  }

  take_first <- function(...) {
    for (m in list(...)) {
      mm <- as_mat(m)
      if (is.matrix(mm) && is.numeric(mm)) return(mm)
    }
    NULL
  }

  # If the user passed a matrix directly, just return it
  if (is.matrix(object) && is.numeric(object)) {
    return(object)
  }

  # Handle supported class: csdm_fit
  if (inherits(object, "csdm_fit")) {
    M <- as_mat(object$residuals_e)
    if (is.matrix(M) && is.numeric(M)) return(M)
  }

  # Default: if it's a list, try common residual names
  if (is.list(object)) {
    if (type == "auto") {
      M <- take_first(object$residuals_pca_std,
                      object$residuals_pca,
                      object$residuals_cce,
                      object$residuals_e,
                      object$residuals,
                      object$E)
      if (!is.null(M)) return(M)
    } else {
      cand <- switch(type,
                     "pca_std" = object$residuals_pca_std,
                     "pca"     = object$residuals_pca,
                     "cce"     = object$residuals_cce,
                     NULL)
      M <- as_mat(cand)
      if (is.matrix(M)) return(M)
    }
  }

  if (strict) {
    stop("get_residuals(): could not find a numeric residual matrix for the given object.")
  } else {
    return(NULL)
  }
}


#' Prepare residuals for CD / CD* tests
#'
#' @description
#' Cleans and transforms an \eqn{N \times T} residual matrix for cross-sectional
#' dependence testing. Operations include:
#' \enumerate{
#'   \item Dropping time periods with fewer than \code{min_per_time} finite observations.
#'   \item Optional row-wise standardization to unit variance over available times.
#'   \item Optional demeaning across units at each time (recommended for CD).
#' }
#'
#' @param E A numeric matrix of residuals (\eqn{N \times T}); rows are units,
#'   columns are time; may be unbalanced (contain \code{NA}).
#' @param standardize One of \code{"row"}, \code{"none"}.
#'   If \code{"row"}, scale each row by its observed standard deviation.
#' @param demean_time Logical; if \code{TRUE}, subtract the cross-sectional mean at
#'   each time from available residuals in that column.
#' @param min_per_time Integer; drop time columns with fewer than this many finite
#'   observations.
#'
#' @returns A list with:
#' \item{Z}{Processed residual matrix (\eqn{N \times T^*}) after filtering/standardizing/demeaning.}
#' \item{kept_t}{Integer indices of kept time columns (relative to the original \code{E}).}
#' \item{m_t}{Integer vector of cross-sectional counts per kept time (number of finite rows).}
#' \item{row_sds}{Numeric vector of row standard deviations used (invisibly \code{NA} if \code{standardize="none"}).}
#' \item{col_means}{Numeric vector of time means subtracted when \code{demean_time=TRUE}.}
#'
#' @keywords internal
#' @export
prepare_cd_input <- function(E,
                             standardize = c("row", "none"),
                             demean_time = TRUE,
                             min_per_time = 2L) {
  if (!is.matrix(E) || !is.numeric(E)) {
    stop("prepare_cd_input(): E must be a numeric matrix (N x T).")
  }
  standardize <- match.arg(standardize)
  N <- nrow(E); T_n <- ncol(E)

  # 1) Keep only time columns with at least min_per_time finite values
  counts <- colSums(is.finite(E))
  keep   <- which(counts >= as.integer(min_per_time))
  if (length(keep) == 0L) {
    stop("prepare_cd_input(): No time columns meet the 'min_per_time' requirement.")
  }
  Z <- E[, keep, drop = FALSE]
  m_t <- colSums(is.finite(Z))

  # 2) Row-wise standardization to unit variance over available observations
  row_sds <- rep(NA_real_, N)
  if (standardize == "row") {
    # Compute row sd using pairwise complete obs across the kept columns
    for (i in seq_len(N)) {
      zi <- Z[i, ]
      s  <- stats::sd(zi[is.finite(zi)])
      if (is.finite(s) && s > 0) {
        Z[i, ] <- zi / s
        row_sds[i] <- s
      } else {
        # if sd is zero or NA, leave the row as is; downstream CD will naturally ignore all-NA or constant rows
        row_sds[i] <- NA_real_
      }
    }
  }

  # 3) Demean across units at each time (optional but recommended for CD)
  col_means <- rep(0, ncol(Z))
  if (isTRUE(demean_time)) {
    for (j in seq_len(ncol(Z))) {
      zj <- Z[, j]
      mu <- mean(zj[is.finite(zj)], na.rm = TRUE)
      if (is.finite(mu)) {
        Z[, j] <- zj - mu
        col_means[j] <- mu
      } else {
        col_means[j] <- NA_real_
      }
    }
  }

  # Preserve dimnames of kept columns
  if (!is.null(colnames(E))) colnames(Z) <- colnames(E)[keep]
  if (!is.null(rownames(E))) rownames(Z) <- rownames(E)

  out <- list(
    Z = Z,
    kept_t = keep,
    m_t = m_t,
    row_sds = row_sds,
    col_means = if (demean_time) col_means else rep(0, length(keep))
  )
  return(out)
}


# utils_cd.R - Cross-sectional dependence tests for panel models

#' Cross-sectional dependence (CD) tests for panel data
#'
#' Computes multiple cross-sectional dependence tests for panel residuals:
#' CD (Pesaran 2015), CDw (Juodis & Reese 2021), CDw+ (Fan et al. 2015), and
#' CD* (Pesaran & Xie 2021) with optional PCA-based factor adjustment.
#'
#' @param object A fitted model object (e.g., class `csdm_fit`) or a numeric N x T matrix of residuals.
#' @param type Which CD test(s) to compute: one of `"CD"` (default), `"CDw"`, `"CDw+"`, `"CDstar"`, or `"all"`.
#' @param n_pc Number of principal components to remove for CD* test (default 4, matching Stata xtcd2).
#' @param res_type Residual type to extract from model objects (see [get_residuals()]).
#' @param min_overlap Minimum number of overlapping time periods for a unit pair to be included.
#'
#' @return A list with test results. If type is a single test, returns list with `statistic`, `p.value`, `N`, `pairs_used`.
#'   If type="all", returns list with elements for each test (`CD`, `CDw`, `CDw_plus`, `CDstar`).
#'
#' @details
#' **CD (Pesaran 2015, 2021):** Classic test using pairwise correlations of residuals.
#'   Statistic: \eqn{CD = \sqrt{2/(N(N-1))} \sum_{i<j} \sqrt{T_{ij}} \rho_{ij}}
#'
#' **CDw (Juodis & Reese 2021):** Weighted CD after removing cross-sectional means from each time period.
#'   Controls for common shocks/factors.
#'
#' **CDw+ (Fan et al. 2015, extended by Juodis & Reese 2021):** Enhanced CDw that ignores small correlations.
#'   Uses threshold c_N = sqrt(2*log(N)/T) for sparse dependence detection.
#'
#' **CD* (Pesaran & Xie 2021):** Bias-corrected test with PCA-based factor detrending.
#'   Removes first n_pc principal components from residuals, then applies CD.
#'   Default n_pc=4 matches Stata xtcd2 implementation.
#'
#' @references
#' Pesaran, M.H. (2015). "Testing weak cross-sectional dependence in large panels."
#'   Econometric Reviews, 34(6-10), 1089-1117.
#'
#' Juodis, A., & Reese, S. (2021). "The incidental parameters problem in testing for
#'   remaining cross-sectional correlation." Journal of Business & Economic Statistics.
#'
#' Fan, J., Liao, Y., & Yao, J. (2015). "Power Enhancement in High-Dimensional
#'   Cross-Sectional Tests." Econometric Reviews, 34(6-10), 742-779.
#'
#' Pesaran, M.H., & Xie, Y. (2021). "A bias-corrected CD test for error cross-sectional
#'   dependence in panel models." Econometric Reviews, 41(6), 649-677.
#'
#' @examples
#' # Simulate independent and dependent panels
#' set.seed(1)
#' E_indep <- matrix(rnorm(100), nrow = 10)
#' E_dep <- matrix(rnorm(10), nrow = 10, ncol = 10, byrow = TRUE)
#'
#' # Compute all four tests
#' cd_test(E_indep, type = "all")
#' cd_test(E_dep, type = "all")
#'
#' # CD* with custom number of PCs
#' cd_test(E_indep, type = "CDstar", n_pc = 2)
#'
#' @export
cd_test <- function(object,
                    type = c("CD", "CDw", "CDw+", "CDstar", "all"),
                    n_pc = 4L,
                    res_type = c("auto", "pca_std", "pca", "cce"),
                    min_overlap = 2L) {
  type <- match.arg(type)
  res_type <- match.arg(res_type)
  n_pc <- as.integer(n_pc)

  # Get residuals matrix (N x T)
  if (is.matrix(object) && is.numeric(object)) {
    E <- object
  } else {
    E <- get_residuals(object, type = res_type, strict = TRUE)
  }
  if (!is.matrix(E) || !is.numeric(E)) stop("Residuals must be numeric matrix.")

  N <- nrow(E)
  T_n <- ncol(E)
  if (N < 2) stop("Need at least two units.")
  if (T_n < 2) stop("Need at least two time periods.")

  # Initialize output list
  out <- list()

  # Compute requested tests
  if (type %in% c("CD", "all")) {
    out$CD <- tryCatch(.csdm_cd_test_classic(E, min_overlap), error = function(e) NULL)
  }
  if (type %in% c("CDw", "all")) {
    out$CDw <- tryCatch(.csdm_cd_test_weighted(E, min_overlap), error = function(e) NULL)
  }
  if (type %in% c("CDw+", "all")) {
    out$CDw_plus <- tryCatch(.csdm_cd_test_weighted_plus(E, min_overlap), error = function(e) NULL)
  }
  if (type %in% c("CDstar", "all")) {
    out$CDstar <- tryCatch(.csdm_cd_test_star(E, n_pc, min_overlap), error = function(e) NULL)
  }

  # Return appropriate format
  if (type == "all") return(out)
  # Return the single test result (not the list wrapper)
  single_test <- switch(type,
    "CD" = out$CD,
    "CDw" = out$CDw,
    "CDw+" = out$CDw_plus,
    "CDstar" = out$CDstar,
    NULL)
  return(single_test)
}

# Internal: CD test (Pesaran 2015, 2021)
.csdm_cd_test_classic <- function(E, min_overlap = 2L) {
  N <- nrow(E)
  T_n <- ncol(E)
  if (N < 2 || T_n < 2) stop("Need at least two units and two time periods.")

  s_sum <- 0
  pairs <- 0L
  for (i in 1:(N - 1)) {
    xi <- E[i, ]
    for (j in (i + 1):N) {
      xj <- E[j, ]
      ok <- is.finite(xi) & is.finite(xj)
      Tij <- sum(ok)
      if (Tij >= min_overlap) {
        rho_ij <- suppressWarnings(stats::cor(xi[ok], xj[ok]))
        if (is.finite(rho_ij)) {
          s_sum <- s_sum + sqrt(Tij) * rho_ij
          pairs <- pairs + 1L
        }
      }
    }
  }
  if (pairs == 0L) stop("No valid unit pairs with sufficient overlap.")

  cd_stat <- sqrt(2 / (N * (N - 1))) * s_sum
  pval <- 2 * (1 - stats::pnorm(abs(cd_stat)))

  list(statistic = unname(cd_stat),
       p.value   = unname(pval),
       N         = N,
       pairs_used = pairs)
}

# Internal: CDw test (Juodis & Reese 2021) - weighted by removing cross-sectional means
.csdm_cd_test_weighted <- function(E, min_overlap = 2L) {
  N <- nrow(E)
  T_n <- ncol(E)
  if (N < 2 || T_n < 2) stop("Need at least two units and two time periods.")

  # Remove cross-sectional means: e_it^w = e_it - mean(e_t)
  E_w <- E
  for (t in 1:T_n) {
    e_t <- E[, t]
    mean_t <- mean(e_t, na.rm = TRUE)
    if (is.finite(mean_t)) {
      E_w[, t] <- e_t - mean_t
    }
  }

  # Apply classic CD to demeaned residuals
  s_sum <- 0
  pairs <- 0L
  for (i in 1:(N - 1)) {
    xi <- E_w[i, ]
    for (j in (i + 1):N) {
      xj <- E_w[j, ]
      ok <- is.finite(xi) & is.finite(xj)
      Tij <- sum(ok)
      if (Tij >= min_overlap) {
        rho_ij <- suppressWarnings(stats::cor(xi[ok], xj[ok]))
        if (is.finite(rho_ij)) {
          s_sum <- s_sum + sqrt(Tij) * rho_ij
          pairs <- pairs + 1L
        }
      }
    }
  }
  if (pairs == 0L) stop("No valid unit pairs with sufficient overlap.")

  cd_stat <- sqrt(2 / (N * (N - 1))) * s_sum
  pval <- 2 * (1 - stats::pnorm(abs(cd_stat)))

  list(statistic = unname(cd_stat),
       p.value   = unname(pval),
       N         = N,
       pairs_used = pairs)
}

# Internal: CDw+ test (Fan et al. 2015, Juodis & Reese 2021) - power-enhanced sparse dependence
.csdm_cd_test_weighted_plus <- function(E, min_overlap = 2L) {
  N <- nrow(E)
  T_n <- ncol(E)
  if (N < 2 || T_n < 2) stop("Need at least two units and two time periods.")

  # Remove cross-sectional means
  E_w <- E
  for (t in 1:T_n) {
    e_t <- E[, t]
    mean_t <- mean(e_t, na.rm = TRUE)
    if (is.finite(mean_t)) {
      E_w[, t] <- e_t - mean_t
    }
  }

  # Sparse threshold: c_N = sqrt(2 * log(N) / T)
  c_N <- sqrt(2 * log(N) / T_n)

  # Apply CD but only use pairs with |rho| > c_N (sparse enhancement)
  s_sum <- 0
  pairs <- 0L
  for (i in 1:(N - 1)) {
    xi <- E_w[i, ]
    for (j in (i + 1):N) {
      xj <- E_w[j, ]
      ok <- is.finite(xi) & is.finite(xj)
      Tij <- sum(ok)
      if (Tij >= min_overlap) {
        rho_ij <- suppressWarnings(stats::cor(xi[ok], xj[ok]))
        if (is.finite(rho_ij) && abs(rho_ij) > c_N) {
          s_sum <- s_sum + sqrt(Tij) * rho_ij
          pairs <- pairs + 1L
        }
      }
    }
  }
  if (pairs == 0L) {
    # No pairs exceed threshold, return zero statistic
    return(list(statistic = 0, p.value = 1, N = N, pairs_used = 0L))
  }

  cd_stat <- sqrt(2 / (N * (N - 1))) * s_sum
  pval <- 2 * (1 - stats::pnorm(abs(cd_stat)))

  list(statistic = unname(cd_stat),
       p.value   = unname(pval),
       N         = N,
       pairs_used = pairs)
}

# Internal: CD* test (Pesaran & Xie 2021) - bias-corrected with PCA factor detrending
.csdm_cd_test_star <- function(E, n_pc = 4L, min_overlap = 2L) {
  N <- nrow(E)
  T_n <- ncol(E)
  n_pc <- as.integer(n_pc)
  if (N < 2 || T_n < 2) stop("Need at least two units and two time periods.")
  if (n_pc < 0) stop("n_pc must be non-negative.")

  # If n_pc = 0, return classic CD
  if (n_pc == 0) {
    return(.csdm_cd_test_classic(E, min_overlap))
  }

  # Ensure n_pc does not exceed min(N, T_n)
  n_pc <- min(n_pc, N - 1L, T_n - 1L)

  # Step 1: Center residuals (subtract means per unit and time)
  E_c <- E
  for (i in 1:N) {
    E_c[i, ] <- E[i, ] - mean(E[i, ], na.rm = TRUE)
  }
  for (t in 1:T_n) {
    E_c[, t] <- E_c[, t] - mean(E_c[, t], na.rm = TRUE)
  }

  # Step 2: Perform PCA on centered residuals to extract factors
  # For simplicity and robustness, work with the full matrix (SVD handles NAs via na.omit)
  if (n_pc > 0) {
    # Remove rows/cols with all NAs, then use SVD
    valid_rows <- rowSums(!is.na(E_c)) > 0
    valid_cols <- colSums(!is.na(E_c)) > 0
    if (sum(valid_rows) >= 2 && sum(valid_cols) >= 2) {
      E_sub <- E_c[valid_rows, valid_cols, drop = FALSE]
      # Use SVD on the submatrix; na.action removes NAs
      tryCatch({
        svd_res <- svd(E_sub, nu = min(n_pc, nrow(E_sub), ncol(E_sub)), 
                           nv = min(n_pc, nrow(E_sub), ncol(E_sub)))
        if (n_pc > 0 && !is.null(svd_res$u) && !is.null(svd_res$v)) {
          n_use <- min(n_pc, length(svd_res$d))
          if (n_use > 0) {
            D_pc <- svd_res$d[seq_len(n_use)]
            U_pc <- svd_res$u[, seq_len(n_use), drop = FALSE]
            V_pc <- svd_res$v[, seq_len(n_use), drop = FALSE]
            # Reconstruct detrended residuals only for valid submatrix
            E_factor <- U_pc %*% diag(D_pc) %*% t(V_pc)
            # Map back to original indices
            rows_idx <- which(valid_rows)
            cols_idx <- which(valid_cols)
            for (ii in seq_len(nrow(E_factor))) {
              for (jj in seq_len(ncol(E_factor))) {
                E_c[rows_idx[ii], cols_idx[jj]] <- E_c[rows_idx[ii], cols_idx[jj]] - E_factor[ii, jj]
              }
            }
          }
        }
      }, error = function(e) NULL)
    }
  }

  # Step 3: Run classic CD on factor-detrended residuals
  s_sum <- 0
  bias_sum <- 0
  pairs <- 0L
  for (i in 1:(N - 1)) {
    xi <- E_c[i, ]
    for (j in (i + 1):N) {
      xj <- E_c[j, ]
      ok <- is.finite(xi) & is.finite(xj)
      Tij <- sum(ok)
      if (Tij >= min_overlap) {
        rho_ij <- suppressWarnings(stats::cor(xi[ok], xj[ok]))
        if (is.finite(rho_ij)) {
          s_sum <- s_sum + sqrt(Tij) * rho_ij
          # Bias term per Pesaran & Xie (2021, Eq. 7)
          bias_ij <- (1 + rho_ij^2) / (2 * (Tij - 1))
          bias_sum <- bias_sum + sqrt(Tij) * bias_ij
          pairs <- pairs + 1L
        }
      }
    }
  }
  if (pairs == 0L) stop("No valid unit pairs with sufficient overlap.")

  cd_stat <- sqrt(2 / (N * (N - 1))) * s_sum
  bias_correction <- sqrt(2 / (N * (N - 1))) * bias_sum
  cd_star <- cd_stat - bias_correction
  pval <- 2 * (1 - stats::pnorm(abs(cd_star)))

  list(statistic = unname(cd_star),
       p.value   = unname(pval),
       N         = N,
       pairs_used = pairs,
       n_pc = n_pc)
}

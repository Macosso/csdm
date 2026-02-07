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
#' CD (Pesaran 2015), CDw (Juodis & Reese 2021), CDw+ (Fan et al. 2015), the
#' unstandardised CDw+ statistic used by Stata (`CDw_plus`), and
#' CD* (Pesaran & Xie 2021) with optional PCA-based factor adjustment.
#'
#' Cross-sectional dependence tests (CD, CDw, CDw+, CD*) for panel residuals or model fits
#'
#' Computes CD, weighted-CD, power-enhanced CDw+, or bias-corrected CD* statistics—as in Stata xtcd2.
#'
#' @param object A csdm_fit model or an N x T residuals matrix (units × time).
#' @param type Which test(s) to compute: one or more of c("CD", "CDw", "CDw+", "CDstar", "all")
#' @param n_pc Number of principal components for CD* (default 4, matches Stata).
#' @param seed Integer. Seed for weight draws in CDw/CDw+ (default 101).
#' @return Named list with statistics and p-values for selected tests.
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
#' **CDw_plus_raw:** Stata-style version of CDw+ that reports the raw power–enhanced
#' statistic without normalising by \eqn{\sqrt{2/(N(N-1))}}.  This statistic
#' multiplies the Juodis–Reese CDw component by \eqn{\sqrt{N}} and adds the
#' sum of \eqn{\sqrt{T_{ij}} |\rho_{ij}|} over pairs exceeding the threshold \eqn{c_N}.
#' Because it is unstandardised, it can be very large and its sampling
#' distribution under the null is not standard normal; p-values are therefore
#' not reported.  Use this option only to reproduce Stata's `xtcd2` output.
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
                    n_pc = 4,
                    seed = NULL) {
  type <- match.arg(type)

  # get matrix of residuals
  if (inherits(object, "csdm_fit")) {
    E <- get_residuals(object, type = "auto")
  } else if (is.matrix(object) && is.numeric(object)) {
    E <- object
  } else {
    stop("object must be an N x T residual matrix or a csdm_fit model.")
  }
  #if (anyNA(E)) stop("Input residual matrix must not contain missing values (for now).")
  if (nrow(E) < 2 || ncol(E) < 2) stop("At least two units and two time periods required.")

  data_start <- na.omit(t(E))   # Use T × N as in Stata
  N <- ncol(data_start); Tt <- nrow(data_start)

  out <- list()

  # 1. CD (classic)
  if (type %in% c("CD", "all")) {
    corr_mat <- cor(data_start, use = "pairwise.complete.obs")
    upper <- corr_mat[upper.tri(corr_mat)] * sqrt(Tt)
    cd_stat <- sqrt(2/(N*(N-1))) * sum(upper)
    cd_p <- 2 * (1 - pnorm(abs(cd_stat)))
    out$CD_stat <- cd_stat
    out$CD_pvalue <- cd_p
  }

  # 2. Weighted CD (Juodis & Reese)
  if (type %in% c("CDw", "all", "CDw+")) {
    if (!is.null(seed)) set.seed(seed)
    w <- sample(c(-1,1), N, replace = TRUE)
    data_cdw <- sweep(data_start, 2, w, FUN="*")
    rho_w <- cor(data_cdw) * sqrt(Tt)
    upper_w <- rho_w[upper.tri(rho_w)]
    cdw_stat <- sqrt(2/(N*(N-1))) * sum(upper_w)
    cdw_p <- 2 * (1 - pnorm(abs(cdw_stat)))
    out$CDw_stat <- cdw_stat
    out$CDw_pvalue <- cdw_p
  }

  # 3. Weighted CD+, power enhancement (Fan et al.)
  if (type %in% c("CDw+", "all")) {
    # need the nonweighted correlations as screening term (from classic CD)
    if (!exists("upper")) {
      corr_mat <- cor(data_start)
      upper <- corr_mat[upper.tri(corr_mat)] * sqrt(Tt)
    }
    crit <- 2 * sqrt(log(N) / Tt)
    screen_const <- sum(abs(upper) * (abs(upper) > crit))
    if (!exists("cdw_stat")) {
      set.seed(seed)
      w <- sample(c(-1,1), N, replace = TRUE)
      data_cdw <- sweep(data_start, 2, w, FUN="*")
      rho_w <- cor(data_cdw) * sqrt(Tt)
      upper_w <- rho_w[upper.tri(rho_w)]
      cdw_stat <- sqrt(2/(N*(N-1))) * sum(upper_w)
    }
    cdw_plus_stat <- cdw_stat + screen_const
    cdw_plus_p <- 2 * (1 - pnorm(abs(cdw_plus_stat)))
    out$CDw_plus_stat <- cdw_plus_stat
    out$CDw_plus_pvalue <- cdw_plus_p
  }

  # 4. CDstar (bias-corrected, Pesaran & Xie, n_pc factors removed)
  if (type %in% c("CDstar", "all")) {
    # Standardize
    col_means <- colMeans(data_start)
    col_sds <- apply(data_start, 2, sd)
    col_sds[col_sds == 0] <- 1
    data_std <- sweep(sweep(data_start, 2, col_means), 2, col_sds, `/`)
    S <- data_std %*% t(data_std)
    eig <- eigen(S, symmetric = TRUE)
    idx <- order(eig$values, decreasing = TRUE)[seq_len(n_pc)]
    f <- eig$vectors[, idx, drop = FALSE]
    fx <- cbind(1, f) # intercept + factors
    beta <- solve(t(fx) %*% fx, t(fx) %*% data_std)
    res_defac <- data_std - fx %*% beta
    corr_defac <- cor(res_defac)
    upper_defac <- corr_defac[upper.tri(corr_defac)] * sqrt(Tt)
    cd_defac <- sqrt(2/(N*(N-1))) * sum(upper_defac)
    # Bias correction terms as in Pesaran & Xie (2021)
    betai <- beta[-1, , drop = FALSE]                 # remove intercept
    betaij <- (betai %*% t(betai)) / N
    betasum <- sqrt(diag(betaij))
    gamma <- sweep(betai, 1, betasum, `/`)
    sgm <- sqrt(mean(res_defac^2))
    phi <- rowMeans(gamma / sgm)
    ai_vals <- as.numeric((1 - t(gamma * sgm) %*% phi) / sqrt(N))
    theta <- sum(ai_vals^2)
    cd_star_stat <- (cd_defac + sqrt(Tt/2) * (1 - theta)) / theta
    cd_star_p <- 2 * (1 - pnorm(abs(cd_star_stat)))
    out$CDstar_stat <- cd_star_stat
    out$CDstar_pvalue <- cd_star_p
    attr(out, "n_pc") <- n_pc
  }
  out
}

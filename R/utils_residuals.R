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
  N <- nrow(E)
  T_n <- ncol(E)

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

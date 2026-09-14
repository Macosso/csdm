# utils_vcov.R

#' Cluster-robust variance-covariance for OLS
#'
#' @description
#' Computes one- or two-way cluster-robust vcov for an OLS design using the
#' Liang-Zeger "meat" and Cameron-Gelbach-Miller inclusion-exclusion for two-way clustering.
#'
#' @param X Numeric design matrix (n x k) used in OLS.
#' @param u Numeric residual vector (length n).
#' @param cluster One of:
#'   \itemize{
#'     \item a vector (length n) of cluster ids for one-way clustering; or
#'     \item a data.frame/list with two vectors (each length n) for two-way clustering.
#'   }
#' @param df_correction Logical; apply small-sample corrections. Default \code{TRUE}.
#' @param type Character, one of \code{"oneway"} or \code{"twoway"}.
#'
#' @returns A \code{k x k} variance-covariance matrix.
#' @keywords internal
#' @export
cluster_vcov <- function(X, u, cluster, df_correction = TRUE,
                         type = c("oneway", "twoway")) {
  type <- match.arg(type)
  X <- as.matrix(X)
  u <- as.numeric(u)
  n <- nrow(X); k <- ncol(X)

  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))

  # helper: meat for one set of clusters
  .meat_oneway <- function(gids) {
    gids <- as.vector(gids)
    # compute S_g = X_g' u_g, accumulate S_g S_g'
    ug <- split(u, gids)
    Xg <- split.data.frame(as.data.frame(X), gids)
    # accumulate efficiently
    acc <- matrix(0, k, k)
    for (g in intersect(names(ug), names(Xg))) {
      Xg_mat <- as.matrix(Xg[[g]])
      ug_vec <- as.numeric(ug[[g]])
      Sg <- crossprod(Xg_mat, ug_vec)      # k x 1
      acc <- acc + tcrossprod(Sg)          # k x k
    }
    # small-sample correction (Liang-Zeger + Bell-McCaffrey style)
    if (df_correction) {
      G <- length(unique(gids))
      if (G <= 1L) warning("Only one cluster found; cluster correction not meaningful.")
      c1 <- G/(G - 1)
      c2 <- (n - 1)/(n - k)
      acc <- acc * (c1 * c2)
    }
    acc
  }

  if (type == "oneway") {
    if (is.data.frame(cluster) || is.list(cluster)) {
      cluster <- unlist(cluster, use.names = FALSE)
    }
    meat <- .meat_oneway(cluster)
  } else { # twoway
    if (!(is.list(cluster) || is.data.frame(cluster)) || length(cluster) != 2L) {
      stop("For type = 'twoway', 'cluster' must be a list/data.frame of two cluster vectors.")
    }
    g1 <- as.vector(cluster[[1]])
    g2 <- as.vector(cluster[[2]])
    g12 <- interaction(g1, g2, drop = TRUE)

    meat1  <- .meat_oneway(g1)
    meat2  <- .meat_oneway(g2)
    meat12 <- .meat_oneway(g12)

    # Inclusion-exclusion
    meat <- meat1 + meat2 - meat12
  }

  vc <- XtX_inv %*% meat %*% XtX_inv
  dimnames(vc) <- list(colnames(X), colnames(X))
  vc
}


#' Heteroskedasticity-robust (HC) sandwich variance-covariance for OLS
#'
#' @description
#' Computes White/Huber HC0-HC3 sandwich vcov for an OLS design.
#'
#' @param X Numeric design matrix (n x k) used in OLS.
#' @param u Numeric residual vector (length n).
#' @param type Character; one of \code{"HC0"}, \code{"HC1"}, \code{"HC2"}, \code{"HC3"}.
#'
#' @returns A \code{k x k} variance-covariance matrix.
#' @keywords internal
#' @export
sandwich_vcov <- function(X, u, type = c("HC0", "HC1", "HC2", "HC3")) {
  type <- match.arg(type)
  design <- .csdm_ols_design(X, u)
  X <- design$X; u <- design$u
  n <- nrow(X); k <- ncol(X)
  h <- rowSums(qr.Q(design$qr)^2)
  if (type %in% c("HC2", "HC3") && any(1 - h <= .Machine$double.eps^0.5)) {
    stop("HC2/HC3 are undefined for leverage-one observations.")
  }
  adjusted <- switch(type, HC0 = u, HC1 = u * sqrt(n / (n - k)),
    HC2 = u / sqrt(1 - h), HC3 = u / (1 - h))
  meat <- crossprod(X * adjusted)
  vc <- design$bread %*% meat %*% design$bread
  dimnames(vc) <- list(colnames(X), colnames(X))
  vc
}


#' Variance-covariance of Mean-Group (MG) averages
#'
#' @description
#' Computes the covariance matrix of the MG estimator \eqn{\bar{\beta} = N^{-1} \sum_i \hat\beta_i},
#' using the cross-sectional covariance of unit-specific slopes and dividing by the
#' effective sample sizes per coefficient (handles missing entries per unit).
#'
#' @param beta_i Numeric matrix of unit-specific coefficients (\eqn{N x K});
#'   rows = units, columns = coefficients. May contain \code{NA}s.
#' @param weights Optional numeric vector of length N with nonnegative weights
#'   summing to 1. If \code{NULL}, uses equal weights.
#' @param pairwise Logical; use pairwise-complete covariances across units (default \code{TRUE}).
#'
#' @returns A \code{K x K} covariance matrix for the MG mean, with
#'   \code{dimnames} inherited from \code{colnames(beta_i)}.
#' @details
#' Weights are fixed relative weights, normalized on the retained units.
#' The calculation assumes independent unit estimates with a common covariance:
#' the weighted sample covariance is divided by \eqn{1-\sum_i w_i^2}, then
#' multiplied by \eqn{\sum_i w_i^2}. Equal weights give sample covariance divided by N.
#' With missing coefficients, use pairwise=FALSE to select one complete-unit
#' sample. Pairwise covariance with missing coefficients is not implemented.
#' This is not an inverse-variance pooled estimator or a general covariance
#' estimator for arbitrary unit-specific covariance matrices.
#'
#' @keywords internal
#' @export
pooled_vcov <- function(beta_i, weights = NULL, pairwise = TRUE) {
  B <- as.matrix(beta_i)
  .csdm_flag(pairwise, "pairwise")
  if (!is.numeric(B) || !ncol(B)) stop("'beta_i' must be a numeric matrix with columns.")
  w <- if (is.null(weights)) rep(1, nrow(B)) else weights
  if (!is.numeric(w) || length(w) != nrow(B) || any(!is.finite(w) | w < 0) || sum(w) <= 0) {
    stop("Weights must be finite, nonnegative, aligned with units, and have positive sum.")
  }
  complete <- rowSums(!is.finite(B)) == 0
  if (pairwise && any(!complete & w > 0)) stop("Pairwise covariance with missing coefficients is not implemented; use pairwise=FALSE for a common complete-unit sample.")
  keep <- complete & w > 0
  B <- B[keep, , drop = FALSE]
  w <- w[keep]
  if (length(w) < 2L) stop("At least two complete units with positive weights are required.")
  w <- w / sum(w)
  concentration <- sum(w^2)
  dev <- sweep(B, 2L, colSums(B * w), "-")
  # Unbiased common-covariance estimate, then Var(sum(w_i * beta_i)).
  V <- crossprod(dev * sqrt(w)) * concentration / (1 - concentration)
  dimnames(V) <- list(colnames(B), colnames(B))
  V
}


.csdm_ols_design <- function(X, u) {
  X <- as.matrix(X)
  if (!is.numeric(X) || !ncol(X) || !is.numeric(u) || length(u) != nrow(X) ||
      any(!is.finite(X)) || any(!is.finite(u))) stop("Supply a finite numeric design and aligned residual vector.")
  q <- qr(X)
  if (q$rank < ncol(X) || nrow(X) <= ncol(X)) stop("A full-rank design with positive residual degrees of freedom is required.")
  inverse <- chol2inv(qr.R(q))
  bread <- inverse[order(q$pivot), order(q$pivot), drop = FALSE]
  list(X = X, u = as.numeric(u), qr = q, bread = bread)
}

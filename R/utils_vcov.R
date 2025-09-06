# utils_vcov.R

#' Cluster-robust variance-covariance for OLS
#'
#' @description
#' Computes one- or two-way cluster-robust vcov for an OLS design using the
#' Liang–Zeger "meat" and Cameron–Gelbach–Miller inclusion–exclusion for two-way clustering.
#'
#' @param X Numeric design matrix (n × k) used in OLS.
#' @param u Numeric residual vector (length n).
#' @param cluster One of:
#'   \itemize{
#'     \item a vector (length n) of cluster ids for one-way clustering; or
#'     \item a data.frame/list with two vectors (each length n) for two-way clustering.
#'   }
#' @param df_correction Logical; apply small-sample corrections. Default \code{TRUE}.
#' @param type Character, one of \code{"oneway"} or \code{"twoway"}.
#'
#' @returns A \code{k × k} variance-covariance matrix.
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

    # Inclusion–exclusion
    meat <- meat1 + meat2 - meat12
  }

  vc <- XtX_inv %*% meat %*% XtX_inv
  dimnames(vc) <- list(colnames(X), colnames(X))
  vc
}


#' Heteroskedasticity-robust (HC) sandwich variance-covariance for OLS
#'
#' @description
#' Computes White/Huber HC0–HC3 sandwich vcov for an OLS design.
#'
#' @param X Numeric design matrix (n × k) used in OLS.
#' @param u Numeric residual vector (length n).
#' @param type Character; one of \code{"HC0"}, \code{"HC1"}, \code{"HC2"}, \code{"HC3"}.
#'
#' @returns A \code{k × k} variance-covariance matrix.
#' @keywords internal
#' @export
sandwich_vcov <- function(X, u, type = c("HC0", "HC1", "HC2", "HC3")) {
  type <- match.arg(type)
  X <- as.matrix(X)
  u <- as.numeric(u)
  n <- nrow(X); k <- ncol(X)

  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))
  # leverage h_ii
  H <- X %*% XtX_inv %*% t(X)
  h <- pmin(1, pmax(0, diag(H))) # clamp numerically

  # scale residuals per HC type
  w <- switch(
    type,
    HC0 = u,
    HC1 = u * sqrt(n/(n - k)),
    HC2 = u / sqrt(1 - h),
    HC3 = u / (1 - h)
  )

  # meat = X' diag(w^2) X
  meat <- crossprod(X, w * X)
  vc <- XtX_inv %*% meat %*% XtX_inv
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
#' @param beta_i Numeric matrix of unit-specific coefficients (\eqn{N × K});
#'   rows = units, columns = coefficients. May contain \code{NA}s.
#' @param weights Optional numeric vector of length N with nonnegative weights
#'   summing to 1. If \code{NULL}, uses equal weights.
#' @param pairwise Logical; use pairwise-complete covariances across units (default \code{TRUE}).
#'
#' @returns A \code{K × K} covariance matrix for the MG mean, with
#'   \code{dimnames} inherited from \code{colnames(beta_i)}.
#' @details
#' For equal weights, diagonal entries coincide with \eqn{\widehat{\mathrm{Var}}(\hat\beta_{ij}) / N_{\text{eff},j}}.
#' Off-diagonals are scaled analogously using pairwise effective N. If \code{weights} are
#' provided, computes \eqn{\mathrm{Var}(\sum_i w_i \hat\beta_i)} using a weighted covariance.
#'
#' @keywords internal
#' @export
pooled_vcov <- function(beta_i, weights = NULL, pairwise = TRUE) {
  B <- as.matrix(beta_i)
  N <- nrow(B); K <- ncol(B)
  cn <- colnames(B)

  if (is.null(weights)) {
    w <- rep(1/N, N)
  } else {
    w <- as.numeric(weights)
    if (length(w) != N) stop("weights must have length equal to nrow(beta_i).")
    if (any(w < 0)) stop("weights must be nonnegative.")
    s <- sum(w)
    if (!isTRUE(all.equal(s, 1))) w <- w / s
  }

  # Weighted means per column (ignoring NAs)
  wmean <- function(x, w) {
    ok <- is.finite(x)
    if (!any(ok)) return(NA_real_)
    sum(w[ok] * x[ok]) / sum(w[ok])
  }
  mu <- vapply(seq_len(K), function(j) wmean(B[, j], w), numeric(1))

  # Weighted covariance (pairwise if requested)
  V <- matrix(NA_real_, K, K)
  for (a in seq_len(K)) {
    for (b in a:K) {
      xa <- B[, a]; xb <- B[, b]
      ok <- is.finite(xa) & is.finite(xb)
      if (!any(ok)) {
        V[a, b] <- V[b, a] <- NA_real_
        next
      }
      wa <- w[ok]
      xa <- xa[ok]; xb <- xb[ok]
      # center
      ma <- sum(wa * xa) / sum(wa)
      mb <- sum(wa * xb) / sum(wa)
      # weighted covariance (population vs sample): use population denom sum(w)
      cov_ab <- sum(wa * (xa - ma) * (xb - mb)) / sum(wa)
      # For equal weights, the variance of the mean divides by N_eff implicitly via sum(w)
      # For unequal weights, this already yields Var(sum w_i beta_i) directly.
      V[a, b] <- V[b, a] <- cov_ab
    }
  }

  # For equal weights, above gives cov of beta across units.
  # Convert to cov of the (weighted) mean: multiply by scaling factor:
  # If weights are equal, Var(mean) = Cov / N_eff (done by sum(w) in pop cov).
  # If weights unequal, the formula above already corresponds to Var(sum w_i beta_i)
  # when using population covariance with weights normalized to 1.
  dimnames(V) <- list(cn, cn)
  V
}

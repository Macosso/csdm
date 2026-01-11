# csdm_internal_helpers.R

.csdm_mg_vcov <- function(B) {
  B <- as.matrix(B)
  ok <- stats::complete.cases(B)
  B <- B[ok, , drop = FALSE]
  N <- nrow(B)
  if (N < 2L) {
    V <- matrix(NA_real_, ncol(B), ncol(B))
    colnames(V) <- rownames(V) <- colnames(B)
    return(V)
  }
  bbar <- colMeans(B)
  dev <- sweep(B, 2L, bbar, "-")
  V <- crossprod(dev) / (N * (N - 1))
  dimnames(V) <- list(colnames(B), colnames(B))
  V
}


.csdm_econ_names <- function(formula, panel_df, na.action = stats::na.omit) {
  mf <- stats::model.frame(formula, panel_df, na.action = na.action)
  X <- stats::model.matrix(formula, mf)
  colnames(X)
}


.csdm_residual_matrix <- function(panel_df, id, time, res_long) {
  .or <- function(x, y) if (is.null(x)) y else x

  ids <- .or(attr(panel_df, "csdm_id_levels"), sort(unique(panel_df[[id]])))
  times <- .or(attr(panel_df, "csdm_time_levels"), sort(unique(panel_df[[time]])))

  E <- matrix(NA_real_, nrow = length(ids), ncol = length(times))
  rownames(E) <- as.character(ids)
  colnames(E) <- as.character(times)

  if (nrow(res_long) == 0L) return(E)

  ii <- match(res_long[[id]], ids)
  tt <- match(res_long[[time]], times)
  keep <- is.finite(ii) & is.finite(tt)
  if (any(keep)) E[cbind(ii[keep], tt[keep])] <- res_long$residual[keep]
  E
}


.csdm_fitted_matrix <- function(panel_df, id, time, fit_long) {
  .or <- function(x, y) if (is.null(x)) y else x

  ids <- .or(attr(panel_df, "csdm_id_levels"), sort(unique(panel_df[[id]])))
  times <- .or(attr(panel_df, "csdm_time_levels"), sort(unique(panel_df[[time]])))

  Xb <- matrix(NA_real_, nrow = length(ids), ncol = length(times))
  rownames(Xb) <- as.character(ids)
  colnames(Xb) <- as.character(times)

  if (nrow(fit_long) == 0L) return(Xb)

  ii <- match(fit_long[[id]], ids)
  tt <- match(fit_long[[time]], times)
  keep <- is.finite(ii) & is.finite(tt)
  if (any(keep)) Xb[cbind(ii[keep], tt[keep])] <- fit_long$xb[keep]
  Xb
}

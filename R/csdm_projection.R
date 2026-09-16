.csdm_unit_regression <- function(X, y, economic, tol = 1e-7) {
  nuisance <- setdiff(colnames(X), economic)
  Z <- X[, economic, drop = FALSE]
  H <- X[, nuisance, drop = FALSE]
  if (ncol(H)) {
    norms <- sqrt(colSums(H^2))
    H <- H[, norms > 0, drop = FALSE]
    if (ncol(H)) {
      H <- sweep(H, 2L, sqrt(colSums(H^2)), "/")
      qh <- qr(H, tol = tol)
      H <- qr.Q(qh)[, seq_len(qh$rank), drop = FALSE]
    }
  }
  projected <- if (ncol(H)) Z - H %*% crossprod(H, Z) else Z
  norms <- sqrt(colSums(Z^2))
  scaled <- sweep(projected, 2L, ifelse(norms > 0, norms, 1), "/")
  zero <- sqrt(colSums(scaled^2)) <= tol
  qx <- qr(scaled, tol = tol)
  aliased <- economic[zero]
  if (qx$rank < ncol(Z)) aliased <- union(aliased, economic[qx$pivot[seq.int(qx$rank + 1L, ncol(Z))]])
  df <- nrow(X) - ncol(H) - ncol(Z)
  if (length(aliased) || df <= 0L) {
    return(list(reason = if (length(aliased)) "unidentified economic terms" else "insufficient residual degrees of freedom",
      aliased = aliased, rank = ncol(H) + qx$rank, df.residual = df))
  }
  if (ncol(H)) colnames(H) <- paste0(".csa_basis", seq_len(ncol(H)))
  fit <- stats::lm.fit(cbind(Z, H), y, tol = tol)
  fit$reason <- NA_character_
  fit$aliased <- character()
  fit
}

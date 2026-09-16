test_that("CDstar preserves unit-specific scales", {
  set.seed(482)
  E <- outer(seq(.1, 3, length.out = 20), rnorm(100)) + matrix(rnorm(2000), 20)
  A <- scale(t(E)); N <- ncol(A); TT <- nrow(A)
  F <- cbind(1, eigen(tcrossprod(A), symmetric = TRUE)$vectors[, 1, drop = FALSE])
  B <- solve(crossprod(F), crossprod(F, A)); U <- A - F %*% B
  G <- B[-1, , drop = FALSE] / sqrt(mean(B[-1, ]^2))
  sigma <- sqrt(colMeans(U^2))
  phi <- mean(G / rep(sigma, each = nrow(G)))
  a <- 1 - as.numeric(G) * sigma * phi
  correction <- mean(a^2)
  cd <- sqrt(2 * TT / (N * (N - 1))) * sum(cor(U)[upper.tri(cor(U))])
  expected <- (cd + sqrt(TT / 2) * (1 - correction)) / correction
  expect_equal(cd_test(E, type = "CDstar", n_pc = 1)$tests$CDstar$statistic, expected, tolerance = 1e-9)
  expect_error(cd_test(E, type = "CDstar", n_pc = 20), "n_pc")
  expect_error(cd_test(E, type = "CDstar", n_pc = 1.5), "integers")
})

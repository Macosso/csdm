test_that("weighted CD matches the paper-defined cross-product calculation", {
  set.seed(729)
  E <- matrix(rnorm(1200), 12) * seq(.5, 2, length.out = 12)
  set.seed(84); w <- sample(c(-1, 1), nrow(E), replace = TRUE)
  U <- E - rowMeans(E)
  products <- crossprod(t(U * w))
  expected <- sqrt(2 / (ncol(E) * nrow(E) * (nrow(E) - 1))) *
    sum(products[upper.tri(products)]) / mean(U^2)
  a <- cd_test(E, type = "CDw+", seed = 84)
  expect_equal(a$tests$CDw$statistic, expected)
  rho <- abs(cor(t(E))[upper.tri(cor(t(E)))])
  enhancement <- sum(rho[rho > 2 * sqrt(log(nrow(E)) / ncol(E))])
  expect_equal(a$tests$CDw_plus$statistic, expected + enhancement)
  E[1, 1] <- NA
  expect_error(cd_test(E, type = "CDw+"), "balanced")
})

# Unit tests for cd_test (classic and bias-corrected)

test_that("cd_test returns classic and CDstar stats for synthetic data", {
  set.seed(123)
  N <- 10; T <- 20
  # Independent panels
  E_indep <- matrix(rnorm(N * T), nrow = N)
  res1 <- cd_test(E_indep, type = "all")
  expect_true(!is.null(res1$classic))
  expect_true(!is.null(res1$CDstar))
  expect_true(abs(res1$classic$statistic) < 2)
  expect_true(abs(res1$CDstar$statistic) < 2)

  # Perfect dependence
  E_dep <- matrix(rnorm(T), nrow = N, ncol = T, byrow = TRUE)
  res2 <- cd_test(E_dep, type = "all")
  expect_true(res2$classic$statistic > 5)
  expect_true(res2$CDstar$statistic > 5)

  # Small N, small T
  E_small <- matrix(rnorm(2 * 3), nrow = 2)
  res3 <- cd_test(E_small, type = "all")
  expect_true(!is.null(res3$classic))
  expect_true(!is.null(res3$CDstar))

  # Missing data
  E_miss <- E_indep
  E_miss[1, 1:5] <- NA
  res4 <- cd_test(E_miss, type = "all")
  expect_true(!is.null(res4$classic))
  expect_true(!is.null(res4$CDstar))
})

test_that("cd_test default is classic", {
  N <- 5; T <- 10
  E <- matrix(rnorm(N * T), nrow = N)
  res <- cd_test(E)
  expect_named(res, c("statistic", "p.value", "N", "pairs_used"))
})

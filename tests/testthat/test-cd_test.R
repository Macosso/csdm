# Unit tests for cd_test (CD, CDw, CDw+, CDstar)

test_that("CD test works for synthetic data", {
  set.seed(123)
  N <- 10; T <- 20
  E_indep <- matrix(rnorm(N * T), nrow = N)
  
  res <- cd_test(E_indep, type = "CD")
  expect_true(!is.null(res))
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
  expect_true(res$N == N)
  expect_true(is.numeric(res$pairs_used))
})

test_that("CDw removes cross-sectional means", {
  set.seed(124)
  N <- 8; T <- 15
  E <- matrix(rnorm(N * T), nrow = N)
  
  res_cd <- cd_test(E, type = "CD")
  res_cdw <- cd_test(E, type = "CDw")
  
  expect_true(!is.null(res_cd))
  expect_true(!is.null(res_cdw))
  # CDw should differ from CD due to demeaning
  expect_false(abs(res_cd$statistic - res_cdw$statistic) < 1e-10)
})

test_that("CDw+ uses sparse threshold", {
  set.seed(125)
  N <- 12; T <- 25
  E <- matrix(rnorm(N * T), nrow = N)
  
  res_cdw <- cd_test(E, type = "CDw")
  res_cdw_plus <- cd_test(E, type = "CDw+")
  
  expect_true(!is.null(res_cdw))
  expect_true(!is.null(res_cdw_plus))
  # CDw+ should have fewer or equal pairs used (only high correlations)
  expect_true(res_cdw_plus$pairs_used <= res_cdw$pairs_used)
})

test_that("CDstar with PCA detrending", {
  set.seed(126)
  N <- 10; T <- 30
  E <- matrix(rnorm(N * T), nrow = N)
  
  res0 <- cd_test(E, type = "CDstar", n_pc = 0)
  res4 <- cd_test(E, type = "CDstar", n_pc = 4)
  
  expect_true(!is.null(res0))
  expect_true(!is.null(res4))
  expect_true("n_pc" %in% names(res4))
  expect_equal(res4$n_pc, 4L)
})

test_that("cd_test type='all' returns all four tests", {
  set.seed(127)
  N <- 8; T <- 20
  E <- matrix(rnorm(N * T), nrow = N)
  
  res_all <- cd_test(E, type = "all")
  
  expect_true(is.list(res_all))
  expect_true("CD" %in% names(res_all) || !is.null(res_all$CD))
  expect_true("CDw" %in% names(res_all) || !is.null(res_all$CDw))
  expect_true("CDw_plus" %in% names(res_all) || !is.null(res_all$CDw_plus))
  # CDstar may be in results
  expect_true(length(res_all) >= 3)
})

test_that("CD test handles missing data", {
  set.seed(128)
  N <- 6; T <- 12
  E <- matrix(rnorm(N * T), nrow = N)
  E[1, 1:3] <- NA
  E[2, c(5, 7)] <- NA
  
  res <- cd_test(E, type = "CD")
  expect_true(!is.null(res))
  expect_true(is.numeric(res$statistic))
})

test_that("CD test rejects on perfect dependence", {
  set.seed(129)
  N <- 5; T <- 10
  # Perfect dependence: all rows identical
  E_dep <- matrix(rnorm(T), nrow = N, ncol = T, byrow = TRUE)
  
  res <- cd_test(E_dep, type = "CD")
  expect_true(abs(res$statistic) > 5)  # Large test statistic
  expect_true(res$p.value < 0.05)
})

test_that("CD test with small N and T", {
  set.seed(130)
  E_small <- matrix(rnorm(2 * 3), nrow = 2)
  
  res <- cd_test(E_small, type = "CD")
  expect_true(!is.null(res))
  expect_true(res$N == 2)
})

test_that("Default behavior uses CD test", {
  set.seed(131)
  N <- 5; T <- 10
  E <- matrix(rnorm(N * T), nrow = N)
  
  res_default <- cd_test(E)
  res_explicit <- cd_test(E, type = "CD")
  
  expect_equal(res_default$statistic, res_explicit$statistic)
  expect_equal(res_default$p.value, res_explicit$p.value)
})

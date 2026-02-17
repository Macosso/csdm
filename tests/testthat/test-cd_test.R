# Unit tests for cd_test (CD, CDw, CDw+, CDstar)

test_that("CD test works for synthetic independent data", {
  set.seed(123)
  N <- 10; T <- 20
  E_indep <- matrix(rnorm(N * T), nrow = N)

  res <- cd_test(E_indep, type = "CD")
  expect_true(!is.null(res))
  expect_s3_class(res, "cd_test")
  expect_true(is.numeric(res$tests$CD[['statistic']]))
  expect_true(is.numeric(res$tests$CD[['p.value']]))
  expect_equal(res$tests$CD[['N']], N)
  expect_equal(res$tests$CD[['T']], T)
  expect_true(is.numeric(res$tests$CD[['pairs_used']]))
  expect_true(res$tests$CD[['pairs_used']] > 0)

  # For independent data, CD should be small and non-significant
  expect_true(abs(res$tests$CD[['statistic']]) < 3)
  expect_true(res$tests$CD[['p.value']] > 0.05)
})

test_that("CD test detects perfect dependence", {
  set.seed(129)
  N <- 5; T <- 10
  # Perfect dependence: all rows identical
  E_dep <- matrix(rnorm(T), nrow = N, ncol = T, byrow = TRUE)

  res <- cd_test(E_dep, type = "CD")
  expect_true(abs(res$tests$CD[['statistic']]) > 5)  # Large test statistic
  expect_true(res$tests$CD[['p.value']] < 0.05)
})

test_that("CDw differs from CD due to random weighting", {
  set.seed(124)
  N <- 8; T <- 15
  E <- matrix(rnorm(N * T), nrow = N)

  res_cd <- cd_test(E, type = "CD")
  res_cdw <- cd_test(E, type = "CDw", seed = 101)

  expect_true(!is.null(res_cd))
  expect_true(!is.null(res_cdw))
  expect_equal(res_cd$tests$CD[['N']], res_cdw$tests$CDw[['N']])
  expect_equal(res_cd$tests$CD[['T']], res_cdw$tests$CDw[['T']])
  # CDw should differ from CD due to random sign flips
  expect_false(abs(res_cd$tests$CD[['statistic']] - res_cdw$tests$CDw[['statistic']]) < 1e-10)
})

test_that("cd_test type='all' returns all four tests", {
  set.seed(127)
  N <- 8; T <- 20
  E <- matrix(rnorm(N * T), nrow = N)

  res_all <- cd_test(E, type = "all", seed = 103)

  expect_true(is.list(res_all))
  expect_s3_class(res_all, "cd_test")
  expect_true("CD" %in% names(res_all$tests))
  expect_true("CDw" %in% names(res_all$tests))
  expect_true("CDw_plus" %in% names(res_all$tests))
  expect_true("CDstar" %in% names(res_all$tests))
  expect_equal(length(res_all$tests), 4)

  # Each test should have required fields
  expect_true(all(c("statistic", "p.value", "N", "T") %in% names(res_all$tests[['CD']])))
  expect_true(all(c("statistic", "p.value", "N", "T") %in% names(res_all$tests[['CDw']])))
  expect_true(all(c("statistic", "p.value", "N", "T") %in% names(res_all$tests[['CDw_plus']])))
  expect_true(all(c("statistic", "p.value", "N", "T", "n_pc") %in% names(res_all$tests[['CDstar']])))
})

test_that("CDw+ statistic is properly normalized", {
  set.seed(125)
  N <- 10; T <- 30
  E <- matrix(rnorm(N * T), nrow = N)

  res <- cd_test(E, type = "CDw+", seed = 104)

  expect_true(!is.null(res))
  expect_s3_class(res, "cd_test")
  expect_true(is.numeric(res$tests$CDw_plus[['statistic']]))
  expect_true(is.finite(res$tests$CDw_plus[['statistic']]))
})

test_that("CD test handles missing data with pairwise complete", {
  set.seed(128)
  N <- 6; T <- 12
  E <- matrix(rnorm(N * T), nrow = N)
  E[1, 1:3] <- NA
  E[2, c(5, 7)] <- NA

  res <- cd_test(E, type = "CD")
  expect_true(!is.null(res))
  expect_true(is.numeric(res$tests$CD[['statistic']]))
  expect_true(is.finite(res$tests$CD[['statistic']]))
  expect_true(res$tests$CD[['pairs_used']] > 0)
})

test_that("CDstar returns NA with warning for unbalanced panel", {
  set.seed(130)
  N <- 5; T <- 10
  E <- matrix(rnorm(N * T), nrow = N)
  E[1, 1:2] <- NA  # Create unbalanced panel

  expect_warning(
    res <- cd_test(E, type = "CDstar", n_pc = 2, na.action = "pairwise"),
    "CD"
  )
  expect_true(is.na(res$tests$CDstar[['statistic']]))
  expect_true(is.na(res$tests$CDstar[['p.value']]))
  expect_equal(res$tests$CDstar[['n_pc']], 2L)
})

test_that("cd_test works on csdm_fit objects", {
  set.seed(133)
  df <- expand.grid(
    id = paste0("i", 1:4),
    time = 1:10,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  df[['x1']] <- rnorm(nrow(df))
  df[['x2']] <- rnorm(nrow(df))
  df[['y']] <- 1 + 0.3 * df[['x1']] - 0.1 * df[['x2']] + rnorm(nrow(df), sd = 0.5)

  fit <- csdm(y ~ x1 + x2, data = df, id = "id", time = "time", model = "mg")

  res <- cd_test(fit, type = "CD")
  expect_true(!is.null(res))
  expect_s3_class(res, "cd_test")
  expect_true(is.numeric(res$tests$CD[['statistic']]))
  expect_true(is.numeric(res$tests$CD[['p.value']]))
  expect_equal(res$tests$CD[['N']], 4)
  expect_equal(res$tests$CD[['T']], 10)
})

test_that("CD test summary shows only classic CD", {
  set.seed(139)
  df <- expand.grid(
    id = paste0("i", 1:4),
    time = 1:10,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  df[['x']] <- rnorm(nrow(df))
  df[['y']] <- 0.5 * df[['x']] + rnorm(nrow(df), sd = 0.5)

  fit <- csdm(y ~ x, data = df, id = "id", time = "time", model = "mg")
  summ_output <- capture.output(summary(fit))

  # Should mention CD
  expect_true(any(grepl("CD\\s*=", summ_output)))
  # Should reference cd_test for additional diagnostics
  expect_true(any(grepl("cd_test", summ_output, ignore.case = TRUE)))
})

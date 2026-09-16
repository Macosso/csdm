test_that("fitting does not consume random numbers", {
  d <- panel_fixture()
  set.seed(729); before <- .Random.seed
  csdm(y ~ x, d, "id", "time")
  expect_identical(.Random.seed, before)
})

test_that("seeded weighted diagnostics preserve the caller RNG", {
  set.seed(44); E <- matrix(rnorm(600), 10); before <- .Random.seed
  a <- cd_test(E, type = "CDw", seed = 2)
  expect_identical(.Random.seed, before)
  expect_equal(a$tests, cd_test(E, type = "CDw", seed = 2)$tests)
})

test_that("classical CD retains pairwise observations and reports exclusions", {
  set.seed(2); E <- matrix(rnorm(300), 5)
  E[1, ] <- NA; E[2, 1:3] <- NA
  a <- cd_test(E)
  expect_equal(a$N, 4)
  expect_equal(a$T, 60)
  expect_identical(a$excluded_units, 1L)
  expect_true(is.finite(a$tests$CD$statistic))
  expect_error(cd_test(E, min_overlap = 1.5), "integers")
})

test_that("empty periods are removed before residual balance is assessed", {
  set.seed(25)
  E <- matrix(rnorm(8 * 20), 8, dimnames = list(NULL, 2001:2020))
  with_empty_period <- cbind(`2000` = NA_real_, E)

  expected <- cd_test(E, type = "all", seed = 7)
  actual <- cd_test(with_empty_period, type = "all", seed = 7)

  expect_equal(actual$tests, expected$tests)
  expect_equal(actual$T, 20L)
  expect_identical(actual$excluded_times, "2000")
  expect_identical(actual$kept_times, colnames(E))
})

test_that("partially observed periods retain the requested missing-data policy", {
  set.seed(26)
  E <- matrix(rnorm(8 * 20), 8)
  E[1, 1] <- NA_real_

  expect_error(cd_test(E, type = "CDw", seed = 7), "balanced sample")
  expect_equal(cd_test(E, type = "CD")$T, 20L)
})

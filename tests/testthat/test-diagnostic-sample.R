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

test_that("complete-time selection records all excluded periods in input order", {
  set.seed(27)
  E <- matrix(rnorm(8 * 20), 8L, dimnames = list(NULL, 2001:2020))
  E[1L, 1L] <- NA_real_
  E[, 2L] <- NA_real_
  E[2L, 3L] <- NA_real_
  expect_message(a <- cd_test(E, type = "CDw", seed = 7,
                              na.action = "drop.incomplete.times"), "Dropped 2")
  expect_identical(a$excluded_times, colnames(E)[1:3])
  expect_identical(a$kept_times, colnames(E)[-(1:3)])
  expect_equal(a$tests, cd_test(E[, -(1:3)], type = "CDw", seed = 7)$tests)
  colnames(E) <- NULL
  expect_message(b <- cd_test(E, na.action = "drop.incomplete.times"), "Dropped 2")
  expect_identical(b$excluded_times, 1:3)
  colnames(E) <- c(rep("same", 4L), as.character(5:20))
  expect_message(c <- cd_test(E, na.action = "drop.incomplete.times"), "Dropped 2")
  expect_identical(c$excluded_times, rep("same", 3L))
  expect_identical(c$kept_times[1L], "same")
})

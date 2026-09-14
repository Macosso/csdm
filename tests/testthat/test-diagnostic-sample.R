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

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

test_that("repeated CDw matches equation 33 using independent cross-products", {
  set.seed(731)
  E <- matrix(rnorm(1200), 12L) * seq(.5, 2, length.out = 12L)
  U <- E - rowMeans(E)
  reps <- 30L
  set.seed(84)
  draws <- replicate(reps, {
    w <- sample(c(-1, 1), nrow(E), replace = TRUE)
    products <- crossprod(t(U * w))
    sqrt(2 / (ncol(E) * nrow(E) * (nrow(E) - 1))) *
      sum(products[upper.tri(products)]) / mean(U^2)
  })
  expected <- sum(draws) / sqrt(reps)
  actual <- cd_test(E, type = "CDw", seed = 84, reps = reps)$tests$CDw
  expect_equal(actual$statistic, expected, tolerance = 1e-12)
  expect_equal(actual$p.value, 2 * pnorm(abs(expected), lower.tail = FALSE))
  expect_identical(actual$reps, reps)
  expect_equal(cd_test(E, type = "CDw", seed = 84, reps = 1L)$tests$CDw$statistic,
               draws[1L], tolerance = 1e-12)
  expect_identical(cd_test(E, type = "CDw", seed = 84)$tests,
                   cd_test(E, type = "CDw", seed = 84, reps = 1L)$tests)
})

test_that("repeated CDw-plus adds its screening statistic exactly once", {
  set.seed(732)
  E <- matrix(rnorm(12 * 100), 12L) + rep(rnorm(100), each = 12L)
  rho <- abs(cor(t(E))[upper.tri(cor(t(E)))])
  enhancement <- sum(rho[rho > 2 * sqrt(log(nrow(E)) / ncol(E))])
  expect_gt(enhancement, 0)
  a <- cd_test(E, type = "CDw+", seed = 18, reps = 30L)
  expected <- cd_test(E, type = "CDw", seed = 18, reps = 30L)$tests$CDw$statistic +
    enhancement
  expect_equal(a$tests$CDw_plus$statistic, expected)
  expect_equal(a$tests$CDw_plus$p.value, 2 * pnorm(abs(expected), lower.tail = FALSE))
  expect_identical(a$tests$CDw_plus$reps, 30L)
  expect_equal(a$tests$CDw_plus$enhancement, enhancement)
  all <- cd_test(E, type = "all", seed = 18, reps = 30L, n_pc = 1L)
  expect_equal(all$tests$CDw, a$tests$CDw)
  expect_equal(all$tests$CDw_plus, a$tests$CDw_plus)
  expect_equal(all$tests$CD, cd_test(E)$tests$CD)
  expect_equal(all$tests$CDstar,
               cd_test(E, type = "CDstar", n_pc = 1L)$tests$CDstar)
})

test_that("repetitions are validated and forwarded through all methods", {
  d <- panel_fixture()
  fit <- csdm(y ~ x, d, "id", "time")
  E <- residuals(fit)
  expected <- cd_test(E, type = "CDw+", reps = 7L, seed = 18)$tests
  expect_equal(cd_test(fit, type = "CDw+", reps = 7L, seed = 18)$tests, expected)
  d$e <- as.numeric(t(E))
  expect_equal(cd_test(d, e, id = "id", time = "time", type = "CDw+",
                       reps = 7L, seed = 18)$e$tests, expected)
  for (invalid in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "30", TRUE, NULL)) {
    expect_error(cd_test(E, type = "CDw", reps = invalid), "reps")
  }
  skip_if_not_installed("plm")
  p <- plm::pdata.frame(d, index = c("id", "time"))
  expect_equal(cd_test(p, e, type = "CDw+", reps = 7L, seed = 18)$e$tests, expected)
})

test_that("repeated weighted tests preserve seeded RNG and consume unseeded draws", {
  set.seed(733)
  E <- matrix(rnorm(600), 10L)
  before <- .Random.seed
  cd_test(E, type = "CDw+", reps = 7L, seed = 9)
  expect_identical(.Random.seed, before)
  expect_error(cd_test(E, type = "all", reps = 7L, seed = 9, n_pc = 100L), "n_pc")
  expect_identical(.Random.seed, before)
  cd_test(E, type = "CDw", reps = 7L)
  after <- .Random.seed
  .Random.seed <<- before
  replicate(7L, sample(c(-1, 1), nrow(E), replace = TRUE))
  expect_identical(.Random.seed, after)
})

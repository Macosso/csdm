test_that("data-frame variables match independently shaped matrices", {
  d <- panel_fixture()
  E <- t(matrix(d$x, nrow = 60L))
  dimnames(E) <- list(as.character(1:8), as.character(1:60))
  expected <- cd_test(E, type = "all", n_pc = 2L, seed = 19)
  shuffled <- d[sample(nrow(d)), ]
  actual <- cd_test(shuffled, x, "y", id = "id", time = "time",
                    type = "all", n_pc = 2L, seed = 19)
  expect_s3_class(actual, "cd_test_list")
  expect_named(actual, c("x", "y"))
  expect_s3_class(actual$x, "cd_test")
  expect_equal(actual$x$tests, expected$tests)
  expect_identical(actual$x$units, rownames(E))
  expect_identical(actual$x$kept_times, colnames(E))
  expect_identical(actual$x$variable, "x")
  expect_identical(attr(actual, "id"), "id")
  expect_identical(attr(actual, "time"), "time")
  expect_true(any(grepl("shuffled", deparse(attr(actual, "call")))))
})

test_that("variable selection supports bare names, strings, and caller vectors", {
  d <- panel_fixture()
  vars <- c("x", "y")
  expected <- cd_test(d, x, y, id = "id", time = "time")
  expect_equal(cd_test(d, "x", "y", id = "id", time = "time")$x$tests,
               expected$x$tests)
  expect_equal(cd_test(d, vars, id = "id", time = "time")$y$tests,
               expected$y$tests)
  expect_named(cd_test(d, c("x", "y"), id = "id", time = "time"), vars)
  # Data column names take precedence over objects in the caller.
  x <- "z"
  expect_named(cd_test(d, x, id = "id", time = "time"), "x")
  d$type <- d$x
  expect_named(cd_test(d, type, id = "id", time = "time"), "type")
})

test_that("data inputs retain variable-specific samples and structural gaps", {
  d <- panel_fixture()
  d$x[d$id == 1L] <- 1
  d$x[d$time == 1L] <- NA_real_
  d <- d[!(d$id == 2L & d$time == 2L), ]
  actual <- cd_test(d, x, y, id = "id", time = "time")
  expect_equal(actual$x$N, 7L)
  expect_equal(actual$x$T, 59L)
  expect_identical(actual$x$excluded_unit_ids, "1")
  expect_identical(actual$x$excluded_times, "1")
  expect_equal(actual$y$N, 8L)
  expect_equal(actual$y$T, 60L)
  E <- matrix(NA_real_, 8L, 60L)
  E[cbind(d$id, d$time)] <- d$y
  expect_equal(actual$y$tests, cd_test(E)$tests)
  expect_error(cd_test(d, y, id = "id", time = "time", type = "CDw"), "balanced")
  expect_message(balanced <- cd_test(d, y, id = "id", time = "time", type = "CDw",
                                    seed = 7, na.action = "drop.incomplete.times"),
                 "Dropped 1")
  expect_equal(balanced$y$T, 59L)
})

test_that("time indexes can be labels or dates without an estimation grid", {
  d <- panel_fixture(periods = 6L)
  expected <- cd_test(d, x, id = "id", time = "time")$x$tests
  d$time <- paste0("period-", d$time)
  expect_equal(cd_test(d, x, id = "id", time = "time")$x$tests, expected)
  d$time <- as.Date("2020-01-01") + rep(c(0, 2, 7, 10, 21, 40), 8L)
  expect_equal(cd_test(d, x, id = "id", time = "time")$x$tests, expected)
})

test_that("pdata.frame methods use stored indexes including removed columns", {
  skip_if_not_installed("plm")
  d <- panel_fixture()
  expected <- cd_test(d, x, y, id = "id", time = "time", type = "CDw", seed = 8)
  for (drop_index in c(FALSE, TRUE)) {
    p <- plm::pdata.frame(d, index = c("id", "time"), drop.index = drop_index)
    actual <- cd_test(p, x, "y", type = "CDw", seed = 8)
    expect_equal(actual$x$tests, expected$x$tests)
    expect_equal(actual$y$tests, expected$y$tests)
    expect_identical(actual$x$kept_times, expected$x$kept_times)
    expect_identical(attr(actual, "id"), "id")
    expect_equal(cd_test(p, x, id = "ignored", time = "ignored")$x$tests,
                 cd_test(d, x, id = "id", time = "time")$x$tests)
  }
})

test_that("invalid selections and panel keys fail explicitly", {
  d <- panel_fixture()
  expect_error(cd_test(d, id = "id", time = "time"), "explicitly")
  expect_error(cd_test(d, x), "distinct columns")
  expect_error(cd_test(d, x, id = "id", time = "id"), "distinct columns")
  expect_error(cd_test(d, "unknown", id = "id", time = "time"), "Unknown")
  expect_error(cd_test(d, x, x, id = "id", time = "time"), "only once")
  expect_error(cd_test(d, x, typo = "CD", id = "id", time = "time"), "unnamed")
  expect_error(cd_test(d, log(x), id = "id", time = "time"), "column names")
  expect_error(cd_test(d, id, id = "id", time = "time"), "indexes")
  d$label <- "a"
  expect_error(cd_test(d, label, id = "id", time = "time"), "numeric vectors")
  expect_error(cd_test(rbind(d, d[1L, ]), x, id = "id", time = "time"), "Duplicate")
  d$id[1L] <- NA
  expect_error(cd_test(d, x, id = "id", time = "time"), "keys")
  d <- panel_fixture()
  d$time[1L] <- Inf
  expect_error(cd_test(d, x, id = "id", time = "time"), "keys")
  d <- panel_fixture()
  d$x[1L] <- Inf
  expect_error(cd_test(d, x, id = "id", time = "time"), "infinite")
  expect_error(cd_test(d[FALSE, ], x, id = "id", time = "time"), "nonempty")
})

test_that("combined printing shows only the requested tests and each sample", {
  d <- panel_fixture()
  d$x[d$id == 1L] <- 0
  one <- cd_test(d, x, id = "id", time = "time")
  expect_named(one, "x")
  output <- capture.output(returned <- print(one))
  expect_identical(returned, one)
  expect_true(any(grepl("variable", output)))
  expect_true(any(grepl("7", output)))
  enhanced <- cd_test(d, x, y, id = "id", time = "time", type = "CDw+", seed = 4)
  output <- capture.output(print(enhanced))
  expect_equal(sum(grepl("CDw+", output, fixed = TRUE)), 2L)
  expect_false(any(grepl("CDstar", output, fixed = TRUE)))
  all <- cd_test(d, x, y, id = "id", time = "time", type = "all", seed = 4)
  expect_equal(sum(grepl("CDstar", capture.output(print(all)), fixed = TRUE)), 2L)
})

test_that("seeded data-frame diagnostics preserve caller RNG and column order", {
  d <- panel_fixture()
  before <- .Random.seed
  a <- cd_test(d, x, y, id = "id", time = "time", type = "CDw", seed = 19)
  expect_identical(.Random.seed, before)
  b <- cd_test(d, y, x, id = "id", time = "time", type = "CDw", seed = 19)
  expect_equal(a$x$tests, b$x$tests)
  expect_equal(a$y$tests, b$y$tests)
})

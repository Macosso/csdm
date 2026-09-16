test_that("lags use time keys instead of row positions", {
  d <- panel_fixture()
  g <- d[!(d$id == 1 & d$time == 20), ]
  a <- suppressMessages(csdm(y ~ x, g, "id", "time", model = "dcce",
    csa = csdm_csa("_none"), lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 1)))
  expect_true(is.na(residuals(a)["1", "21"]))
  d[d$id == 1 & d$time == 20, c("y", "x")] <- NA
  b <- suppressMessages(csdm(y ~ x, d, "id", "time", model = "dcce",
    csa = csdm_csa("_none"), lr = csdm_lr(type = "ardl", ylags = 1, xdlags = 1)))
  expect_equal(coef(a), coef(b))
  expect_equal(residuals(a), residuals(b))
})

test_that("CSA lags preserve gaps shared by all units", {
  d <- panel_fixture()
  d <- d[d$time != 20, ]
  a <- suppressMessages(csdm(y ~ x, d, "id", "time", model = "dcce", csa = csdm_csa(lags = 1)))
  expect_true(all(is.na(residuals(a)[, "21"])))
})

test_that("non-unit grids require an explicit time step", {
  d <- panel_fixture(); d$time <- d$time / 4
  expect_error(csdm(y ~ x, d, "id", "time"), "time_step")
  expect_no_error(suppressMessages(csdm(y ~ x, d, "id", "time", time_step = 0.25)))
  expect_equal(csdm:::.csdm_time_index(c(1, 3, 4)), c(1, 3, 4))
})

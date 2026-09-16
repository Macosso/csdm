test_that("panel keys are validated before fitting", {
  d <- panel_fixture()
  expect_error(csdm(y ~ x, rbind(d, d[1, ]), "id", "time"), "Duplicate")
  d$id[1] <- NA
  expect_error(csdm(y ~ x, d, "id", "time"), "keys")
  d <- panel_fixture(); d$time[1] <- Inf
  expect_error(csdm(y ~ x, d, "id", "time"), "keys")
  d <- panel_fixture(); d$.csdm_rowid__ <- 1
  expect_error(csdm(y ~ x, d, "id", "time"), "reserved")
})

test_that("sorting retains original observation identity", {
  d <- panel_fixture()
  d <- d[sample(nrow(d)), ]
  p <- csdm:::.csdm_prepare_panel_df(d, "id", "time")
  expect_equal(p$y, d$y[p$.csdm_rowid__])
  a <- suppressMessages(csdm(y ~ x, d, "id", "time"))
  b <- suppressMessages(csdm(y ~ x, d[order(d$id, d$time), ], "id", "time"))
  expect_equal(coef(a), coef(b))
  expect_equal(residuals(a), residuals(b))
})

test_that("numeric pdata.frame indexes are normalized safely", {
  skip_if_not_installed("plm")
  d <- panel_fixture()
  p <- plm::pdata.frame(d, index = c("id", "time"))
  a <- suppressMessages(csdm(y ~ x, p))
  b <- suppressMessages(csdm(y ~ x, d, "id", "time"))
  expect_equal(coef(a), coef(b))
})

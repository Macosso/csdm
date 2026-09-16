test_that("CSA-spanned regressors are not reported as identified", {
  d <- panel_fixture(); d$x <- sin(d$time)
  expect_error(suppressWarnings(csdm(y ~ x + z, d, "id", "time", model = "cce")), "No identified")
})

test_that("redundant nuisance columns preserve identified slopes", {
  d <- panel_fixture(); d$xx <- d$x
  a <- csdm(y ~ x + z, d, "id", "time", model = "cce", csa = csdm_csa(c("y", "x", "z")))
  b <- csdm(y ~ x + z, d, "id", "time", model = "cce", csa = csdm_csa(c("y", "x", "z", "xx")))
  expect_equal(coef(a), coef(b), tolerance = 1e-9)
  expect_equal(residuals(a), residuals(b), tolerance = 1e-9)
})

test_that("unidentified units are explicitly excluded", {
  d <- panel_fixture(); d$z[d$id == 1] <- 0
  expect_warning(a <- suppressMessages(csdm(y ~ x + z, d, "id", "time")), "Excluded")
  expect_identical(a$meta$dropped_units, "1")
  expect_true("z" %in% a$units[["1"]]$aliased)
  expect_true(all(is.na(residuals(a)["1", ])))
})

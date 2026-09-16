test_that("tidy inference and augmentation use public contracts", {
  d <- panel_fixture(); d$x[7] <- NA
  a <- csdm(y ~ x + z, d, "id", "time")
  td <- generics::tidy(a, conf.int = TRUE)
  expect_s3_class(td, "tbl_df")
  expect_equal(td$estimate, unname(coef(a)))
  expect_equal(td$std.error, unname(sqrt(diag(vcov(a)))))
  expect_equal(td$conf.low, unname(confint(a)[, 1]))
  expect_equal(generics::glance(a)$nobs, nobs(a))
  aug <- generics::augment(a)
  expect_equal(nrow(aug), nrow(d))
  expect_false(aug$.used[7])
  expect_equal(aug$.fitted + aug$.resid, replace(d$y, 7, NA))
})

test_that("lmtest and broom use the same inference", {
  skip_if_not_installed("lmtest"); skip_if_not_installed("broom")
  d <- panel_fixture()
  a <- csdm(y ~ x, d, "id", "time")
  td <- broom::tidy(a)
  ct <- lmtest::coeftest(a)
  expect_equal(td$p.value, unname(ct[, 4]))
  expect_equal(td$std.error, unname(ct[, 2]))
})

test_that("modelsummary extracts coefficients", {
  skip_if_not_installed("modelsummary")
  d <- panel_fixture()
  a <- csdm(y ~ x, d, "id", "time")
  out <- modelsummary::modelsummary(list(MG = a), output = "data.frame")
  expect_s3_class(out, "data.frame")
  expect_true("x" %in% out$term)
})

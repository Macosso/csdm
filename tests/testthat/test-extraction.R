test_that("standard generics expose the fitted sample and original row mapping", {
  d <- panel_fixture(); d <- d[sample(nrow(d)), ]; d$x[7] <- NA
  a <- csdm(y ~ x + z, d, "id", "time", subset = time > 5, na.action = na.exclude)
  expect_s3_class(model.frame(a), "data.frame")
  expect_equal(nrow(model.frame(a)), nobs(a))
  expect_equal(nrow(model.matrix(a)), nobs(a))
  expect_s3_class(terms(a), "terms")
  expect_equal(df.residual(a), Inf)
  f <- fitted(a, format = "vector"); e <- residuals(a, format = "vector")
  expect_length(f, nrow(d))
  ok <- is.finite(f)
  expect_equal(unname(f[ok] + e[ok]), d$y[ok])
  expect_true(all(is.na(f[d$time <= 5])))
  expect_equal(fitted(a, format = "long")$.row, seq_len(nrow(d)))
  expect_equal(fitted(a), predict(a))
})

test_that("update works with stored data", {
  a <- local({
    local_data <- panel_fixture()
    csdm(y ~ x, local_data, "id", "time")
  })
  b <- update(a, . ~ . + z)
  reference <- csdm(y ~ x + z, a$data, "id", "time")
  expect_equal(coef(b), coef(reference))
  expect_true(is.call(update(a, evaluate = FALSE)))
})

test_that("formula exposes constructed economic lags", {
  d <- panel_fixture()
  a <- csdm(y ~ x, d, "id", "time", model = "dcce", csa = csdm_csa("_none"),
    lr = csdm_lr(type = "ardl", ylags = 1))
  expect_true("lag1_y" %in% all.vars(formula(a)))
})

test_that("update retains resolved local options and pdata time indexes", {
  skip_if_not_installed("plm")
  a <- local({
    d <- plm::pdata.frame(panel_fixture(), index = c("id", "time"))
    cutoff <- 5
    spacing <- 1
    action <- na.exclude
    estimator <- "mg"
    csdm(y ~ x, d, model = estimator, subset = as.numeric(as.character(time)) > cutoff,
      na.action = action, time_step = spacing)
  })
  b <- update(a)
  expect_equal(coef(b), coef(a))
  expect_equal(nobs(b), nobs(a))
  expect_equal(b$meta$selected_rows, a$meta$selected_rows)
})

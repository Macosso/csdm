test_that("cluster covariance matches sandwich on supported designs", {
  skip_if_not_installed("sandwich")
  d <- panel_fixture()
  m <- lm(y ~ x + z, d)
  X <- model.matrix(m); u <- residuals(m)
  expect_equal(suppressWarnings(cluster_vcov(X, u, d$id)), sandwich::vcovCL(m, cluster = d$id, type = "HC1"), tolerance = 1e-9)
  expect_equal(suppressWarnings(cluster_vcov(X, u, d[c("id", "time")], type = "twoway")),
    sandwich::vcovCL(m, cluster = d[c("id", "time")], type = "HC1"), tolerance = 1e-9)
  expect_error(suppressWarnings(cluster_vcov(X, u, rep(1, nrow(d)))), "two clusters")
  expect_error(suppressWarnings(cluster_vcov(X, u, c(NA, d$id[-1]))), "nonmissing")
  expect_error(suppressWarnings(cluster_vcov(X, u, d$id[-1])), "aligned")
})

test_that("orphaned covariance utilities are deprecated", {
  d <- panel_fixture()
  m <- lm(y ~ x + z, d)
  X <- model.matrix(m); u <- residuals(m)

  expect_warning(sandwich_vcov(X, u), class = "deprecatedWarning")
  expect_warning(cluster_vcov(X, u, d$id), class = "deprecatedWarning")
  expect_warning(pooled_vcov(cbind(a = 1:4, b = 4:1)), class = "deprecatedWarning")
})

test_that("residual type requests cannot silently return another type", {
  d <- panel_fixture(); m <- csdm(y ~ x, d, "id", "time")
  expect_error(suppressWarnings(get_residuals(m, "pca")), "PCA")
  expect_null(suppressWarnings(get_residuals(m, "pca", strict = FALSE)))
  expect_warning(get_residuals(m), class = "deprecatedWarning")
  E <- residuals(m)
  expect_equal(suppressWarnings(prepare_cd_input(E, standardize = "none"))$Z, E)
  expect_warning(prepare_cd_input(E), class = "deprecatedWarning")
})

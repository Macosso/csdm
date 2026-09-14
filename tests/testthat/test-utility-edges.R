test_that("cluster covariance matches sandwich on supported designs", {
  skip_if_not_installed("sandwich")
  d <- panel_fixture()
  m <- lm(y ~ x + z, d)
  X <- model.matrix(m); u <- residuals(m)
  expect_equal(cluster_vcov(X, u, d$id), sandwich::vcovCL(m, cluster = d$id, type = "HC1"), tolerance = 1e-9)
  expect_equal(cluster_vcov(X, u, d[c("id", "time")], type = "twoway"),
    sandwich::vcovCL(m, cluster = d[c("id", "time")], type = "HC1"), tolerance = 1e-9)
  expect_error(cluster_vcov(X, u, rep(1, nrow(d))), "two clusters")
  expect_error(cluster_vcov(X, u, c(NA, d$id[-1])), "nonmissing")
  expect_error(cluster_vcov(X, u, d$id[-1]), "aligned")
})

test_that("residual type requests cannot silently return another type", {
  d <- panel_fixture(); m <- csdm(y ~ x, d, "id", "time")
  expect_error(get_residuals(m, "pca"), "PCA")
  expect_null(get_residuals(m, "pca", strict = FALSE))
  E <- residuals(m)
  expect_equal(prepare_cd_input(E, standardize = "none")$Z, E)
})

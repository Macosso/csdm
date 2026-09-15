test_that("CSA lag names survive validation and support partial specification", {
  expect_identical(csdm_csa(lags = c(y = 2, x = 1))$lags, c(y = 2L, x = 1L))
  expect_identical(csdm:::.csdm_csa_lags(c(y = 2L), c("y", "x")), c(y = 2L, x = 0L))
  expect_error(csdm_csa(lags = c(y = 1, y = 2)), "unique")
  for (v in list(0.5, -1, Inf, NA_real_, numeric())) expect_error(csdm_csa(lags = v), "integers")
  expect_error(csdm_lr(ylags = 1.5), "integers")
  d <- panel_fixture()
  expect_no_error(suppressMessages(csdm(y ~ x, d, "id", "time", model = "dcce",
    csa = csdm_csa(lags = c(y = 2, x = 1)))))
  expect_error(csdm(y ~ x, d, "id", "time", model = "dcce", csa = csdm_csa(lags = c(typo = 1))), "Unknown CSA")
})

test_that("unsupported settings cannot silently succeed", {
  d <- panel_fixture()
  expect_error(csdm(y ~ x, d, "id", "time", vcov = csdm_vcov("nw")), "Only vcov")
  expect_error(csdm(y ~ x, d, "id", "time", weights = rep(1, nrow(d))), "Unused")
  expect_error(csdm(y ~ x, d, "id", "time", pooled = suppressWarnings(csdm_pooled("x"))), "Pooled")
  expect_error(csdm(y ~ x, d, "id", "time", lr = csdm_lr(type = "ardl", ylags = 1)), "Use model")
  expect_error(csdm(y ~ x, d, "id", "time", csa = list()), "construct")
})

test_that("the unused pooled specification is deprecated without warning on default fits", {
  d <- panel_fixture()
  expect_warning(csdm_pooled(), class = "deprecatedWarning")
  expect_no_warning(csdm(y ~ x, d, "id", "time"))
})

test_that("inactive CSA requests cannot silently succeed", {
  d <- panel_fixture()
  expect_error(csdm(y ~ x, d, "id", "time", csa = csdm_csa("x")), "MG does not")
  expect_error(csdm(y ~ x, d, "id", "time", csa = csdm_csa(lags = 1)), "MG does not")
  expect_error(csdm(y ~ x, d, "id", "time", model = "dcce",
    csa = csdm_csa("_none", lags = 1)), "CSA lags require")
})

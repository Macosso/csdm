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
  expect_error(csdm_vcov("nw"), "Only type='mg'")
  expect_error(csdm_vcov(adjust = TRUE), "options are not implemented")
  expect_error(csdm_csa(scope = "global"), "scope='estimation'")
  expect_error(csdm_csa(cluster = "id"), "cluster.*not implemented")
  expect_error(csdm_lr(type = "ecm"), "one of")
  expect_error(csdm_lr(vars = "x"), "vars.*not implemented")
  expect_error(csdm_lr(options = list(foo = TRUE)), "options.*not implemented")
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

test_that("fullsample averages each variable over its available observations", {
  d <- data.frame(
    id = rep(1:3, each = 2),
    time = rep(1:2, 3),
    y = c(1, 2, 30, 4, 5, 6),
    x = c(1, 2, NA, 4, 5, 6),
    .csdm_rowid__ = 1:6
  )
  attr(d, "csdm_time_step") <- 1

  estimation <- csdm:::.csdm_augment(
    d, y ~ x, y ~ x, "id", "time", csdm_csa(), fullsample = FALSE
  )
  full <- csdm:::.csdm_augment(
    d, y ~ x, y ~ x, "id", "time", csdm_csa(), fullsample = TRUE
  )

  expect_equal(estimation$data$.csdm_csa1_lag0[d$time == 1], rep(3, 3))
  expect_equal(full$data$.csdm_csa1_lag0[d$time == 1], rep(12, 3))
  expect_equal(estimation$data$.csdm_csa2_lag0[d$time == 1], rep(3, 3))
  expect_equal(full$data$.csdm_csa2_lag0[d$time == 1], rep(3, 3))
  expect_identical(estimation$csa$source_n, 5L)
  expect_identical(full$csa$source_n_by_variable, c(y = 6L, x = 5L))
  expect_true(full$csa$fullsample)
})

test_that("fullsample is available only when averages are active", {
  d <- panel_fixture()
  expect_error(
    csdm(y ~ x, d, "id", "time", model = "mg", fullsample = TRUE),
    "models with cross-sectional averages"
  )
  expect_error(
    csdm(
      y ~ x, d, "id", "time", model = "dcce",
      csa = csdm_csa("_none"), fullsample = TRUE
    ),
    "requires active cross-sectional averages"
  )
  expect_no_error(
    suppressMessages(csdm(
      y ~ x, d, "id", "time", model = "cce", fullsample = TRUE
    ))
  )
})

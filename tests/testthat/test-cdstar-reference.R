test_that("CDstar preserves unit-specific scales", {
  set.seed(482)
  E <- outer(seq(.1, 3, length.out = 20), rnorm(100)) + matrix(rnorm(2000), 20)
  A <- scale(t(E)); N <- ncol(A); TT <- nrow(A)
  F <- cbind(1, eigen(tcrossprod(A), symmetric = TRUE)$vectors[, 1, drop = FALSE])
  B <- solve(crossprod(F), crossprod(F, A)); U <- A - F %*% B
  G <- B[-1, , drop = FALSE] / sqrt(mean(B[-1, ]^2))
  sigma <- sqrt(colMeans(U^2))
  phi <- mean(G / rep(sigma, each = nrow(G)))
  a <- 1 - as.numeric(G) * sigma * phi
  correction <- mean(a^2)
  cd <- sqrt(2 * TT / (N * (N - 1))) * sum(cor(U)[upper.tri(cor(U))])
  expected <- (cd + sqrt(TT / 2) * (1 - correction)) / correction
  expect_equal(cd_test(E, type = "CDstar", n_pc = 1)$tests$CDstar$statistic, expected, tolerance = 1e-9)
  expect_error(cd_test(E, type = "CDstar", n_pc = 20), "n_pc")
  expect_error(cd_test(E, type = "CDstar", n_pc = 1.5), "integers")
})

test_that("standardized CDstar agrees with the loading-normalized eigenvector reference", {
  set.seed(485)
  N <- 18L; TT <- 90L; r <- 2L
  E <- matrix(rnorm(N * TT), N) * seq(.3, 2.5, length.out = N) +
    outer(seq(-1, 2, length.out = N), rnorm(TT)) +
    outer(rep(c(-.7, .8, 1.4), length.out = N), rnorm(TT))
  A <- scale(t(E))
  # Pesaran-Xie loading normalization, applied to the chosen standardized input.
  Q <- eigen(crossprod(A), symmetric = TRUE)$vectors[, seq_len(r), drop = FALSE]
  Gamma <- sqrt(N) * Q
  factors <- A %*% Gamma / N
  U <- A - factors %*% t(Gamma)
  sigma <- sqrt(colMeans(U^2))
  phi <- colMeans(Gamma / sigma)
  a <- 1 - sigma * as.numeric(Gamma %*% phi)
  q <- mean(a^2)
  epsilon <- sweep(U, 2L, sigma, "/")
  rho <- crossprod(epsilon) / TT
  cd <- sqrt(2 * TT / (N * (N - 1))) * sum(rho[upper.tri(rho)])
  expected <- (cd + sqrt(TT / 2) * (1 - q)) / q
  actual <- cd_test(E, type = "CDstar", n_pc = r)$tests$CDstar
  expect_equal(crossprod(Gamma) / N, diag(r), ignore_attr = TRUE)
  expect_equal(actual$statistic, expected, tolerance = 1e-10)
  expect_equal(actual$p.value, 2 * pnorm(abs(expected), lower.tail = FALSE))
  expect_equal(cd_test(E * seq(.5, 3, length.out = N), type = "CDstar", n_pc = r)$tests,
               cd_test(E, type = "CDstar", n_pc = r)$tests, tolerance = 1e-10)
  expect_equal(cd_test(E, type = "CDstar", n_pc = 0L)$tests$CDstar$statistic,
               cd_test(E)$tests$CD$statistic)
})

test_that("fitted CDstar retains the factor component while other tests use full residuals", {
  d <- panel_fixture(n = 16L, periods = 90L)
  fit <- csdm(y ~ x, d, "id", "time", model = "cce")
  # Construct the paper's regression input directly from the observations.
  V <- residuals(fit)
  V[] <- NA_real_
  for (uid in fit$meta$included_units) {
    beta <- fit$coef_i[uid, ]
    rows <- which(d$id == uid)
    V[uid, as.character(d$time[rows])] <- d$y[rows] - beta["(Intercept)"] - beta["x"] * d$x[rows]
  }
  expected <- cd_test(V, type = "CDstar", n_pc = 1L)
  actual <- cd_test(fit, type = "CDstar", n_pc = 1L)
  expect_equal(actual$tests$CDstar$statistic, expected$tests$CDstar$statistic)
  expect_equal(actual$tests$CDstar$p.value, expected$tests$CDstar$p.value)
  expect_identical(actual$tests$CDstar$input, "partial_residuals")
  expect_false(isTRUE(all.equal(V, residuals(fit))))
  all <- cd_test(fit, type = "all", n_pc = 1L, seed = 19, reps = 7L)
  full <- cd_test(residuals(fit), type = "all", n_pc = 1L, seed = 19, reps = 7L)
  expect_equal(all$tests[c("CD", "CDw", "CDw_plus")],
               full$tests[c("CD", "CDw", "CDw_plus")])
  expect_equal(all$tests$CDstar, actual$tests$CDstar)
  expect_false(isTRUE(all.equal(actual$tests$CDstar$statistic, full$tests$CDstar$statistic)))
})

test_that("partial CDstar input respects dynamic terms and original sample rows", {
  set.seed(486)
  d <- panel_fixture(n = 12L, periods = 90L)
  d <- d[sample(nrow(d)), ]
  d$y[d$id == 1L & d$time == 25] <- NA_real_
  fit <- csdm(y ~ x, d, "id", "time", model = "cs_ardl",
               lr = csdm_lr(type = "ardl", ylags = 1L, xdlags = 1L),
               subset = time > 2L)
  V <- csdm:::.cdstar_fit_input(fit, residuals(fit))
  used <- fit$sample$row[fit$sample$used]
  keys <- d[used, c("id", "time")]
  # Add back the CSA fitted component as an independent reconstruction.
  W <- residuals(fit)
  nuisance <- setdiff(colnames(fit$augmented_matrix), colnames(fit$model_matrix))
  for (uid in fit$meta$included_units) {
    rows <- which(keys$id == uid)
    unit_fit <- lm.fit(fit$augmented_matrix[rows, , drop = FALSE],
                        model.response(model.frame(fit))[rows])
    beta <- unit_fit$coefficients[nuisance]
    beta[is.na(beta)] <- 0
    W[uid, as.character(keys$time[rows])] <-
      W[uid, as.character(keys$time[rows])] +
      as.numeric(fit$augmented_matrix[rows, nuisance, drop = FALSE] %*% beta)
  }
  expect_equal(V, W, tolerance = 1e-10)
  expect_identical(is.na(V), is.na(residuals(fit)))
  expected <- suppressMessages(cd_test(V, type = "CDstar", n_pc = 1L,
                         na.action = "drop.incomplete.times"))
  actual <- suppressMessages(cd_test(fit, type = "CDstar", n_pc = 1L,
                         na.action = "drop.incomplete.times"))
  expect_equal(actual$tests$CDstar$statistic, expected$tests$CDstar$statistic)
  expect_equal(actual$excluded_times, expected$excluded_times)
  mg <- csdm(y ~ x, d, "id", "time", model = "mg")
  expect_equal(csdm:::.cdstar_fit_input(mg, residuals(mg)), residuals(mg))
  trend <- csdm(y ~ x, d, "id", "time", model = "cce", trend = "unit")
  V_trend <- csdm:::.cdstar_fit_input(trend, residuals(trend))
  for (uid in trend$meta$included_units) {
    rows <- which(d$id == uid & is.finite(d$y))
    beta <- trend$coef_i[uid, ]
    expected <- d$y[rows] - beta["(Intercept)"] - beta["x"] * d$x[rows] -
      beta[".csdm_trend__"] * d$time[rows]
    expect_equal(as.numeric(V_trend[uid, as.character(d$time[rows])]), expected)
  }
  old <- fit
  old$model_matrix <- NULL
  expect_error(cd_test(old, type = "CDstar", n_pc = 1L), "refit")
})

test_that("fitted CDstar maps pdata.frame factor time labels to the fitted numeric grid", {
  skip_if_not_installed("plm")
  d <- panel_fixture(n = 12L, periods = 90L)
  numeric_fit <- csdm(y ~ x, d, "id", "time", model = "cce")
  d$time <- sprintf("%03d", d$time)
  panel <- plm::pdata.frame(d, index = c("id", "time"))
  fit <- csdm(y ~ x, panel, "id", "time", model = "cce")
  expect_equal(cd_test(fit, type = "CDstar", n_pc = 1L)$tests,
               cd_test(numeric_fit, type = "CDstar", n_pc = 1L)$tests)
})

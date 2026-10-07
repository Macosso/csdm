test_that("unbalanced CD uses overlap means and Pesaran's unit normalization", {
  E <- rbind(c(1, 4, 2, 7, NA, 5, 8, 3),
             c(NA, 3, 6, 2, 9, 4, NA, 7),
             c(2, NA, 1, 8, 3, NA, 6, 4),
             c(3, 5, NA, 1, 8, 7, 2, NA))
  pairs <- combn(nrow(E), 2L)
  terms <- apply(pairs, 2L, function(pair) {
    observed <- is.finite(E[pair[1L], ]) & is.finite(E[pair[2L], ])
    x <- E[pair[1L], observed]; y <- E[pair[2L], observed]
    x <- x - mean(x); y <- y - mean(y)
    sqrt(sum(observed)) * sum(x * y) / sqrt(sum(x^2) * sum(y^2))
  })
  expected <- sqrt(2 / (nrow(E) * (nrow(E) - 1))) * sum(terms)
  actual <- cd_test(E)$tests$CD
  expect_equal(actual$statistic, expected)
  expect_equal(actual$p.value, 2 * pnorm(abs(expected), lower.tail = FALSE))
  expect_equal(actual$pairs_used, ncol(pairs))
  balanced <- E
  balanced[is.na(balanced)] <- seq_len(sum(is.na(balanced)))
  expected <- sqrt(2 * ncol(balanced) / (nrow(E) * (nrow(E) - 1))) *
    sum(cor(t(balanced))[upper.tri(cor(t(balanced)))])
  expect_equal(cd_test(balanced)$tests$CD$statistic, expected)
})

test_that("CDw uses equation 30 pooled variance rather than pair-specific correlations", {
  E <- rbind(c(1, 3, 2, 5, 4, 8), c(4, 1, 8, 3, 12, 6),
             c(20, 8, 4, 16, 28, 12), c(6, 12, 3, 9, 18, 15))
  N <- nrow(E); TT <- ncol(E)
  U <- E - rowMeans(E)
  set.seed(42)
  w <- sample(c(-1, 1), N, replace = TRUE)
  pooled <- sum(w^2 * U^2) / (N * TT)
  numerator <- 0
  for (t in seq_len(TT)) for (i in 2:N) for (j in seq_len(i - 1L)) {
    numerator <- numerator + w[i] * U[i, t] * w[j] * U[j, t]
  }
  expected <- pooled^(-1) * sqrt(2 / (TT * N * (N - 1))) * numerator
  actual <- cd_test(E, type = "CDw", seed = 42)$tests$CDw$statistic
  expect_equal(actual, expected)
  correlations <- cor(t(E)) * outer(w, w)
  correlation_variant <- sqrt(2 * TT / (N * (N - 1))) *
    sum(correlations[upper.tri(correlations)])
  expect_gt(abs(actual - correlation_variant), .1)
})

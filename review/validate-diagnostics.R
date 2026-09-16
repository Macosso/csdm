.libPaths(c(normalizePath(".dev-library"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
set.seed(9142026)
reps <- 500L
results <- list()
for (N in c(20L, 50L)) for (TT in c(50L, 100L)) {
  for (alternative in c(FALSE, TRUE)) {
    p <- replicate(reps, {
      E <- matrix(rnorm(N * TT), N)
      E <- E * seq(.5, 1.5, length.out = N)
      if (alternative) E <- E + rep(rnorm(TT), each = N)
      a <- cd_test(E, type = "CDw+")
      c(CDw = a$tests$CDw$p.value, CDw_plus = a$tests$CDw_plus$p.value)
    })
    for (method in rownames(p)) {
      rejects <- sum(p[method, ] < .05)
      interval <- binom.test(rejects, reps)$conf.int
      results[[length(results) + 1L]] <- data.frame(N, T = TT, alternative, method,
        rejection = rejects / reps, lower = interval[1], upper = interval[2])
    }
  }
}
out <- do.call(rbind, results)
print(out, row.names = FALSE)
write.csv(out, "review/diagnostic-calibration.csv", row.names = FALSE)
# A broad predeclared smoke bound: reject severe null miscalibration.
stopifnot(all(out$rejection[!out$alternative] < .12))

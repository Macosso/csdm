# Reproducible size/power smoke checks for one and thirty CDw weight draws.
dir.create(".dev-check/config", recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CONFIG_DIR = normalizePath(".dev-check/config"))
.libPaths(c(normalizePath(".dev-library"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
set.seed(10072026)
simulations <- 500L
designs <- data.frame(N = c(50L, 100L), T = c(100L, 50L))
results <- list()
for (design in seq_len(nrow(designs))) {
  N <- designs$N[design]
  TT <- designs$T[design]
  for (alternative in c(FALSE, TRUE)) for (weight_reps in c(1L, 30L)) {
    p <- replicate(simulations, {
      E <- matrix(rnorm(N * TT), N) * seq(.5, 1.5, length.out = N)
      if (alternative) E <- E + rep(rnorm(TT), each = N)
      tests <- cd_test(E, type = "CDw+", reps = weight_reps)$tests
      c(CDw = tests$CDw$p.value, CDw_plus = tests$CDw_plus$p.value)
    })
    for (method in rownames(p)) {
      rejects <- sum(p[method, ] < .05)
      interval <- binom.test(rejects, simulations)$conf.int
      results[[length(results) + 1L]] <- data.frame(
        N, T = TT, alternative, reps = weight_reps, method, simulations,
        rejection = rejects / simulations, lower = interval[1L], upper = interval[2L]
      )
    }
  }
}
out <- do.call(rbind, results)
print(out, row.names = FALSE)
write.csv(out, "review/cd-reps-calibration.csv", row.names = FALSE)
# Broad predeclared bounds detect severe failures in these selected designs.
stopifnot(all(out$rejection[!out$alternative] < .12))
stopifnot(all(out$rejection[out$alternative & out$method == "CDw_plus"] > .8))

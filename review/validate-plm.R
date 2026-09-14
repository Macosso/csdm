.libPaths(c(normalizePath(".dev-library"), .libPaths()))
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-panel.R")
library(plm)
d <- panel_fixture()
a <- csdm(y ~ x + z, d, "id", "time", model = "cce")
b <- plm::pcce(y ~ x + z, d, index = c("id", "time"), model = "mg")
terms <- c("x", "z")
errors <- c(coef = max(abs(coef(a)[terms] - coef(b)[terms])),
  vcov = max(abs(vcov(a)[terms, terms] - vcov(b)[terms, terms])))
print(errors)
stopifnot(all(errors < 1e-8))
writeLines(c(paste(names(errors), format(errors, digits = 16)),
  paste("plm", packageVersion("plm")), capture.output(sessionInfo())),
  "review/plm-reference.txt")

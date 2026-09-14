.libPaths(c(normalizePath('.dev-library'), .libPaths()))
args <- commandArgs(TRUE)
pkgload::load_all('.', quiet = TRUE)
testthat::test_dir('tests/testthat', filter = if (length(args)) args[1] else NULL,
                   reporter = 'summary', stop_on_failure = TRUE)

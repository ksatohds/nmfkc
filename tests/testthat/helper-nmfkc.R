# CRAN caps the whole check at 10 minutes.  This package's regression suite is
# built around fitting the model many times -- bootstrap replicates, CV folds,
# multistart restarts -- so it cannot be made cheap by shrinking data alone:
# 0.9.4 spent 34 minutes here, and 0.9.5, with the worst blocks already reduced
# and skipped, still spent 442 seconds.
#
# Following CRAN's own suggestion (Uwe Ligges, 2026-08-23: "running less
# important tests only conditionally if some environment variable is set that
# you only define on your machine"), the full suite now runs only when
# NMFKC_FULL_TESTS is set.  Everywhere else -- CRAN, and any user who simply
# installs from source -- only tests/testthat/test-cran-smoke.R runs: it
# exercises every exported fitter and its S3 methods on toy data in a few
# seconds, which is what a check on someone else's machine is for.
#
# To run everything (this is what CI and the maintainer do):
#
#   Sys.setenv(NMFKC_FULL_TESTS = "true"); devtools::test()
#   NMFKC_FULL_TESTS=true R CMD check nmfkc_*.tar.gz
#
skip_unless_full <- function() {
  if (!nzchar(Sys.getenv("NMFKC_FULL_TESTS"))) {
    testthat::skip("full suite: set NMFKC_FULL_TESTS=true to run")
  }
  invisible(TRUE)
}

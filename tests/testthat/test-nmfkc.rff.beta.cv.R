## nmfkc.rff.beta.cv(): bandwidth / feature-count selection for RFF and
## positive RFF covariates by column-wise CV, mirroring nmfkc.kernel.beta.cv().

make_case <- function(seed = 4, p = 3, N = 120, P = 3) {
  set.seed(seed)
  centers <- matrix(stats::rnorm(p * P, sd = 3), p, P)
  lab <- sample.int(P, N, replace = TRUE)
  U <- centers[, lab] + matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(0, P, N); Y[cbind(lab, seq_len(N))] <- 1
  list(U = U, Y = Y)
}

test_that("signed RFF: selects a beta from the default grid and returns the CV matrix", {
  skip_unless_full()
  d <- make_case()
  sel <- suppressMessages(nmfkc.rff.beta.cv(d$Y, rank = 3, d$U, D = 40, plot = FALSE, nfolds = 3))
  expect_identical(sel$type, "signed")
  expect_length(sel$beta.candidates, 7L)
  expect_equal(dim(sel$objfunc), c(7L, 1L))
  expect_true(all(is.finite(sel$objfunc)))
  expect_true(sel$beta %in% sel$beta.candidates)
  expect_equal(sel$D, 40L)
  expect_equal(sel$n.used, 120L)
  ## the selected candidate is the minimizer
  expect_equal(sel$objfunc[which(sel$beta.candidates == sel$beta), 1], min(sel$objfunc))
  ## reproducible: same seed -> same objective values
  sel2 <- suppressMessages(nmfkc.rff.beta.cv(d$Y, rank = 3, d$U, D = 40, plot = FALSE, nfolds = 3))
  expect_equal(sel$objfunc, sel2$objfunc)
  ## the caller's random stream is untouched
  set.seed(77); r0 <- stats::runif(1)
  set.seed(77); invisible(suppressMessages(nmfkc.rff.beta.cv(d$Y, rank = 3, d$U, D = 10, beta = 0.2, plot = FALSE, nfolds = 2)))
  expect_identical(stats::runif(1), r0)
})

test_that("positive RFF over a (beta, D) grid, with a subsample", {
  skip_unless_full()
  d <- make_case()
  sel <- suppressWarnings(suppressMessages(
    nmfkc.rff.beta.cv(d$Y, rank = 3, d$U, beta = c(0.05, 0.2, 0.8), D = c(50, 200),
                      type = "positive", plot = FALSE, nfolds = 3, sample.size = 90,
                      verbose = FALSE)))
  expect_identical(sel$type, "positive")
  expect_equal(dim(sel$objfunc), c(3L, 2L))
  expect_equal(sel$n.used, 90L)
  expect_true(sel$beta %in% c(0.05, 0.2, 0.8)); expect_true(sel$D %in% c(50L, 200L))
  expect_equal(min(sel$objfunc), sel$objfunc[match(sel$beta, sel$beta.candidates), match(sel$D, sel$D.candidates)])
  ## the plot path runs
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_silent(suppressWarnings(suppressMessages(
    nmfkc.rff.beta.cv(d$Y, rank = 3, d$U, beta = c(0.1, 0.5), D = c(20, 40), type = "positive",
                      nfolds = 2, verbose = FALSE))))
})

test_that("the selected beta feeds the Gram route", {
  skip_unless_full()
  d <- make_case()
  sel <- suppressMessages(nmfkc.rff.beta.cv(d$Y, rank = 3, d$U, D = 30, plot = FALSE, nfolds = 3))
  g <- nmfkc.signed.rff.gram(d$Y, d$U, beta = sel$beta, D = sel$D, seed = 1)
  fit <- suppressMessages(nmfkc.signed(d$Y, A = g, rank = 3, maxit = 200, verbose = FALSE))
  expect_equal(dim(fit$C), c(3L, 30L))
  expect_error(nmfkc.rff.beta.cv(d$Y[, 1:50], rank = 3, d$U, D = 10), "ncol")
})

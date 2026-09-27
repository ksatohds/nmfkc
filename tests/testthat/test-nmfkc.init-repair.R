## A start that the column normalization cannot handle, and a NaN objective.
##
## With X.init = "kmeans" (the default) or "kmeans++", data with two or more
## all-zero observation columns could make k-means return the zero vector as a
## cluster centre.  Normalizing that basis column gave 0/0, the NaN spread
## through every update, and the fit stopped at the first convergence check
## with "missing value where TRUE/FALSE needed" -- or, with maxit < 10,
## returned an all-NaN result as if it were a fit.
##
## The repair acts only on a start that normalizes to a non-finite value, so
## every fit that could already succeed is bit-identical.  These tests pin both
## halves: the failing cases now fit, and the cases that already worked --
## including degenerate ones -- do not see the repair at all.

zero_cols <- function(z, seed = 1, P = 100, N = 10) {
  set.seed(seed)
  Y <- matrix(stats::rexp(P * N), P, N)
  if (z > 0) Y[, seq_len(z)] <- 0
  Y
}

test_that("zero observation columns no longer stop kmeans / kmeans++ starts", {
  skip_unless_full()
  ## the data of the original report: every rank failed here with "kmeans"
  for (z in 2:3) for (Q in 2:5) for (ini in c("kmeans", "kmeans++")) {
    Y <- zero_cols(z)
    f <- NULL
    msgs <- character(0)
    withCallingHandlers(
      f <- nmfkc(Y, rank = Q, X.init = ini, epsilon = 1e-6, maxit = 3000,
                 verbose = FALSE),
      message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                              invokeRestart("muffleMessage") })
    expect_true(is.finite(f$objfunc))
    expect_false(any(colSums(f$X) == 0))
    ## the fit still uses every column: B keeps all N of them
    expect_identical(ncol(f$B), ncol(Y))
    ## a message is given only when the repair ran (kmeans++ at rank 2 needed
    ## none), and it names the cause
    if (length(msgs)) expect_match(msgs, "all-zero column", all = FALSE)
  }
})

test_that("the repair says what it did", {
  skip_unless_full()
  expect_message(nmfkc(zero_cols(2), rank = 4, epsilon = 1e-6, maxit = 3000,
                       verbose = FALSE),
                 "rerun on the other 8")
})

test_that("fits that already worked do not see the repair", {
  skip_unless_full()
  ## one zero column: k-means does not isolate it, so no zero basis arises
  for (Q in 2:5)
    expect_no_message(nmfkc(zero_cols(1), rank = Q, epsilon = 1e-6,
                            maxit = 3000, verbose = FALSE))
  ## inits that never produce a zero basis
  for (ini in c("kmeansar", "nndsvd", "runif"))
    expect_no_message(nmfkc(zero_cols(2), rank = 4, X.init = ini,
                            epsilon = 1e-6, maxit = 3000, verbose = FALSE))
  ## X.restriction = "none" does not divide by the column sum, so a zero
  ## basis column is not a failure there: it stays zero, as it always did
  ## (a dead basis column keeps the fit from converging; the maxit warning is
  ## the behaviour it always had)
  expect_no_message(f <- suppressWarnings(
    nmfkc(zero_cols(2), rank = 4, X.restriction = "none",
          epsilon = 1e-6, maxit = 3000, verbose = FALSE)))
  expect_true(is.finite(f$objfunc))
})

test_that("Y with no positive entries stops with that reason", {
  skip_unless_full()
  Y0 <- matrix(0, 30, 10)
  for (ini in c("kmeans", "kmeans++", "kmeansar", "nndsvd", "runif"))
    expect_error(nmfkc(Y0, rank = 2, X.init = ini, verbose = FALSE),
                 "no positive entries")
  ## with maxit < 10 the per-iteration check never runs; this used to return
  ## an all-NaN object
  expect_error(nmfkc(Y0, rank = 2, maxit = 5, verbose = FALSE),
               "no positive entries")
  ## ...but X.restriction = "none" never divides, and still returns the
  ## all-zero fit it returned before
  f <- suppressWarnings(nmfkc(Y0, rank = 2, X.restriction = "none",
                              verbose = FALSE))
  expect_equal(f$objfunc, 0)
})

test_that("a user-supplied X.init is reported, not altered", {
  skip_unless_full()
  set.seed(2)
  Y  <- matrix(stats::rexp(300), 30, 10)
  X0 <- cbind(matrix(stats::runif(60) + 0.1, 30, 2), 0)
  expect_error(nmfkc(Y, rank = 3, X.init = X0, verbose = FALSE),
               "X.init has zero column\\(s\\) 3")
  ## the same matrix is fine where no normalization divides by the sum
  expect_no_error(suppressWarnings(
    nmfkc(Y, rank = 3, X.init = X0, X.restriction = "none",
          epsilon = 1e-6, maxit = 3000, verbose = FALSE)))
})

test_that("a repaired start leaves the caller's random stream alone", {
  skip_unless_full()
  ## the repair may draw random numbers (the refit, or the positive fill);
  ## nmfkc() restores the caller's RNG state on exit, repair or not
  set.seed(42); u1 <- stats::runif(3)
  set.seed(42)
  suppressMessages(nmfkc(zero_cols(2), rank = 4, epsilon = 1e-6, maxit = 3000,
                         verbose = FALSE))
  u2 <- stats::runif(3)
  expect_identical(u1, u2)
})

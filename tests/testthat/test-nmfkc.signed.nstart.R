## Multi-start for nmfkc.signed().  Signed models have many local minima
## (Theta = Cp - Cn takes both signs), so nstart > 1 repeats the fit from
## consecutive seeds and keeps the smallest objective.  These tests pin down
## (i) that the returned fit really is the best of the starts and equals the
## single-start fit from that seed, (ii) that $restarts records every start,
## (iii) that nstart = 1 is unchanged, and (iv) that Gram input works too.

make_case <- function(seed = 11, p = 4, N = 60, Q_obs = 3, D = 10) {
  set.seed(seed)
  U <- matrix(stats::rnorm(p * N), p, N)
  lab <- sample.int(Q_obs, N, replace = TRUE)
  Y <- matrix(0, Q_obs, N); Y[cbind(lab, seq_len(N))] <- 1
  rownames(Y) <- paste0("c", seq_len(Q_obs))
  list(U = U, Y = Y, D = D, beta = 0.5)
}

test_that("nstart > 1 returns the start with the smallest objective", {
  skip_unless_full()
  d <- make_case()
  Z <- nmfkc.signed.rff(d$U, beta = d$beta, D = d$D, seed = 1)$Z
  fit <- suppressMessages(
    nmfkc.signed(d$Y, A = Z, rank = 3, maxit = 60, seed = 100, nstart.signed = 4,
                 warm.start = FALSE, verbose = FALSE))
  expect_true(is.data.frame(fit$restarts))
  expect_equal(nrow(fit$restarts), 4L)
  expect_equal(fit$restarts$seed, 100:103)
  expect_equal(fit$objfunc, min(fit$restarts$objfunc))
  expect_identical(fit$nstart.signed, 4L)
  ## the winner equals the single-start fit run from that same seed
  best_seed <- fit$restarts$seed[which.min(fit$restarts$objfunc)]
  solo <- suppressMessages(
    nmfkc.signed(d$Y, A = Z, rank = 3, maxit = 60, seed = best_seed,
                 warm.start = FALSE, verbose = FALSE))
  expect_equal(fit$objfunc, solo$objfunc, tolerance = 1e-12)
  expect_equal(fit$X, solo$X, tolerance = 1e-10)
  expect_equal(fit$C, solo$C, tolerance = 1e-10)
})

test_that("nstart = 1 (the default) is unchanged and carries no restarts", {
  skip_unless_full()
  d <- make_case()
  Z <- nmfkc.signed.rff(d$U, beta = d$beta, D = d$D, seed = 1)$Z
  a <- suppressMessages(nmfkc.signed(d$Y, A = Z, rank = 3, maxit = 60, seed = 7,
                                     warm.start = FALSE, verbose = FALSE))
  b <- suppressMessages(nmfkc.signed(d$Y, A = Z, rank = 3, maxit = 60, seed = 7,
                                     warm.start = FALSE, verbose = FALSE, nstart.signed = 1))
  expect_null(a$restarts)
  expect_null(b$restarts)
  expect_equal(a$objfunc, b$objfunc, tolerance = 1e-12)
})

test_that("multi-start also works on Gram input", {
  skip_unless_full()
  d <- make_case()
  g <- nmfkc.signed.rff.gram(d$Y, d$U, beta = d$beta, D = d$D, seed = 1)
  fit <- suppressMessages(
    nmfkc.signed(d$Y, A = g, rank = 3, maxit = 60, seed = 5, nstart.signed = 3, verbose = FALSE))
  expect_equal(nrow(fit$restarts), 3L)
  expect_equal(fit$objfunc, min(fit$restarts$objfunc))
  expect_equal(dim(fit$B), c(3L, ncol(d$Y)))
})

test_that("nstart.signed does not leak into nstart through partial matching", {
  skip_unless_full()
  ## `$` on a list partial-matches, so a naive extra_args$nstart would read
  ## nstart.signed.  nstart (initialization of X in the non-negative warm
  ## start) and nstart.signed (restarts of the signed fit) are different
  ## things and must stay independent.
  d <- make_case()
  Z <- nmfkc.signed.rff(d$U, beta = d$beta, D = d$D, seed = 1)$Z
  ## warm.start = TRUE so that nstart would reach the warm-start fit if leaked
  a <- suppressMessages(nmfkc.signed(d$Y, A = Z, rank = 3, maxit = 40, seed = 3,
                                     verbose = FALSE))
  b <- suppressMessages(nmfkc.signed(d$Y, A = Z, rank = 3, maxit = 40, seed = 3,
                                     verbose = FALSE, nstart.signed = 2))
  ## the first restart of b uses seed 3, i.e. exactly the fit a computed
  expect_equal(min(b$restarts$objfunc), min(a$objfunc, b$objfunc), tolerance = 1e-12)
  expect_equal(b$restarts$objfunc[1], a$objfunc, tolerance = 1e-12)
})

test_that("an invalid nstart is rejected", {
  skip_unless_full()
  d <- make_case()
  Z <- nmfkc.signed.rff(d$U, beta = d$beta, D = d$D, seed = 1)$Z
  expect_error(nmfkc.signed(d$Y, A = Z, rank = 3, nstart.signed = 0), "nstart")
  expect_error(nmfkc.signed(d$Y, A = Z, rank = 3, nstart.signed = c(2, 3)), "nstart")
})

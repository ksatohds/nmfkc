## Regression tests for the WEIGHTED multiplicative-update path of
## nmfkc.signed().  Historically this path used an approximate (non-majorizing)
## sign-split for signed A, so on high-magnitude weighted data the Cp/Cn and X
## updates overshot and the objective diverged to Inf/NaN within a few
## iterations (found via Rprof profiling, 2026-07-20).  The path now damps each
## multiplicative factor and backtracks the damping exponent so the (penalized)
## objective is monotone non-increasing, which guarantees finiteness.  These
## tests lock in that behaviour.

make_signed_weighted <- function(scale_mag, seed = 42,
                                  Q_obs = 10, N = 60, D = 8, rank = 3,
                                  holdout = 0.2) {
  set.seed(seed)
  Xt <- matrix(abs(stats::rnorm(Q_obs * rank)), Q_obs, rank)
  Ct <- matrix(stats::rnorm(rank * D), rank, D)
  A  <- matrix(stats::rnorm(D * N) * scale_mag, D, N)         # signed, large
  Y  <- Xt %*% Ct %*% A + matrix(stats::rnorm(Q_obs * N) * scale_mag, Q_obs, N)
  W  <- matrix(1, Q_obs, N)
  set.seed(seed + 7L)
  W[sample.int(length(W), floor(holdout * length(W)))] <- 0   # binary mask
  list(Y = Y, A = A, W = W, rank = rank)
}

test_that("weighted nmfkc.signed stays finite on high-magnitude data (no NaN)", {
  d <- make_signed_weighted(scale_mag = 1000)
  fit <- suppressWarnings(suppressMessages(
    nmfkc.signed(d$Y, d$A, rank = d$rank, Y.weights = d$W,
                 maxit = 300, verbose = FALSE)))

  ## The whole solution and its objective trace must be finite (the bug
  ## produced objfunc.iter -> ...e+68 -> NaN here).
  expect_true(all(is.finite(fit$objfunc.iter)))
  expect_true(is.finite(fit$objfunc))
  expect_true(all(is.finite(fit$X)))
  expect_true(all(is.finite(fit$Cp)) && all(is.finite(fit$Cn)))
  expect_true(all(is.finite(fit$B)))
  ## X remains non-negative (Semi-NMF constraint) and Cp/Cn are non-negative.
  expect_true(all(fit$X >= 0))
  expect_true(all(fit$Cp >= 0) && all(fit$Cn >= 0))
  ## The fit is genuinely useful, not a degenerate stall.
  expect_gt(fit$r.squared.uncentered, 0.5)
})

test_that("weighted nmfkc.signed objective is monotone non-increasing", {
  ## The backtracking guard accepts a damped step only when it does not raise
  ## the objective, so the recorded trace must never increase.
  for (mag in c(1, 100, 1e4, 1e6)) {
    d <- make_signed_weighted(scale_mag = mag, seed = 99)
    fit <- suppressWarnings(suppressMessages(
      nmfkc.signed(d$Y, d$A, rank = d$rank, Y.weights = d$W,
                   maxit = 400, verbose = FALSE)))
    oi <- fit$objfunc.iter
    expect_true(all(is.finite(oi)), info = paste("mag", mag))
    ## Allow a negligible relative slack for floating-point noise.
    if (length(oi) > 1) {
      incr <- diff(oi)
      tol  <- 1e-8 * abs(utils::head(oi, -1))
      expect_true(all(incr <= tol), info = paste("mag", mag))
    }
  }
})

test_that("weighted nmfkc.signed fit quality is magnitude-independent", {
  ## After internal init-alignment the same problem scaled up by 10^k should
  ## reach essentially the same fit; the pre-fix code diverged for large k.
  r2 <- vapply(c(1, 1e3, 1e6, 1e9), function(mag) {
    d <- make_signed_weighted(scale_mag = mag, seed = 7)
    fit <- suppressWarnings(suppressMessages(
      nmfkc.signed(d$Y, d$A, rank = d$rank, Y.weights = d$W,
                   maxit = 500, verbose = FALSE)))
    fit$r.squared.uncentered
  }, numeric(1))
  expect_true(all(is.finite(r2)))
  expect_true(all(r2 > 0.5))
  expect_lt(diff(range(r2)), 0.2)   # tightly clustered across 9 decades
})

test_that("weighted element-wise CV runs to completion on high-magnitude data", {
  ## nmfkc.signed.ecv() drives the weighted path via a 0/1 hold-out mask; it
  ## must produce finite CV errors rather than NaN.
  d <- make_signed_weighted(scale_mag = 1000, N = 40)
  cv <- suppressWarnings(suppressMessages(
    nmfkc.signed.ecv(d$Y, d$A, rank = 2:3, nfolds = 3, maxit = 200)))
  expect_true(all(is.finite(cv$objfunc)))
  expect_true(all(is.finite(cv$sigma)))
})

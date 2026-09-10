## update.power in nmfkc.signed(): the square-root (Ding) form of the updates.

make_signed <- function(seed = 1L, P = 6L, N = 90L, D = 40L) {
  set.seed(seed)
  lab <- rep(seq_len(P), length.out = N)
  Y <- matrix(0, P, N); Y[cbind(lab, seq_len(N))] <- 1
  U <- matrix(rnorm(2 * N), 2, N) + 1.5 * rbind(cos(lab), sin(lab))
  list(Y = Y, Z = nmfkc.signed.rff(U, beta = 1, D = D, seed = seed)$Z, P = P)
}

test_that("update.power = 1 is the default and unchanged", {
  skip_unless_full()
  d <- make_signed()
  f0 <- nmfkc.signed(d$Y, d$Z, rank = 4, epsilon = 1e-6, maxit = 300, verbose = FALSE, seed = 1)
  f1 <- nmfkc.signed(d$Y, d$Z, rank = 4, epsilon = 1e-6, maxit = 300, verbose = FALSE, seed = 1,
                     update.power = 1)
  expect_equal(f0$X, f1$X); expect_equal(f0$objfunc, f1$objfunc)
})

test_that("the square-root form is monotone and reaches the same objective", {
  skip_unless_full()
  d <- make_signed()
  f1 <- nmfkc.signed(d$Y, d$Z, rank = 4, epsilon = 1e-7, maxit = 5000, verbose = FALSE, seed = 1)
  fr <- nmfkc.signed(d$Y, d$Z, rank = 4, epsilon = 1e-7, maxit = 5000, verbose = FALSE, seed = 1,
                     update.power = 0.5)
  o <- fr$objfunc.iter
  expect_true(all(diff(o) <= 1e-9 * abs(o[-length(o)])))
  expect_true(all(is.finite(fr$X)) && all(fr$X >= 0))
  expect_equal(unname(colSums(fr$X)), rep(1, 4), tolerance = 1e-10)
  expect_lt(abs(fr$objfunc - f1$objfunc) / f1$objfunc, 0.02)
})

test_that("the square-root form works with the row-sum floor", {
  skip_unless_full()
  d <- make_signed()
  tau <- 0.5 * 4 / d$P
  fr <- nmfkc.signed(d$Y, d$Z, rank = 4, epsilon = 1e-6, maxit = 2000, verbose = FALSE, seed = 1,
                     update.power = 0.5, X.rowSums.min = tau)
  expect_true(all(rowSums(fr$X) >= tau * (1 - 1e-6)))
})

test_that("bad update.power is refused", {
  skip_unless_full()
  d <- make_signed()
  expect_error(nmfkc.signed(d$Y, d$Z, rank = 4, maxit = 5, verbose = FALSE, update.power = 0), "update.power")
  expect_error(nmfkc.signed(d$Y, d$Z, rank = 4, maxit = 5, verbose = FALSE, update.power = 2), "update.power")
})

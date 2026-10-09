## nmfkc.ecv() keeps a user's Y.weights (2026-10-09).  Before, the fold mask
## was spliced in ahead of `...`, so a Y.weights given here was dropped.

ecv_data <- function(P = 12, N = 30, seed = 1) {
  set.seed(seed)
  matrix(stats::rexp(P * N), P, N)
}
q <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("nmfkc.ecv: Y.weights is used, and zero-weight cells are never held out", {
  skip_unless_full()
  Y <- ecv_data(); P <- nrow(Y); N <- ncol(Y)
  W <- matrix(1, P, N); W[, 1:10] <- 0
  a <- q(nmfkc.ecv(Y, rank = 1:2, nfolds = 3))
  b <- q(nmfkc.ecv(Y, rank = 1:2, nfolds = 3, Y.weights = W))
  expect_false(identical(a$objfunc, b$objfunc))
  expect_true(all(W[unlist(b$folds)] > 0))
  expect_setequal(unlist(b$folds), which(W > 0))
  ## the per-column vector form nmfkc() accepts gives the same result
  v <- q(nmfkc.ecv(Y, rank = 1:2, nfolds = 3, Y.weights = rep(c(0, 1), c(10, N - 10))))
  expect_identical(v$objfunc, b$objfunc)
  ## all-one weights: same folds and fits as no weights
  o <- q(nmfkc.ecv(Y, rank = 1:2, nfolds = 3, Y.weights = matrix(1, P, N)))
  expect_identical(o$folds, a$folds)
  expect_equal(o$objfunc, a$objfunc, tolerance = 1e-12)
})

test_that("nmfkc.ecv: bad Y.weights are refused", {
  skip_unless_full()
  Y <- ecv_data()
  expect_error(nmfkc.ecv(Y, rank = 1, nfolds = 3, Y.weights = matrix(1, 3, 3)),
               "same dimensions")
  expect_error(nmfkc.ecv(Y, rank = 1, nfolds = 3, Y.weights = rep(1, 7)),
               "length ncol")
  expect_error(nmfkc.ecv(Y, rank = 1, nfolds = 3, Y.weights = -matrix(1, nrow(Y), ncol(Y))),
               "non-negative")
})

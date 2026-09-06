## X.rowSums.min in nmfkc(): a row-sum floor on the basis matrix, the same
## constraint as in nmfkc.signed().

make_data <- function(seed = 1L, P = 6L, N = 90L) {
  set.seed(seed)
  lab <- rep(seq_len(P), length.out = N)
  Y <- matrix(0, P, N); Y[cbind(lab, seq_len(N))] <- 1
  U <- matrix(rnorm(2 * N), 2, N) + 1.5 * rbind(cos(lab), sin(lab))
  A <- nmfkc.kernel(U, beta = 1)
  list(Y = Y, A = A, U = U, P = P)
}

test_that("tau = 0 leaves nmfkc() unchanged", {
  d <- make_data()
  f0 <- nmfkc(d$Y, d$A, Q = 4, epsilon = 1e-6, maxit = 300, print.dims = FALSE)
  f1 <- nmfkc(d$Y, d$A, Q = 4, epsilon = 1e-6, maxit = 300, print.dims = FALSE,
              X.rowSums.min = 0)
  expect_equal(f0$X, f1$X)
  expect_equal(f0$objfunc, f1$objfunc)
})

test_that("the floor is met and the column sums are kept, at Q < P", {
  d <- make_data()
  tau <- 0.5 * 4 / d$P
  f <- nmfkc(d$Y, d$A, Q = 4, epsilon = 1e-6, maxit = 300, print.dims = FALSE,
             X.rowSums.min = tau)
  expect_true(all(rowSums(f$X) >= tau * (1 - 1e-6)))
  expect_equal(unname(colSums(f$X)), rep(1, 4), tolerance = 1e-10)
  expect_true(all(is.finite(f$X)))
  expect_true(all(is.finite(f$objfunc)))
})

test_that("the floor is met under the other column restrictions", {
  d <- make_data()
  for (r in c("colSqSums", "totalSum")) {
    f <- nmfkc(d$Y, d$A, Q = 4, epsilon = 1e-6, maxit = 100, print.dims = FALSE,
               X.restriction = r, X.rowSums.min = 0.05)
    expect_true(all(rowSums(f$X) >= 0.05 * (1 - 1e-6)), info = r)
  }
})

test_that("infeasible or incompatible settings are refused", {
  d <- make_data()
  expect_error(nmfkc(d$Y, d$A, Q = 4, maxit = 5, print.dims = FALSE,
                     X.rowSums.min = 4 / d$P + 0.1), "infeasible")
  expect_error(nmfkc(d$Y, d$A, Q = 4, maxit = 5, print.dims = FALSE,
                     X.rowSums.min = -1), "single number")
  expect_error(nmfkc(d$Y, d$A, Q = d$P, maxit = 5, print.dims = FALSE,
                     X.init = diag(d$P), X.restriction = "fixed",
                     X.rowSums.min = 0.1), "fixed")
})

test_that("the floor also holds on the Gram route and with weights", {
  d <- make_data()
  tau <- 0.5 * 4 / d$P
  g <- nmfkc.kernel.gram(d$Y, d$U, d$U, beta = 1)
  fg <- nmfkc(d$Y, g, Q = 4, epsilon = 1e-6, maxit = 300, print.dims = FALSE,
              X.rowSums.min = tau)
  expect_true(all(rowSums(fg$X) >= tau * (1 - 1e-6)))
  W <- matrix(runif(length(d$Y), 0.5, 1), nrow(d$Y))
  fw <- nmfkc(d$Y, d$A, Q = 4, epsilon = 1e-6, maxit = 300, print.dims = FALSE,
              Y.weights = W, X.rowSums.min = tau)
  expect_true(all(rowSums(fw$X) >= tau * (1 - 1e-6)))
})

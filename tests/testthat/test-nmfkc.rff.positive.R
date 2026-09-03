## Positive random features (Choromanski et al. 2021, transported to the
## Gaussian kernel): phi(u) = exp(sqrt(2b) w'u - 2b||u||^2 + kappa)/sqrt(D),
## E[phi(u)'phi(u')] = exp(2 kappa) exp(-b||u-u'||^2).  Non-negative, so the
## features and their Gram object go to nmfkc() (standard NMF-LAB).

make_prf_case <- function(seed = 11, p = 3, N = 40) {
  set.seed(seed)
  U <- scale(matrix(stats::rnorm(p * N), p, N))          # centred, scaled columns... (rows are features)
  U <- matrix(as.numeric(U), p, N)
  lab <- sample.int(3, N, replace = TRUE)
  Y <- matrix(0, 3, N); Y[cbind(lab, seq_len(N))] <- 1
  list(U = U, Y = Y, beta = 0.4)
}

test_that("positive random features are non-negative and unbiased for the Gaussian kernel", {
  skip_unless_full()
  d <- make_prf_case()
  ## The estimator is log-normal-tailed: its variance grows with beta*||u||^2,
  ## so the Monte-Carlo check uses a moderate bandwidth (beta = 0.1 on scaled
  ## inputs) where D = 40000 pins the kernel down, and a scaling check at the
  ## working bandwidth (error must shrink with D).
  beta0 <- 0.1
  K0 <- unclass(nmfkc.kernel(d$U, beta = beta0)); attributes(K0) <- list(dim = dim(K0))
  for (hyp in c(FALSE, TRUE)) {
    pr <- nmfkc.rff.positive(d$U, beta = beta0, D = 40000, seed = 2, hyperbolic = hyp)
    expect_equal(dim(pr$Z), c(40000L, 40L))
    expect_true(all(pr$Z >= 0))
    expect_identical(pr$pars$type, "positive")
    Kp <- crossprod(pr$Z) / exp(2 * pr$pars$kappa)
    expect_lt(max(abs(Kp - K0)), 0.05)
    expect_lt(mean(abs(Kp - K0)), 0.01)
    expect_equal(diag(Kp), rep(1, 40), tolerance = 0.05)
  }
  ## At the working bandwidth the error decreases with D (unbiased, finite variance)
  K <- unclass(nmfkc.kernel(d$U, beta = d$beta)); attributes(K) <- list(dim = dim(K))
  rmse <- function(D, seed) {
    e <- numeric(length(seed))
    for (i in seq_along(seed)) {
      pr <- nmfkc.rff.positive(d$U, beta = d$beta, D = D, seed = seed[i])
      e[i] <- sqrt(mean((crossprod(pr$Z) / exp(2 * pr$pars$kappa) - K)^2))
    }
    mean(e)
  }
  expect_lt(rmse(32000, 1:3), rmse(2000, 1:3) / 2)
})

test_that("pars reuse reproduces the map; kappa is fixed by the training data", {
  skip_unless_full()
  d <- make_prf_case()
  pr <- nmfkc.rff.positive(d$U, beta = d$beta, D = 30, seed = 2)
  pr2 <- nmfkc.rff.positive(d$U[, 1:7], pars = pr$pars)
  expect_equal(pr2$Z, pr$Z[, 1:7])
  expect_identical(pr2$pars, pr$pars)
  ## stabilizer: training features are bounded by 1/sqrt(D) (largest exponent is 0)
  expect_true(is.finite(pr$pars$kappa))
  expect_equal(max(pr$Z) * sqrt(30), 1, tolerance = 1e-12)
  ## Wrong input dimension is caught; seed does not disturb the caller's stream
  expect_error(nmfkc.rff.positive(d$U[1:2, ], pars = pr$pars), "p = 3")
  set.seed(77); r0 <- stats::runif(1)
  set.seed(77); invisible(nmfkc.rff.positive(d$U, beta = 1, D = 5, seed = 9)); expect_identical(stats::runif(1), r0)
  expect_error(nmfkc.rff.positive(d$U, D = 5), "'beta'")
})

test_that("nmfkc.rff.positive.gram() equals the one-shot products and nmfkc() accepts it", {
  skip_unless_full()
  d <- make_prf_case(N = 90)
  g <- nmfkc.rff.positive.gram(d$Y, d$U, beta = d$beta, D = 24, seed = 3, block.size = 17)
  expect_s3_class(g, "nmfkc.gram")
  expect_identical(g$type, "prf"); expect_false(g$signed)
  expect_equal(g$n.blocks, 6L)
  Z <- nmfkc.rff.positive(d$U, pars = g$pars)$Z
  expect_true(all(Z >= 0))
  expect_equal(g$S, tcrossprod(Z), tolerance = 1e-12)
  expect_equal(unname(g$G0), unname(tcrossprod(d$Y, Z)), tolerance = 1e-12)
  expect_equal(g$A.block(c(2L, 90L)), Z[, c(2L, 90L)])
  expect_output(print(g), "positive random features")
  ## Same fit as the explicit non-negative matrix, in nmfkc()
  fg <- suppressWarnings(nmfkc(d$Y, A = g, rank = 3, verbose = FALSE, maxit = 300))
  fm <- suppressWarnings(nmfkc(d$Y, A = Z, rank = 3, verbose = FALSE, maxit = 300))
  expect_equal(fg$iter, fm$iter)
  expect_equal(fg$C, fm$C, tolerance = 1e-8)
  expect_equal(fg$objfunc.iter, fm$objfunc.iter, tolerance = 1e-8)
  expect_true(all(fg$B >= 0))                       # memberships stay non-negative
  Znew <- nmfkc.rff.positive(d$U[, 1:10], pars = g$pars)$Z
  expect_equal(predict(fg, newA = Znew), fm$XB[, 1:10], tolerance = 1e-8)
  ## nmfkc.signed() accepts it too; the hyperbolic variant builds a 2*ceiling(D/2)-feature object
  fs <- suppressMessages(nmfkc.signed(d$Y, A = g, rank = 3, maxit = 100, verbose = FALSE))
  expect_equal(dim(fs$C), c(3L, 24L))
  gh <- nmfkc.rff.positive.gram(d$Y, d$U, beta = d$beta, D = 25, seed = 3, hyperbolic = TRUE)
  expect_equal(gh$D, 26L); expect_equal(dim(gh$S), c(26L, 26L))
  expect_error(nmfkc.rff.positive.gram(d$Y, d$U, beta = d$beta), "'D'")
})

## Stationarity + latent transition matrices G_d = Theta_d %*% X of an NMF-VAR fit.
## The reference numbers are the published Canada fit (paper eq. 19), for which
## the DQ x DQ latent route and the PD x PD companion route must agree to machine
## precision (the two companion matrices share their non-zero eigenvalues).

make_canada_fit <- function() {
  X  <- matrix(c(0.55, 0.40, 0.06, 0.00,
                 0.00, 0.15, 0.41, 0.44), nrow = 4, ncol = 2)
  Th <- matrix(c(1.01, 0.00,
                 0.53, 0.00,
                 0.00, 0.54,
                 0.00, 0.88), nrow = 2, ncol = 4)   # Q x P, D = 1, no intercept
  dimnames(X) <- list(paste0("v", 1:4), c("Condition1", "Condition2"))
  structure(list(X = X, C = Th), class = "nmfkc")
}

test_that("the legacy interface is preserved", {
  st <- nmfkc.ar.stationarity(make_canada_fit())
  ## the two components the function has always returned
  expect_true(is.numeric(st$spectral.radius) && length(st$spectral.radius) == 1L)
  expect_true(is.logical(st$stationary) && length(st$stationary) == 1L)
  expect_identical(st$stationary, st$spectral.radius < 1)
  expect_s3_class(st, "nmfkc.ar.stationarity")
})

test_that("latent transition matrix reproduces the published Canada fit", {
  st <- nmfkc.ar.stationarity(make_canada_fit())
  expect_equal(unname(st$G$lag1),
               matrix(c(0.7675, 0.0324, 0.0795, 0.6086), 2, 2), tolerance = 1e-9)
  ## row = effect at t, column = cause at t-d
  expect_equal(dimnames(st$G$lag1),
               list(c("Condition1", "Condition2"), c("Condition1", "Condition2")))
  expect_equal(st$dims, c(P = 4L, Q = 2L, D = 1L))
  expect_equal(st$G.sum, st$G$lag1)          # D = 1: the sum is the single block
})

test_that("both methods give the same spectral radius", {
  fit <- make_canada_fit()
  lat <- nmfkc.ar.stationarity(fit)                        # default
  cmp <- nmfkc.ar.stationarity(fit, method = "companion")
  expect_identical(lat$method, "latent")
  expect_identical(cmp$method, "companion")
  expect_equal(lat$spectral.radius, cmp$spectral.radius, tolerance = 1e-12)
  expect_equal(lat$spectral.radius, 0.7823267336, tolerance = 1e-9)
  expect_identical(lat$stationary, cmp$stationary)
  expect_true(lat$stationary)
})

test_that("column sums bracket rho and can be inconclusive", {
  st <- nmfkc.ar.stationarity(make_canada_fit())
  expect_equal(unname(st$colsum), c(1.01, 0.53, 0.54, 0.88), tolerance = 1e-12)
  expect_equal(st$colsum.max, 1.01, tolerance = 1e-12)
  ## max_j c_j >= 1 is inconclusive, yet the exact criterion says stationary
  expect_gt(st$colsum.max, 1)
  expect_true(st$stationary)
  ## min_j c_j <= rho(sum_d G_d) <= max_j c_j
  expect_gte(st$spectral.radius.sum, min(st$colsum) - 1e-12)
  expect_lte(st$spectral.radius.sum, st$colsum.max + 1e-12)
})

test_that("the two methods agree for multi-lag models", {
  set.seed(7)
  P <- 6; Q <- 3; D <- 3
  for (i in 1:10) {
    X <- matrix(runif(P * Q), P, Q); X <- sweep(X, 2, colSums(X), "/")
    Th <- matrix(runif(Q * P * D) * 0.15, Q, P * D)
    fit <- structure(list(X = X, C = cbind(Th, runif(Q))), class = "nmfkc")
    lat <- nmfkc.ar.stationarity(fit)
    cmp <- nmfkc.ar.stationarity(fit, method = "companion")
    expect_equal(lat$spectral.radius, cmp$spectral.radius, tolerance = 1e-8)
    expect_identical(lat$stationary, cmp$stationary)
  }
})

test_that("stationary means solve the VAR fixed point and df is the effective dimension", {
  set.seed(11)
  P <- 5; Q <- 2; D <- 2
  X <- matrix(runif(P * Q), P, Q); X <- sweep(X, 2, colSums(X), "/")
  Th <- matrix(runif(Q * P * D) * 0.2, Q, P * D); th0 <- runif(Q)
  fit <- structure(list(X = X, C = cbind(Th, th0)), class = "nmfkc")
  st <- nmfkc.ar.stationarity(fit)
  expect_true(st$stationary)
  ## mu_y must satisfy mu = sum_d Xi_d mu + X theta
  Xi <- lapply(1:D, function(d) X %*% Th[, ((d - 1) * P + 1):(d * P), drop = FALSE])
  rhs <- as.numeric(Reduce(`+`, lapply(Xi, function(M) M %*% st$mu.y)) + X %*% th0)
  expect_equal(unname(st$mu.y), rhs, tolerance = 1e-10)
  expect_equal(unname(st$mu.y), as.numeric(X %*% st$mu.b), tolerance = 1e-12)
  ## effective dimension Q(P + PD + 1 - Q), not the naive Q(P + PD + 1)
  expect_equal(st$df, Q * (P + P * D + 1 - Q))
})

test_that("the lag/intercept split comes from nmfkc.ar metadata, not dimensions", {
  ## P = 1 is the ambiguous case: P*D and P*D+1 are both multiples of P, so a
  ## dimension-only rule always reports an intercept.  nmfkc.ar() records the
  ## truth in attributes of A, which nmfkc() keeps in x$A.attr.
  d <- log(AirPassengers)
  for (ic in c(TRUE, FALSE)) {
    a <- nmfkc.ar(d, degree = 12, intercept = ic)
    f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 1, seed = 1, verbose = FALSE))
    st <- suppressWarnings(nmfkc.ar.stationarity(f))
    expect_length(st$G, 12L)                       # D = 12 either way
    expect_equal(unname(st$dims["D"]), 12L)
    expect_identical(is.null(st$theta0), !ic)      # intercept present iff asked for
  }
})

test_that("df is the effective dimension with and without an intercept", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  for (ic in c(TRUE, FALSE)) {
    a <- nmfkc.ar(Y, degree = 2, intercept = ic)
    f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, seed = 1, verbose = FALSE))
    st <- suppressWarnings(nmfkc.ar.stationarity(f))
    m <- 4 * 2 + as.integer(ic)                    # columns of Theta
    expect_equal(st$df, 2 * (4 + m - 2))           # Q(P + m - Q)
  }
})

test_that("a fit from a non-AR design is rejected", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  U <- matrix(rnorm(3 * 40), 3, 40)
  fk <- suppressWarnings(nmfkc(Y, nmfkc.kernel(U, beta = 0.5), rank = 2,
                               seed = 1, verbose = FALSE))
  expect_error(nmfkc.ar.stationarity(fk), "not fitted on an nmfkc.ar")
})

test_that("negative coefficients warn that the Perron-Frobenius bounds may fail", {
  X  <- matrix(c(0.6, 0.4, 0.3, 0.7), 2, 2)
  Th <- matrix(c(0.3, -0.1, 0.05, 0.2), 2, 2)      # signed Theta
  fit <- structure(list(X = X, C = Th), class = "nmfkc")
  expect_warning(nmfkc.ar.stationarity(fit), "negative entries")
})

test_that("a non-stationary fit warns and returns NA means", {
  X  <- matrix(c(1, 0, 0, 1), 2, 2)
  Th <- matrix(c(1.5, 0, 0, 1.5), 2, 2)          # rho = 1.5
  fit <- structure(list(X = X, C = cbind(Th, c(0.1, 0.1))), class = "nmfkc")
  expect_warning(st <- nmfkc.ar.stationarity(fit), "not stationary")
  expect_false(st$stationary)
  expect_true(all(is.na(st$mu.b)))
})

test_that("nmfkc.ar.stationarity works end to end on a real fit and prints", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  a <- nmfkc.ar(Y, degree = 2, intercept = TRUE)
  f <- nmfkc(a$Y, a$A, rank = 2, verbose = FALSE)
  st <- nmfkc.ar.stationarity(f)
  expect_length(st$G, 2)
  expect_equal(dim(st$G$lag1), c(2L, 2L))
  expect_true(all(st$G$lag1 >= 0))          # G_d = Theta_d X is non-negative
  expect_equal(st$spectral.radius,
               nmfkc.ar.stationarity(f, method = "companion")$spectral.radius,
               tolerance = 1e-8)
  expect_length(st$separability, 2)
  expect_true(all(st$separability >= 0 & st$separability <= 1))
  expect_output(print(st), "Latent transition matrices")
})

test_that("nmfkc.ar.stationarity carries X and Theta for the inference step", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  st <- suppressWarnings(nmfkc.ar.stationarity(f))
  expect_identical(st$X, f$X)
  expect_identical(st$Theta, f$C)          # FULL coefficient matrix, intercept included
})

test_that("nmfkc.ar.stationarity.inference bootstraps the invariants", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  st <- suppressWarnings(nmfkc.ar.stationarity(f))
  inf <- suppressWarnings(nmfkc.ar.stationarity.inference(st, a$Y, a$A, wild.B = 25))
  expect_s3_class(inf, "nmfkc.ar.stationarity.inference")
  ## the point estimate is the one nmfkc.ar.stationarity already reported
  expect_equal(inf$spectral.radius, st$spectral.radius, tolerance = 1e-8)
  ## uncertainty summaries are finite and ordered
  expect_true(is.finite(inf$spectral.radius.se) && inf$spectral.radius.se >= 0)
  expect_lte(inf$spectral.radius.ci.lower, inf$spectral.radius.ci.upper)
  expect_true(inf$p.nonstationary >= 0 && inf$p.nonstationary <= 1)
  expect_true(inf$truncate.rate >= 0 && inf$truncate.rate <= 1)
  expect_equal(inf$wild.B + inf$n.fail, 25L)
  expect_output(print(inf), "Bootstrap inference for NMF-VAR invariants")
})

test_that("the bootstrap is reproducible and leaves the RNG stream alone", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  st <- suppressWarnings(nmfkc.ar.stationarity(f))
  run <- function() suppressWarnings(
    nmfkc.ar.stationarity.inference(st, a$Y, a$A, wild.B = 15, wild.seed = 7))
  expect_equal(run()$spectral.radius.se, run()$spectral.radius.se)
  ## resetting the stream would make what came before the call irrelevant
  set.seed(1); invisible(run()); u1 <- runif(1)
  set.seed(2); invisible(run()); u2 <- runif(1)
  expect_false(isTRUE(all.equal(u1, u2)))
})

test_that("a fit whose object predates X/Theta can still be passed via fit=", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  st <- suppressWarnings(nmfkc.ar.stationarity(f))
  st$X <- NULL; st$Theta <- NULL                     # emulate an older object
  expect_error(nmfkc.ar.stationarity.inference(st, a$Y, a$A, wild.B = 5),
               "pass the fitted nmfkc object")
  inf <- suppressWarnings(nmfkc.ar.stationarity.inference(st, a$Y, a$A,
                                                          wild.B = 10, fit = f))
  expect_s3_class(inf, "nmfkc.ar.stationarity.inference")
})

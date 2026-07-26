## Latent VAR representation G_d = Theta_d %*% X of an NMF-VAR fit, its
## stationarity verdict, its plots and its bootstrap.
##
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

# ---- nmfkc.ar.latent -------------------------------------------------------

test_that("latent transition matrix reproduces the published Canada fit", {
  lat <- nmfkc.ar.latent(make_canada_fit())
  expect_s3_class(lat, "nmfkc.ar.latent")
  expect_equal(unname(lat$G$lag1),
               matrix(c(0.7675, 0.0324, 0.0795, 0.6086), 2, 2), tolerance = 1e-9)
  ## row = effect at t, column = cause at t-d
  expect_equal(dimnames(lat$G$lag1),
               list(c("Condition1", "Condition2"), c("Condition1", "Condition2")))
  expect_equal(lat$dims, c(P = 4L, Q = 2L, D = 1L))
  expect_equal(lat$G.sum, lat$G$lag1)          # D = 1: the sum is the single block
})

test_that("the eigenvalues are all of the model's roots", {
  lat <- nmfkc.ar.latent(make_canada_fit())
  expect_length(lat$eigenvalues, 2L)           # D * Q
  expect_equal(lat$spectral.radius, max(Mod(lat$eigenvalues)))
  expect_equal(lat$spectral.radius, 0.7823267336, tolerance = 1e-9)
  ## all real here, so no cycle is reported
  expect_true(all(abs(Im(lat$eigenvalues)) < 1e-8))
  expect_true(is.na(lat$cycle.period))
})

test_that("a seasonal fit reports the period of its dominant complex root", {
  a <- nmfkc.ar(log(AirPassengers), degree = 12, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 1, seed = 1, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  expect_length(lat$eigenvalues, 12L)
  expect_true(is.finite(lat$cycle.period))
  ## monthly data with a 12-lag design: the dominant cycle is the year
  expect_gt(lat$cycle.period, 9)
  expect_lt(lat$cycle.period, 15)
})

test_that("stationary means solve the VAR fixed point and df is the effective dimension", {
  set.seed(11)
  P <- 5; Q <- 2; D <- 2
  X <- matrix(runif(P * Q), P, Q); X <- sweep(X, 2, colSums(X), "/")
  Th <- matrix(runif(Q * P * D) * 0.2, Q, P * D); th0 <- runif(Q)
  fit <- structure(list(X = X, C = cbind(Th, th0)), class = "nmfkc")
  lat <- nmfkc.ar.latent(fit)
  expect_true(lat$stationary)
  ## mu_y must satisfy mu = sum_d Xi_d mu + X theta
  Xi <- lapply(1:D, function(d) X %*% Th[, ((d - 1) * P + 1):(d * P), drop = FALSE])
  rhs <- as.numeric(Reduce(`+`, lapply(Xi, function(M) M %*% lat$mu.y)) + X %*% th0)
  expect_equal(unname(lat$mu.y), rhs, tolerance = 1e-10)
  expect_equal(unname(lat$mu.y), as.numeric(X %*% lat$mu.b), tolerance = 1e-12)
  ## p* is mu.b normalized to sum 1
  expect_equal(sum(lat$p.star), 1, tolerance = 1e-12)
  expect_equal(unname(lat$p.star), unname(lat$mu.b / sum(lat$mu.b)), tolerance = 1e-14)
  ## effective dimension Q(P + PD + 1 - Q), not the naive Q(P + PD + 1)
  expect_equal(lat$df, Q * (P + P * D + 1 - Q))
})

test_that("the lag/intercept split comes from nmfkc.ar metadata, not dimensions", {
  ## P = 1 is the ambiguous case: P*D and P*D+1 are both multiples of P, so a
  ## dimension-only rule always reports an intercept.  nmfkc.ar() records the
  ## truth in attributes of A, which nmfkc() keeps in x$A.attr.
  d <- log(AirPassengers)
  for (ic in c(TRUE, FALSE)) {
    a <- nmfkc.ar(d, degree = 12, intercept = ic)
    f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 1, seed = 1, verbose = FALSE))
    lat <- suppressWarnings(nmfkc.ar.latent(f))
    expect_length(lat$G, 12L)                      # D = 12 either way
    expect_equal(unname(lat$dims["D"]), 12L)
    expect_identical(is.null(lat$theta0), !ic)     # intercept present iff asked for
  }
})

test_that("df is the effective dimension with and without an intercept", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  for (ic in c(TRUE, FALSE)) {
    a <- nmfkc.ar(Y, degree = 2, intercept = ic)
    f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, seed = 1, verbose = FALSE))
    lat <- suppressWarnings(nmfkc.ar.latent(f))
    m <- 4 * 2 + as.integer(ic)                    # columns of Theta
    expect_equal(lat$df, 2 * (4 + m - 2))          # Q(P + m - Q)
  }
})

test_that("nmfkc() records the constraint that decides identification", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, seed = 1, verbose = FALSE))
  expect_identical(f$X.restriction, "colSums")     # the default pins the scale
  expect_identical(nmfkc.ar.latent(f)$X.restriction, "colSums")
  f2 <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, seed = 1, verbose = FALSE,
                               X.restriction = "none"))
  expect_identical(f2$X.restriction, "none")
})

test_that("a fit from a non-AR design is rejected", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  U <- matrix(rnorm(3 * 40), 3, 40)
  fk <- suppressWarnings(nmfkc(Y, nmfkc.kernel(U, beta = 0.5), rank = 2,
                               seed = 1, verbose = FALSE))
  expect_error(nmfkc.ar.latent(fk), "not fitted on an nmfkc.ar")
  expect_error(nmfkc.ar.stationarity(fk), "not fitted on an nmfkc.ar")
})

test_that("negative coefficients warn that the Perron-Frobenius bounds may fail", {
  X  <- matrix(c(0.6, 0.4, 0.3, 0.7), 2, 2)
  Th <- matrix(c(0.3, -0.1, 0.05, 0.2), 2, 2)      # signed Theta
  fit <- structure(list(X = X, C = Th), class = "nmfkc")
  expect_warning(nmfkc.ar.latent(fit), "negative entries")
})

test_that("a non-stationary fit warns and returns NA means", {
  X  <- matrix(c(1, 0, 0, 1), 2, 2)
  Th <- matrix(c(1.5, 0, 0, 1.5), 2, 2)          # rho = 1.5
  fit <- structure(list(X = X, C = cbind(Th, c(0.1, 0.1))), class = "nmfkc")
  expect_warning(lat <- nmfkc.ar.latent(fit), "not stationary")
  expect_false(lat$stationary)
  expect_true(all(is.na(lat$mu.b)))
  expect_true(all(is.na(lat$p.star)))
})

test_that("nmfkc.ar.latent works end to end on a real fit and prints", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  a <- nmfkc.ar(Y, degree = 2, intercept = TRUE)
  f <- nmfkc(a$Y, a$A, rank = 2, verbose = FALSE)
  lat <- nmfkc.ar.latent(f)
  expect_length(lat$G, 2)
  expect_equal(dim(lat$G$lag1), c(2L, 2L))
  expect_true(all(lat$G$lag1 >= 0))          # G_d = Theta_d X is non-negative
  expect_length(lat$separability, 2)
  expect_true(all(lat$separability >= 0 & lat$separability <= 1))
  expect_output(print(lat), "Latent transition matrices")
  ## the pieces the inference step needs, so fit= is unnecessary
  expect_identical(lat$X, f$X)
  expect_identical(lat$Theta, f$C)          # FULL coefficient matrix, intercept included
})

# ---- nmfkc.ar.stationarity -------------------------------------------------

test_that("nmfkc.ar.stationarity answers the yes/no question only", {
  st <- nmfkc.ar.stationarity(make_canada_fit())
  expect_s3_class(st, "nmfkc.ar.stationarity")
  expect_true(is.numeric(st$spectral.radius) && length(st$spectral.radius) == 1L)
  expect_true(is.logical(st$stationary) && length(st$stationary) == 1L)
  expect_identical(st$stationary, st$spectral.radius < 1)
  ## the dynamics live in nmfkc.ar.latent(), not here
  expect_null(st$G)
  expect_null(st$mu.b)
  expect_output(print(st), "nmfkc.ar.latent")
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

test_that("stationarity accepts a fit or an already-computed latent object", {
  fit <- make_canada_fit()
  expect_identical(nmfkc.ar.stationarity(nmfkc.ar.latent(fit)),
                   nmfkc.ar.stationarity(fit))
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

# ---- plotting --------------------------------------------------------------

test_that("plot.nmfkc.ar.latent draws every panel and restores par()", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 60)) + 1, 4, 60)
  a <- nmfkc.ar(Y, degree = 3, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, seed = 1, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  before <- par("mfrow")
  for (ty in c("roots", "graph", "irf", "lag")) {
    expect_silent(plot(lat, type = ty))
    expect_identical(par("mfrow"), before)    # multi-panel types must restore it
  }
  expect_silent(plot(lat, type = "graph", lag = 2))
  expect_silent(plot(lat, type = "irf", horizon = 4))
  ## `...` reaches the underlying plot rather than colliding with our defaults
  expect_silent(plot(lat, type = "roots", main = "custom", col = "black"))
  ## D = 3 has no phase portrait, and the bracket lives on the other class
  expect_error(plot(lat, type = "phase"), "single lag")
  expect_error(plot(lat, type = "bracket"))
})

test_that("plot.nmfkc.ar.stationarity draws the bracket", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  st <- nmfkc.ar.stationarity(make_canada_fit())
  expect_silent(plot(st))
  expect_silent(plot(st, main = "custom"))
})

test_that("the phase portrait needs a stationary D = 1, Q = 2 fit with an intercept", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, seed = 1, verbose = FALSE))
  expect_silent(plot(suppressWarnings(nmfkc.ar.latent(f)), type = "phase"))
  ## no intercept -> no fixed point to draw
  a0 <- nmfkc.ar(Y, degree = 1, intercept = FALSE)
  f0 <- suppressWarnings(nmfkc(a0$Y, a0$A, rank = 2, seed = 1, verbose = FALSE))
  expect_error(plot(suppressWarnings(nmfkc.ar.latent(f0)), type = "phase"),
               "intercept")
  ## non-stationary -> mu.b is NA
  X  <- matrix(c(1, 0, 0, 1), 2, 2)
  Th <- matrix(c(1.5, 0, 0, 1.5), 2, 2)
  ns <- suppressWarnings(nmfkc.ar.latent(
    structure(list(X = X, C = cbind(Th, c(0.1, 0.1))), class = "nmfkc")))
  expect_error(plot(ns, type = "phase"), "stationary")
})

test_that("the rank-1 panels work (AirPassengers-style Q = 1)", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  a <- nmfkc.ar(log(AirPassengers), degree = 12, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 1, seed = 1, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  before <- par("mfrow")
  for (ty in c("roots", "graph", "lag", "irf"))
    expect_silent(plot(lat, type = ty))
  expect_identical(par("mfrow"), before)      # Q = 1 "lag" must not set mfrow
  expect_silent(plot(nmfkc.ar.stationarity(f)))
})

# ---- nmfkc.ar.latent.inference ---------------------------------------------

test_that("nmfkc.ar.latent.inference bootstraps the latent VAR", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  inf <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 25))
  expect_s3_class(inf, "nmfkc.ar.latent.inference")
  ## the point estimates are the ones nmfkc.ar.latent already reported
  expect_equal(inf$G$lag1, lat$G$lag1, tolerance = 0)
  expect_equal(inf$spectral.radius, lat$spectral.radius, tolerance = 1e-8)
  ## uncertainty summaries are finite and ordered
  expect_true(is.finite(inf$spectral.radius.se) && inf$spectral.radius.se >= 0)
  expect_lte(inf$spectral.radius.ci.lower, inf$spectral.radius.ci.upper)
  expect_true(inf$p.nonstationary >= 0 && inf$p.nonstationary <= 1)
  expect_true(inf$truncate.rate >= 0 && inf$truncate.rate <= 1)
  expect_equal(inf$wild.B + inf$n.fail, 25L)
  expect_output(print(inf), "Bootstrap inference for the latent VAR")
})

test_that("the coefficient table covers every entry of every G_d", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
  a <- nmfkc.ar(Y, degree = 2, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  inf <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 20))
  ct <- inf$coefficients
  expect_equal(nrow(ct), 2 * 2 * 2)                       # Q * Q * D
  expect_named(ct, c("Lag", "Basis", "Covariate", "Estimate", "BSE",
                     "CI_low", "CI_high", "p_value"))
  expect_setequal(ct$Lag, 1:2)
  ## the table and the matrix form must agree, entry by entry
  for (d in 1:2) {
    sub <- ct[ct$Lag == d, ]
    m <- matrix(sub$Estimate, 2, 2,
                dimnames = list(unique(sub$Basis), unique(sub$Covariate)))
    expect_equal(m, inf$G[[d]], tolerance = 0)
  }
  ## p-values invert the centred distribution, so they are never all zero
  expect_true(all(ct$p_value >= 0 & ct$p_value <= 1))
  expect_true(any(ct$p_value > 0))
  expect_true(all(is.finite(ct$BSE) & ct$BSE >= 0))
})

test_that("entry-level inference warns when the column scale of X is free", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE,
                              X.restriction = "none"))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  ## the KKT diagnostic may fire too on this small design; only the
  ## identification warning is under test here
  inf <- NULL
  withCallingHandlers(
    expect_warning(inf <- nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 10),
                   "not identified"),
    warning = function(w) if (grepl("complementarity", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  expect_false(inf$identified)
  expect_output(print(inf), "only the diagonal")
  ## the default constraint does pin it down, so no warning there
  f2 <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  lat2 <- suppressWarnings(nmfkc.ar.latent(f2))
  i2 <- suppressWarnings(nmfkc.ar.latent.inference(lat2, a$Y, a$A, wild.B = 10))
  expect_true(i2$identified)
})

test_that("the bootstrap is reproducible, parallel-safe and leaves the RNG alone", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  run <- function(cores = 1L) suppressWarnings(
    nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 15, wild.seed = 7,
                              cores = cores))
  r1 <- run()
  expect_identical(run(), r1)
  ## the multipliers are pre-drawn, so more workers cannot change the answer
  expect_identical(run(cores = 2L), r1)
  ## resetting the stream would make what came before the call irrelevant
  set.seed(1); invisible(run()); u1 <- runif(1)
  set.seed(2); invisible(run()); u2 <- runif(1)
  expect_false(isTRUE(all.equal(u1, u2)))
})

test_that("an object without X/Theta can still be given the fit", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  lat$X <- NULL; lat$Theta <- NULL
  expect_error(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 5),
               "pass the fitted nmfkc object")
  inf <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A,
                                                    wild.B = 10, fit = f))
  expect_s3_class(inf, "nmfkc.ar.latent.inference")
})

test_that("the wrong class is rejected by each entry point", {
  fit <- make_canada_fit()
  expect_error(nmfkc.ar.latent.inference(nmfkc.ar.stationarity(fit),
                                         matrix(1, 4, 4), matrix(1, 4, 4)),
               "nmfkc.ar.latent")
})

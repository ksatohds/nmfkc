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
  st <- nmfkc.ar.stationarity(fit <- make_canada_fit())
  ## The c_j are the column sums of X Lambda, i.e. colSums(X) %*% Lambda.  The
  ## paper's printed X is rounded, so colSums(X) = (1.01, 1.00) rather than
  ## (1, 1) and the weighted values differ slightly from colSums(Lambda) =
  ## (1.01, 0.53, 0.54, 0.88) -- which is exactly why the weights are needed.
  expect_equal(colSums(fit$X), c(Condition1 = 1.01, Condition2 = 1.00),
               tolerance = 1e-12)
  expect_equal(unname(st$colsum), c(1.0201, 0.5353, 0.54, 0.88), tolerance = 1e-12)
  expect_equal(st$colsum.max, 1.0201, tolerance = 1e-12)
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

## Assert one warning while letting the unrelated ones (e.g. the KKT diagnostic
## on these small designs) pass without failing the run.
quietly <- function(expr, keep) {
  val <- NULL
  withCallingHandlers(
    expect_warning(val <- expr, keep),
    warning = function(w) if (!grepl(keep, conditionMessage(w)))
      invokeRestart("muffleWarning"))
  val
}

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
  expect_true(inf$prob.nonstationary >= 0 && inf$prob.nonstationary <= 1)
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

test_that("a separable, column-normalised X does NOT certify uniqueness", {
  ## The counterexample that demotes `identified`: X = I is perfectly separable
  ## with unit column sums, yet a non-permutation T keeps both factors
  ## non-negative and the product unchanged -- so G is not pinned down.
  X <- diag(2); Th <- matrix(c(2, 1, 3, 1), 2, 2)     # [[2,3],[1,1]]
  T <- matrix(c(.9, .1, .1, .9), 2, 2)
  expect_equal(colSums(X), c(1, 1))                   # scale is fixed
  expect_equal(unname(apply(X / rowSums(X), 2, max)), c(1, 1))  # separable
  expect_true(all(X %*% T >= 0))
  expect_equal(colSums(X %*% T), c(1, 1))             # still admissible
  expect_true(all(solve(T) %*% Th >= -1e-12))
  expect_equal(X %*% Th, (X %*% T) %*% (solve(T) %*% Th))
  expect_false(all(T %in% c(0, 1)))                   # not a permutation
  ## and the latent transition matrix genuinely differs between the two
  expect_false(isTRUE(all.equal(Th %*% X, (solve(T) %*% Th) %*% (X %*% T))))
})

test_that("rank 1 has no rotation, and Q > 1 is only ever APPROXIMATE", {
  a <- nmfkc.ar(log(AirPassengers), degree = 2, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 1, seed = 1, verbose = FALSE))
  inf <- suppressWarnings(nmfkc.ar.latent.inference(
    suppressWarnings(nmfkc.ar.latent(f)), a$Y, a$A, wild.B = 5))
  expect_true(inf$identified.approx)                 # no rotation at Q = 1

  ## Q = 2: the anchor / pure tests can pass, but only up to the thresholds, so
  ## the function still warns that this is not a certificate.
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a2 <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f2 <- suppressWarnings(nmfkc(a2$Y, a2$A, rank = 2, verbose = FALSE))
  lat2 <- suppressWarnings(nmfkc.ar.latent(f2))
  expect_identical(f2$X.restriction, "colSums")
  i2 <- quietly(nmfkc.ar.latent.inference(lat2, a2$Y, a2$A, wild.B = 10),
                "APPROXIMATELY|not even approximately")
  if (i2$identified.approx) {
    ## passing means it is approximate, and the print-out says so
    expect_output(print(i2), "APPROXIMATELY")
    expect_true(is.finite(i2$identification$anchor.margin))
  } else {
    expect_output(print(i2), "DESCRIPTIVE")
  }
  ## a scale that is not per-column normalised is called out by name
  f3 <- suppressWarnings(nmfkc(a2$Y, a2$A, rank = 2, verbose = FALSE,
                               X.restriction = "none"))
  i3 <- quietly(nmfkc.ar.latent.inference(
    suppressWarnings(nmfkc.ar.latent(f3)), a2$Y, a2$A, wild.B = 5),
    "does not fix the per-column scale")
  expect_false(i3$identified.approx)
})

test_that("the KKT diagnostic reports both margins", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, epsilon = 1e-8, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  inf <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 5))
  expect_named(inf$kkt, c("n.boundary", "all.positive", "ratio",
                          "delta.dual", "delta.prim"))
  if (inf$kkt$n.boundary > 0 && inf$kkt$n.boundary < length(lat$Theta)) {
    ## delta.prim is the smallest surviving coefficient, so it clears kkt.tol
    expect_gt(inf$kkt$delta.prim, 1e-3)
    expect_true(is.finite(inf$kkt$delta.dual))
  }
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

test_that("the colsum bracket is weighted by colSums(X)", {
  ## Counterexample: the raw column sums of Lambda would certify stationarity
  ## for a fit whose spectral radius is 1.6.
  X  <- matrix(c(20, 20), 2, 1)
  Th <- matrix(c(0.04, 0.04), 1, 2)
  lat <- suppressWarnings(nmfkc.ar.latent(
    structure(list(X = X, C = Th), class = "nmfkc")))
  expect_equal(unname(lat$colsum), c(1.6, 1.6), tolerance = 1e-12)
  expect_equal(lat$colsum.max, 1.6, tolerance = 1e-12)
  expect_equal(lat$spectral.radius, 1.6, tolerance = 1e-12)
  expect_false(lat$stationary)
  expect_gt(lat$colsum.max, 1)                 # correctly inconclusive
  expect_lt(max(colSums(Th)), 1)               # the unweighted rule would not be
  ## the bracket must actually bracket
  expect_gte(lat$spectral.radius.sum, min(lat$colsum) - 1e-12)
  expect_lte(lat$spectral.radius.sum, lat$colsum.max + 1e-12)
})

test_that("the weighted and unweighted rules agree under column normalisation", {
  set.seed(5)
  for (i in 1:5) {
    P <- 5; Q <- 2; D <- 2
    X <- matrix(runif(P * Q), P, Q); X <- sweep(X, 2, colSums(X), "/")
    Th <- matrix(runif(Q * P * D) * 0.15, Q, P * D)
    lat <- nmfkc.ar.latent(structure(list(X = X, C = cbind(Th, runif(Q))),
                                     class = "nmfkc"))
    Lam <- Th[, 1:P, drop = FALSE] + Th[, (P + 1):(2 * P), drop = FALSE]
    expect_equal(unname(lat$colsum), unname(colSums(Lam)), tolerance = 1e-12)
  }
})

test_that("the bootstrapped radius is the companion radius, not rho(sum G_d)", {
  ## For D > 1 the two differ: with G_1 = 0.2, G_2 = 0.1 the companion radius is
  ## 0.4317 while rho(G_1 + G_2) = 0.3.  The inference must report the former,
  ## which is what nmfkc.ar.latent() reports.
  ## P = 2, Q = 1, D = 2 (ncol(Theta) = 4 = P*D, so the split is unambiguous;
  ## P = 1 would be the ambiguous case documented for .nmfvar.parts)
  X  <- matrix(c(0.5, 0.5), 2, 1)
  Th <- matrix(c(0.2, 0.2, 0.1, 0.1), 1, 4)       # G_1 = 0.2, G_2 = 0.1
  lat <- nmfkc.ar.latent(structure(list(X = X, C = Th), class = "nmfkc"))
  expect_equal(unname(unlist(lat$G)), c(0.2, 0.1), tolerance = 1e-12)
  expect_equal(lat$spectral.radius, 0.4316625, tolerance = 1e-6)
  expect_equal(lat$spectral.radius.sum, 0.3, tolerance = 1e-12)
  expect_true(lat$alternating)                  # the second root is negative
  expect_false(lat$non.oscillatory)
  expect_output(print(lat), "alternating")

  set.seed(3)
  Y <- matrix(abs(rnorm(3 * 40)) + 1, 3, 40)
  a <- nmfkc.ar(Y, degree = 3, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, epsilon = 1e-8, verbose = FALSE))
  l <- suppressWarnings(nmfkc.ar.latent(f))
  inf <- suppressWarnings(nmfkc.ar.latent.inference(l, a$Y, a$A, wild.B = 8))
  ## the two must now agree; before the fix inf used rho(sum_d G_d)
  expect_equal(inf$spectral.radius, l$spectral.radius, tolerance = 1e-10)
  expect_true(abs(l$spectral.radius - l$spectral.radius.sum) > 1e-6)  # D > 1
  ## and the stationarity verdict is the same either way under non-negativity
  expect_identical(l$spectral.radius < 1, l$spectral.radius.sum < 1)
})

test_that("Q = 2, D = 1 with G >= 0 can never produce a complex pair", {
  ## discriminant = (a - d)^2 + 4bc >= 0, so the phase portrait cannot spiral;
  ## the absence of a cycle there is a property of the design, not a finding.
  set.seed(11)
  for (i in 1:200) {
    X <- matrix(runif(4), 2, 2); X <- sweep(X, 2, colSums(X), "/")
    Th <- matrix(runif(4) * 0.4, 2, 2)
    lat <- suppressWarnings(nmfkc.ar.latent(
      structure(list(X = X, C = Th), class = "nmfkc")))
    expect_true(all(abs(Im(lat$eigenvalues)) < 1e-8))
    expect_true(is.na(lat$cycle.period))
  }
  ## the print-out says so
  lat <- nmfkc.ar.latent(make_canada_fit())
  expect_output(print(lat), "not a finding")
})

test_that("the latent recursion is VARMA: dropping the MA term is not exact", {
  set.seed(1)
  P <- 3; Q <- 2; n <- 40
  X <- matrix(runif(P * Q), P, Q); X <- sweep(X, 2, colSums(X), "/")
  Th <- matrix(runif(Q * P) * 0.3, Q, P)
  G <- Th %*% X
  Y <- matrix(abs(rnorm(P * (n + 1))) + 1, P, n + 1)
  b <- Th %*% Y[, 1:n]                       # b_t = Theta y_{t-1}
  e <- Y[, 2:(n + 1)] - X %*% b              # y_t = X b_t + e_t
  lhs <- b[, 2:n]
  ## the VAR form is NOT an identity ...
  expect_gt(max(abs(lhs - G %*% b[, 1:(n - 1)])), 1e-3)
  ## ... but the VARMA form is, to machine precision
  expect_equal(lhs, G %*% b[, 1:(n - 1)] + Th %*% e[, 1:(n - 1)],
               tolerance = 1e-12)
})

test_that("identification needs BOTH anchor rows in X and pure columns in Theta", {
  ## X separable but Theta with no pure column -> the counterexample family
  X  <- diag(2); Th <- matrix(c(2, 1, 3, 1), 2, 2)
  id <- nmfkc:::.nmfvar.identified(X, Th)
  expect_true(all(id$anchor))          # X = I has an anchor row per column
  expect_false(all(id$pure))           # no column of Theta is pure
  expect_false(id$ok)

  ## Theta with pure columns but X not separable -> T >= 0 not forced
  X2 <- matrix(c(0.6, 0.4, 0.4, 0.6), 2, 2)
  Th2 <- matrix(c(0.5, 0, 0, 0.5), 2, 2)
  id2 <- nmfkc:::.nmfvar.identified(X2, Th2)
  expect_false(all(id2$anchor))
  expect_true(all(id2$pure))
  expect_false(id2$ok)

  ## both sides -> identified (anchor rows for BOTH columns of X this time)
  X3 <- matrix(c(1, 0.4, 0, 0, 0.6, 1), 3, 2)
  id3 <- nmfkc:::.nmfvar.identified(X3, Th2)
  expect_true(all(id3$anchor))
  expect_true(all(id3$pure))
  expect_true(id3$ok)
})

test_that("the monomial step is what makes T a permutation", {
  ## T >= 0 and T^-1 >= 0 together force T to be a scaled permutation, and the
  ## column-sum constraint removes the scale.  Positive test:
  for (p in list(c(1, 2), c(2, 1))) {
    T <- diag(2)[p, , drop = FALSE]
    expect_true(all(T >= 0) && all(solve(T) >= 0))
  }
  ## and the counterexample's T is non-negative but its inverse is not, which is
  ## exactly why the anchor condition alone does not suffice
  T <- matrix(c(.9, .1, .1, .9), 2, 2)
  expect_true(all(T >= 0))
  expect_true(any(solve(T) < 0))
  ## a non-monomial non-negative T can never have a non-negative inverse
  set.seed(2)
  for (i in 1:2000) {
    M <- matrix(runif(4), 2, 2)          # dense, so not monomial
    expect_true(any(solve(M) < 0))
  }
})

test_that("Canada is identified, and the diagnostic says which side would fail", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  ## a fit with a non-separable basis reports the anchor side as the culprit
  X  <- matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2)
  Th <- matrix(c(0.4, 0, 0, 0.4), 2, 2)
  lat <- nmfkc.ar.latent(structure(list(X = X, C = Th, X.restriction = "colSums"),
                                   class = "nmfkc"))
  inf <- quietly(nmfkc.ar.latent.inference(lat, matrix(1, 2, 8), matrix(1, 2, 8),
                                           wild.B = 5),
                 "no anchor row")
  expect_false(inf$identified.approx)
  expect_false(all(inf$identification$anchor))
  expect_true(all(inf$identification$pure))
})

test_that("the bootstrap re-fits inherit the original estimator's settings", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE,
                              method = "KL", X.restriction = "colSqSums"))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  expect_identical(lat$fit.args$method, "KL")
  expect_identical(lat$fit.args$X.restriction, "colSqSums")
  inf <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 5))
  expect_identical(inf$refit.args$method, "KL")
  expect_identical(inf$refit.args$X.restriction, "colSqSums")
  ## refit.args= overrides
  i2 <- suppressWarnings(nmfkc.ar.latent.inference(
    lat, a$Y, a$A, wild.B = 5, refit.args = list(method = "EU")))
  expect_identical(i2$refit.args$method, "EU")
})

test_that("only per-column normalisation fixes the scale, and 'fixed' is special", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  ## totalSum constrains the grand total, so columns can still trade scale
  ftot <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE,
                                 X.restriction = "totalSum"))
  ltot <- suppressWarnings(nmfkc.ar.latent(ftot))
  itot <- quietly(nmfkc.ar.latent.inference(ltot, a$Y, a$A, wild.B = 5),
                  "does not fix the per-column scale")
  expect_false(itot$identified.approx)
  ## X.restriction = "fixed" does not estimate X at all, so T = I
  X0 <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))$X
  ffix <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE,
                                 X.init = X0, X.restriction = "fixed"))
  lfix <- suppressWarnings(nmfkc.ar.latent(ffix))
  ifix <- suppressWarnings(nmfkc.ar.latent.inference(lfix, a$Y, a$A, wild.B = 5,
                                                     refit.args = list(X.init = X0)))
  expect_true(ifix$identified.approx)
})

test_that("the identification diagnostic reports how far it is from exact", {
  X <- diag(2); Th <- matrix(c(0.5, 0, 0, 0.5), 2, 2)
  id <- nmfkc:::.nmfvar.identified(X, Th)
  expect_true(id$ok)
  expect_equal(id$anchor.margin, 0)          # exact anchors
  expect_equal(id$pure.margin, 0)            # exact zeros
  ## leak a little and the margins pick it up
  X2 <- matrix(c(1, 1e-4, 1e-4, 1), 2, 2)
  Th2 <- matrix(c(0.5, 1e-5, 1e-5, 0.5), 2, 2)
  id2 <- nmfkc:::.nmfvar.identified(X2, Th2)
  expect_true(id2$ok)
  expect_gt(id2$anchor.margin, 0)
  expect_gt(id2$pure.margin, 0)
})

test_that("a fixed basis is carried into the re-fits, not re-initialised", {
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  X0 <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))$X
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE,
                              X.init = X0, X.restriction = "fixed"))
  ## with "fixed" the returned X *is* the matrix that was held
  expect_equal(f$X, X0, tolerance = 0)
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  expect_equal(lat$fit.args$X.init, X0, tolerance = 0)
  inf <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 5))
  expect_equal(inf$refit.args$X.init, X0, tolerance = 0)
  ## without it the re-fits would hold a completely different basis
  g <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE,
                              X.restriction = "fixed"))
  expect_gt(max(abs(g$X - X0)), 1e-3)
  ## and no warning about an unrecoverable argument, since X.init is inherited
  expect_silent(suppressWarnings(
    nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 5)))
})

test_that("a legitimate rho = 0 replicate is not counted as a failure", {
  ## rho = 0 is a valid, stationary value: it happens when Theta collapses.
  Xz  <- matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2)
  latz <- suppressWarnings(nmfkc.ar.latent(
    structure(list(X = Xz, C = matrix(0, 2, 2)), class = "nmfkc")))
  expect_equal(latz$spectral.radius, 0)
  expect_true(latz$stationary)

  ## every usable replicate must be counted, so B + n.fail = B requested and the
  ## count comes from the failure flag rather than from rho > 0
  set.seed(3)
  Y <- matrix(abs(rnorm(4 * 30)) + 1, 4, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  inf <- suppressWarnings(nmfkc.ar.latent.inference(
    suppressWarnings(nmfkc.ar.latent(f)), a$Y, a$A, wild.B = 12))
  expect_equal(inf$wild.B + inf$n.fail, inf$wild.B.requested)
  expect_equal(inf$wild.B.requested, 12L)
})

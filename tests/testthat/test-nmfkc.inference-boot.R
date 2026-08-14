## Regression tests for two defects found by review (2026-07-26):
##  B1  method = "refit" reported p = 0 for every coefficient, because the
##      bootstrap p-value compared the RAW replicates with zero.  Under the
##      non-negativity constraint every replicate is >= 0, so mean(v < 0) was
##      always 0 and the two-sided formula collapsed to 0.
##  B2  nmfkc() seeds itself (default seed = 123) and left .Random.seed at that
##      state, so a user loop drawing random numbers after a fit silently
##      repeated itself and an outer set.seed() had no effect.

## A design with a hard zero pattern in C, so the fit has BOTH interior
## coefficients and boundary (essentially zero) ones -- the mix the p-value
## must be able to tell apart.  A fixture whose coefficients are all large
## would legitimately give p = 0 everywhere and could not detect the defect.
## A design with a hard zero pattern in C, fitted with the basis held fixed so
## that pattern stays identifiable.  The fit then has BOTH interior
## coefficients and near-zero ones -- the mix the p-value must tell apart.  A
## fixture whose coefficients are all large would legitimately give p = 0
## everywhere and could not detect the defect.
##
## The fit uses a tight tolerance on purpose: at nmfkc()'s default 1e-4 the two
## structural zeros stall around 0.05 instead of descending to ~0.002, i.e. the
## fixture would not contain a boundary coefficient at all.
make_var_fit <- function(seed = 11, P = 8, N = 80) {
  set.seed(seed)
  X <- matrix(abs(rnorm(P * 2)) + 0.2, P, 2)
  C <- rbind(c(0.5, 3, 0),                       # basis 1 responds to a1 only
             c(0.5, 0, 3))                       # basis 2 responds to a2 only
  A <- rbind(1, runif(N), runif(N))
  Y <- X %*% C %*% A + matrix(abs(rnorm(P * N)) * 0.05, P, N)
  list(Y = Y, A = A,
       fit = suppressWarnings(nmfkc(Y, A = A, rank = 2, X.init = X,
                                    X.restriction = "fixed", epsilon = 1e-8,
                                    maxit = 100000, verbose = FALSE)))
}

test_that("refit bootstrap p-values are not all zero and separate signal from noise", {
  d <- make_var_fit()
  inf <- suppressWarnings(nmfkc.inference(d$fit, d$Y, d$A, method = "refit",
                                          wild.unit = "column", wild.B = 100))
  p <- inf$coefficients$p_value
  est <- inf$coefficients$Estimate
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
  ## the defect: every coefficient got exactly the same (zero) p-value
  expect_gt(length(unique(p)), 1L)
  expect_false(all(p == 0))
  ## the two smallest coefficients (the planted zeros) must be less significant
  ## than the two largest ones
  ord <- order(est)
  expect_gt(mean(p[utils::head(ord, 2)]), mean(p[utils::tail(ord, 2)]))
})

test_that("the bootstrap p-value inverts the CENTRED replicate distribution", {
  ## direct check of the helper: replicates centred on est, so a coefficient
  ## whose replicates sit well below 2*est gets a small one-sided p, and one
  ## whose estimate is at the boundary gets a large (conservative) p.
  set.seed(1)
  draws <- rbind(rnorm(500, mean = 1, sd = 0.1),   # est = 1, clearly positive
                 abs(rnorm(500, mean = 0, sd = 0.1)))  # est ~ 0, at the boundary
  est <- c(1, 0)
  bs <- nmfkc:::.boot.summarize(draws, est = est, level = 0.95, p.side = "one.sided")
  expect_lt(bs$p.boot[1], 0.05)
  expect_gt(bs$p.boot[2], 0.5)
  expect_length(bs$se, 2L)
  expect_true(all(bs$ci.lower <= bs$ci.upper))
})

test_that("nmfkc() does not leave the caller's RNG stream at a fixed state", {
  d <- make_var_fit()
  draw_after_fit <- function() {
    invisible(suppressWarnings(nmfkc(d$Y, rank = 2, epsilon = 1e-4,
                                     maxit = 50, verbose = FALSE)))
    runif(1)
  }
  set.seed(99); got <- c(draw_after_fit(), draw_after_fit(), draw_after_fit())
  ## the defect: all three draws were identical
  expect_equal(length(unique(got)), 3L)
  ## and the stream must be exactly the one the caller would have had
  set.seed(99); want <- c(runif(1), runif(1), runif(1))
  expect_equal(got, want)
})

test_that("the refit bootstrap converges tightly enough for boundary coefficients", {
  ## Under multiplicative updates a coefficient heading for the non-negativity
  ## boundary approaches 0 slowly.  At nmfkc()'s own default tolerance (1e-4)
  ## every replicate stops while it is still spuriously positive, the
  ## replicates pile up above the point estimate, and the percentile interval
  ## can exclude the estimate.  The refit now defaults to 1e-8.
  d <- make_var_fit()
  tight <- suppressWarnings(nmfkc.inference(d$fit, d$Y, d$A, method = "refit",
                                            wild.unit = "column", wild.B = 60))
  cf <- tight$coefficients
  expect_true(all(cf$Estimate >= cf$CI_low - 1e-8))
  expect_true(all(cf$Estimate <= cf$CI_high + 1e-8))
  ## the tolerance is reachable from the user side
  loose <- suppressWarnings(nmfkc.inference(d$fit, d$Y, d$A, method = "refit",
                                            wild.unit = "column", wild.B = 60,
                                            refit.epsilon = 1e-4))
  ## Compare BSE, not SE: SE is now the sandwich standard error, which does not
  ## depend on how tightly the re-fits converged.  BSE is the bootstrap one.
  expect_false(isTRUE(all.equal(tight$coefficients$BSE, loose$coefficients$BSE)))
  ## the sandwich SE is indeed unaffected by refit.epsilon
  expect_equal(tight$coefficients$SE, loose$coefficients$SE, tolerance = 1e-12)
})

test_that("no inference or CV wrapper resets the caller's RNG stream", {
  ## Consuming randomness is fine; RESETTING the stream to a fixed state is the
  ## bug -- it makes every call leave .Random.seed in the same place, so a user
  ## loop repeats itself.  Detect it by checking that what came before the call
  ## still matters.
  set.seed(11); P <- 8; N <- 30
  Y  <- matrix(abs(rnorm(P * 2)), P, 2) %*% matrix(abs(rnorm(2 * N)), 2, N) +
        matrix(abs(rnorm(P * N)) * 0.1, P, N)
  A  <- rbind(1, runif(N))
  Y1 <- matrix(abs(rnorm(6 * N)), 6, N); Y2 <- matrix(abs(rnorm(5 * N)), 5, N)
  qq <- function(e) suppressWarnings(suppressMessages(e))
  fk  <- qq(nmfkc(Y, A = A, rank = 2, verbose = FALSE))
  ff  <- qq(nmf.ffb(Y1, Y2, rank = 2, verbose = FALSE))
  calls <- list(
    nmfkc.inference   = function() qq(nmfkc.inference(fk, Y, A, wild.B = 20)),
    nmf.ffb.inference = function() qq(nmf.ffb.inference(ff, Y1, Y2, B = 10,
                                                        ncores = 1, print.trace = FALSE)),
    nmfkc.cv          = function() qq(nmfkc.cv(Y, rank = 2, verbose = FALSE)),
    nmfkc.ecv         = function() qq(nmfkc.ecv(Y, rank = 2, verbose = FALSE)),
    nmf.ffb.cv        = function() qq(nmf.ffb.cv(Y1, Y2, rank = 2, seed = 7,
                                                 verbose = FALSE))
  )
  for (nm in names(calls)) {
    f <- calls[[nm]]
    set.seed(1); invisible(f()); a <- runif(1)
    set.seed(2); invisible(f()); b <- runif(1)
    expect_false(isTRUE(all.equal(a, b)), info = nm)
  }
})

test_that("nmf.ffb() with seed = NULL does not re-initialize the generator", {
  ## set.seed(NULL) re-seeds from the clock rather than doing nothing, which
  ## would make a seed = NULL fit irreproducible; nmf.ffb.cv() passes NULL to
  ## mean "use the stream as it is".
  set.seed(21); N <- 30
  Y1 <- matrix(abs(rnorm(6 * N)), 6, N); Y2 <- matrix(abs(rnorm(5 * N)), 5, N)
  qq <- function(e) suppressWarnings(suppressMessages(e))
  set.seed(5); a <- qq(nmf.ffb(Y1, Y2, rank = 2, seed = NULL, verbose = FALSE))
  set.seed(5); b <- qq(nmf.ffb(Y1, Y2, rank = 2, seed = NULL, verbose = FALSE))
  expect_equal(a$X, b$X)          # same ambient stream => same fit
})

test_that("no fitter leaves the caller's RNG stream at a fixed state", {
  set.seed(11); P <- 8; N <- 40
  Y  <- matrix(abs(rnorm(P * 2)), P, 2) %*% matrix(abs(rnorm(2 * N)), 2, N) +
        matrix(abs(rnorm(P * N)) * 0.1, P, N)
  A  <- rbind(1, runif(N))
  Ys <- matrix(rnorm(P * N), P, N); As <- matrix(rnorm(3 * N), 3, N)
  Y1 <- matrix(abs(rnorm(6 * N)), 6, N); Y2 <- matrix(abs(rnorm(5 * N)), 5, N)
  Yn <- crossprod(matrix(abs(rnorm(4 * 16)), 4, 16)); diag(Yn) <- 0
  qq <- function(e) suppressWarnings(suppressMessages(e))
  fitters <- list(
    nmfkc        = function() qq(nmfkc(Y, rank = 2, verbose = FALSE)),
    nmfre        = function() qq(nmfre(Y, A = A, rank = 2, verbose = FALSE)),
    nmf.ffb      = function() qq(nmf.ffb(Y1, Y2, rank = 2, verbose = FALSE)),
    nmfkc.net    = function() qq(nmfkc.net(Yn, rank = 2, verbose = FALSE)),
    nmfkc.signed = function() qq(nmfkc.signed(Ys, As, rank = 2, maxit = 50,
                                              verbose = FALSE)),
    nmfkc.ard    = function() qq(nmfkc.ard(Y, rank = 3, nrun = 2, plot = FALSE))
  )
  for (nm in names(fitters)) {
    f <- fitters[[nm]]
    set.seed(42); invisible(f()); after <- runif(1)
    set.seed(42); clean <- runif(1)
    expect_equal(after, clean, info = nm)
  }
})

test_that("seeding behaviour of the fit itself is unchanged", {
  d <- make_var_fit()
  ## a seeded fit stays reproducible ...
  a <- suppressWarnings(nmfkc(d$Y, rank = 2, seed = 7, verbose = FALSE))
  b <- suppressWarnings(nmfkc(d$Y, rank = 2, seed = 7, verbose = FALSE))
  expect_identical(a$X, b$X); expect_identical(a$C, b$C)
  ## ... and different seeds really do give different fits (they could not,
  ## while the internal default seed overrode everything)
  cc <- suppressWarnings(nmfkc(d$Y, rank = 2, seed = 7, X.init = "runif",
                               verbose = FALSE))
  dd <- suppressWarnings(nmfkc(d$Y, rank = 2, seed = 4321, X.init = "runif",
                               verbose = FALSE))
  expect_false(isTRUE(all.equal(cc$X, dd$X)))
  ## seed = NULL is left completely alone: we neither set nor restore, so that
  ## path keeps whatever RNG behaviour the initializer had (the default kmeans
  ## initializer happens to draw nothing).
  set.seed(5); invisible(suppressWarnings(nmfkc(d$Y, rank = 2, seed = NULL,
                                                epsilon = 1e-4, maxit = 50,
                                                verbose = FALSE)))
  after <- runif(1)
  set.seed(5); before <- runif(1)
  expect_equal(after, before)
})

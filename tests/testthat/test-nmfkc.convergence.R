## Convergence diagnostics of nmfkc() and nmfkc.signed(): `epsilon.iter`,
## `objfunc.increases`, and the `converged` verdict built from them.
##
## `converged` used to be `iter < maxit` in nmfkc.signed(), which called a run
## converged whenever it stopped for any reason at all -- including the early
## break on a non-finite objective.  It is now the same test the loop applies,
## so these tests pin the two together.

make_case <- function(seed = 5, P = 10, N = 30, D = 3) {
  set.seed(seed)
  list(Y = matrix(abs(stats::rnorm(P * N)), P, N),
       A = matrix(abs(stats::rnorm(D * N)), D, N),   # non-negative: usable by nmfkc()
       S = matrix(stats::rnorm(D * N), D, N))        # signed: for nmfkc.signed()
}

test_that("nmfkc() reports the quantity its own stopping rule compared", {
  skip_unless_full()
  d <- make_case()
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-6, maxit = 20000, verbose = FALSE)
  expect_true(all(c("epsilon.iter", "objfunc.increases", "converged", "iter",
                    "maxit", "epsilon") %in% names(f)))
  ## the denominator is the CURRENT objective, floored at 1 (see ?nmfkc); this
  ## floor is what makes the test absolute once the objective drops below 1
  tr <- f$objfunc.iter
  n <- length(tr)
  expect_gt(n, 1)
  expect_equal(f$epsilon.iter,
               abs(tr[n] - tr[n - 1]) / max(abs(tr[n]), 1),
               tolerance = 1e-12)
  expect_true(f$converged)
  expect_lte(f$epsilon.iter, f$epsilon)
})

test_that("nmfkc.signed() reports the quantity its own stopping rule compared", {
  skip_unless_full()
  d <- make_case()
  f <- nmfkc.signed(d$Y, d$S, rank = 2, epsilon = 1e-6, maxit = 20000,
                    verbose = FALSE, seed = 1)
  expect_true(all(c("epsilon.iter", "objfunc.increases", "converged") %in% names(f)))
  ## here the denominator is the PREVIOUS objective, floored at 1e-12 -- a
  ## different quantity from nmfkc()'s, deliberately: the rule stays relative
  ## all the way down.  The two epsilon.iter are comparable only above 1.
  tr <- f$objfunc.iter
  n <- length(tr)
  expect_identical(n, as.integer(f$iter))   # this trace is NOT trimmed
  expect_equal(f$epsilon.iter,
               abs(tr[n] - tr[n - 1]) / max(abs(tr[n - 1]), 1e-12),
               tolerance = 1e-12)
  expect_true(f$converged)
  expect_lt(f$epsilon.iter, f$epsilon)
})

test_that("a run stopped by maxit is not reported as converged", {
  skip_unless_full()
  d <- make_case()
  ## nmfkc.signed(): this is the case the old `iter < maxit` test got wrong
  s <- suppressWarnings(nmfkc.signed(d$Y, d$S, rank = 2, epsilon = 1e-14,
                                     maxit = 6, verbose = FALSE, seed = 1))
  expect_identical(as.integer(s$iter), 6L)
  expect_false(s$converged)
  expect_gt(s$epsilon.iter, s$epsilon)
  expect_match(paste(utils::capture.output(print(s)), collapse = " "),
               "NOT converged")
  ## nmfkc() behaved correctly already; keep it pinned beside it
  g <- suppressWarnings(nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-14, maxit = 12,
                              verbose = FALSE))
  expect_false(g$converged)
  expect_gt(g$epsilon.iter, abs(g$epsilon))
})

test_that("objfunc.increases counts the rises in the recorded trace", {
  skip_unless_full()
  d <- make_case()
  for (f in list(nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-6, maxit = 20000,
                       verbose = FALSE),
                 nmfkc.signed(d$Y, d$S, rank = 2, epsilon = 1e-6, maxit = 20000,
                              verbose = FALSE, seed = 1))) {
    expect_identical(as.integer(f$objfunc.increases),
                     as.integer(sum(diff(f$objfunc.iter) > 0)))
  }
})

test_that("the gauge-fixing paths stay monotone, so the count is zero", {
  skip_unless_full()
  ## A watchdog, not a property of the data.  A column restriction is a gauge
  ## fix -- the scale it removes from X is handed to C, so the model is
  ## unchanged and the multiplicative form is intact -- and the weighted sweep
  ## has its own backtracking guard.  Any non-zero count here means a change
  ## has broken that.  The one known offender of this kind, the X.rowSums.min
  ## floor, was removed in 0.9.8 because it oscillated to maxit on Covertype.
  d <- make_case()
  W <- matrix(stats::runif(length(d$Y), 0.2, 2), nrow(d$Y), ncol(d$Y))
  expect_identical(as.integer(nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-7,
                                    maxit = 20000, verbose = FALSE)$objfunc.increases), 0L)
  for (args in list(list(), list(update.power = 0.5), list(Y.weights = W))) {
    f <- suppressWarnings(do.call(nmfkc.signed,
      c(list(d$Y, d$S, rank = 2, epsilon = 1e-7, maxit = 5000,
             verbose = FALSE, seed = 1), args)))
    expect_identical(as.integer(f$objfunc.increases), 0L)
  }
})

test_that("X.restriction = 'rowSums' is refused, and says why", {
  skip_unless_full()
  ## Removed in 0.9.8.  It scaled the ROWS of X, which changes X %*% C %*% A,
  ## so it was a restriction acting rather than a reparametrization and did not
  ## preserve monotonicity: the diagnostic above counted 45 rises on this very
  ## data, alternating every second step -- an oscillation, the same defect
  ## that removed X.rowSums.min.  It is refused rather than silently re-mapped
  ## to "colSums", which would change results without saying so.
  d <- make_case()
  expect_error(nmfkc.signed(d$Y, d$S, rank = 2, X.restriction = "rowSums",
                            verbose = FALSE),
               "removed in 0.9.8")
  ## the gauge fixes still on the menu are all accepted
  for (r in c("colSums", "colSqSums", "totalSum", "none")) {
    f <- suppressWarnings(nmfkc.signed(d$Y, d$S, rank = 2, X.restriction = r,
                                       epsilon = 1e-6, maxit = 3000,
                                       verbose = FALSE, seed = 1))
    expect_identical(f$X.restriction, r)
    expect_identical(as.integer(f$objfunc.increases), 0L)
  }
})

test_that("print() and summary() of nmfkc.signed give the same convergence line", {
  skip_unless_full()
  ## the same contract test-nmfkc.inference.R applies to nmfkc(): the shared
  ## .print.convergence() helper must be fed the same fields by both paths
  d <- make_case()
  f <- nmfkc.signed(d$Y, d$S, rank = 2, epsilon = 1e-6, maxit = 20000,
                    verbose = FALSE, seed = 1)
  grab <- function(x, tag) {
    ln <- grep(tag, utils::capture.output(print(x)), value = TRUE)
    expect_length(ln, 1)
    trimws(sub(paste0("^\\s*", tag), "", ln))
  }
  expect_identical(grab(f, "Convergence: "),
                   grab(summary(f), "Iterations: {2,}"))
})


test_that("SE is the sandwich SE and BSE the bootstrap one, not the same number", {
  set.seed(4)
  Y <- matrix(abs(rnorm(5 * 40)) + 1, 5, 40)
  A <- rbind(1, abs(rnorm(40)))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE))
  inf <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit",
                                          wild.B = 40, wild.seed = 1))
  ct <- inf$coefficients
  expect_true("on.boundary" %in% names(ct))
  expect_false(isTRUE(all.equal(ct$SE, ct$BSE)))
  expect_equal(ct$SE, as.vector(inf$C.se.sandwich), tolerance = 1e-12)
  expect_equal(ct$BSE, as.vector(inf$C.se.boot), tolerance = 1e-12)
  ## z still uses the primary SE, which is the bootstrap one under "refit"
  expect_equal(ct$z_value, ifelse(ct$BSE > 0, ct$Estimate / ct$BSE, NA_real_),
               tolerance = 1e-12)
  ## and the boundary flag matches the estimates
  expect_identical(ct$on.boundary, ct$Estimate <= 1e-3)
})

test_that("boundary.tol controls the on.boundary flag", {
  set.seed(4)
  Y <- matrix(abs(rnorm(5 * 40)) + 1, 5, 40)
  A <- rbind(1, abs(rnorm(40)))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE))
  inf <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.B = 20,
                                          wild.seed = 1, boundary.tol = 1))
  expect_true(all(inf$coefficients$on.boundary[inf$coefficients$Estimate <= 1]))
})

test_that("print.nmfkc is an lm-style one-screen display, not a list dump", {
  Y <- matrix(cars$dist, nrow = 1); A <- rbind(1, cars$speed)
  f <- suppressWarnings(nmfkc(Y, A, rank = 1, verbose = FALSE))
  out <- capture.output(print(f))
  expect_lt(length(out), 15)                     # unclassed it is 350+ lines
  expect_true(any(grepl("^Call:", out)))
  expect_true(any(grepl("^Model:", out)))
  expect_true(any(grepl("^Convergence:", out)))
  ## C can be Q x N, so it is deliberately NOT printed
  expect_false(any(grepl("Coefficients", out)))
  expect_true(any(grepl("coef()", out, fixed = TRUE)))
  expect_identical(withVisible(print(f))$visible, FALSE)   # returns invisibly
})

test_that("print.nmfkc says whether the fit converged", {
  set.seed(2)
  Y <- matrix(abs(rnorm(5 * 60)) + 1, 5, 60); A <- rbind(1, abs(rnorm(60)))
  ok <- suppressWarnings(nmfkc(Y, A, rank = 2, verbose = FALSE))
  expect_true(ok$converged)
  ## print() and summary() go through the same .print.convergence() helper, so
  ## they must render the fact identically -- "<iter> / <maxit> (verdict)".
  pr <- grep("^Convergence:", capture.output(print(ok)), value = TRUE)
  sm <- grep("^Iterations:", capture.output(print(summary(ok))), value = TRUE)
  expect_match(pr, "68 / 5000 \\(converged\\)", fixed = FALSE)
  expect_identical(sub("^Convergence: ", "", pr), sub("^Iterations: ", "", sm))
  expect_false(any(grepl("NOT converged", capture.output(print(ok)))))
  ## a run capped by maxit used to be indistinguishable from a converged one
  cap <- suppressWarnings(nmfkc(Y, A, rank = 2, maxit = 3, verbose = FALSE))
  expect_false(cap$converged)
  expect_equal(cap$maxit, 3)
  expect_true(any(grepl("3 / 3 \\(NOT converged\\)", capture.output(print(cap)))))
})

test_that("fitted() works even when XB was not stored", {
  set.seed(2)
  Y <- matrix(abs(rnorm(5 * 40)) + 1, 5, 40); A <- rbind(1, abs(rnorm(40)))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, verbose = FALSE))
  m <- suppressWarnings(nmfkc(Y, A, rank = 2, verbose = FALSE, detail = "minimal"))
  expect_true(length(m$XB) == 1L && is.na(m$XB))   # not stored
  expect_equal(fitted(m), f$X %*% f$B, tolerance = 1e-12)
  expect_equal(fitted(f), f$XB, tolerance = 0)
})

test_that("the unweighted short-circuit is bit-identical", {
  ## Y.weights is all ones unless ECV / BiCV / an NA mask / explicit weights
  ## put zeros in it, and multiplying by exactly 1 is the identity -- so the
  ## skip must not move a single bit on either path.
  set.seed(9)
  Y <- matrix(abs(rnorm(8 * 80)) + 1, 8, 80); A <- rbind(1, abs(rnorm(80)))
  w1 <- matrix(1, 8, 80)                     # explicitly all ones
  a <- suppressWarnings(nmfkc(Y, A, rank = 3, maxit = 300, verbose = FALSE))
  b <- suppressWarnings(nmfkc(Y, A, rank = 3, maxit = 300, verbose = FALSE,
                              Y.weights = w1))
  expect_identical(a$X, b$X)
  expect_identical(a$C, b$C)
  expect_identical(a$objfunc, b$objfunc)
  ## and a genuinely weighted fit still takes the weighted path
  w2 <- w1; w2[1, 1:10] <- 0
  d <- suppressWarnings(nmfkc(Y, A, rank = 3, maxit = 300, verbose = FALSE,
                              Y.weights = w2))
  expect_false(isTRUE(all.equal(a$objfunc, d$objfunc)))
})

test_that("detail defaults to fast and drops only the O(N^2) criteria", {
  set.seed(2)
  Y <- matrix(abs(rnorm(5 * 60)) + 1, 5, 60); A <- rbind(1, abs(rnorm(60)))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, verbose = FALSE))
  g <- suppressWarnings(nmfkc(Y, A, rank = 2, verbose = FALSE, detail = "full"))
  ## default no longer carries them at all (absent, not NA)
  expect_null(f$criterion$silhouette)
  expect_null(f$criterion$CPCC)
  expect_null(f$criterion$dist.cor)
  expect_false(is.null(g$criterion$silhouette))
  ## B.prob.max.mean is gone in both: nothing consumed it
  expect_null(f$criterion$B.prob.max.mean)
  expect_null(g$criterion$B.prob.max.mean)
  ## and the fit itself is untouched by the switch
  expect_equal(f$X, g$X, tolerance = 0)
  expect_equal(f$C, g$C, tolerance = 0)
  expect_equal(f$r.squared, g$r.squared, tolerance = 0)
})

test_that("effective.rank.index is the [0,1] broken-stick index nmfkc.rank plots", {
  set.seed(2)
  Y <- matrix(abs(rnorm(5 * 60)) + 1, 5, 60); A <- rbind(1, abs(rnorm(60)))
  fits <- lapply(2:3, function(q)
    suppressWarnings(nmfkc(Y, A, rank = q, verbose = FALSE)))
  rk <- suppressWarnings(nmfkc.rank(Y, A, rank = 2:3, detail = "fast", plot = FALSE))
  ## one helper, so the two must agree exactly
  expect_equal(as.numeric(rk$criteria$effective.rank.index),
               vapply(fits, function(f) f$criterion$effective.rank.index,
                      numeric(1)), tolerance = 0)
  ## it is in [0, 1], unlike the raw effective.rank / Q
  idx <- vapply(fits, function(f) f$criterion$effective.rank.index, numeric(1))
  expect_true(all(idx >= 0 & idx <= 1))
  ## summary prints it in place of the old crispness
  out <- capture.output(print(summary(fits[[1]])))
  expect_true(any(grepl("Factor variance share", out)))
  expect_false(any(grepl("Clustering Crispness", out)))
})

test_that("the broken-stick index measures variance share, not usefulness", {
  eri <- nmfkc:::.effective.rank.index
  er  <- nmfkc:::.effective.rank
  set.seed(1); N <- 200; b <- runif(N)
  ## a large but CONSTANT factor carries no variance: the raw index goes
  ## negative (below the random null) and is clamped at 0
  B.const <- rbind(rep(10, N), runif(N), runif(N))
  expect_equal(eri(er(B.const), 3), 0)
  ## two DUPLICATED factors split the variance evenly and score near 1,
  ## which is why the label says "variance share", not "all factors useful"
  expect_gt(eri(er(rbind(b, b, runif(N))), 3), 0.9)
  ## undefined at Q = 1, where the broken-stick expectation equals Q
  expect_true(is.na(eri(er(matrix(runif(N), 1, N)), 1)))
})

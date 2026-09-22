## predict.nmfkc.signed(): probabilities from signed least-squares scores.
## With signed covariates b_n = C a_n can be negative, so type = "prob" maps
## yhat_n = X b_n to the probability simplex by Euclidean projection (Duchi
## et al. 2008); the old clip-and-renormalize rule stays as prob.method =
## "clip"; type = "class" is the argmax of yhat and is the same under both.

test_that(".proj_simplex_cols() is the Euclidean projection onto the simplex", {
  skip_unless_full()
  set.seed(3)
  Y <- matrix(stats::rnorm(6 * 200, sd = 0.6), 6, 200)
  Y[, 1] <- c(0.2, 0.3, 0.5, 0, 0, 0)          # already on the simplex
  Y[, 2] <- c(2, -1, 0, 0, 0, 0)               # one dominant coordinate
  Y[, 3] <- rep(1 / 6, 6)                      # uniform
  P <- nmfkc:::.proj_simplex_cols(Y)
  expect_equal(dim(P), dim(Y))
  expect_true(all(P >= 0))
  expect_equal(colSums(P), rep(1, ncol(Y)), tolerance = 1e-12)
  expect_equal(P[, 1], Y[, 1], tolerance = 1e-12)      # fixed point on the simplex
  expect_equal(P[, 2], c(1, 0, 0, 0, 0, 0))
  expect_equal(P[, 3], Y[, 3], tolerance = 1e-12)
  ## KKT characterization: on the support y - p is a constant tau; off the
  ## support y <= tau.
  for (n in 4:ncol(Y)) {
    y <- Y[, n]; p <- P[, n]; on <- p > 0
    tau <- unique(round(y[on] - p[on], 10))
    expect_length(tau, 1L)
    expect_true(all(y[!on] <= tau + 1e-10))
  }
  ## Optimality against random feasible points: no p' on the simplex is closer.
  for (n in 4:8) {
    y <- Y[, n]; p <- P[, n]
    for (k in 1:20) {
      q <- stats::runif(6); q <- q / sum(q)
      expect_lte(sum((p - y)^2), sum((q - y)^2) + 1e-12)
    }
  }
  ## Argmax is preserved
  expect_equal(apply(P, 2, which.max), apply(Y, 2, which.max))
  ## P = 1 degenerate case
  expect_equal(as.vector(nmfkc:::.proj_simplex_cols(matrix(c(-2, 5), 1, 2))), c(1, 1))
})

test_that("predict.nmfkc.signed(type = 'prob') projects onto the simplex; 'clip' is kept", {
  skip_unless_full()
  set.seed(7)
  p <- 5; N <- 90; P <- 4; D <- 12
  U <- matrix(stats::rnorm(p * N), p, N)
  lab <- sample.int(P, N, replace = TRUE)
  Y <- matrix(0, P, N); Y[cbind(lab, seq_len(N))] <- 1
  rownames(Y) <- paste0("c", seq_len(P))
  rff <- nmfkc.signed.rff(U, beta = 0.3, D = D, seed = 3)
  fit <- suppressMessages(nmfkc.signed(Y, A = rff$Z, rank = P, maxit = 300,
                                       warm.start = FALSE, verbose = FALSE))
  Ynew <- nmfkc.signed.rff(U[, 1:30], pars = rff$pars)$Z
  yhat <- predict(fit, newA = Ynew)                         # signed response
  expect_true(any(yhat < 0))                                # the case that motivates this
  pr <- predict(fit, newA = Ynew, type = "prob")
  expect_equal(dim(pr), c(P, 30L))
  expect_identical(rownames(pr), rownames(Y))
  expect_true(all(pr >= 0)); expect_equal(colSums(pr), rep(1, 30), tolerance = 1e-12)
  expect_equal(pr, nmfkc:::.proj_simplex_cols(yhat))
  ## Legacy rule still available and different where yhat has negatives
  pc <- predict(fit, newA = Ynew, type = "prob", prob.method = "clip")
  expect_true(all(pc >= 0)); expect_equal(colSums(pc), rep(1, 30), tolerance = 1e-9)
  expect_false(isTRUE(all.equal(pr, pc)))
  ## Class = argmax of the response, identical under both probability rules
  cl <- predict(fit, newA = Ynew, type = "class")
  expect_identical(cl, rownames(Y)[apply(yhat, 2, which.max)])
  expect_identical(cl, rownames(Y)[apply(pr, 2, which.max)])
  expect_identical(cl, rownames(Y)[apply(pc, 2, which.max)])
  expect_error(predict(fit, newA = Ynew, type = "prob", prob.method = "softmax"))
})

## Gram-object input to nmfkc.signed(): the large-N route for Random
## Fourier Features.  nmfkc.signed.rff.gram() accumulates S = Z Z^T and
## G0 = Y Z^T over column blocks so the D x N feature matrix never exists;
## nmfkc.signed() accepts the object in place of A.  These tests pin down
## (i) that the accumulated S, G0 equal the one-shot products, (ii) that the
## fit equals the explicit-matrix fit with warm.start = FALSE, (iii) that
## the documented restrictions (warm.start, Y.weights, NA, fold helpers)
## are enforced, and (iv) that a matrix A still takes the old path exactly.

make_gram_case <- function(seed = 7, p = 5, N = 90, Q_obs = 4, D = 12) {
  set.seed(seed)
  U <- matrix(stats::rnorm(p * N), p, N)
  lab <- sample.int(Q_obs, N, replace = TRUE)
  Y <- matrix(0, Q_obs, N); Y[cbind(lab, seq_len(N))] <- 1
  rownames(Y) <- paste0("c", seq_len(Q_obs))
  list(U = U, Y = Y, D = D, beta = 0.3)
}

test_that("block-wise S and G0 equal the one-shot products", {
  skip_unless_full()
  d <- make_gram_case()
  ## block.size = 17 with N = 90 -> 6 blocks, the last one ragged (5 columns)
  g <- nmfkc.signed.rff.gram(d$Y, d$U, beta = d$beta, D = d$D, seed = 3,
                             block.size = 17)
  expect_s3_class(g, "nmfkc.signed.gram")
  expect_equal(g$n.blocks, 6L)
  expect_equal(g$N, 90L); expect_equal(g$D, 12L)
  Z <- nmfkc.signed.rff(d$U, pars = g$pars)$Z
  expect_equal(dim(Z), c(12L, 90L))
  expect_equal(g$S,  tcrossprod(Z),      tolerance = 1e-12)
  expect_equal(g$G0, tcrossprod(d$Y, Z), tolerance = 1e-12)
  ## The block generator reproduces columns of Z
  expect_equal(g$A.block(c(3L, 50L, 90L)), Z[, c(3L, 50L, 90L)])
  expect_output(print(g), "Gram summary")
})

test_that("pars reuse gives the same Gram object without drawing new omega, b", {
  skip_unless_full()
  d <- make_gram_case()
  g1 <- nmfkc.signed.rff.gram(d$Y, d$U, beta = d$beta, D = d$D, seed = 3)
  g2 <- nmfkc.signed.rff.gram(d$Y, d$U, pars = g1$pars, block.size = 7)
  expect_identical(g1$pars, g2$pars)
  expect_equal(g1$S, g2$S, tolerance = 1e-12)
  expect_equal(g1$G0, g2$G0, tolerance = 1e-12)
  ## D is required when pars are not given
  expect_error(nmfkc.signed.rff.gram(d$Y, d$U, beta = d$beta), "'D'")
  expect_error(nmfkc.signed.rff.gram(d$Y, d$U, D = d$D), "'beta'")
})

test_that("Gram-object fit equals the explicit-matrix fit (warm.start = FALSE)", {
  skip_unless_full()
  d <- make_gram_case()
  g <- nmfkc.signed.rff.gram(d$Y, d$U, beta = d$beta, D = d$D, seed = 3,
                             block.size = 20)
  Z <- nmfkc.signed.rff(d$U, pars = g$pars)$Z
  expect_message(
    fg <- nmfkc.signed(d$Y, A = g, rank = 4, maxit = 300, seed = 11),
    "warm.start is not available")
  fm <- suppressMessages(
    nmfkc.signed(d$Y, A = Z, rank = 4, maxit = 300, seed = 11,
                 warm.start = FALSE))
  expect_s3_class(fg, "nmfkc.signed")
  expect_equal(fg$iter, fm$iter)
  expect_equal(fg$X,  fm$X,  tolerance = 1e-8)
  expect_equal(fg$C,  fm$C,  tolerance = 1e-8)
  ## `%*%` leaves an all-NULL dimnames list on B in the matrix path; the
  ## values are what matters here.
  expect_equal(unname(fg$B),  unname(fm$B),  tolerance = 1e-8)
  expect_equal(unname(fg$XB), unname(fm$XB), tolerance = 1e-8)
  expect_equal(fg$objfunc, fm$objfunc, tolerance = 1e-8)
  expect_equal(fg$objfunc.iter, fm$objfunc.iter, tolerance = 1e-8)
  expect_equal(fg$r.squared, fm$r.squared, tolerance = 1e-8)
  expect_equal(fg$sigma, fm$sigma, tolerance = 1e-8)
  expect_equal(fg$D, 12L)
  ## RFF pars are inherited from the Gram object for summary()/predict()
  expect_identical(fg$pars, g$pars)
  expect_output(print(summary(fg)), "beta")
  ## predict() is unchanged: features for new data come from the stored pars
  Znew <- nmfkc.signed.rff(d$U[, 1:10], pars = g$pars)$Z
  expect_equal(predict(fg, newA = Znew), fm$XB[, 1:10], tolerance = 1e-8)
  expect_length(predict(fg, newA = Znew, type = "class"), 10L)
})

test_that("Gram input without a block generator still returns the closed-form objective", {
  skip_unless_full()
  d <- make_gram_case()
  g <- nmfkc.signed.rff.gram(d$Y, d$U, beta = d$beta, D = d$D, seed = 3)
  g0 <- g; g0$A.block <- NULL
  fg  <- suppressMessages(nmfkc.signed(d$Y, A = g,  rank = 4, maxit = 200, seed = 11))
  fg0 <- suppressMessages(nmfkc.signed(d$Y, A = g0, rank = 4, maxit = 200, seed = 11))
  expect_null(fg0$B); expect_null(fg0$XB)
  expect_equal(fg0$C, fg$C, tolerance = 1e-10)
  expect_equal(fg0$objfunc, fg$objfunc, tolerance = 1e-8)
  expect_equal(fg0$r.squared.uncentered, fg$r.squared.uncentered, tolerance = 1e-8)
  expect_equal(fg0$sigma, fg$sigma, tolerance = 1e-8)
  expect_true(is.na(fg0$r.squared)); expect_true(is.na(fg0$mae))
})

test_that("Gram input enforces its documented restrictions", {
  skip_unless_full()
  d <- make_gram_case()
  g <- nmfkc.signed.rff.gram(d$Y, d$U, beta = d$beta, D = d$D, seed = 3)
  W <- matrix(1, nrow(d$Y), ncol(d$Y)); W[1, 1] <- 0
  expect_error(nmfkc.signed(d$Y, A = g, rank = 4, Y.weights = W, verbose = FALSE),
               "Y.weights is not supported")
  Yna <- d$Y; Yna[1, 1] <- NA
  expect_error(nmfkc.signed(Yna, A = g, rank = 4, verbose = FALSE), "NA")
  ## Built for a different N / a different Y
  expect_error(nmfkc.signed(d$Y[, 1:80], A = g, rank = 4, verbose = FALSE), "N = 90")
  expect_error(nmfkc.signed(d$Y[1:3, ], A = g, rank = 3, verbose = FALSE), "different Y")
  ## Fold-based helpers refuse a Gram object with a pointer to the alternative
  expect_error(nmfkc.signed.cv(d$Y, g, rank = 4), "cannot be split into folds")
  expect_error(nmfkc.signed.ecv(d$Y, g, rank = 2:3), "cannot be split into folds")
  expect_error(nmfkc.signed.rank(d$Y, g, rank = 2:3, plot = FALSE),
               "cannot be split into folds")
  ## An explicit warm.start = FALSE is silent (nothing was overridden)
  expect_silent(nmfkc.signed(d$Y, A = g, rank = 4, maxit = 50, verbose = FALSE))
})

test_that("a matrix A still takes the original path (regression guard)", {
  skip_unless_full()
  d <- make_gram_case()
  Z <- nmfkc.signed.rff(d$U, beta = d$beta, D = d$D, seed = 3)$Z
  rownames(Z) <- paste0("f", seq_len(nrow(Z)))
  fit <- suppressMessages(nmfkc.signed(d$Y, A = Z, rank = 4, maxit = 200))
  ## warm-start ran (default TRUE, Y >= 0) and the covariate names survive
  expect_identical(colnames(fit$C), rownames(Z))
  expect_equal(dim(fit$B), c(4L, 90L))
  expect_equal(fit$objfunc, sum((d$Y - fit$XB)^2), tolerance = 1e-10)
  expect_false(is.na(fit$mae))
})

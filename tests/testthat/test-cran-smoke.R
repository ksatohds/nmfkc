# The suite that runs everywhere, including on CRAN.  Every other test file is
# behind skip_unless_full() (see helper-nmfkc.R): the regression suite fits the
# model thousands of times across bootstrap replicates, CV folds and restarts,
# which does not fit CRAN's 10-minute budget for the whole check.
#
# What belongs here: does each exported fitter run on toy data, return its
# documented class and shapes, and do its S3 verbs work?  That is what a check
# on a machine that is not the maintainer's can usefully establish.  Numerical
# regressions -- the ones about p-values, RNG streams, convergence tolerances,
# identification -- live in the full suite.
#
# Keep it fast.  Every fit below is a handful of rows on a handful of columns.

set.seed(1)
P <- 6; N <- 20; Q <- 2
Y  <- matrix(abs(rnorm(P * N)) + 0.5, P, N)
A  <- rbind(1, runif(N))
Y1 <- matrix(abs(rnorm(4 * N)) + 0.5, 4, N)
Y2 <- matrix(abs(rnorm(3 * N)) + 0.5, 3, N)
SM <- 0.5 * (function(M) M + t(M))(matrix(abs(rnorm(6 * 6)), 6, 6)); diag(SM) <- 0

test_that("nmfkc() fits, and its S3 verbs work", {
  f <- nmfkc(Y, A = A, rank = Q, verbose = FALSE)
  expect_s3_class(f, "nmfkc")
  expect_equal(dim(f$X), c(P, Q))
  expect_equal(dim(fitted(f)), dim(Y))
  expect_equal(dim(residuals(f, Y)), dim(Y))
  expect_true(all(f$X >= 0) && all(f$B >= 0))
  expect_true(is.data.frame(coef(f)) || is.matrix(coef(f)))
  expect_output(print(f))
  expect_output(print(summary(f)))
})

test_that("nmfkc() reduces to plain NMF without covariates", {
  f <- nmfkc(Y, rank = Q, verbose = FALSE)
  expect_s3_class(f, "nmfkc")
  expect_equal(dim(f$X), c(P, Q))
  expect_true(all(is.finite(f$objfunc)))
  expect_lte(utils::tail(f$objfunc, 1), utils::head(f$objfunc, 1) + 1e-8)
})

test_that("the kernel path runs", {
  U <- matrix(rnorm(2 * N), 2, N)
  K <- nmfkc.kernel(U, beta = 0.5)
  expect_equal(dim(K), c(N, N))
  f <- nmfkc(Y, A = K, rank = Q, verbose = FALSE)
  expect_equal(dim(f$X), c(P, Q))
})

test_that("nmfkc.inference() returns a coefficient table", {
  f <- nmfkc(Y, A = A, rank = Q, verbose = FALSE)
  inf <- suppressWarnings(nmfkc.inference(f, Y, A, method = "onestep"))
  cf <- inf$coefficients
  expect_true(all(c("Basis", "Covariate", "Estimate", "SE") %in% names(cf)))
  expect_true(all(is.finite(cf$Estimate)))
})

test_that("nmfre() fits and infers", {
  f <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE)
  expect_s3_class(f, "nmfre")
  expect_equal(dim(f$X), c(P, Q))
  expect_equal(dim(fitted(f)), dim(Y))
})

test_that("nmf.rrr() returns both score matrices", {
  f <- suppressWarnings(nmf.rrr(Y1, Y2, rank1 = Q, maxit = 200, verbose = FALSE))
  expect_s3_class(f, "nmfae")
  expect_equal(nrow(f$B1), Q)
  expect_true(!is.null(f$B2))
  expect_equal(dim(fitted(f)), dim(Y1))
})

test_that("nmf.ffb() fits its two blocks", {
  f <- suppressWarnings(nmf.ffb(Y1, Y2, rank = Q, verbose = FALSE))
  expect_equal(nrow(f$X), nrow(Y1))
  expect_output(print(f))
})

test_that("nmfkc.net() fits a symmetric matrix", {
  f <- suppressWarnings(nmfkc.net(SM, rank = Q, verbose = FALSE))
  expect_equal(nrow(f$X), nrow(SM))
  expect_true(all(f$X >= 0))
})

test_that("nmfkc.signed() accepts signed covariates", {
  As <- matrix(rnorm(2 * N), 2, N)
  f <- suppressWarnings(nmfkc.signed(matrix(rnorm(P * N), P, N), As,
                                     rank = Q, maxit = 50, verbose = FALSE))
  expect_equal(dim(f$X), c(P, Q))
  expect_true(all(is.finite(f$X)))
})

test_that("nmfkc.signed() accepts a block-wise RFF Gram object in place of A", {
  U <- matrix(rnorm(3 * N), 3, N)
  g <- nmfkc.signed.rff.gram(Y, U, beta = 0.5, D = 6, seed = 1, block.size = 8)
  expect_s3_class(g, "nmfkc.gram")
  Z <- nmfkc.signed.rff(U, pars = g$pars)$Z
  expect_equal(g$S, tcrossprod(Z), tolerance = 1e-12)
  fg <- suppressMessages(nmfkc.signed(Y, A = g, rank = Q, maxit = 30, verbose = FALSE))
  fm <- suppressMessages(nmfkc.signed(Y, A = Z, rank = Q, maxit = 30, verbose = FALSE,
                                      warm.start = FALSE))
  expect_equal(fg$C, fm$C, tolerance = 1e-8)
  expect_equal(dim(fg$B), c(Q, N))
})

test_that("nmfkc() accepts a block-wise Nystroem Gram object in place of A", {
  U <- matrix(rnorm(3 * N), 3, N)
  g <- nmfkc.kernel.gram(Y, U, V = U[, 1:4], beta = 0.5, block.size = 7)
  expect_s3_class(g, "nmfkc.gram")
  A <- nmfkc.kernel(U[, 1:4], U, beta = 0.5)
  expect_equal(g$S, unname(tcrossprod(unclass(A))), tolerance = 1e-12)
  fg <- suppressWarnings(nmfkc(Y, A = g, rank = Q, maxit = 30, verbose = FALSE))
  fm <- suppressWarnings(nmfkc(Y, A = A, rank = Q, maxit = 30, verbose = FALSE))
  expect_equal(fg$C, fm$C, tolerance = 1e-8)
  expect_equal(dim(fg$B), c(Q, N))
})

test_that("nmfkc.ar() builds a lagged design that nmfkc() can fit", {
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = Q, verbose = FALSE))
  lat <- suppressWarnings(nmfkc.ar.latent(f))
  expect_s3_class(lat, "nmfkc.ar.latent")
  expect_equal(dim(lat$G$lag1), c(Q, Q))
})

test_that("rank selection and cross-validation run", {
  expect_s3_class(suppressWarnings(nmfkc.cv(Y, A = A, rank = Q, div = 2)), "nmfkc.cv")
  ec <- suppressWarnings(nmfkc.ecv(Y, A = A, rank = Q, div = 2))
  expect_true(is.finite(ec$objfunc[1]) || is.finite(ec$sigma[1]))
})

test_that("the DOT writers produce graph source", {
  f <- nmfkc(Y, A = A, rank = Q, verbose = FALSE)
  d <- nmfkc.DOT(f)
  expect_type(as.character(d), "character")
  expect_match(paste(as.character(d), collapse = " "), "digraph|graph")
})

test_that("nmf.gmm and its two-stage baseline fit on toy data", {
  Afit <- rbind(1, rnorm(ncol(Y)))
  g <- suppressWarnings(nmf.gmm(Y, Afit, rank = Q, K = 2, nstart = 1,
                                maxit = 50, seed = 1))
  expect_s3_class(g, "nmf.gmm")
  expect_length(g$cluster, ncol(Y))
  expect_true(is.finite(g$BIC))
  expect_length(predict(g), ncol(Y))
  ts <- suppressWarnings(nmf.gmm.twostage(Y, Afit, rank = Q, K = 2, nstart = 1,
                                          maxit = 50, seed = 1))
  expect_s3_class(ts, "nmf.gmm.twostage")
  expect_gte(ts$twostage$shift, 0)
  df <- data.frame(a = Afit[2, ])
  gf <- suppressWarnings(nmf.gmm(Y, ~ a, rank = Q, K = 2, nstart = 1,
                                 maxit = 50, seed = 1, data = df))
  expect_equal(nrow(gf$A), 2L)
})

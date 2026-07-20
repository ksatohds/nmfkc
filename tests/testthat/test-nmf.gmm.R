## small synthetic NMF-GMM data: K=2 well-separated classes on the scores
make_gmm_data <- function(seed = 1, P = 10, N = 80, Q = 2) {
  set.seed(seed)
  X <- matrix(abs(rnorm(P * Q)), P, Q); X <- sweep(X, 2, colSums(X), "/")
  A <- rbind(1, rnorm(N))                                   # intercept + covariate
  C <- matrix(c(2, 1.5, 0.8, -0.6), Q, 2)
  z <- rep(1:2, length.out = N)
  muk <- cbind(c(3, -3), c(-3, 3))                          # two separated class means
  B <- C %*% A + muk[, z] + matrix(rnorm(Q * N) * 0.3, Q, N)
  Y <- X %*% B + matrix(rnorm(P * N) * 0.2, P, N)
  list(Y = Y, A = A, X = X, z = z, Q = Q)
}

test_that("nmf.gmm fits, uses house-style returns, and recovers clusters", {
  d <- make_gmm_data()
  fit <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, X.init = d$X, seed = 1)
  expect_s3_class(fit, "nmf.gmm")
  ## house-style fields; parameter matrix is C, not Theta
  expect_true(all(c("X", "C", "mu", "tau2", "sigma2", "xi", "gamma", "cluster",
                    "loglik", "BIC", "ICL", "objfunc", "objfunc.iter", "iter",
                    "converged", "dims", "runtime") %in% names(fit)))
  expect_null(fit$Theta)
  ## shapes & constraints
  expect_equal(dim(fit$C), c(d$Q, 2L))
  expect_equal(dim(fit$mu), c(d$Q, 2L))
  expect_equal(dim(fit$gamma), c(ncol(d$Y), 2L))
  expect_true(all(fit$X >= 0))
  expect_equal(unname(colSums(fit$X)), rep(1, d$Q), tolerance = 1e-6)  # column-normalized
  expect_equal(sum(fit$xi), 1, tolerance = 1e-8)
  expect_equal(rowSums(fit$gamma), rep(1, ncol(d$Y)), tolerance = 1e-8)
  ## clustering recovers the two planted classes
  expect_gt(nmfkc:::.nmfgmm.ARI(fit$cluster, d$z), 0.8)
  ## S3
  expect_equal(coef(fit), fit$C)
  expect_equal(dim(fitted(fit)), dim(d$Y))
  expect_length(predict(fit), ncol(d$Y))
  expect_equal(dim(predict(fit, type = "responsibility")), dim(fit$gamma))
  expect_s3_class(summary(fit), "summary.nmf.gmm")
})

test_that("K=1 nmf.gmm degenerates to a single-class (NMF-RE-type) fit", {
  d <- make_gmm_data()
  g1 <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 1, X.init = d$X, seed = 1)
  expect_equal(g1$K, 1L)
  expect_true(all(g1$cluster == 1L))
  expect_equal(max(abs(g1$mu)), 0, tolerance = 1e-6)   # sum-to-zero => mu==0 at K=1
  expect_equal(unname(g1$xi), 1)
  expect_equal(g1$BIC, g1$ICL, tolerance = 1e-8)       # zero entropy at K=1
  expect_true(is.finite(g1$loglik))
})

test_that("nmf.gmm.inference builds the C coefficients table (tied)", {
  d <- make_gmm_data()
  fit <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, X.init = d$X, seed = 1)
  inf <- nmf.gmm.inference(fit, d$Y, d$A, boot = 200, seed = 1)
  expect_s3_class(inf, "nmf.gmm.inference")
  expect_true(all(c("Basis", "Covariate", "Estimate", "SE", "BSE", "z_value",
                    "p_value", "CI_low", "CI_high") %in% names(inf$coefficients)))
  expect_equal(nrow(inf$coefficients), prod(dim(fit$C)))
  expect_true(all(is.finite(inf$coefficients$SE)))
  expect_true(all(is.finite(inf$coefficients$BSE)))
  expect_error(nmf.gmm.inference(nmf.gmm(d$Y, d$A, rank = d$Q, K = 2,
                                         X.init = d$X, cov = "free"), d$Y, d$A),
               "tied")
})

test_that("nmf.gmm.select sweeps K with BIC/ICL/ARI", {
  d <- make_gmm_data()
  sel <- nmf.gmm.select(d$Y, d$A, rank = d$Q, K = 1:3, X.init = d$X,
                        truth = d$z, verbose = FALSE)
  expect_s3_class(sel, "nmf.gmm.select")
  expect_true(all(c("K", "logLik", "BIC", "ICL", "ARI") %in% names(sel$table)))
  expect_equal(nrow(sel$table), 3L)
  expect_true(sel$K.best %in% 1:3)
  ## the planted K=2 should be recovered by BIC on clean synthetic data
  expect_equal(sel$K.best, 2L)
})

test_that("cov='free' (per-class variances) also runs", {
  d <- make_gmm_data()
  ff <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, cov = "free", X.init = d$X, seed = 1)
  expect_equal(dim(ff$tau2), c(d$Q, 2L))               # per-class diagonal
  expect_true(is.finite(ff$loglik))
})

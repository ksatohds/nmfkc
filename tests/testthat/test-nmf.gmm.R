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
  skip_unless_full()
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
  ## score plot support
  expect_equal(dim(fit$scores), c(d$Q, ncol(d$Y)))
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  expect_null(plot(fit, type = "scores"))
  expect_null(plot(fit, type = "scores", group = d$z))
  expect_null(plot(fit, type = "convergence"))
  dev.off(); unlink(tmp)
})

test_that("K=1 nmf.gmm degenerates to a single-class (NMF-RE-type) fit", {
  skip_unless_full()
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
  skip_unless_full()
  d <- make_gmm_data()
  fit <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, X.init = d$X, seed = 1)
  inf <- nmf.gmm.inference(fit, d$Y, d$A, wild.B = 200, seed = 1)
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
  skip_unless_full()
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
  skip_unless_full()
  d <- make_gmm_data()
  ff <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, cov = "free", X.init = d$X, seed = 1)
  expect_equal(dim(ff$tau2), c(d$Q, 2L))               # per-class diagonal
  expect_true(is.finite(ff$loglik))
})

test_that("cov='scalar' is the isotropic (single-variance) variant", {
  skip_unless_full()
  d <- make_gmm_data()
  sc <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, cov = "scalar", X.init = d$X, seed = 1)
  ti <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, cov = "tied",   X.init = d$X, seed = 1)
  expect_length(sc$tau2, 1L)                            # a single isotropic variance
  expect_length(ti$tau2, d$Q)
  expect_equal(sc$n.params, ti$n.params - (d$Q - 1))    # Q-1 fewer variance params
  ## EM stays monotone (valid GEM)
  ll <- -sc$objfunc.iter
  expect_true(all(diff(ll) >= -1e-6))
  ## K=1 scalar: single cluster, one variance (the nmfre model)
  g1 <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 1, cov = "scalar", X.init = d$X, seed = 1)
  expect_length(g1$tau2, 1L)
  expect_true(all(g1$cluster == 1L))
})

test_that("the scatterplot draws the ADJUSTED scores, and says so", {
  skip_unless_full()
  set.seed(1)
  Y <- matrix(abs(rnorm(6 * 60)) + 1, 6, 60)
  A <- rbind(1, rnorm(60))
  f <- nmf.gmm(Y, A, rank = 3, K = 2, seed = 1, verbose = FALSE)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ## "scores" is kept as an alias for the same panel
  expect_silent(plot(f, type = "adjusted.scores"))
  expect_silent(plot(f, type = "scores"))
  ## the alias really lands on the same panel
  plot(f, type = "scores");          u1 <- par("usr")
  plot(f, type = "adjusted.scores"); u2 <- par("usr")
  expect_identical(u1, u2)
  ## and that panel plots b - C a, not b
  adj <- f$scores - f$C %*% f$A
  expect_false(isTRUE(all.equal(as.numeric(adj), as.numeric(f$scores))))
  ## the drawn x range covers the ADJUSTED first PC, not the raw one
  pc.adj <- prcomp(t(adj))$x[, 1]
  expect_lte(u2[1], min(pc.adj))
  expect_gte(u2[2], max(pc.adj))
})

test_that("plotting never advances the caller's random stream", {
  skip_unless_full()
  set.seed(1)
  Y <- matrix(abs(rnorm(6 * 40)) + 1, 6, 40)
  A <- rbind(1, rnorm(40))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  for (rk in 1:3) {                      # rank 1 used to call runif()
    f <- nmf.gmm(Y, A, rank = rk, K = 2, seed = 1, verbose = FALSE)
    for (ty in c("convergence", "adjusted.scores")) {
      set.seed(99); before <- runif(3)
      set.seed(99); plot(f, type = ty); after <- runif(3)
      expect_identical(before, after)
    }
  }
})

test_that("the rank-1 strip plot is deterministic", {
  skip_unless_full()
  set.seed(1)
  Y <- matrix(abs(rnorm(6 * 40)) + 1, 6, 40)
  A <- rbind(1, rnorm(40))
  f <- nmf.gmm(Y, A, rank = 1, K = 2, seed = 1, verbose = FALSE)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ## two draws must land the points in exactly the same places
  grab <- function() { plot(f, type = "adjusted.scores"); par("usr") }
  expect_identical(grab(), grab())
})

test_that("a rank-2 fit is drawn in its own coordinates, not rotated", {
  skip_unless_full()
  set.seed(1)
  Y <- matrix(abs(rnorm(6 * 40)) + 1, 6, 40)
  A <- rbind(1, rnorm(40))
  f <- nmf.gmm(Y, A, rank = 2, K = 2, seed = 1, verbose = FALSE)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  plot(f, type = "adjusted.scores")
  adj <- f$scores - f$C %*% f$A
  ## the x range is that of adjusted basis 1 itself; a PCA rotation would
  ## generally give a different one
  expect_equal(par("usr")[1] <= min(adj[1, ]), TRUE)
  expect_equal(par("usr")[2] >= max(adj[1, ]), TRUE)
})

test_that("nmf.gmm accepts a one-sided formula covariate with data", {
  skip_unless_full()
  d <- make_gmm_data()
  df <- data.frame(a = d$A[2, ])
  fit_mat <- nmf.gmm(d$Y, d$A, rank = d$Q, K = 2, X.init = d$X, seed = 1)
  ## standardize = FALSE reproduces the numeric-A fit exactly
  fit_fml <- nmf.gmm(d$Y, ~ a, rank = d$Q, K = 2, X.init = d$X, seed = 1,
                     data = df, standardize = FALSE)
  expect_equal(unname(fit_fml$A[2, ]), unname(d$A[2, ]))
  expect_equal(fit_fml$cluster, fit_mat$cluster)
  expect_equal(fit_fml$loglik, fit_mat$loglik, tolerance = 1e-8)
  expect_s3_class(fit_fml$A.formula, "formula")
  ## standardized default still recovers the classes and records the transform
  fit_std <- nmf.gmm(d$Y, ~ a, rank = d$Q, K = 2, X.init = d$X, seed = 1,
                     data = df)
  expect_gt(nmfkc:::.nmfgmm.ARI(fit_std$cluster, d$z), 0.8)
  expect_length(fit_std$A.center, 1L)
  expect_length(fit_std$A.scale, 1L)
  ## a factor covariate expands to indicator rows
  df$g <- factor(rep(1:4, length.out = ncol(d$Y)))
  fit_fac <- nmf.gmm(d$Y, ~ a + g, rank = d$Q, K = 2, X.init = d$X, seed = 1,
                     data = df)
  expect_equal(nrow(fit_fac$A), 1L + 1L + 3L)   # intercept + a + 3 indicators
  ## formula without data is an error
  expect_error(nmf.gmm(d$Y, ~ a, rank = d$Q, K = 2), "data")
})

test_that("nmf.gmm.twostage runs the matched adjust-then-cluster baseline", {
  skip_unless_full()
  d <- make_gmm_data()
  ts <- nmf.gmm.twostage(d$Y, d$A, rank = d$Q, K = 2, X.init = d$X, seed = 1)
  expect_s3_class(ts, "nmf.gmm.twostage")
  expect_s3_class(ts, "nmf.gmm")
  expect_true(is.list(ts$twostage))
  expect_gte(ts$twostage$shift, 0)
  expect_equal(ts$twostage$X0, d$X / rep(colSums(d$X), each = nrow(d$X)))
  ## on well-separated data with a mild covariate both routes recover the classes
  expect_gt(nmfkc:::.nmfgmm.ARI(ts$cluster, d$z), 0.5)
  ## all nmf.gmm S3 methods apply to the returned fit
  expect_equal(dim(fitted(ts)), dim(d$Y))
  expect_length(predict(ts), ncol(d$Y))
  ## intercept-only A is refused
  expect_error(nmf.gmm.twostage(d$Y, matrix(1, 1, ncol(d$Y)), rank = d$Q, K = 2),
               "intercept")
})

test_that("class.effects gives the chosen covariates class-specific coefficients", {
  skip_unless_full()
  ## two classes whose covariate slope differs in sign on part 1
  set.seed(3)
  P <- 10; N <- 120; Q <- 2
  X <- matrix(abs(rnorm(P * Q)), P, Q); X <- sweep(X, 2, colSums(X), "/")
  a <- rnorm(N); A <- rbind(1, a)
  z <- rep(1:2, length.out = N)
  slope <- cbind(c(2, 0.5), c(-2, 0.5))                     # Q x K
  muk <- cbind(c(2, -1), c(-2, 1))
  B <- 3 + sweep(slope[, z], 2, a, "*") + muk[, z] + matrix(rnorm(Q * N) * 0.3, Q, N)
  Y <- X %*% B + matrix(rnorm(P * N) * 0.1, P, N)
  f0 <- nmf.gmm(Y, A, rank = Q, K = 2, X.init = X, seed = 1)
  f1 <- nmf.gmm(Y, A, rank = Q, K = 2, X.init = X, seed = 1, class.effects = 2)
  expect_null(f0$C.class)
  expect_equal(f1$class.rows, 2L)
  expect_equal(dim(f1$C.class), c(Q, 1L, 2L))
  expect_equal(f1$n.params, f0$n.params + (2 - 1) * Q * 1)
  expect_gt(f1$loglik, f0$loglik)
  expect_gt(nmfkc:::.nmfgmm.ARI(f1$cluster, z), 0.8)
  ## C[, S] is the average of the class-specific coefficients weighted by the class
  ## shares of the final responsibilities (fit$xi lags them by one EM step)
  xi <- colMeans(f1$gamma)
  avg <- apply(f1$C.class[, 1, , drop = FALSE], 1, function(v) sum(v * xi))
  expect_equal(unname(avg), unname(f1$C[, 2]), tolerance = 1e-8)
  ## the two classes' slopes on part 1 have opposite signs, as planted
  expect_lt(prod(f1$C.class[1, 1, ]), 0)
  ## the EM stays monotone
  expect_true(all(diff(-f1$objfunc.iter) >= -1e-6))
  ## the reconstruction includes the class-specific terms
  expect_equal(dim(fitted(f1)), dim(Y))
  expect_lt(mean((Y - fitted(f1))^2), mean((Y - fitted(f0))^2))
  ## by name, and at K = 1 the model is the common-effect one
  rownames(A) <- c("Intercept", "a")
  f1n <- nmf.gmm(Y, A, rank = Q, K = 2, X.init = X, seed = 1, class.effects = "a")
  expect_equal(f1n$loglik, f1$loglik, tolerance = 1e-10)
  g0 <- nmf.gmm(Y, A, rank = Q, K = 1, X.init = X, seed = 1)
  g1 <- nmf.gmm(Y, A, rank = Q, K = 1, X.init = X, seed = 1, class.effects = "a")
  expect_equal(g1$loglik, g0$loglik, tolerance = 1e-6)
  ## a formula term name frees every indicator row of a factor
  df <- data.frame(a = a, g = factor(rep(1:3, length.out = N)))
  ff <- nmf.gmm(Y, ~ a + g, rank = Q, K = 2, X.init = X, seed = 1, data = df,
                class.effects = "g", nstart = 2, maxit = 200)
  expect_equal(ff$class.rows, 3:4)
  ## refusals
  expect_error(nmf.gmm(Y, A, rank = Q, K = 2, class.effects = 1), "intercept")
  expect_error(nmf.gmm(Y, A, rank = Q, K = 2, class.effects = "b"), "named")
  expect_error(nmf.gmm.inference(f1, Y, A), "class-specific")
  expect_error(nmf.gmm.twostage(Y, A, rank = Q, K = 2, class.effects = 2), "class-specific")
  ## the adjusted-scores plot removes each unit's own class effect
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_null(plot(f1, type = "adjusted.scores"))
})

test_that("X.init = 'nmf' starts from the ordinary NMF basis, and X0 is returned", {
  skip_unless_full()
  d <- make_gmm_data()
  Y <- d$Y - min(0, min(d$Y))                  # nmfkc() needs non-negative data
  norm <- function(X) X / rep(colSums(X), each = nrow(X))
  Xnmf <- norm(pmax(nmfkc(Y, rank = d$Q, seed = 7, verbose = FALSE)$X, 1e-8))
  fit <- nmf.gmm(Y, d$A, rank = d$Q, K = 2, X.init = "nmf", seed = 7,
                 nstart = 2, maxit = 200)
  expect_equal(unname(fit$X0), unname(Xnmf))
  expect_equal(dimnames(fit$X0), dimnames(fit$X))
  ## a matrix X.init is returned column-normalized, and the default is NNDSVD
  fit2 <- nmf.gmm(Y, d$A, rank = d$Q, K = 2, X.init = d$X, nstart = 2, maxit = 200)
  expect_equal(unname(fit2$X0), unname(norm(d$X)))
  fit3 <- nmf.gmm(Y, d$A, rank = d$Q, K = 2, nstart = 2, maxit = 200)
  expect_equal(unname(fit3$X0),
               unname(norm(nmfkc:::.init_X_method("nndsvd", Y, d$Q, seed = 1))))
  ## the two-stage baseline starts from the same basis as the joint fit
  ts <- nmf.gmm.twostage(Y, d$A, rank = d$Q, K = 2, X.init = "nmf", seed = 7,
                         nstart = 2, maxit = 200)
  expect_equal(unname(ts$twostage$X0), unname(Xnmf))
})

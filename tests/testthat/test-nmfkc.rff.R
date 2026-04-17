## Tests for nmfkc.rff.R (Kernel-faithful NMF with Random Fourier Features)

test_that("nmfkc.rff.random returns expected structure", {
  set.seed(1)
  p <- 5; N <- 12
  U <- matrix(stats::rnorm(p * N), p, N)

  ## default D = ceiling(N/2), nonneg = TRUE
  Z <- nmfkc.rff.random(U, beta = 0.5, seed = 1)
  expect_true(is.list(Z))
  expect_named(Z, c("Zp", "Zn", "pars"))
  expect_equal(dim(Z$Zp), c(ceiling(N / 2), N))
  expect_equal(dim(Z$Zn), c(ceiling(N / 2), N))
  expect_true(all(Z$Zp >= 0))
  expect_true(all(Z$Zn >= 0))
  ## Complementary support
  expect_true(all(Z$Zp * Z$Zn == 0))
  ## pars carries generating parameters
  expect_named(Z$pars, c("omega", "b", "D", "beta"))
  expect_equal(dim(Z$pars$omega), c(ceiling(N / 2), p))
  expect_equal(length(Z$pars$b), ceiling(N / 2))
  expect_equal(Z$pars$D, ceiling(N / 2))
  expect_equal(Z$pars$beta, 0.5)

  ## Explicit D
  Z2 <- nmfkc.rff.random(U, beta = 0.5, D = 30, seed = 1)
  expect_equal(dim(Z2$Zp), c(30, N))

  ## nonneg = FALSE returns signed matrix with pars attribute
  Zraw <- nmfkc.rff.random(U, beta = 0.5, D = 30, seed = 1, nonneg = FALSE)
  expect_true(is.matrix(Zraw))
  expect_equal(dim(Zraw), c(30, N))
  expect_true(!is.null(attr(Zraw, "pars")))
  ## Compare values only (Zraw has $pars attribute, Z2$Zp - Z2$Zn does not)
  expect_equal(unname(as.matrix(Zraw)), Z2$Zp - Z2$Zn, tolerance = 1e-12,
               ignore_attr = TRUE)
})


test_that("nmfkc.rff.random: pars reuse on new data", {
  set.seed(1)
  p <- 5
  U.tr <- matrix(stats::rnorm(p * 40), p, 40)
  U.te <- matrix(stats::rnorm(p * 10), p, 10)

  Ztr <- nmfkc.rff.random(U.tr, beta = 0.5, seed = 1)
  ## Reuse the same omega, b on test data
  Zte <- nmfkc.rff.random(U.te, pars = Ztr$pars)

  expect_equal(dim(Zte$Zp), c(Ztr$pars$D, 10))
  expect_identical(Zte$pars$omega, Ztr$pars$omega)
  expect_identical(Zte$pars$b, Ztr$pars$b)

  ## Error: neither beta nor pars supplied
  expect_error(nmfkc.rff.random(U.te), "beta.*pars")
})


test_that("nmfkc.rff converges and preserves non-negativity", {
  set.seed(1)
  p <- 5; N <- 40; Q_obs <- 6; Q <- 3
  U <- matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(abs(stats::rnorm(Q_obs * N)), Q_obs, N)

  Z <- nmfkc.rff.random(U, beta = 0.5, D = 30, seed = 1)
  res <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = Q, maxit = 200, epsilon = 1e-5)

  expect_s3_class(res, "nmfkc.rff")
  expect_s3_class(res, "nmfkc")
  expect_equal(dim(res$X),  c(Q_obs, Q))
  expect_equal(dim(res$Hp), c(Q, 30))
  expect_equal(dim(res$Hn), c(Q, 30))
  expect_true(all(res$X  >= 0))
  expect_true(all(res$Hp >= 0))
  expect_true(all(res$Hn >= 0))
  expect_true(res$iter >= 1 && res$iter <= 200)
  expect_length(res$objfunc.iter, res$iter)

  diffs <- diff(res$objfunc.iter)
  expect_true(all(diffs <= 1e-8))
})


test_that("predict.nmfkc.rff returns correct shapes and types", {
  set.seed(1)
  p <- 4; N <- 30; Q_obs <- 5; Q <- 2
  U <- matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(abs(stats::rnorm(Q_obs * N)), Q_obs, N)
  rownames(Y) <- paste0("class", 1:Q_obs)

  Z <- nmfkc.rff.random(U, beta = 0.5, D = 20, seed = 1)
  res <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = Q, maxit = 100)

  Yhat <- predict(res, newZp = Z$Zp, newZn = Z$Zn, type = "response")
  expect_equal(dim(Yhat), c(Q_obs, N))
  expect_true(all(Yhat >= 0))

  probs <- predict(res, newZp = Z$Zp, newZn = Z$Zn, type = "prob")
  expect_equal(dim(probs), c(Q_obs, N))
  expect_true(all(abs(colSums(probs) - 1) < 1e-6))

  cls <- predict(res, newZp = Z$Zp, newZn = Z$Zn, type = "class")
  expect_length(cls, N)
  expect_true(all(cls %in% rownames(Y)))

  expect_error(predict(res, newZp = Z$Zp), "newZp.*newZn")
})


test_that("nmfkc.kernel.beta.nearest.med supports candidates argument", {
  set.seed(1)
  U <- matrix(stats::rnorm(5 * 50), 5, 50)

  info7 <- nmfkc.kernel.beta.nearest.med(U)
  expect_length(info7$beta_candidates, 7)
  expect_equal(info7$beta_candidates[4], info7$beta, tolerance = 1e-12)

  info4 <- nmfkc.kernel.beta.nearest.med(U, candidates = "4points")
  expect_length(info4$beta_candidates, 4)

  t_custom <- c(-0.5, 0, 0.5)
  info_c <- nmfkc.kernel.beta.nearest.med(U, candidates = t_custom)
  expect_length(info_c$beta_candidates, 3)
  expect_equal(info_c$beta_candidates,
               info_c$beta * 10^(-2 * t_custom), tolerance = 1e-12)

  expect_error(nmfkc.kernel.beta.nearest.med(U, candidates = "bogus"),
               "7points|4points|numeric")

  Uk <- U[, 1:5, drop = FALSE]
  infoUk <- nmfkc.kernel.beta.nearest.med(U, Uk = Uk)
  expect_length(infoUk$beta_candidates, 7)
})


test_that("warm.start default (TRUE) triggers internal nmfkc() warm-start", {
  set.seed(1)
  p <- 4; N <- 30; Q_obs <- 5; Q <- 2
  U <- matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(abs(stats::rnorm(Q_obs * N)), Q_obs, N)
  Z <- nmfkc.rff.random(U, beta = 0.5, D = 20, seed = 1)

  ## default: warm.start = TRUE, no explicit X.init/C.init
  res_warm <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = Q, maxit = 100)
  expect_s3_class(res_warm, "nmfkc.rff")
  expect_true(all(res_warm$X >= 0))

  ## warm.start = FALSE: random init
  res_rand <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = Q, maxit = 100,
                         warm.start = FALSE)
  expect_s3_class(res_rand, "nmfkc.rff")
  expect_true(all(res_rand$X >= 0))
})


test_that("explicit X.init / C.init override warm.start", {
  set.seed(1)
  p <- 4; N <- 30; Q_obs <- 5; Q <- 2
  U <- matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(abs(stats::rnorm(Q_obs * N)), Q_obs, N)
  Z <- nmfkc.rff.random(U, beta = 0.5, D = 20, seed = 1)
  D_rff <- nrow(Z$Zp)

  ## Externally compute a nmfkc warm-start and pass pieces explicitly.
  res_posneg <- nmfkc(Y, A = rbind(Z$Zp, Z$Zn),
                       rank = Q, maxit = 200, epsilon = 1e-5,
                       print.dims = FALSE)

  res_explicit <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = Q, maxit = 100,
                             X.init = res_posneg$X, C.init = res_posneg$C)
  expect_s3_class(res_explicit, "nmfkc.rff")

  ## Dimension checks
  expect_error(
    nmfkc.rff(Y, Z$Zp, Z$Zn, rank = Q,
              X.init = matrix(0, 2, 2), C.init = res_posneg$C),
    "X.init.*dimensions"
  )
  expect_error(
    nmfkc.rff(Y, Z$Zp, Z$Zn, rank = Q,
              X.init = res_posneg$X, C.init = matrix(0, 1, 1)),
    "C.init.*dimensions"
  )
})

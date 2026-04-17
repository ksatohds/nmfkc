## Tests for nmfkc.rff.R (Kernel-faithful NMF with Random Fourier Features)

test_that("RFF helpers produce expected shapes", {
  set.seed(1)
  p <- 5; D <- 20; N <- 12
  U <- matrix(stats::rnorm(p * N), p, N)

  pars <- nmfkc.rff.params(p = p, D = D, beta = 0.5, seed = 1)
  expect_equal(dim(pars$omega), c(D, p))
  expect_equal(length(pars$b), D)

  Z <- nmfkc.rff.apply(U, pars)
  expect_equal(dim(Z), c(D, N))

  Zaug <- nmfkc.rff.nonneg(Z, posneg.split = TRUE)
  expect_equal(dim(Zaug), c(2 * D, N))
  expect_true(all(Zaug >= 0))

  Zpos_only <- nmfkc.rff.nonneg(Z, posneg.split = FALSE)
  expect_equal(dim(Zpos_only), c(D, N))
  expect_true(all(Zpos_only >= 0))
})


test_that("nmfkc.rff.direct converges and preserves non-negativity", {
  set.seed(1)
  p <- 5; D <- 30; N <- 40; Q_obs <- 6; Q <- 3
  U <- matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(abs(stats::rnorm(Q_obs * N)), Q_obs, N)

  pars <- nmfkc.rff.params(p = p, D = D, beta = 0.5, seed = 1)
  Zraw <- nmfkc.rff.apply(U, pars)
  Zp   <- pmax(Zraw, 0); Zn <- pmax(-Zraw, 0)

  res <- nmfkc.rff.direct(Y, Zp, Zn, rank = Q, maxit = 200, epsilon = 1e-5)

  expect_s3_class(res, "nmfkc.rff")
  expect_s3_class(res, "nmfkc")
  expect_equal(dim(res$X),  c(Q_obs, Q))
  expect_equal(dim(res$Hp), c(Q, D))
  expect_equal(dim(res$Hn), c(Q, D))
  expect_true(all(res$X  >= 0))
  expect_true(all(res$Hp >= 0))
  expect_true(all(res$Hn >= 0))
  expect_true(res$iter >= 1 && res$iter <= 200)
  expect_length(res$objfunc.iter, res$iter)

  ## Objective must be (weakly) monotonically non-increasing
  diffs <- diff(res$objfunc.iter)
  expect_true(all(diffs <= 1e-8))
})


test_that("predict.nmfkc.rff returns correct shapes and types", {
  set.seed(1)
  p <- 4; D <- 20; N <- 30; Q_obs <- 5; Q <- 2
  U <- matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(abs(stats::rnorm(Q_obs * N)), Q_obs, N)
  rownames(Y) <- paste0("class", 1:Q_obs)

  pars <- nmfkc.rff.params(p = p, D = D, beta = 0.5, seed = 1)
  Zraw <- nmfkc.rff.apply(U, pars)
  Zp <- pmax(Zraw, 0); Zn <- pmax(-Zraw, 0)

  res <- nmfkc.rff.direct(Y, Zp, Zn, rank = Q, maxit = 100)

  Yhat <- predict(res, newZp = Zp, newZn = Zn, type = "response")
  expect_equal(dim(Yhat), c(Q_obs, N))
  expect_true(all(Yhat >= 0))

  probs <- predict(res, newZp = Zp, newZn = Zn, type = "prob")
  expect_equal(dim(probs), c(Q_obs, N))
  expect_true(all(abs(colSums(probs) - 1) < 1e-6))

  cls <- predict(res, newZp = Zp, newZn = Zn, type = "class")
  expect_length(cls, N)
  expect_true(all(cls %in% rownames(Y)))

  ## Missing Zp or Zn raises an error
  expect_error(predict(res, newZp = Zp), "newZp.*newZn")
})


test_that("warm-start from nmfkc() works", {
  set.seed(1)
  p <- 4; D <- 20; N <- 30; Q_obs <- 5; Q <- 2
  U <- matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(abs(stats::rnorm(Q_obs * N)), Q_obs, N)

  pars <- nmfkc.rff.params(p = p, D = D, beta = 0.5, seed = 1)
  Zraw <- nmfkc.rff.apply(U, pars)
  Zp <- pmax(Zraw, 0); Zn <- pmax(-Zraw, 0)
  Zaug <- rbind(Zp, Zn)

  res_posneg <- nmfkc(Y, A = Zaug, rank = Q, maxit = 200, epsilon = 1e-5,
                       print.dims = FALSE)
  res_warm <- nmfkc.rff.direct(Y, Zp, Zn, rank = Q, maxit = 100,
                                res.init = res_posneg)

  expect_s3_class(res_warm, "nmfkc.rff")
  expect_true(all(res_warm$X >= 0))
})

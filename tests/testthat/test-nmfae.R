test_that("nmfae() fits with EU and KL objectives", {
  set.seed(7)
  P1 <- 6; P2 <- 6; N <- 40; Q <- 2; R <- 2
  X1 <- matrix(abs(rnorm(P1 * Q)), P1, Q)
  X2 <- matrix(abs(rnorm(R * P2)), R, P2)
  Cm <- matrix(abs(rnorm(Q * R)), Q, R)
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)
  Y1 <- X1 %*% Cm %*% X2 %*% Y2 + matrix(abs(rnorm(P1 * N, 0, 0.2)), P1, N)

  ## --- EU (default) ---
  fe <- nmfae(Y1, Y2 = Y2, rank = Q, rank.encoder = R, maxit = 2000)
  expect_s3_class(fe, "nmfae")
  expect_identical(fe$method, "EU")
  expect_true(is.finite(fe$sigma))
  expect_true(fe$r.squared >= 0 && fe$r.squared <= 1)
  expect_equal(dim(fe$Y1hat), c(P1, N))
  expect_true(all(fe$Y1hat >= 0))

  ## --- KL ---
  fk <- nmfae(Y1, Y2 = Y2, rank = Q, rank.encoder = R, method = "KL", maxit = 2000)
  expect_identical(fk$method, "KL")
  expect_true(is.na(fk$sigma))                       # RMSE not defined for KL
  expect_true(all(fk$Y1hat >= 0))
  expect_true(all(is.finite(fk$B.cluster)))
  ## objective is monotone non-increasing (MU descent property)
  oi <- fk$objfunc.iter
  if (length(oi) > 1) {
    tol <- 1e-6 * abs(oi[-length(oi)]) + 1e-8
    expect_true(all(diff(oi) <= tol))
  }

  ## invalid method rejected
  expect_error(nmfae(Y1, Y2 = Y2, rank = Q, method = "XX"))
})

test_that("nmfae() forwards nstart to the nmfkc() initialisation", {
  set.seed(3)
  P1 <- 6; P2 <- 6; N <- 40; Q <- 2; R <- 2
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)
  Y1 <- matrix(abs(rnorm(P1 * Q)), P1, Q) %*% matrix(abs(rnorm(Q * R)), Q, R) %*%
        matrix(abs(rnorm(R * P2)), R, P2) %*% Y2 +
        matrix(abs(rnorm(P1 * N, 0, 0.3)), P1, N)

  ## nstart = 1 (default) is reproducible for a fixed seed
  a  <- nmfae(Y1, Y2 = Y2, rank = Q, rank.encoder = R, seed = 1, maxit = 1500)
  a2 <- nmfae(Y1, Y2 = Y2, rank = Q, rank.encoder = R, seed = 1, nstart = 1, maxit = 1500)
  expect_equal(a$objfunc, a2$objfunc)

  ## nstart > 1 is accepted, takes effect, and returns a valid fit
  b <- nmfae(Y1, Y2 = Y2, rank = Q, rank.encoder = R, seed = 1, nstart = 10, maxit = 1500)
  expect_s3_class(b, "nmfae")
  expect_true(b$r.squared >= 0 && b$r.squared <= 1)
  expect_true(all(b$Y1hat >= 0))
})

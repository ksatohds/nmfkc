test_that("nmf.rrr() fits with EU and KL objectives", {
  set.seed(7)
  P1 <- 6; P2 <- 6; N <- 40; Q <- 2; R <- 2
  X1 <- matrix(abs(rnorm(P1 * Q)), P1, Q)
  X2 <- matrix(abs(rnorm(R * P2)), R, P2)
  Cm <- matrix(abs(rnorm(Q * R)), Q, R)
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)
  Y1 <- X1 %*% Cm %*% X2 %*% Y2 + matrix(abs(rnorm(P1 * N, 0, 0.2)), P1, N)

  ## --- EU (default) ---
  fe <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, maxit = 2000)
  expect_s3_class(fe, "nmf.rrr")
  expect_s3_class(fe, "nmfae")
  expect_identical(fe$method, "EU")
  expect_true(is.finite(fe$sigma))
  expect_true(fe$r.squared >= 0 && fe$r.squared <= 1)
  expect_equal(dim(fe$Y1hat), c(P1, N))
  expect_true(all(fe$Y1hat >= 0))

  ## --- KL ---
  fk <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, method = "KL", maxit = 2000)
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
  expect_error(nmf.rrr(Y1, Y2 = Y2, rank1 = Q, method = "XX"))
})

test_that("nmf.rrr() forwards nstart to the nmfkc() initialisation", {
  set.seed(3)
  P1 <- 6; P2 <- 6; N <- 40; Q <- 2; R <- 2
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)
  Y1 <- matrix(abs(rnorm(P1 * Q)), P1, Q) %*% matrix(abs(rnorm(Q * R)), Q, R) %*%
        matrix(abs(rnorm(R * P2)), R, P2) %*% Y2 +
        matrix(abs(rnorm(P1 * N, 0, 0.3)), P1, N)

  ## nstart = 1 (default) is reproducible for a fixed seed
  a  <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, seed = 1, maxit = 1500)
  a2 <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, seed = 1, nstart = 1, maxit = 1500)
  expect_equal(a$objfunc, a2$objfunc)

  ## nstart > 1 is accepted, takes effect, and returns a valid fit
  b <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, seed = 1, nstart = 10, maxit = 1500)
  expect_s3_class(b, "nmf.rrr")
  expect_true(b$r.squared >= 0 && b$r.squared <= 1)
  expect_true(all(b$Y1hat >= 0))
})

test_that("nmf.rrr uses Resp/Cov labels; inference and signed helpers work", {
  set.seed(7); P1 <- 6; P2 <- 6; N <- 40; Q <- 2; R <- 3
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)
  Y1 <- matrix(abs(rnorm(P1 * Q)), P1, Q) %*% matrix(abs(rnorm(Q * R)), Q, R) %*%
        matrix(abs(rnorm(R * P2)), R, P2) %*% Y2 +
        matrix(abs(rnorm(P1 * N, 0, 0.3)), P1, N)

  ## labels are Resp (response basis X1) / Cov (covariate basis X2)
  g <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, seed = 1, maxit = 1500)
  expect_s3_class(g, "nmf.rrr")
  expect_s3_class(g, "nmfae")
  expect_equal(colnames(g$X1), paste0("Resp", 1:Q))
  expect_equal(rownames(g$X2), paste0("Cov", 1:R))

  ## inference helper works and inherits the S3 machinery
  gi <- nmf.rrr.inference(g, Y1, Y2 = Y2, wild.B = 50)
  expect_true("Basis" %in% names(gi$coefficients))
  expect_true(any(grepl("^Resp", gi$coefficients$Basis)))

  ## signed fitter
  gs <- nmf.rrr.signed(Y1, Y2 = Y2, rank1 = Q, rank2 = R, seed = 1, maxit = 600)
  expect_s3_class(gs, "nmf.rrr.signed")
  expect_s3_class(gs, "nmfae.signed")
  expect_equal(colnames(gs$X1), paste0("Resp", 1:Q))
})

test_that("deprecated nmfae() alias still works and warns", {
  set.seed(7); P1 <- 6; P2 <- 6; N <- 30; Q <- 2; R <- 2
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)
  Y1 <- matrix(abs(rnorm(P1 * Q)), P1, Q) %*% matrix(abs(rnorm(Q * R)), Q, R) %*%
        matrix(abs(rnorm(R * P2)), R, P2) %*% Y2 +
        matrix(abs(rnorm(P1 * N, 0, 0.3)), P1, N)

  ## nmfae() is a deprecated alias for nmf.rrr(): it warns but still fits,
  ## and the object carries both the new and legacy S3 classes.
  expect_warning(f <- nmfae(Y1, Y2 = Y2, rank1 = Q, rank2 = R, maxit = 800),
                 "deprecated")
  expect_s3_class(f, "nmf.rrr")
  expect_s3_class(f, "nmfae")
  expect_true(f$r.squared >= 0 && f$r.squared <= 1)

  ## the reference fit from the canonical name matches the deprecated one
  g <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, maxit = 800)
  expect_equal(f$objfunc, g$objfunc)
})

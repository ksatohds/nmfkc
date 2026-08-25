test_that("nmf.rrr() fits with EU and KL objectives", {
  skip_unless_full()
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

  ## B1 is the decoder-side score (Q x N), B2 the encoder-side one (R x N);
  ## the old single `B` and its .prob/.cluster are gone.
  expect_equal(dim(fe$B1), c(Q, N))
  expect_equal(dim(fe$B2), c(R, N))
  expect_identical(fe$H, fe$B1)                      # $H is a deprecated alias
  expect_equal(fe$B1, fe$C %*% fe$B2, tolerance = 1e-12)
  expect_null(fe$B); expect_null(fe$B.prob); expect_null(fe$B.cluster)
  ## Columns of each .prob are memberships: non-negative and summing to one.
  expect_equal(unname(colSums(fe$B1.prob)), rep(1, N), tolerance = 1e-10)
  expect_equal(unname(colSums(fe$B2.prob)), rep(1, N), tolerance = 1e-10)
  expect_true(all(fe$B1.prob >= 0) && all(fe$B2.prob >= 0))
  expect_identical(fe$B1.cluster, apply(fe$B1.prob, 2, which.max))
  expect_identical(fe$B2.cluster, apply(fe$B2.prob, 2, which.max))
  expect_true(all(fe$B2.cluster %in% seq_len(R)))

  ## --- KL ---
  fk <- nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, method = "KL", maxit = 2000)
  expect_identical(fk$method, "KL")
  expect_true(is.na(fk$sigma))                       # RMSE not defined for KL
  expect_true(all(fk$Y1hat >= 0))
  expect_true(all(is.finite(fk$B1.cluster)))
  expect_true(all(is.finite(fk$B2.cluster)))
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
  skip_unless_full()
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
  skip_unless_full()
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
  skip_unless_full()
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

test_that("the KL convergence test survives a negative objective", {
  skip_unless_full()
  ## KL objective sum(-Y log Yhat + Yhat) goes negative once Y exceeds e over
  ## most cells.  Without abs() on the denominator the relative change came out
  ## NEGATIVE, which is < any epsilon, so the fit "converged" after two
  ## iterations no matter how tight the tolerance was.
  set.seed(2)
  Y1 <- matrix(abs(rnorm(4 * 30)) * 50 + 50, 4, 30)
  Y2 <- matrix(abs(rnorm(4 * 30)) * 50 + 50, 4, 30)
  f <- suppressWarnings(nmf.rrr(Y1, Y2, rank1 = 2, rank2 = 2, method = "KL",
                                epsilon = 1e-4, maxit = 2000, print.trace = FALSE))
  expect_true(any(f$objfunc.iter < 0))       # the trace really is negative
  ## tightening epsilon must now change the run; before the fix it could not
  g <- suppressWarnings(nmf.rrr(Y1, Y2, rank1 = 2, rank2 = 2, method = "KL",
                                epsilon = 1e-9, maxit = 2000, print.trace = FALSE))
  expect_gt(g$iter, f$iter)
})

test_that("nmf.rrr(maxit = 1) does not error", {
  skip_unless_full()
  ## rel_change is only bound from the second iteration, and `&&` does not
  ## short-circuit past an unbound name.
  set.seed(3)
  Y <- matrix(abs(rnorm(6 * 20)) + 1, 6, 20)
  f <- suppressWarnings(nmf.rrr(Y[1:3, ], Y[4:6, ], rank1 = 2, rank2 = 2,
                                maxit = 1, print.trace = FALSE))
  expect_false(f$converged)
  expect_equal(f$iter, 1)
})

test_that("nmf.rrr.signed returns B1/B2 but no memberships for them", {
  skip_unless_full()
  set.seed(11)
  P1 <- 5; P2 <- 4; N <- 12; Q <- 2; R <- 3
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)
  Y1 <- matrix(runif(P1 * N) + 1, P1, N)
  g <- suppressWarnings(nmf.rrr.signed(Y1, Y2 = Y2, rank1 = Q, rank2 = R, maxit = 100))

  expect_equal(dim(g$B1), c(Q, N))
  expect_equal(dim(g$B2), c(R, N))
  expect_equal(g$B1, g$C %*% g$B2, tolerance = 1e-12)
  expect_equal(g$Y1hat, g$X1 %*% g$B1, tolerance = 1e-12)
  expect_identical(g$H, g$B1)             # deprecated alias
  expect_true(all(g$B2 >= 0))             # X2, Y2 both non-negative

  ## No .prob / .cluster for B1 or B2: C is signed, so a column-normalized
  ## membership would not be a distribution.  The clipped B.prob stays.
  expect_null(g$B1.prob); expect_null(g$B1.cluster)
  expect_null(g$B2.prob); expect_null(g$B2.cluster)
  expect_false(is.null(g$B.prob))
})

test_that("sample names survive when only Y1 is labelled", {
  skip_unless_full()
  set.seed(12)
  P1 <- 5; P2 <- 4; N <- 9; Q <- 2; R <- 3
  Y2 <- matrix(rpois(P2 * N, 4) + 0.1, P2, N)   # deliberately unlabelled
  Y1 <- matrix(runif(P1 * N) + 1, P1, N)
  rownames(Y1) <- paste0("resp", 1:P1); colnames(Y1) <- paste0("s", 1:N)

  ## Everything downstream of Y2 used to inherit its names from Y2 alone, so
  ## labelling only Y1 left B1/B2/Y1hat unnamed even though the columns are
  ## the same N samples.
  for (f in list(suppressWarnings(nmf.rrr(Y1, Y2 = Y2, rank1 = Q, rank2 = R, maxit = 100)),
                 suppressWarnings(nmf.rrr.signed(Y1, Y2 = Y2, rank1 = Q, rank2 = R,
                                                 maxit = 100)))) {
    expect_identical(colnames(f$B1), colnames(Y1))
    expect_identical(colnames(f$B2), colnames(Y1))
    expect_identical(colnames(f$Y1hat), colnames(Y1))
    expect_identical(rownames(f$Y1hat), rownames(Y1))
    expect_identical(rownames(f$B1), colnames(f$X1))   # Resp1..Q
    expect_identical(rownames(f$B2), rownames(f$X2))   # Cov1..R
  }
})

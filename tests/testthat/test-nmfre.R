test_that("nmfre() X.init paths run (CRAN-sized)", {
  skip_unless_full()
  # The full block below fits nmfre() eight times (~6s locally, x17 on CRAN's
  # regression covered there: character X.init used to error in
  # .nmfre.normalize.X ("'x' must be an array of at least two dimensions").
  set.seed(5); P <- 5; N <- 30; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  Y  <- pmax(Xt %*% (matrix(abs(rnorm(Q * 2)), Q, 2) %*% A +
                       matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)
  d <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE)
  expect_s3_class(d, "nmfre")
  fm <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = "nndsvd", verbose = FALSE)
  expect_true(all(is.finite(fm$X)))
})

test_that("nmfre() accepts X.init as NULL, a method name, or a matrix", {
  skip_unless_full()
  set.seed(5); P <- 5; N <- 50; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  Th <- matrix(abs(rnorm(Q * 2)), Q, 2)
  Y  <- pmax(Xt %*% (Th %*% A + matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)

  ## default (NULL) works
  d <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE)
  expect_s3_class(d, "nmfre")
  expect_equal(dim(d$X), c(P, Q))

  ## character init methods no longer crash (regression: used to error in
  ## .nmfre.normalize.X with "'x' must be an array of at least two dimensions")
  for (m in c("runif", "nndsvd", "kmeans")) {
    fm <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = m, verbose = FALSE)
    expect_s3_class(fm, "nmfre")
    expect_equal(dim(fm$X), c(P, Q))
    expect_true(all(is.finite(fm$X)))
  }

  ## a supplied basis matrix works
  Xm <- matrix(runif(P * Q), P, Q)
  fmat <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = Xm, verbose = FALSE)
  expect_equal(dim(fmat$X), c(P, Q))

  ## nstart is accepted and the default is reproducible
  a  <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE)
  a1 <- nmfre(Y, A = A, rank = Q, seed = 1, nstart = 1, verbose = FALSE)
  expect_equal(a$objfunc, a1$objfunc)
  expect_s3_class(nmfre(Y, A = A, rank = Q, seed = 1, nstart = 5, verbose = FALSE), "nmfre")
})

test_that("nmfre.inference: identification is NA for a signed C, by the data for C >= 0", {
  skip_unless_full()
  set.seed(5); P <- 5; N <- 50; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  Th <- matrix(abs(rnorm(Q * 2)), Q, 2)
  Y  <- pmax(Xt %*% (Th %*% A + matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)
  ## C.signed = TRUE (default): non-negativity cannot pin T^-1
  f <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE)
  inf <- suppressWarnings(nmfre.inference(f, Y, A, wild.B = 5))
  expect_true(is.na(inf$identified.approx))
  expect_match(inf$identification$note, "C is signed")
  ## Q = 2 with tau2 > 0: the isotropic random effects allow only a permutation
  expect_identical(inf$identification$by.assumption, isTRUE(f$tau2 > 0))
  expect_output(print(summary(inf)), "Identification:")
  ## C >= 0: the data diagnostic, exactly the shared check on (X, C)
  f2 <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE, C.signed = FALSE)
  inf2 <- suppressWarnings(nmfre.inference(f2, Y, A, wild.B = 5))
  expect_false(is.na(inf2$identified.approx))
  expect_identical(inf2$identified.approx, nmfkc:::.ident.verdict(f2$X, f2$C, A)$identified.approx)
})

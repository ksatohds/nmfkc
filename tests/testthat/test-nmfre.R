test_that("nmfre() accepts X.init as NULL, a method name, or a matrix", {
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

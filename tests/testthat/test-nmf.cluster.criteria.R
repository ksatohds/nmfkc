## Tests for nmf.cluster.criteria() sample-clustering diagnostics

test_that("nmf.cluster.criteria computes silhouette/CPCC/dist.cor for non-negative B", {
  Y <- t(as.matrix(iris[, 1:4]))
  fit <- nmfkc(Y, Q = 3, print.dims = FALSE)
  cl <- nmf.cluster.criteria(fit, Y)

  expect_s3_class(cl, "nmf.cluster.criteria")
  expect_true(cl$hard)
  expect_true(is.finite(cl$silhouette))
  expect_true(is.finite(cl$CPCC))
  expect_true(is.finite(cl$dist.cor))
  expect_length(cl$cluster, ncol(Y))
  expect_equal(cl$rank, 3)
})

test_that("nmf.cluster.criteria sets silhouette NA when the coefficient is signed", {
  set.seed(7)
  Y <- matrix(rnorm(8 * 20), 8, 20)
  A <- rbind(intercept = 1, x = rnorm(20))
  fit <- suppressWarnings(nmfkc.signed(Y, A, rank = 3, maxit = 500))
  skip_if_not(any(fit$B < 0), "coefficient happened to be non-negative")

  cl <- nmf.cluster.criteria(fit, Y)
  expect_false(cl$hard)
  expect_true(is.na(cl$silhouette))
  ## distance-based criteria are still defined for a signed B
  expect_true(is.finite(cl$CPCC))
  expect_true(is.finite(cl$dist.cor))
})

test_that("nmf.cluster.criteria matches nmfkc.criterion's values (shared helper)", {
  Y <- t(as.matrix(iris[, 1:4]))
  fit <- nmfkc(Y, Q = 3, print.dims = FALSE)
  cl  <- nmf.cluster.criteria(fit, Y)
  crit <- fit$criterion
  expect_equal(cl$CPCC,     crit$CPCC,     tolerance = 1e-8)
  expect_equal(cl$dist.cor, crit$dist.cor, tolerance = 1e-8)
  expect_equal(cl$silhouette, crit$silhouette, tolerance = 1e-8)
})

test_that("nmf.cluster.criteria requires Y", {
  Y <- t(as.matrix(iris[, 1:4]))
  fit <- nmfkc(Y, Q = 2, print.dims = FALSE)
  expect_error(nmf.cluster.criteria(fit), "requires the original data")
})

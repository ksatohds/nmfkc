## Tests for nmf.cluster.criteria() sample-clustering quality across ranks

test_that("nmf.cluster.criteria builds a per-rank criteria table from a fit list", {
  Y <- t(as.matrix(iris[, 1:4]))
  fits <- lapply(2:4, function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  cc <- nmf.cluster.criteria(fits, Y, plot = FALSE)

  expect_s3_class(cc, "nmf.cluster.criteria")
  expect_equal(cc$criteria$rank, 2:4)
  expect_setequal(names(cc$criteria),
                  c("rank", "silhouette", "CPCC", "dist.cor", "hard"))
  expect_true(all(cc$criteria$hard))                 # non-negative B
  expect_true(all(is.finite(cc$criteria$silhouette)))
  expect_true(all(is.finite(cc$criteria$CPCC)))
  expect_true(all(is.finite(cc$criteria$dist.cor)))
})

test_that("nmf.cluster.criteria accepts a single fit (wrapped)", {
  Y <- t(as.matrix(iris[, 1:4]))
  fit <- nmfkc(Y, Q = 3, print.dims = FALSE)
  cc <- nmf.cluster.criteria(fit, Y, plot = FALSE)
  expect_equal(nrow(cc$criteria), 1L)
  expect_equal(cc$criteria$rank, 3L)
})

test_that("nmf.cluster.criteria sets silhouette NA when the coefficient is signed", {
  set.seed(7)
  Y <- matrix(rnorm(8 * 20), 8, 20)
  A <- rbind(intercept = 1, x = rnorm(20))
  fit <- suppressWarnings(nmfkc.signed(Y, A, rank = 3, maxit = 500))
  skip_if_not(any(fit$B < 0), "coefficient happened to be non-negative")

  cc <- nmf.cluster.criteria(fit, Y, plot = FALSE)
  expect_false(cc$criteria$hard[1])
  expect_true(is.na(cc$criteria$silhouette[1]))
  ## distance-based criteria are still defined for a signed B
  expect_true(is.finite(cc$criteria$CPCC[1]))
  expect_true(is.finite(cc$criteria$dist.cor[1]))
})

test_that("nmf.cluster.criteria matches nmfkc.criterion's values (shared helper)", {
  Y <- t(as.matrix(iris[, 1:4]))
  fit <- nmfkc(Y, Q = 3, print.dims = FALSE)
  cc  <- nmf.cluster.criteria(fit, Y, plot = FALSE)
  crit <- fit$criterion
  expect_equal(cc$criteria$CPCC,       crit$CPCC,       tolerance = 1e-8)
  expect_equal(cc$criteria$dist.cor,   crit$dist.cor,   tolerance = 1e-8)
  expect_equal(cc$criteria$silhouette, crit$silhouette, tolerance = 1e-8)
})

test_that("nmf.cluster.criteria requires Y", {
  Y <- t(as.matrix(iris[, 1:4]))
  fit <- nmfkc(Y, Q = 2, print.dims = FALSE)
  expect_error(nmf.cluster.criteria(fit), "requires the original data")
})

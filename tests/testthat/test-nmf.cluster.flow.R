## Tests for nmf.cluster.flow() cluster-flow diagram across ranks

test_that("nmf.cluster.flow returns an N x R cluster table", {
  Y <- t(as.matrix(iris[, 1:4]))
  fits <- lapply(2:5, function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  fl <- nmf.cluster.flow(fits, reference = 3, plot = FALSE)

  expect_equal(dim(fl$clusters), c(ncol(Y), 4L))
  expect_equal(fl$ranks, 2:5)
  expect_equal(fl$reference, 3)
  expect_equal(colnames(fl$clusters), paste0("rank", 2:5))
  ## entries are valid cluster numbers (1..q at each rank)
  for (j in seq_along(fl$ranks)) {
    expect_true(all(fl$clusters[, j] >= 1))
    expect_true(all(fl$clusters[, j] <= fl$ranks[j]))
  }
  ## one colour per individual
  expect_length(fl$colors, ncol(Y))
})

test_that("nmf.cluster.flow sorts fits by rank regardless of input order", {
  Y <- t(as.matrix(iris[, 1:4]))
  fits <- lapply(c(4, 2, 3), function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  fl <- nmf.cluster.flow(fits, plot = FALSE)
  expect_equal(fl$ranks, c(2L, 3L, 4L))
})

test_that("nmf.cluster.flow rejects a single fit and an unknown reference", {
  Y <- t(as.matrix(iris[, 1:4]))
  one <- list(nmfkc(Y, Q = 2, print.dims = FALSE))
  expect_error(nmf.cluster.flow(one, plot = FALSE), "at least two")

  fits <- lapply(2:4, function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  expect_error(nmf.cluster.flow(fits, reference = 9, plot = FALSE),
               "must be one of the fitted ranks")
})

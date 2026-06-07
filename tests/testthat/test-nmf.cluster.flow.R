## Tests for nmf.cluster.flow() cluster-flow diagram across a sequence of fits

test_that("nmf.cluster.flow returns an N x R cluster table", {
  Y <- t(as.matrix(iris[, 1:4]))
  fits <- lapply(2:5, function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  fl <- nmf.cluster.flow(fits, reference = 2, plot = FALSE)

  expect_equal(dim(fl$clusters), c(ncol(Y), 4L))
  expect_equal(fl$ranks, 2:5)
  expect_equal(fl$labels, as.character(2:5))
  expect_equal(fl$reference, 2L)                 # reference is an index
  expect_equal(colnames(fl$clusters), as.character(2:5))
  ## entries are valid cluster numbers (1..rank at each result)
  for (j in seq_along(fl$ranks)) {
    expect_true(all(fl$clusters[, j] >= 1))
    expect_true(all(fl$clusters[, j] <= fl$ranks[j]))
  }
  ## one colour per individual
  expect_length(fl$colors, ncol(Y))
})

test_that("nmf.cluster.flow keeps the given order (no rank sorting)", {
  Y <- t(as.matrix(iris[, 1:4]))
  fits <- lapply(c(4, 2, 3), function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  fl <- nmf.cluster.flow(fits, plot = FALSE)
  expect_equal(fl$ranks, c(4L, 2L, 3L))          # order preserved
  expect_equal(fl$labels, c("4", "2", "3"))
})

test_that("nmf.cluster.flow honours a custom `names` for the labels", {
  Y <- t(as.matrix(iris[, 1:4]))
  fits <- lapply(c(3, 3, 3), function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  fl <- nmf.cluster.flow(fits, names = c("modelA", "modelB", "modelC"),
                         plot = FALSE)
  expect_equal(fl$labels, c("modelA", "modelB", "modelC"))
  expect_equal(colnames(fl$clusters), c("modelA", "modelB", "modelC"))
  expect_error(nmf.cluster.flow(fits, names = c("a", "b"), plot = FALSE),
               "length\\(fits\\)")
})

test_that("nmf.cluster.flow rejects a single fit and an out-of-range reference", {
  Y <- t(as.matrix(iris[, 1:4]))
  one <- list(nmfkc(Y, Q = 2, print.dims = FALSE))
  expect_error(nmf.cluster.flow(one, plot = FALSE), "at least two")

  fits <- lapply(2:4, function(q) nmfkc(Y, Q = q, print.dims = FALSE))
  expect_error(nmf.cluster.flow(fits, reference = 9, plot = FALSE),
               "result index in 1")
})

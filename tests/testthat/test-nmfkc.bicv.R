## Tests for nmfkc.bicv() (bi-cross-validation; lightweight, like nmfkc.ecv)

test_that("nmfkc.bicv returns sigma per rank and recovers a rank-3 signal", {
  set.seed(1)
  X <- matrix(abs(rnorm(30 * 3)), 30, 3)
  B <- matrix(abs(rnorm(3 * 40)), 3, 40)
  Y <- X %*% B                         # clean rank-3 non-negative matrix
  bv <- suppressMessages(nmfkc.bicv(Y, rank = 1:6, nfolds = 2, seed = 7))

  ## plain list, same spirit as nmfkc.ecv (objfunc / sigma)
  expect_type(bv, "list")
  expect_equal(bv$rank, 1:6)
  expect_length(bv$sigma, 6)
  expect_equal(bv$sigma, sqrt(bv$objfunc), tolerance = 1e-8)
  ## clean rank-3 data: error drops sharply from rank 2 to rank 3
  expect_lt(bv$sigma[3], bv$sigma[2])
  expect_true(bv$rank[which.min(bv$sigma)] >= 3)
})

test_that("nmfkc.bicv warns and returns NA for ranks too large for the block", {
  set.seed(2)
  Y <- matrix(abs(rnorm(8 * 50)), 8, 50)   # only 8 rows -> nfolds=2 keeps ~4
  ## rank 6 leaves <=4 kept rows per fold -> infeasible -> warning + NA
  expect_warning(
    bv <- suppressMessages(nmfkc.bicv(Y, rank = c(2, 6), nfolds = 2)),
    "infeasible")
  expect_false(is.na(bv$sigma[1]))
  expect_true(is.na(bv$sigma[2]))
})

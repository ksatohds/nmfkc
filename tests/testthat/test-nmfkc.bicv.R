## Tests for nmfkc.bicv() (bi-cross-validation prototype)

test_that("nmfkc.bicv returns a sigma per rank and recovers a rank-3 signal", {
  set.seed(1)
  X <- matrix(abs(rnorm(30 * 3)), 30, 3)
  B <- matrix(abs(rnorm(3 * 40)), 3, 40)
  Y <- X %*% B                         # clean rank-3 non-negative matrix
  bv <- suppressMessages(nmfkc.bicv(Y, rank = 1:6, nfolds = 2, seed = 7))

  expect_s3_class(bv, "nmfkc.bicv")
  expect_equal(bv$rank, 1:6)
  expect_length(bv$sigma, 6)
  ## clean rank-3 data: held-out error should be (near) minimal at >= 3
  expect_true(bv$rank.best >= 3)
  ## error should drop sharply from rank 2 to rank 3
  expect_lt(bv$sigma[3], bv$sigma[2])
})

test_that("nmfkc.bicv skips ranks too large for the retained block", {
  set.seed(2)
  Y <- matrix(abs(rnorm(8 * 50)), 8, 50)   # only 8 rows -> nfolds=2 keeps ~4
  bv <- suppressMessages(nmfkc.bicv(Y, rank = c(2, 6), nfolds = 2))
  ## rank 6 leaves <=4 kept rows per fold -> NA; rank 2 is fine
  expect_false(is.na(bv$sigma[1]))
  expect_true(is.na(bv$sigma[2]))
})

## Tests for nmfkc.bicv() (bi-cross-validation, returns an nmf.rank object)

test_that("nmfkc.bicv returns an nmf.rank object with R2 / eff.rank / sigma.bicv", {
  set.seed(1)
  X <- matrix(abs(rnorm(30 * 3)), 30, 3)
  B <- matrix(abs(rnorm(3 * 40)), 3, 40)
  Y <- X %*% B                         # clean rank-3 non-negative matrix
  bv <- suppressMessages(nmfkc.bicv(Y, rank = 1:6, nfolds = 2, seed = 7,
                                    plot = FALSE))

  ## same shape / methods as nmfkc.rank
  expect_s3_class(bv, "nmf.rank")
  expect_equal(bv$rank, 1:6)
  expect_length(bv$sigma, 6)
  ## criteria table carries the requested diagnostics
  expect_true(all(c("rank", "r.squared", "effective.rank",
                    "effective.rank.index", "sigma.bicv") %in%
                  names(bv$criteria)))
  expect_equal(bv$cv.label, "bi-CV")
  ## clean rank-3 data: bi-CV best at >= 3, error drops sharply 2 -> 3
  expect_true(bv$rank.best >= 3)
  expect_lt(bv$criteria$sigma.bicv[3], bv$criteria$sigma.bicv[2])
  ## print / plot run without error
  expect_output(print(bv), "bi-CV min")
})

test_that("nmfkc.bicv skips ranks too large for the retained block", {
  set.seed(2)
  Y <- matrix(abs(rnorm(8 * 50)), 8, 50)   # only 8 rows -> nfolds=2 keeps ~4
  bv <- suppressMessages(nmfkc.bicv(Y, rank = c(2, 6), nfolds = 2, plot = FALSE))
  ## rank 6 leaves <=4 kept rows per fold -> NA; rank 2 is fine
  expect_false(is.na(bv$criteria$sigma.bicv[1]))
  expect_true(is.na(bv$criteria$sigma.bicv[2]))
})

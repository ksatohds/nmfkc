## Tests for nmfkc.ard() (ARD-NMF rank determination prototype)

test_that("nmfkc.ard prunes an over-complete fit toward the true rank", {
  set.seed(1)
  X <- matrix(abs(rnorm(40 * 3)), 40, 3)
  B <- matrix(abs(rnorm(3 * 60)), 3, 60)
  ar <- nmfkc.ard(X %*% B, rank = 10, seed = 1)   # clean rank-3 signal

  expect_s3_class(ar, "nmfkc.ard")
  expect_equal(ar$rank.init, 10)
  expect_length(ar$relevance, 10)
  ## relevance is in [0, 1], descending, first = 1
  expect_true(all(ar$relevance >= 0 & ar$relevance <= 1))
  expect_equal(ar$relevance[1], 1)
  expect_false(is.unsorted(rev(ar$relevance)))
  ## an over-complete fit should prune: estimated rank well below the start
  expect_lt(ar$rank, 10)
  expect_true(ar$rank >= 2 && ar$rank <= 5)   # near the true 3
  ## W, H are ordered/sized to the starting rank
  expect_equal(ncol(ar$W), 10)
  expect_equal(nrow(ar$H), 10)
})

test_that("nmfkc.ard print and plot run without error", {
  set.seed(2)
  Y <- matrix(abs(rnorm(30 * 40)), 30, 40)
  ar <- nmfkc.ard(Y, rank = 8, seed = 2)
  expect_output(print(ar), "estimated rank")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ar))
})

test_that("L1 prior also returns a valid pruned rank", {
  set.seed(3)
  X <- matrix(abs(rnorm(40 * 3)), 40, 3)
  B <- matrix(abs(rnorm(3 * 60)), 3, 60)
  ar <- nmfkc.ard(X %*% B, rank = 10, prior = "L1", seed = 3)
  expect_equal(ar$prior, "L1")
  expect_true(ar$rank >= 1 && ar$rank <= 10)
})

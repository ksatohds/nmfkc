## Tests for nmfkc.consensus() (Brunet 2004 consensus-clustering rank selection)

make_blocks <- function(seed = 1) {
  set.seed(seed)
  X <- matrix(0, 30, 3); X[1:10, 1] <- 1; X[11:20, 2] <- 1; X[21:30, 3] <- 1
  B <- matrix(0, 3, 60); for (j in 1:60) B[((j - 1) %/% 20) + 1, j] <- 1
  X %*% B + matrix(abs(rnorm(30 * 60, 0, 0.1)), 30, 60)   # 3 clear clusters
}

test_that("nmfkc.consensus returns stability scores per rank", {
  Y <- make_blocks()
  cs <- suppressMessages(nmfkc.consensus(Y, rank = 2:5, nrun = 15, seed = 1))

  expect_type(cs, "list")
  expect_equal(cs$rank, 2:5)
  expect_equal(cs$nrun, 15)
  expect_length(cs$cophenetic, 4)
  expect_length(cs$dispersion, 4)
  ## valid ranges
  expect_true(all(cs$dispersion >= 0 & cs$dispersion <= 1, na.rm = TRUE))
  expect_true(all(abs(cs$cophenetic) <= 1, na.rm = TRUE))
  expect_null(cs$consensus)
})

test_that("nmfkc.consensus is most stable at the true rank (3 clusters)", {
  Y <- make_blocks()
  cs <- suppressMessages(nmfkc.consensus(Y, rank = 2:5, nrun = 15, seed = 1))
  ## true rank 3: dispersion ~1 and clearly above rank 2 (under-split)
  k3 <- which(cs$rank == 3)
  k2 <- which(cs$rank == 2)
  expect_gt(cs$dispersion[k3], 0.95)
  expect_gt(cs$dispersion[k3], cs$dispersion[k2])
})

test_that("keep.consensus returns N x N consensus matrices", {
  Y <- make_blocks()
  cs <- suppressMessages(nmfkc.consensus(Y, rank = 2:3, nrun = 10,
                                         keep.consensus = TRUE))
  expect_length(cs$consensus, 2)
  expect_equal(dim(cs$consensus[[1]]), c(ncol(Y), ncol(Y)))
  ## consensus entries are in [0, 1]
  expect_true(all(cs$consensus[[2]] >= 0 & cs$consensus[[2]] <= 1))
})

test_that("workwell", {
  X <- cbind(c(1,0,1),c(0,1,0))
  B <- cbind(c(1,0),c(0,1),c(1,1))
  Y <- X %*% B
  res <- nmfkc::nmfkc(Y,rank=2,epsilon=1e-12,maxit=50000,print.dims = FALSE)
  tolerance <- 1e-3
  expect_equal(Y,res$XB, tolerance = tolerance)
})

test_that("X.init = kmeans++ works and is reproducible", {
  Y <- t(iris[, -5])
  r1 <- nmfkc(Y, rank = 3, X.init = "kmeans++", seed = 1, verbose = FALSE, print.dims = FALSE)
  r2 <- nmfkc(Y, rank = 3, X.init = "kmeans++", seed = 1, verbose = FALSE, print.dims = FALSE)
  expect_true(r1$r.squared > 0.9)
  expect_equal(r1$XB, r2$XB)                       # reproducible with fixed seed
  expect_true(all(r1$X >= 0))                      # basis stays non-negative
  # alias
  r3 <- nmfkc(Y, rank = 3, X.init = "kmeanspp", seed = 1, verbose = FALSE, print.dims = FALSE)
  expect_equal(r1$XB, r3$XB)
})

test_that(".kmeanspp.seed returns k distinct-ish centres of right shape", {
  set.seed(1)
  pts <- matrix(rnorm(200), nrow = 40, ncol = 5)   # 40 points in R^5
  s <- nmfkc:::.kmeanspp.seed(pts, 4)
  expect_equal(dim(s), c(4L, 5L))
  expect_true(all(is.finite(s)))
})

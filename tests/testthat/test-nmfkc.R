test_that("workwell", {
  X <- cbind(c(1,0,1),c(0,1,0))
  B <- cbind(c(1,0),c(0,1),c(1,1))
  Y <- X %*% B
  res <- nmfkc::nmfkc(Y,Q=2,epsilon=1e-12,maxit=50000,print.dims = FALSE)
  tolerance <- 1e-3
  expect_equal(Y,res$XB, tolerance = tolerance)
})

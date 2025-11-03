test_that("workwell", {
  (X <- cbind(c(1,0,1),c(0,1,0)))
  (B <- cbind(c(1,0),c(0,1),c(1,1)))
  (Y <- X %*% B)
  res <- nmfkc::nmfkc(Y,Q=2,epsilon=1e-12,maxit=50000)
  tolerance <- 1e-3
  # Y_reconstructed が Y_original に（ほぼ）等しいか？
  # expect_equal は行列の各要素を比較します
  expect_equal(Y,res$XB, tolerance = tolerance)
})

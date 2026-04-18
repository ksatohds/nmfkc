## Tests for nmfkc.net family (symmetric NMF for networks)

## Small symmetric test network (two communities)
make_test_network <- function() {
  Y <- matrix(c(0, 1, 1, 0, 0, 0,
                1, 0, 1, 0, 0, 0,
                1, 1, 0, 1, 0, 0,
                0, 0, 1, 0, 1, 1,
                0, 0, 0, 1, 0, 1,
                0, 0, 0, 1, 1, 0), 6, 6)
  rownames(Y) <- colnames(Y) <- paste0("N", 1:6)
  Y
}


test_that("nmfkc.net (tri) fits a small symmetric network", {
  Y <- make_test_network()
  res <- nmfkc.net(Y, rank = 2, type = "tri", nstart = 5, maxit = 200)

  expect_s3_class(res, "nmfkc.net.tri")
  expect_s3_class(res, "nmfkc.net")
  expect_s3_class(res, "nmfkc")
  expect_equal(dim(res$X), c(6, 2))
  expect_equal(dim(res$C), c(2, 2))
  expect_true(all(res$X >= 0))
  expect_true(all(res$C >= 0))

  ## C is symmetric by design
  expect_true(isSymmetric(res$C, tol = 1e-8))

  ## X.prob rows sum to 1
  expect_true(all(abs(rowSums(res$X.prob) - 1) < 1e-8))
  expect_length(res$X.cluster, 6)

  expect_true(is.finite(res$r.squared))
})


test_that("nmfkc.net (bi) fits with C = I_Q fixed", {
  Y <- make_test_network()
  res <- nmfkc.net(Y, rank = 2, type = "bi", nstart = 5, maxit = 200)
  expect_s3_class(res, "nmfkc.net.bi")
  expect_equal(unname(res$C), diag(2))
  expect_true(all(res$X >= 0))
})


test_that("nmfkc.net.signed produces signed C", {
  Y <- make_test_network()
  res <- nmfkc.net.signed(Y, rank = 2, nstart = 5, maxit = 200)

  expect_s3_class(res, "nmfkc.net.signed")
  expect_s3_class(res, "nmfkc.net")
  expect_true(all(res$X >= 0))
  expect_true(all(res$Cp >= 0))
  expect_true(all(res$Cn >= 0))

  ## C = Cp - Cn
  expect_equal(res$C, res$Cp - res$Cn, tolerance = 1e-10)
  ## C is symmetric by design
  expect_true(isSymmetric(res$C, tol = 1e-8))
})


test_that("nmfkc.net.DOT works for tri, bi, and signed", {
  Y <- make_test_network()
  res_tri    <- nmfkc.net(Y, rank = 2, type = "tri", nstart = 3, maxit = 100)
  res_bi     <- nmfkc.net(Y, rank = 2, type = "bi",  nstart = 3, maxit = 100)
  res_signed <- nmfkc.net.signed(Y, rank = 2, nstart = 3, maxit = 100)

  dot_tri <- nmfkc.net.DOT(res_tri)
  expect_s3_class(dot_tri, "nmfkc.DOT")
  expect_true(nchar(dot_tri) > 100)

  dot_bi <- nmfkc.net.DOT(res_bi)
  expect_s3_class(dot_bi, "nmfkc.DOT")

  ## Signed should auto-detect and produce a DOT string that references style="dashed"
  dot_signed <- nmfkc.net.DOT(res_signed)
  expect_s3_class(dot_signed, "nmfkc.DOT")
})


test_that("nmfkc.net.ecv runs on small data", {
  Y <- make_test_network()
  out <- nmfkc.net.ecv(Y, rank = c(1, 2), type = "tri",
                       nfolds = 3, nstart = 3, maxit = 100)
  expect_length(out$objfunc, 2)
  expect_length(out$sigma, 2)
  expect_s3_class(out, "nmfkc.net.ecv")
})


test_that("nmfkc(Y.symmetric = 'tri') emits Deprecated warning", {
  Y <- make_test_network()
  expect_warning(
    res <- nmfkc(Y, rank = 2, Y.symmetric = "tri",
                  nstart = 3, maxit = 100, print.dims = FALSE),
    "deprecated"
  )
  ## Still returns a valid object
  expect_s3_class(res, "nmfkc")
})

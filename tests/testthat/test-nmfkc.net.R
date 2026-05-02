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


test_that("nmfkc.net(type='signed') produces signed C", {
  Y <- make_test_network()
  res <- nmfkc.net(Y, rank = 2, type = "signed", nstart = 5, maxit = 200)

  expect_s3_class(res, "nmfkc.net.signed")
  expect_s3_class(res, "nmfkc.net")
  expect_s3_class(res, "nmfkc")
  expect_true(all(res$X >= 0))
  expect_true(all(res$Cp >= 0))
  expect_true(all(res$Cn >= 0))

  ## C = Cp - Cn
  expect_equal(res$C, res$Cp - res$Cn, tolerance = 1e-10)
  ## C is symmetric by design
  expect_true(isSymmetric(res$C, tol = 1e-8))
})


test_that("nmfkc.net return structure is uniform (Cp/Cn NULL for tri/bi)", {
  Y <- make_test_network()
  res_tri    <- nmfkc.net(Y, rank = 2, type = "tri",    nstart = 3, maxit = 100)
  res_bi     <- nmfkc.net(Y, rank = 2, type = "bi",     nstart = 3, maxit = 100)
  res_signed <- nmfkc.net(Y, rank = 2, type = "signed", nstart = 3, maxit = 100)

  ## Uniform field presence
  for (res in list(res_tri, res_bi, res_signed)) {
    expect_true("X"  %in% names(res))
    expect_true("C"  %in% names(res))
    expect_true("Cp" %in% names(res))
    expect_true("Cn" %in% names(res))
  }
  ## Cp, Cn are NULL for tri/bi
  expect_null(res_tri$Cp); expect_null(res_tri$Cn)
  expect_null(res_bi$Cp);  expect_null(res_bi$Cn)
  ## Cp, Cn are matrices for signed
  expect_true(is.matrix(res_signed$Cp))
  expect_true(is.matrix(res_signed$Cn))
})


test_that("nmfkc.net.DOT works for tri, bi, and signed", {
  Y <- make_test_network()
  res_tri    <- nmfkc.net(Y, rank = 2, type = "tri", nstart = 3, maxit = 100)
  res_bi     <- nmfkc.net(Y, rank = 2, type = "bi",  nstart = 3, maxit = 100)
  res_signed <- nmfkc.net(Y, rank = 2, type = "signed", nstart = 3, maxit = 100)

  dot_tri <- nmfkc.net.DOT(res_tri)
  expect_s3_class(dot_tri, "nmfkc.DOT")
  expect_true(nchar(dot_tri) > 100)

  dot_bi <- nmfkc.net.DOT(res_bi)
  expect_s3_class(dot_bi, "nmfkc.DOT")

  ## Signed should auto-detect and produce a DOT string that references style="dashed"
  dot_signed <- nmfkc.net.DOT(res_signed)
  expect_s3_class(dot_signed, "nmfkc.DOT")
})


test_that("nmfkc.net.ecv supports tri/bi/signed via type argument", {
  Y <- make_test_network()

  ## maxit large enough to avoid the "maximum iterations reached"
  ## warning that nmfkc.net's MU loop now emits (commit 1d93c48 et al.).
  out_tri <- nmfkc.net.ecv(Y, rank = c(1, 2), type = "tri",
                            nfolds = 3, nstart = 3, maxit = 5000)
  expect_length(out_tri$objfunc, 2)
  expect_s3_class(out_tri, "nmfkc.net.ecv")
  expect_equal(out_tri$type, "tri")

  out_bi <- nmfkc.net.ecv(Y, rank = c(1, 2), type = "bi",
                           nfolds = 3, nstart = 3, maxit = 5000)
  expect_equal(out_bi$type, "bi")

  out_signed <- nmfkc.net.ecv(Y, rank = c(1, 2), type = "signed",
                               nfolds = 3, nstart = 3, maxit = 5000)
  expect_s3_class(out_signed, "nmfkc.net.signed.ecv")
  expect_s3_class(out_signed, "nmfkc.net.ecv")
  expect_equal(out_signed$type, "signed")
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

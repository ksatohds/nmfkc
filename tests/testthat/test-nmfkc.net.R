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
  skip_unless_full()
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
  skip_unless_full()
  Y <- make_test_network()
  res <- nmfkc.net(Y, rank = 2, type = "bi", nstart = 5, maxit = 200)
  expect_s3_class(res, "nmfkc.net.bi")
  expect_equal(unname(res$C), diag(2))
  expect_true(all(res$X >= 0))
})


test_that("nmfkc.net(type='signed') produces signed C", {
  skip_unless_full()
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
  skip_unless_full()
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
  skip_unless_full()
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
  skip_unless_full()
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


test_that("nmfkc(Y.symmetric = ...) stops and redirects to nmfkc.net()", {
  skip_unless_full()
  Y <- make_test_network()
  ## Symmetric NMF was removed from nmfkc(); passing Y.symmetric now
  ## errors with a message pointing to nmfkc.net().
  expect_error(
    nmfkc(Y, rank = 2, Y.symmetric = "tri",
          nstart = 3, maxit = 100, print.dims = FALSE),
    "nmfkc.net"
  )
  expect_error(
    nmfkc(Y, rank = 2, Y.symmetric = "bi",
          nstart = 3, maxit = 100, print.dims = FALSE),
    "no longer supported"
  )
})

test_that("nmfkc.ecv(Y.symmetric = ...) stops and redirects to nmfkc.net.ecv()", {
  skip_unless_full()
  Y <- make_test_network()
  expect_error(
    nmfkc.ecv(Y, rank = c(1, 2), Y.symmetric = "tri", nfolds = 3),
    "nmfkc.net.ecv"
  )
})

test_that("type='signed' does not collapse X to zero", {
  skip_unless_full()
  ## The default init drew every entry of C0 from U(-1, 1) and symmetrised it,
  ## so with probability (1/2)^3 at Q = 2 the whole matrix came out negative:
  ## Cp was identically 0, C = Cp - Cn wholly negative against a non-negative
  ## Y, the X-step numerator vanished and X collapsed on the first iteration
  ## (X == 0, r.squared = NA, iter = 2).
  set.seed(1)
  Y <- matrix(abs(rnorm(8 * 40)) + 1, 8, 40)
  S <- Y %*% t(Y)
  for (sd in 1:6) for (q in 2:3) {
    f <- suppressWarnings(nmfkc.net(S, rank = q, verbose = FALSE,
                                    type = "signed", seed = sd))
    expect_false(all(f$X == 0))
    expect_true(all(is.finite(f$X)))
    expect_true(is.finite(f$r.squared))
  }
  ## and it fits about as well as the non-negative variant
  a <- suppressWarnings(nmfkc.net(S, rank = 2, verbose = FALSE, type = "tri"))
  b <- suppressWarnings(nmfkc.net(S, rank = 2, verbose = FALSE, type = "signed"))
  expect_gt(b$r.squared, 0.5 * a$r.squared)
})

test_that("X.restriction='fixed' actually holds X fixed in every type", {
  skip_unless_full()
  set.seed(1)
  Y <- matrix(abs(rnorm(8 * 40)) + 1, 8, 40); S <- Y %*% t(Y)
  for (ty in c("tri", "bi", "signed")) {
    X0 <- suppressWarnings(nmfkc.net(S, rank = 2, verbose = FALSE, type = ty))$X
    g  <- suppressWarnings(nmfkc.net(S, rank = 2, verbose = FALSE, type = ty,
                                     X.init = X0, X.restriction = "fixed"))
    expect_equal(g$X, X0, tolerance = 0, info = ty)
  }
})

test_that("X.L2.ortho reaches the signed path", {
  skip_unless_full()
  set.seed(1)
  Y <- matrix(abs(rnorm(8 * 40)) + 1, 8, 40); S <- Y %*% t(Y)
  a <- suppressWarnings(nmfkc.net(S, rank = 2, verbose = FALSE,
                                  type = "signed", X.L2.ortho = 0))
  b <- suppressWarnings(nmfkc.net(S, rank = 2, verbose = FALSE,
                                  type = "signed", X.L2.ortho = 1e4))
  ## it was simply not forwarded, so the option was a silent no-op
  expect_false(isTRUE(all.equal(a$X, b$X, tolerance = 1e-10)))
})

test_that("nmfkc.net.DOT hides nodes on the membership scale, not the raw X scale", {
  skip_unless_full()
  Y <- make_test_network()
  Q <- 2
  count_nodes <- function(dot)
    sum(grepl("^[[:space:]]+Y_[0-9]+ [[]",
              strsplit(as.character(dot), "\n", fixed = TRUE)[[1]]))

  ## The raw basis X carries the scale of X.restriction, which differs by type
  ## ("none" for bi, "colSums" for tri), so a threshold applied to it meant
  ## different things per type -- for tri it could hide every node at once.
  ## The filter is documented as "no X edge above threshold", and the X edges
  ## are drawn from X.prob, so X.prob is what must be tested.
  for (ty in c("tri", "bi")) {
    fit <- nmfkc.net(Y, rank = Q, type = ty, nstart = 3, maxit = 100)

    ## Row sums of X.prob are one, so a row maximum is at least 1/Q and
    ## nothing can be hidden below that -- correctly, since every node then
    ## does have an edge.
    expect_equal(unname(rowSums(fit$X.prob)), rep(1, nrow(fit$X.prob)),
                 tolerance = 1e-8)
    th <- 1 / Q - 1e-6
    expect_identical(count_nodes(nmfkc.net.DOT(fit, threshold = th,
                                               hide.isolated = TRUE)),
                     count_nodes(nmfkc.net.DOT(fit, threshold = th,
                                               hide.isolated = FALSE)))

    ## Every displayed node has at least one X edge: the filter and the edge
    ## loop now read the same matrix.
    kept <- which(apply(fit$X.prob, 1L, function(r) any(r >= th)))
    expect_identical(length(kept), nrow(fit$X.prob))
  }
})

test_that("summary of nmfkc.net.inference reports the model type, not 'unknown'", {
  skip_unless_full()
  Y <- make_test_network()
  ## $Y.symmetric is the pre-0.9.x name and was removed from nmfkc(); the
  ## variant now lives in $type.  The summary has to carry it across, or the
  ## printer falls through to "unknown" for every nmfkc.net object.
  for (ty in c("tri", "bi", "signed")) {
    fit <- nmfkc.net(Y, rank = 2, type = ty, nstart = 3, maxit = 100)
    inf <- nmfkc.net.inference(fit, Y)
    expect_identical(summary(inf)$type, ty)
    out <- capture.output(print(inf))
    expect_true(any(out == paste0("Symmetric NMF type: ", ty)))
    expect_false(any(grepl("Symmetric NMF type: unknown", out, fixed = TRUE)))
  }
})

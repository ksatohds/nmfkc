## nmfkc.ident(): anchor test, identification region, repair and profile intervals

make_ident_data <- function(anchored = TRUE, seed = 1, sd = 0.1, n = 30) {
  set.seed(seed)
  X <- if (anchored) cbind(c(0.40, 0.30, 0.20, 0.10, 0.00, 0.00),
                           c(0.00, 0.00, 0.10, 0.20, 0.30, 0.40))
       else cbind(c(0.30, 0.25, 0.20, 0.10, 0.10, 0.05),
                  c(0.05, 0.10, 0.10, 0.20, 0.25, 0.30))
  Th <- if (anchored) diag(c(10, 8)) else matrix(c(7, 3, 2, 6), 2, 2)
  A <- rbind(g1 = rep(c(1, 0), each = n), g2 = rep(c(0, 1), each = n))
  ## multiplicative noise keeps the zeros of X Theta A at exactly zero (additive noise clipped at 0 would make
  ## them significantly positive, and the anchor test would rightly reject)
  Y <- (X %*% Th %*% A) * exp(matrix(stats::rnorm(6 * 2 * n, sd = sd), 6, 2 * n))
  rownames(Y) <- paste0("v", 1:6)
  list(Y = Y, A = A)
}

test_that("anchored data: the test does not reject and the region is narrow", {
  skip_unless_full()
  d <- make_ident_data(TRUE)
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(f, d$Y, d$A, region.n = 5000)
  expect_s3_class(id, "nmfkc.ident")
  expect_false(id$test$reject)
  expect_gt(id$test$pi.max, 0.9)
  expect_true(id$identified)                         # exact zeros in Xi.hat: the region is a point
  expect_lt(max(id$C.region.upper - id$C.region.lower), 1e-6 * max(colSums(id$Xi.hat)))
  expect_true(all(id$C >= id$C.region.lower - 1e-6 & id$C <= id$C.region.upper + 1e-6))
})

test_that("mixed data: the region has width and the repaired fit is exact", {
  skip_unless_full()
  d <- make_ident_data(FALSE)
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(f, d$Y, d$A, region.n = 5000)
  expect_false(id$identified)
  expect_gt(max(id$C.region.upper - id$C.region.lower), 0.1)
  ## the repaired fit reproduces the maximum exactly and keeps object's basis names
  expect_lt(id$repair$gap, 1e-6)
  expect_s3_class(id$fit, "nmfkc")
  expect_identical(dimnames(id$C), dimnames(f$C))
  expect_true(all(id$C >= id$C.region.lower - 1e-6 & id$C <= id$C.region.upper + 1e-6))
  ## X stays in its region, columns sum to 1
  expect_equal(unname(colSums(id$X)), c(1, 1), tolerance = 1e-8)
  expect_true(all(id$X >= id$X.region.lower - 1e-6 & id$X <= id$X.region.upper + 1e-6))
  ## start = "center" also lands in the region
  idc <- nmfkc.ident(f, d$Y, d$A, region.n = 5000, start = "center")
  expect_true(all(idc$C >= idc$C.region.lower - 1e-6 & idc$C <= idc$C.region.upper + 1e-6))
})

test_that("profile intervals contain the identification region", {
  skip_unless_full()
  d <- make_ident_data(FALSE)
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(f, d$Y, d$A, region.n = 5000, ci = TRUE)
  expect_true(all(id$C.ci.lower <= id$C.region.lower + 1e-6))
  expect_true(all(id$C.ci.upper >= id$C.region.upper - 1e-6))
  expect_true(all(is.finite(id$C.ci.upper)))
  expect_lt(id$ci.dip, 1e-6)
  s <- summary(id)
  expect_true(all(c("CI_low", "CI_high", "region_low", "region_high") %in% names(s$coefficients)))
})

test_that("more covariates than bases: only the test is returned", {
  skip_unless_full()
  set.seed(2)
  A <- rbind(g1 = rep(c(1, 0, 0), each = 15), g2 = rep(c(0, 1, 0), each = 15), g3 = rep(c(0, 0, 1), each = 15))
  X <- cbind(c(0.4, 0.3, 0.2, 0.1, 0), c(0, 0.1, 0.2, 0.3, 0.4))
  Y <- pmax(X %*% matrix(c(8, 0, 4, 4, 0, 8), 2, 3) %*% A + matrix(rnorm(225, sd = 0.2), 5, 45), 0)
  f <- nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(f, Y, A)
  expect_null(id$Xi.hat)
  expect_true(is.logical(id$test$reject))
  expect_identical(id$fit, f)
})

test_that("the caller's random stream is restored", {
  skip_unless_full()
  d <- make_ident_data(FALSE)
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  set.seed(99); a <- stats::runif(3)
  set.seed(99); invisible(nmfkc.ident(f, d$Y, d$A, region.n = 500, region.nstart = 1)); b <- stats::runif(3)
  expect_identical(a, b)
})

test_that("input checks", {
  skip_unless_full()
  d <- make_ident_data(TRUE)
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  expect_error(nmfkc.ident(list(), d$Y, d$A), "class 'nmfkc'")
  expect_error(nmfkc.ident(f, d$Y, d$A[, -1]), "ncol")
  f3 <- nmfkc(d$Y, rbind(d$A, d$A[1, ] + d$A[2, ]), rank = 2, epsilon = 1e-4, verbose = FALSE)
  expect_error(nmfkc.ident(f3, d$Y, rbind(d$A, d$A[1, ] + d$A[2, ])), "full row rank")
})

test_that("region = FALSE returns the test only", {
  skip_unless_full()
  d <- make_ident_data(TRUE)
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  a <- nmfkc.ident(f, d$Y, d$A, region = FALSE)
  b <- nmfkc.ident(f, d$Y, d$A, region.n = 500, region.nstart = 1)
  expect_null(a$Xi.hat)
  expect_identical(a$test, b$test)                    # the test does not depend on the region
})

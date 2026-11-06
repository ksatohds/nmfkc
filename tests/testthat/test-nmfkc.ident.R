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

test_that("summary(), print() and plot() of nmfkc.ident", {
  skip_unless_full()
  d <- make_ident_data(FALSE)
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(f, d$Y, d$A, region.n = 2000, ci = TRUE)
  s <- summary(id)
  expect_s3_class(s, "summary.nmfkc.ident")
  expect_true(all(c("Basis", "Covariate", "Estimate", "region_low", "region_high", "CI_low", "CI_high", "nonzero")
                  %in% names(s$coefficients)))
  expect_identical(s$coefficients$nonzero, s$coefficients$CI_low > 0)
  expect_true(s$region.share >= 0 && s$region.share <= 1)
  expect_output(print(id), "Anchor test")
  expect_output(print(id), "Identification region of C")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(id))
  expect_silent(plot(id, type = "X"))
  ## test only: summary says so and plot() refuses
  a <- nmfkc.ident(f, d$Y, d$A, region = FALSE)
  expect_output(print(a), "not computed")
  expect_error(plot(a), "no identification region")
})

test_that("H(rho) tests the entries of Xi.hat one by one; pi.max = 1 exactly when H(rho) is not rejected", {
  skip_unless_full()
  for (anch in c(TRUE, FALSE)) {
    d <- make_ident_data(anch, seed = 3)
    f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
    id <- nmfkc.ident(f, d$Y, d$A, region = FALSE)
    ## direct implementation for Q = R = 2: rows j1 != j2, columns in both orders
    A <- d$A; Y <- d$Y; N <- ncol(Y); R <- nrow(A); Q <- 2; rho <- 0.5 / nrow(Y)
    W <- solve(A %*% t(A), A); Xi <- Y %*% t(W); E <- Y - Xi %*% A; hc <- N / (N - R); s <- colSums(Xi)
    se <- sapply(1:R, function(l) sqrt(rowSums(sweep(E, 2, W[l, ], "*")^2) * hc))
    tE <- ifelse(se > 0, Xi / se, ifelse(Xi > 0, Inf, ifelse(Xi < 0, -Inf, 0)))   # exact zeros: t = 0
    tS <- function(j, l) { psi <- E[j, ] * W[l, ] - rho * colSums(E) * W[l, ]
                           g <- Xi[j, l] - rho * s[l]; v <- sum(psi^2) * hc
                           if (v > 0) g / sqrt(v) else if (g >= 0) Inf else -Inf }
    cv <- qt(1 - 0.05 / Q^2, N - R); survive <- FALSE
    for (L in list(c(1, 2), c(2, 1))) for (j1 in 1:nrow(Y)) for (j2 in 1:nrow(Y)) if (j1 != j2)
      if (tE[j1, L[2]] <= cv && tE[j2, L[1]] <= cv && tS(j1, L[1]) >= -cv && tS(j2, L[2]) >= -cv) survive <- TRUE
    expect_identical(id$test$reject, !survive)
    expect_equal(id$test$crit, cv)
    expect_identical(id$test$pi.max == 1, !id$test$reject)
  }
})

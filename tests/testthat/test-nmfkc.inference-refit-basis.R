## nmfkc.inference(method = "refit.basis"): the residual wild bootstrap of
## method = "refit" with the basis re-estimated in every replicate and its bases
## matched to the original ones by permutation.

make_rb_data <- function(anchored = TRUE, seed = 1, n = 30) {
  set.seed(seed)
  X <- if (anchored) cbind(c(0.40, 0.30, 0.20, 0.10, 0.00, 0.00),
                           c(0.00, 0.00, 0.10, 0.20, 0.30, 0.40))
       else cbind(c(0.30, 0.25, 0.20, 0.10, 0.10, 0.05),
                  c(0.05, 0.10, 0.10, 0.20, 0.25, 0.30))
  Th <- if (anchored) diag(c(10, 8)) else matrix(c(7, 3, 2, 6), 2, 2)
  A <- rbind(g1 = rep(c(1, 0), each = n), g2 = rep(c(0, 1), each = n))
  Y <- (X %*% Th %*% A) * exp(matrix(stats::rnorm(6 * 2 * n, sd = 0.1), 6, 2 * n))
  rownames(Y) <- paste0("v", 1:6)
  list(Y = Y, A = A)
}

test_that("refit.basis equals a direct loop: same multipliers, re-fit, permutation, basic p", {
  skip_unless_full()
  d <- make_rb_data()
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  B <- 15
  r <- nmfkc.inference(f, d$Y, d$A, method = "refit.basis", wild.B = B, wild.seed = 7)
  Yhat <- f$X %*% f$C %*% d$A; E <- d$Y - Yhat
  P <- nrow(d$Y); N <- ncol(d$Y)
  perms <- rbind(c(1, 2), c(2, 1))
  set.seed(7)
  D <- sapply(seq_len(B), function(b) {
    W <- matrix(sample(c(-1, 1), N, replace = TRUE), P, N, byrow = TRUE)
    Ys <- pmax(Yhat + W * E, 0)
    g <- nmfkc(Ys, d$A, rank = 2, epsilon = 1e-8, maxit = 1e5, detail = "fast", verbose = FALSE,
               method = f$method, X.restriction = f$X.restriction)
    pm <- perms[which.min(apply(perms, 1, function(p) sum((g$X[, p] - f$X)^2))), ]
    as.vector(g$C[pm, ])
  })
  expect_identical(r$C.boot.draws, D)
  est <- as.vector(f$C)
  p1 <- sapply(seq_along(est), function(i) mean(D[i, ] >= 2 * est[i]))
  expect_identical(r$coefficients$p_value, pmin(1, p1))
  expect_identical(r$coefficients$BSE, apply(D, 1, sd))
  expect_identical(r$boot.method, "refit.basis")
  expect_identical(r$wild.unit, "column")
  expect_identical(r$refit.info$n.fail, 0L)
  expect_equal(sum(r$refit.info$perm), B)
  expect_true(r$refit.info$truncate.rate >= 0 && r$refit.info$truncate.rate <= 1)
  expect_output(print(summary(r)), "basis re-estimated")
  ## the sandwich SE does not depend on the method
  r0 <- nmfkc.inference(f, d$Y, d$A, wild.bootstrap = FALSE)
  expect_identical(r$coefficients$SE, r0$coefficients$SE)
})

test_that("refit.basis: bases given in the other order are matched back to that order", {
  skip_unless_full()
  d <- make_rb_data()
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  fs <- f; fs$X <- f$X[, 2:1]; fs$C <- f$C[2:1, ]; fs$B <- f$B[2:1, ]
  colnames(fs$X) <- rownames(fs$C) <- rownames(fs$B) <- colnames(f$X)
  a <- nmfkc.inference(f, d$Y, d$A, method = "refit.basis", wild.B = 10)
  b <- nmfkc.inference(fs, d$Y, d$A, method = "refit.basis", wild.B = 10)
  ## the same replicates, rows of vec(C) swapped (equal to rounding: X C A is summed in
  ## another order, so Y* differs in the last bit)
  sw <- as.vector(matrix(seq_len(4), 2)[2:1, ])
  expect_equal(b$C.boot.draws, a$C.boot.draws[sw, ], tolerance = 1e-8)
  expect_equal(unname(b$refit.info$perm["2-1"]), unname(a$refit.info$perm["1-2"]))
})

test_that("refit.basis: warnings for arguments of the call and for a non-unique factorization", {
  skip_unless_full()
  d <- make_rb_data()
  f <- nmfkc(d$Y, d$A, rank = 2, X.init = "nndsvd", epsilon = 1e-8, verbose = FALSE)
  expect_warning(nmfkc.inference(f, d$Y, d$A, method = "refit.basis", wild.B = 3), "X.init")
  r <- expect_no_warning(nmfkc.inference(f, d$Y, d$A, method = "refit.basis", wild.B = 3,
                                         refit.args = list(X.init = "nndsvd")))
  expect_identical(r$refit.info$refit.args$X.init, "nndsvd")
  dm <- make_rb_data(anchored = FALSE)
  fm <- nmfkc(dm$Y, dm$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  rm <- nmfkc.inference(fm, dm$Y, dm$A, wild.bootstrap = FALSE)
  if (!isTRUE(rm$identified.approx))
    expect_warning(nmfkc.inference(fm, dm$Y, dm$A, method = "refit.basis", wild.B = 3), "descriptive")
})

test_that("refit.basis without covariates: C = B, one column per sample", {
  skip_unless_full()
  d <- make_rb_data()
  f <- nmfkc(d$Y, rank = 2, epsilon = 1e-8, verbose = FALSE)
  r <- suppressWarnings(nmfkc.inference(f, d$Y, method = "refit.basis", wild.B = 5))
  expect_identical(dim(r$C.boot.draws), c(2L * ncol(d$Y), 5L))
  expect_identical(nrow(r$coefficients), 2L * ncol(d$Y))
})

test_that("cores = 2 gives the same draws as cores = 1 (refit and refit.basis)", {
  skip_unless_full()
  skip_on_os("windows")     # PSOCK workers would load the installed nmfkc
  d <- make_rb_data()
  f <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  for (m in c("refit", "refit.basis")) for (u in c("column", "element", "block")) {
    a <- nmfkc.inference(f, d$Y, d$A, method = m, wild.unit = u, wild.B = 8, cores = 1)
    b <- nmfkc.inference(f, d$Y, d$A, method = m, wild.unit = u, wild.B = 8, cores = 2)
    expect_identical(a$C.boot.draws, b$C.boot.draws)
    expect_identical(a$coefficients, b$coefficients)
  }
})

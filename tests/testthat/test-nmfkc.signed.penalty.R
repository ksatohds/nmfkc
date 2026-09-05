## Lasso penalty on the signed coefficient matrix of nmfkc.signed().
##
## With the split C = Cp - Cn the penalty C.L1 * sum(Cp + Cn) equals
## C.L1 * ||C||_1 at any minimizer, because a cell in which both parts were
## positive can be reduced in both without changing C.  Its gradient is the
## constant C.L1 in each factor, so it enters only the denominators.
##
## The joint problem in (X, C) is not convex, so the path in C.L1 is checked
## with X held fixed, where the problem in C alone IS convex and the solution
## is characterised by the lasso subgradient conditions.

make_signed_case <- function(seed = 1, p = 3, N = 200, P = 4, D = 30) {
  set.seed(seed)
  U <- matrix(stats::rnorm(p * N), p, N)
  lab <- sample.int(P, N, replace = TRUE)
  Y <- matrix(0, P, N); Y[cbind(lab, seq_len(N))] <- 1
  list(Y = Y, U = U, Z = nmfkc.signed.rff(U, beta = 0.5, D = D, seed = seed)$Z, P = P)
}

fit_fixedX <- function(cs, lam, maxit = 200000L) {
  nmfkc.signed(cs$Y, A = cs$Z, rank = cs$P, epsilon = 1e-10, maxit = maxit,
               verbose = FALSE, warm.start = FALSE, X.init = diag(cs$P),
               X.restriction = "fixed", C.L1 = lam)
}

test_that("C.L1 = 0 leaves the fit unchanged", {
  cs <- make_signed_case()
  a <- fit_fixedX(cs, 0, maxit = 3000L)
  b <- nmfkc.signed(cs$Y, A = cs$Z, rank = cs$P, epsilon = 1e-10, maxit = 3000L,
                    verbose = FALSE, warm.start = FALSE, X.init = diag(cs$P),
                    X.restriction = "fixed")
  expect_equal(a$C, b$C)
})

test_that("with X fixed the lasso path is monotone in both directions", {
  cs <- make_signed_case()
  lams <- c(0, 0.03, 0.1, 0.3, 1)
  l1  <- numeric(length(lams)); sse <- numeric(length(lams))
  for (i in seq_along(lams)) {
    f <- fit_fixedX(cs, lams[i])
    l1[i]  <- sum(abs(f$C))
    sse[i] <- sum((cs$Y - f$X %*% f$C %*% cs$Z)^2)
  }
  expect_true(all(diff(l1)  <= 1e-4))   # more shrinkage as lambda grows
  expect_true(all(diff(sse) >= -1e-4))  # and a worse fit: a penalty, not a free lunch
})

test_that("the returned C satisfies the lasso subgradient conditions", {
  cs <- make_signed_case(); lam <- 0.3
  f <- fit_fixedX(cs, lam)
  G <- 2 * (f$C %*% cs$Z - cs$Y) %*% t(cs$Z)   # gradient of the squared loss at X = I
  nz <- abs(f$C) > 1e-6
  expect_lt(max(abs(G[nz] + lam * sign(f$C[nz]))), 1e-2)  # stationary where C != 0
  expect_lt(max(abs(G[!nz])), lam * 1.01)                 # subgradient where C = 0
})

test_that("Cp and Cn end with disjoint supports, so sum(Cp + Cn) is ||C||_1", {
  cs <- make_signed_case()
  f <- fit_fixedX(cs, 0.2, maxit = 20000L)
  expect_equal(sum(pmin(pmax(f$C, 0), pmax(-f$C, 0))), 0)
})

test_that("C.L1 acts identically on the Gram route", {
  cs <- make_signed_case()
  g <- nmfkc.signed.rff.gram(cs$Y, cs$U, beta = 0.5, D = 30L, seed = 1)
  fg <- nmfkc.signed(cs$Y, A = g, rank = cs$P, epsilon = 1e-8, maxit = 5000L,
                     verbose = FALSE, warm.start = FALSE, X.init = diag(cs$P), C.L1 = 0.05)
  fd <- nmfkc.signed(cs$Y, A = cs$Z, rank = cs$P, epsilon = 1e-8, maxit = 5000L,
                     verbose = FALSE, warm.start = FALSE, X.init = diag(cs$P), C.L1 = 0.05)
  expect_lt(max(abs(fg$X %*% fg$C - fd$X %*% fd$C)), 1e-6)
})

test_that("C.L1 combines with C.L2 and with a free X", {
  cs <- make_signed_case()
  f <- nmfkc.signed(cs$Y, A = cs$Z, rank = cs$P, epsilon = 1e-8, maxit = 3000L,
                    verbose = FALSE, warm.start = FALSE, C.L1 = 0.05, C.L2 = 0.05)
  expect_true(all(is.finite(f$C)))
  expect_true(is.finite(f$objfunc))
})

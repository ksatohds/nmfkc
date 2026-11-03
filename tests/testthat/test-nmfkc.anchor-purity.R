## summary.nmfkc(): the anchor purity of X and of C, counted by the tolerances
## of the identification verdict, and that verdict (.ident.verdict, the check
## of nmfkc.inference() and the other inference functions).

anc <- function(...) nmfkc:::.nmfkc.anchor(...)
vrd <- function(X, C, a = anc(X, C), xr = "colSums", ...)
  nmfkc:::.nmfkc.ident.summary(X, C, a, xr, ...)
colnorm <- function(X) sweep(X, 2, colSums(X), "/")

test_that("anchors on both sides: every basis anchored, unique", {
  skip_unless_full()
  X <- colnorm(rbind(diag(3), c(0.2, 0.3, 0.5), c(0.6, 0.1, 0.3)))
  C <- cbind(diag(c(2, 1, 4)), c(1, 1, 1), c(0.5, 2, 1))
  a <- anc(X, C)
  expect_true(all(a$X.anchor)); expect_true(all(a$C.anchor)); expect_true(a$anchor.both)
  expect_equal(unname(a$X.anchor.purity), rep(1, 3))
  expect_equal(unname(a$C.anchor.col), 1:3)
  expect_true(vrd(X, C)$identified.approx)
})

test_that("anchors in X alone are not enough; pinning finds what anchors miss", {
  skip_unless_full()
  ## X = I_2, Theta = [1 1; 1 2]: T_t = [1-t 0; t 1] is admissible for 0 < t < 1/2
  a <- anc(diag(2), matrix(c(1, 1, 1, 2), 2, 2))
  expect_true(all(a$X.anchor)); expect_false(any(a$C.anchor))
  expect_equal(unname(a$C.anchor.purity), c(1 / 2, 2 / 3))
  expect_false(vrd(diag(2), matrix(c(1, 1, 1, 2), 2, 2))$identified.approx)
  ## basis 3 has no anchor column, yet it is pinned by the columns (1,0,1) and
  ## (0,1,1) once bases 1 and 2 are: unique, though not anchored on both sides
  X <- colnorm(rbind(diag(3), c(0.2, 0.3, 0.5)))
  C <- cbind(c(1, 0, 0), c(0, 1, 0), c(1, 0, 1), c(0, 1, 1))
  b <- anc(X, C)
  expect_identical(unname(b$C.anchor), c(TRUE, TRUE, FALSE))
  expect_false(b$anchor.both)
  expect_true(vrd(X, C)$identified.approx)
})

test_that("an all-positive X (Orthodont, AJMMS 2026) has no anchor", {
  skip_unless_full()
  Xh <- matrix(c(0.284, 0.317, 0.227, 0.172, 0.175, 0.160, 0.287, 0.377), 4, 2)
  a <- anc(Xh, matrix(c(49, 41, 0, 8), 2, 2))
  expect_false(any(a$X.anchor))
  expect_identical(unname(a$C.anchor), c(FALSE, TRUE))     # the male column anchors basis 2
  expect_false(vrd(Xh, matrix(c(49, 41, 0, 8), 2, 2), a)$identified.approx)
})

test_that("the counts use the verdict's tolerances and agree with it", {
  skip_unless_full()
  X <- rbind(c(0.9995, 0.0005), c(0.01, 0.99), c(0.5, 0.5))
  expect_identical(unname(anc(X, diag(2))$X.anchor), c(TRUE, FALSE))   # 0.999 boundary
  C <- cbind(c(1, 0.005), c(0.02, 1), c(5e-4, 3))
  expect_identical(unname(anc(diag(2), C)$C.anchor), c(TRUE, TRUE))   # 0.5%, and 5e-4 <= 1e-3
  expect_identical(unname(anc(diag(2), C[, 1:2])$C.anchor), c(TRUE, FALSE))
  ## on random near-sparse fits the counts are those of .nmfvar.identified()
  for (s in 1:20) {
    set.seed(s)
    Q <- sample(2:5, 1); P <- 3 * Q; K <- 2 * Q
    X <- matrix(stats::runif(P * Q), P, Q) * (matrix(stats::runif(P * Q), P, Q) > 0.6)
    X[seq_len(Q), ] <- diag(Q) + 1e-4 * stats::runif(Q)
    X <- X[rowSums(X) > 0, , drop = FALSE]
    C <- matrix(stats::runif(Q * K), Q, K) * (matrix(stats::runif(Q * K), Q, K) > 0.5) +
         5e-4 * (matrix(stats::runif(Q * K), Q, K) > 0.5)
    id <- nmfkc:::.nmfvar.identified(X, C)
    a <- anc(X, C)
    expect_identical(unname(a$X.anchor), id$anchor)
    expect_identical(unname(a$C.anchor), unname(id$pure))
  }
})

test_that("the verdict follows nmfkc.inference(): scale, X held, rank of A, large Q", {
  skip_unless_full()
  X <- colnorm(rbind(diag(2), c(1, 1)))
  C <- cbind(diag(2), c(1, 1))
  expect_true(vrd(X, C)$identified.approx)
  expect_false(vrd(X, C, xr = "none")$identified.approx)       # column scale free
  expect_true(vrd(diag(2), matrix(1, 2, 2), xr = "fixed")$identified.approx)
  a <- anc(X, C, A = rbind(1, 1, stats::runif(5)))             # A rank deficient
  expect_false(a$A.full.rank)
  expect_false(vrd(X, C, a)$identified.approx)
  ## above max.work only anchors on both sides are decided
  expect_true(vrd(X, C, max.work = 1)$identified.approx)
  big <- vrd(diag(2), matrix(c(1, 1, 1, 2), 2, 2), max.work = 1)
  expect_true(is.na(big$identified.approx))
  expect_match(big$identification$note, "not assessed")
})

test_that("negative entries, an all-zero fit, a zero row of X", {
  skip_unless_full()
  s <- anc(diag(2), matrix(c(1, -1, 0, 1), 2, 2))
  expect_true(all(is.na(s$C.anchor)))
  expect_true(is.na(vrd(diag(2), matrix(c(1, -1, 0, 1), 2, 2), s)$identified.approx))
  f <- suppressWarnings(nmfkc(matrix(0, 30, 10), rank = 2,
                              X.restriction = "none", verbose = FALSE))
  expect_false(f$criterion$anchors$anchor.both)
  expect_output(print(summary(f)), "Anchor purity")
  ## a zero row no longer stops the rows of X from being read as shares
  Xz <- rbind(diag(2), 0, c(0.5, 0.5))
  expect_true(all(nmfkc:::.nmfvar.identified(Xz, diag(2))$anchor))
})

test_that("summary() prints the purity lines and one Identification line", {
  skip_unless_full()
  set.seed(7)
  P <- 20; N <- 60
  Xt <- matrix(stats::runif(P * 2, 0.2, 1), P, 2)
  Xt[3, ] <- c(1, 0); Xt[11, ] <- c(0, 1)
  A <- rbind(intercept = 1, dose = stats::runif(N))
  ## Theta with structural zeros (kept exactly by the multiplicative updates)
  Y <- Xt %*% diag(c(2, 3)) %*% A
  fz <- nmfkc(Y, A, rank = 2, X.anchor = c(3, 11), C.init = diag(2),
              epsilon = 1e-10, maxit = 1e5, verbose = FALSE)
  sz <- summary(fz)
  expect_true(sz$anchors$anchor.both)
  expect_true(sz$identified.approx)
  out <- utils::capture.output(print(sz))
  expect_true(any(grepl("Anchor purity \\(X\\): +min 1\\.0000; anchors 2/2", out)))
  expect_true(any(grepl("Anchor purity \\(C\\): +min 1\\.0000; anchors 2/2", out)))
  expect_true(any(grepl("Identification: unique up to a permutation", out)))
  expect_false(any(grepl("Separability|Two-sided", out)))
  expect_null(fz$criterion$separability)                       # removed in 1.0.3
  ## a dense Theta: anchors in X only, not identified
  fa <- nmfkc(Xt %*% matrix(c(2, 1, 0.5, 3), 2, 2) %*% A, A, rank = 2,
              X.anchor = "spa", epsilon = 1e-10, maxit = 1e5, verbose = FALSE)
  sa <- summary(fa)
  expect_identical(sa$anchors, fa$criterion$anchors)
  expect_true(all(sa$anchors$X.anchor)); expect_false(sa$identified.approx)
  expect_true(any(grepl("Identification: not shown", utils::capture.output(print(sa)))))
  ## nmfkc.inference(): the same verdict, printed once, in its own section
  inf <- nmfkc.inference(fa, Xt %*% matrix(c(2, 1, 0.5, 3), 2, 2) %*% A, A,
                         wild.bootstrap = FALSE)
  expect_identical(summary(inf)$identified.approx, sa$identified.approx)
  oi <- utils::capture.output(print(summary(inf)))
  expect_identical(sum(grepl("Identification:", oi)), 1L)
  expect_true(any(grepl("Anchor purity \\(X\\)", oi)))
  ## an object saved before 1.0.3: recomputed, the rank of A unknown
  old <- fa; old$criterion$anchors <- NULL
  so <- summary(old)
  expect_true(is.na(so$anchors$A.full.rank))
  expect_identical(so$anchors$X.anchor, sa$anchors$X.anchor)
  expect_identical(so$identified.approx, sa$identified.approx)
  ## rank 1: nothing printed
  s1 <- summary(nmfkc(Y, A, rank = 1, verbose = FALSE))
  expect_false(any(grepl("Anchor purity|Identification", utils::capture.output(print(s1)))))
})

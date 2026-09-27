## The separability shown by summary.nmfkc(): for each basis the largest share
## it takes in any row of the column-normalized X, minimized over bases.

test_that("separability is 1 for an anchored basis and 1/Q for equal mixing", {
  skip_unless_full()
  sep <- nmfkc:::.nmfkc.separability
  X <- rbind(c(1, 0), c(0, 1), c(0.5, 0.5))
  colnames(X) <- c("Basis1", "Basis2")
  s <- sep(X)
  expect_equal(s$index, 1)
  expect_equal(unname(s$row), c(1L, 2L))
  ## every row mixes both bases in equal shares
  expect_equal(sep(matrix(1, 4, 2))$index, 1 / 2)
  expect_equal(sep(matrix(1, 4, 3))$index, 1 / 3)
})

test_that("separability does not depend on the column scale of X", {
  skip_unless_full()
  ## so it reads the same whatever X.restriction was used
  sep <- nmfkc:::.nmfkc.separability
  set.seed(3)
  X <- matrix(stats::runif(30), 10, 3)
  expect_equal(sep(X)$index, sep(X %*% diag(c(2, 7, 0.1)))$index)
  expect_equal(sep(X)$purity, sep(X %*% diag(c(2, 7, 0.1)))$purity)
})

test_that("a dead basis is not anchored, and rank 1 has nothing to separate", {
  skip_unless_full()
  sep <- nmfkc:::.nmfkc.separability
  X <- cbind(c(1, 2, 3), 0)
  expect_equal(unname(sep(X)$purity), c(1, 0))
  expect_equal(sep(X)$index, 0)
  expect_true(is.na(sep(matrix(1:3, 3, 1))$index))
  ## an all-zero X (all-zero Y under X.restriction = "none" returns one):
  ## no basis appears anywhere, and computing this must not break the fit
  z <- sep(matrix(0, 4, 2))
  expect_equal(unname(z$purity), c(0, 0))
  expect_equal(z$index, 0)
  expect_true(all(is.na(z$row)))
  f <- suppressWarnings(nmfkc(matrix(0, 30, 10), rank = 2,
                              X.restriction = "none", verbose = FALSE))
  expect_equal(f$criterion$separability, 0)
})

test_that("summary() reports it, and an X.anchor fit scores exactly 1", {
  skip_unless_full()
  set.seed(7)
  P <- 20; N <- 60
  Xt <- matrix(stats::runif(P * 2, 0.2, 1), P, 2)
  Xt[3, ] <- c(1, 0); Xt[11, ] <- c(0, 1)
  A <- rbind(intercept = 1, dose = stats::runif(N))
  Y <- Xt %*% matrix(c(2, 1, 0.5, 3), 2, 2) %*% A
  fa <- nmfkc(Y, A, rank = 2, X.anchor = "spa", epsilon = 1e-10, maxit = 1e5,
              verbose = FALSE)
  s <- summary(fa)
  expect_equal(s$separability, 1)
  expect_setequal(unname(s$anchor.row), unname(fa$X.anchor))
  ## the fit carries the same values in its criterion list
  expect_identical(fa$criterion$separability, s$separability)
  expect_identical(fa$criterion$anchor.purity, s$anchor.purity)
  expect_identical(fa$criterion$anchor.row, s$anchor.row)
  ## an object saved before 1.0.0 has none of them: summary() recomputes
  old <- fa
  old$criterion[c("separability", "anchor.purity", "anchor.row")] <- NULL
  expect_equal(summary(old)$separability, s$separability)
  expect_equal(summary(old)$anchor.row, s$anchor.row)
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("Separability:", out)))
  ## rank 1: NA, and no line printed
  s1 <- summary(nmfkc(Y, A, rank = 1, verbose = FALSE))
  expect_true(is.na(s1$separability))
  expect_false(any(grepl("Separability:", utils::capture.output(print(s1)))))
})

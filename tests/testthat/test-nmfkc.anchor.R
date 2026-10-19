## Anchor rows of X: X.init = "spa" (a start) and X.anchor (a constraint).
##
## An anchor row carries one basis alone.  Imposed as zeros, anchors make the
## factorization unique up to column scale -- another solution XR with the
## same zeros needs a diagonal R -- and, with covariates, identify C as well.
## The data below are separable by construction, so the truth is known and
## every claim can be checked against it.

anchor_data <- function(seed = 7, P = 20, N = 60, noise = 0) {
  set.seed(seed)
  Xt <- matrix(stats::runif(P * 2, 0.2, 1), P, 2)
  Xt[3, ] <- c(1, 0); Xt[11, ] <- c(0, 1)             # the anchors
  Xt <- sweep(Xt, 2, colSums(Xt), "/")
  A  <- rbind(intercept = 1, dose = stats::runif(N))  # nrow(A) = 2 = rank
  Th <- matrix(c(2, 1, 0.5, 3), 2, 2)
  Y  <- Xt %*% Th %*% A
  if (noise > 0) Y <- pmax(Y + matrix(stats::rnorm(length(Y), 0, noise), P, N), 0)
  rownames(Y) <- paste0("v", seq_len(P))
  list(Y = Y, A = A, Xt = Xt, Th = Th)
}
best_perm_err <- function(X, Xt) {               # error after the best column order
  min(max(abs(X - Xt)), max(abs(X[, 2:1] - Xt)))
}

test_that("SPA finds the anchor rows of separable data, with or without noise", {
  skip_unless_full()
  d <- anchor_data()
  expect_setequal(nmfkc:::.nmfkc_spa_rows(d$Y, 2), c(3L, 11L))
  d <- anchor_data(noise = 0.002)
  expect_setequal(nmfkc:::.nmfkc_spa_rows(d$Y, 2), c(3L, 11L))
})

test_that("X.anchor = 'spa' recovers X and the covariate effects", {
  skip_unless_full()
  d <- anchor_data()
  f <- nmfkc(d$Y, d$A, rank = 2, X.anchor = "spa", epsilon = 1e-12,
             maxit = 1e5, verbose = FALSE)
  expect_setequal(unname(f$X.anchor), c(3L, 11L))
  expect_identical(names(f$X.anchor), c("Basis1", "Basis2"))
  ## the anchors' zeros are exact and stay so
  for (q in 1:2) expect_true(all(f$X[f$X.anchor[q], -q] == 0))
  expect_lt(best_perm_err(f$X, d$Xt), 1e-4)
  ## C = Theta, up to the order of the bases
  Cp <- if (max(abs(f$C - d$Th)) < max(abs(f$C[2:1, ] - d$Th))) f$C else f$C[2:1, ]
  expect_lt(max(abs(unname(Cp) - d$Th)), 1e-3)
})

test_that("without anchors the exact fit is a rotated one; X.init = 'spa' is only a start", {
  skip_unless_full()
  d <- anchor_data()
  f0 <- nmfkc(d$Y, d$A, rank = 2, epsilon = 1e-12, maxit = 1e5, verbose = FALSE)
  fs <- nmfkc(d$Y, d$A, rank = 2, X.init = "spa", epsilon = 1e-12, maxit = 1e5,
              verbose = FALSE)
  fa <- nmfkc(d$Y, d$A, rank = 2, X.anchor = "spa", epsilon = 1e-12, maxit = 1e5,
              verbose = FALSE)
  ## all three fit the data essentially exactly ...
  for (f in list(f0, fs, fa)) expect_lt(f$objfunc, 1e-6)
  ## ... but only the anchored one is pinned to the truth
  expect_gt(best_perm_err(f0$X, d$Xt), 10 * best_perm_err(fa$X, d$Xt))
  ## X.init = "spa" leaves no structural zero behind, and records no anchors
  expect_false(any(fs$X == 0))
  expect_null(fs$X.anchor)
})

test_that("given rows are used as given, and keep the basis order", {
  skip_unless_full()
  d <- anchor_data()
  fr <- nmfkc(d$Y, d$A, rank = 2, X.anchor = c("v11", "v3"), epsilon = 1e-12,
              maxit = 1e5, verbose = FALSE)
  expect_identical(unname(fr$X.anchor), c(11L, 3L))
  expect_true(fr$X[11, 2] == 0 && fr$X[3, 1] == 0)
  fi <- nmfkc(d$Y, d$A, rank = 2, X.anchor = c(11, 3), epsilon = 1e-12,
              maxit = 1e5, verbose = FALSE)
  expect_equal(fi$X, fr$X)
})

test_that("X.anchor takes precedence over X.init, and says so", {
  skip_unless_full()
  d <- anchor_data()
  expect_message(nmfkc(d$Y, d$A, rank = 2, X.anchor = "spa", X.init = "kmeans",
                       verbose = FALSE),
                 "X.init was not used")
  ## the same start is no conflict
  expect_no_message(nmfkc(d$Y, d$A, rank = 2, X.anchor = "spa", X.init = "spa",
                          verbose = FALSE))
})

test_that("anchors that cannot identify the bases are refused", {
  skip_unless_full()
  d <- anchor_data()
  ## more bases than covariates: B = C A has rank <= nrow(A)
  expect_error(nmfkc(d$Y, d$A, rank = 3, X.anchor = "spa", verbose = FALSE),
               "rank <= nrow\\(A\\)")
  expect_error(nmfkc(d$Y, d$A, rank = 2, X.anchor = 3, verbose = FALSE),
               "one anchor row per basis")
  expect_error(nmfkc(d$Y, d$A, rank = 2, X.anchor = c(3, 3), verbose = FALSE),
               "repeats a row")
  expect_error(nmfkc(d$Y, d$A, rank = 2, X.anchor = c(3, 99), verbose = FALSE),
               "between 1 and nrow")
  expect_error(nmfkc(d$Y, d$A, rank = 2, X.anchor = c("v3", "nope"), verbose = FALSE),
               "not in Y")
  Yz <- d$Y; Yz[5, ] <- 0
  expect_error(nmfkc(Yz, d$A, rank = 2, X.anchor = c(3, 5), verbose = FALSE),
               "no positive entry")
  Yd <- d$Y; Yd[5, ] <- 2 * Yd[3, ]
  expect_error(nmfkc(Yd, d$A, rank = 2, X.anchor = c(3, 5), verbose = FALSE),
               "linearly dependent")
})

test_that("an ordinary fit carries no X.anchor field", {
  skip_unless_full()
  d <- anchor_data()
  f <- nmfkc(d$Y, d$A, rank = 2, verbose = FALSE)
  expect_false("X.anchor" %in% names(f))
})

test_that("an anchored fit goes through nmfkc.inference()", {
  skip_unless_full()
  d <- anchor_data(noise = 0.002)
  f <- nmfkc(d$Y, d$A, rank = 2, X.anchor = "spa", epsilon = 1e-10,
             maxit = 1e5, verbose = FALSE)
  inf <- nmfkc.inference(f, d$Y, d$A, wild.B = 50)
  expect_true(all(is.finite(inf$coefficients$Estimate)))
})

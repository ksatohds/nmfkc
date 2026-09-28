## X.anchor and X.init = "spa" beyond nmfkc(): nmfre, nmf.ffb, nmf.rrr,
## nmf.gmm (and its two-stage route) and nmfkc.signed, which used to drop
## X.anchor without a word, and nmfkc.net, which cannot use it and now says so.

fam_data <- function(seed = 7, P = 20, N = 60) {
  set.seed(seed)
  Xt <- matrix(stats::runif(P * 2, 0.2, 1), P, 2)
  Xt[3, ] <- c(1, 0); Xt[11, ] <- c(0, 1)            # the anchors
  Xt <- sweep(Xt, 2, colSums(Xt), "/")
  A  <- rbind(intercept = 1, dose = stats::runif(N))
  Y  <- Xt %*% matrix(c(2, 1, 0.5, 3), 2, 2) %*% A +
        matrix(stats::rexp(P * N, 50), P, N)
  Y2 <- rbind(A[2, ] * 3 + 0.1, stats::runif(N) + 0.1)
  list(Y = Y, A = A, Y2 = Y2, Xt = Xt)
}
basis_of <- function(fit) if (!is.null(fit$X1)) fit$X1 else fit$X
## every anchor row keeps exact zeros off its own basis, and a positive entry on it
anchors_hold <- function(fit) {
  X <- basis_of(fit); J <- fit$X.anchor
  !is.null(J) && all(vapply(seq_along(J), function(k)
    all(X[J[k], -k] == 0) && X[J[k], k] > 0, logical(1)))
}
q <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("X.anchor = 'spa' is imposed and kept by every family that takes it", {
  skip_unless_full()
  d <- fam_data()
  fits <- list(
    nmfre    = q(nmfre(d$Y, d$A, rank = 2, X.anchor = "spa", verbose = FALSE)),
    ffb      = q(nmf.ffb(d$Y, d$Y2, rank = 2, X.anchor = "spa")),
    rrr      = q(nmf.rrr(d$Y, d$Y2, rank1 = 2, X.anchor = "spa")),
    gmm      = q(nmf.gmm(d$Y, rank = 2, K = 2, X.anchor = "spa")),
    twostage = q(nmf.gmm.twostage(d$Y, d$A, rank = 2, K = 2, X.anchor = "spa")),
    signed   = q(nmfkc.signed(d$Y, d$A, rank = 2, X.anchor = "spa", verbose = FALSE)))
  for (nm in names(fits)) {
    expect_setequal(unname(fits[[nm]]$X.anchor), c(3L, 11L))
    expect_true(anchors_hold(fits[[nm]]), label = nm)
  }
})

test_that("given anchor rows are used as given, in the given order", {
  skip_unless_full()
  d <- fam_data()
  fits <- list(
    nmfre  = q(nmfre(d$Y, d$A, rank = 2, X.anchor = c(11, 3), verbose = FALSE)),
    ffb    = q(nmf.ffb(d$Y, d$Y2, rank = 2, X.anchor = c(11, 3))),
    rrr    = q(nmf.rrr(d$Y, d$Y2, rank1 = 2, X.anchor = c(11, 3))),
    gmm    = q(nmf.gmm(d$Y, rank = 2, K = 2, X.anchor = c(11, 3))),
    signed = q(nmfkc.signed(d$Y, d$A, rank = 2, X.anchor = c(11, 3), verbose = FALSE)))
  for (nm in names(fits)) {
    expect_identical(unname(fits[[nm]]$X.anchor), c(11L, 3L), label = nm)
    expect_true(anchors_hold(fits[[nm]]), label = nm)
  }
})

test_that("X.init = 'spa' is a start only, and adds no X.anchor field", {
  skip_unless_full()
  d <- fam_data()
  fits <- list(
    nmfre  = q(nmfre(d$Y, d$A, rank = 2, X.init = "spa", verbose = FALSE)),
    ffb    = q(nmf.ffb(d$Y, d$Y2, rank = 2, X.init = "spa")),
    rrr    = q(nmf.rrr(d$Y, d$Y2, rank1 = 2, X.init = "spa")),
    gmm    = q(nmf.gmm(d$Y, rank = 2, K = 2, X.init = "spa")),
    signed = q(nmfkc.signed(d$Y, d$A, rank = 2, X.init = "spa", verbose = FALSE)))
  for (nm in names(fits)) {
    expect_null(fits[[nm]]$X.anchor)
    expect_true(all(is.finite(basis_of(fits[[nm]]))), label = nm)
  }
})

test_that("nmfkc.signed: a signed Y takes given anchors, and recovers the basis", {
  skip_unless_full()
  set.seed(7); P <- 20; N <- 60
  Xt <- matrix(stats::runif(P * 2, 0.2, 1), P, 2); Xt[3, ] <- c(1, 0); Xt[11, ] <- c(0, 1)
  As <- rbind(intercept = 1, dose = stats::rnorm(N))
  Ys <- Xt %*% matrix(c(2, 1, 0.5, -0.8), 2, 2) %*% As +
        matrix(stats::rnorm(P * N, 0, 0.01), P, N)
  expect_true(any(Ys < 0))
  Xt <- sweep(Xt, 2, colSums(Xt), "/")
  f <- q(nmfkc.signed(Ys, As, rank = 2, X.anchor = c(11, 3), epsilon = 1e-10,
                      maxit = 1e5, verbose = FALSE))
  expect_true(anchors_hold(f))
  ## basis 1 is anchored at row 11, which is the second true column
  expect_lt(max(abs(f$X - Xt[, c(2, 1)])), 1e-3)
  ## ... but SPA cannot find them in a signed Y
  expect_error(nmfkc.signed(Ys, As, rank = 2, X.anchor = "spa", verbose = FALSE),
               "needs a non-negative Y")
  expect_error(nmfkc.signed(Ys, As, rank = 2, X.init = "spa", verbose = FALSE),
               "needs a non-negative Y")
})

test_that("anchors that cannot identify the bases are refused, family by family", {
  skip_unless_full()
  d <- fam_data()
  ## B = C A in nmfkc.signed: rank <= nrow(A)
  expect_error(nmfkc.signed(d$Y, d$A, rank = 3, X.anchor = "spa", verbose = FALSE),
               "rank <= nrow\\(A\\)")
  ## nmf.rrr: the scores of X1 are C X2 Y2
  expect_error(nmf.rrr(d$Y, d$Y2, rank1 = 3, X.anchor = "spa"),
               "rank1 <= min\\(rank2, nrow\\(Y2\\)\\)")
  ## nmf.ffb: stage 1 is nmfkc(Y1, A = Y2)
  expect_error(nmf.ffb(d$Y, d$Y2, rank = 3, X.anchor = "spa"), "rank <= nrow\\(Y2\\)")
  ## nmf.gmm: SPA needs a non-negative Y
  expect_error(nmf.gmm(d$Y - mean(d$Y), rank = 2, K = 2, X.anchor = "spa"),
               "needs a non-negative Y")
})

test_that("where the options cannot apply, they are refused instead of dropped", {
  skip_unless_full()
  d <- fam_data()
  expect_error(nmfkc.net(tcrossprod(d$Y), rank = 2, X.anchor = "spa"),
               "not available for the symmetric model")
  expect_error(nmfkc.net(tcrossprod(d$Y), rank = 2, X.init = "spa"),
               "not available for the symmetric model")
  expect_error(nmf.ffb(d$Y, d$Y2, rank = 2, method = "mu", X.anchor = "spa"),
               "method = \"fiml\" only")
  expect_error(nmf.ffb(d$Y, d$Y2, rank = 2, method = "mu", X.init = "spa"),
               "method = \"fiml\" only")
  ffb0 <- q(nmf.ffb(d$Y, d$Y2, rank = 2))
  expect_error(nmf.ffb(d$Y, d$Y2, X = ffb0$X, X.anchor = "spa"),
               "cannot be combined with a supplied X")
})

test_that("X.anchor takes precedence over an X.init given with it, and says so", {
  skip_unless_full()
  d <- fam_data()
  expect_message(nmfre(d$Y, d$A, rank = 2, X.anchor = "spa", X.init = "kmeans",
                       verbose = FALSE), "X.init was not used")
  expect_message(nmf.gmm(d$Y, rank = 2, K = 2, X.anchor = "spa", X.init = "nndsvd"),
                 "X.init was not used")
  expect_message(nmfkc.signed(d$Y, d$A, rank = 2, X.anchor = "spa", X.init = "kmeans",
                              verbose = FALSE), "X.init was not used")
})

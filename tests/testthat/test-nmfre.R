test_that("nmfre() X.init paths run (CRAN-sized)", {
  skip_unless_full()
  # The full block below fits nmfre() eight times (~6s locally, x17 on CRAN's
  # regression covered there: character X.init used to error in
  # .nmfre.normalize.X ("'x' must be an array of at least two dimensions").
  set.seed(5); P <- 5; N <- 30; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  Y  <- pmax(Xt %*% (matrix(abs(rnorm(Q * 2)), Q, 2) %*% A +
                       matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)
  d <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE)
  expect_s3_class(d, "nmfre")
  fm <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = "nndsvd", verbose = FALSE)
  expect_true(all(is.finite(fm$X)))
})

test_that("nmfre() accepts X.init as NULL, a method name, or a matrix", {
  skip_unless_full()
  set.seed(5); P <- 5; N <- 50; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  ## one basis per row raised: without it the columns are alike (cosine 0.91), the
  ## fit is flagged degenerate, and X.init = "mix2" restarts 20 times
  ij <- cbind(seq_len(P), (seq_len(P) - 1) %% Q + 1); Xt[ij] <- Xt[ij] + 1
  Th <- matrix(abs(rnorm(Q * 2)), Q, 2)
  Y  <- pmax(Xt %*% (Th %*% A + matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)

  ## default (NULL) works
  d <- nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE)
  expect_s3_class(d, "nmfre")
  expect_equal(dim(d$X), c(P, Q))

  ## character init methods no longer crash (regression: used to error in
  ## .nmfre.normalize.X with "'x' must be an array of at least two dimensions")
  for (m in c("runif", "nndsvd", "kmeans")) {
    fm <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = m, verbose = FALSE)
    expect_s3_class(fm, "nmfre")
    expect_equal(dim(fm$X), c(P, Q))
    expect_true(all(is.finite(fm$X)))
  }

  ## a supplied basis matrix works
  Xm <- matrix(runif(P * Q), P, Q)
  fmat <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = Xm, verbose = FALSE)
  expect_equal(dim(fmat$X), c(P, Q))

  ## nstart is accepted and the default is reproducible
  a  <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = "kmeans", verbose = FALSE)
  a1 <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = "kmeans", nstart = 1, verbose = FALSE)
  expect_equal(a$objfunc, a1$objfunc)
  expect_s3_class(nmfre(Y, A = A, rank = Q, seed = 1, X.init = "kmeans", nstart = 5, verbose = FALSE), "nmfre")
})

test_that("nmfre.inference: identification is NA for a signed C, by the data for C >= 0", {
  skip_unless_full()
  set.seed(5); P <- 5; N <- 50; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  ij <- cbind(seq_len(P), (seq_len(P) - 1) %% Q + 1); Xt[ij] <- Xt[ij] + 1   # distinct bases (see above)
  Th <- matrix(abs(rnorm(Q * 2)), Q, 2)
  Y  <- pmax(Xt %*% (Th %*% A + matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)
  ## C.signed = TRUE (default): non-negativity cannot pin T^-1
  f <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = "kmeans", verbose = FALSE)
  inf <- suppressWarnings(nmfre.inference(f, Y, A, wild.B = 5))
  expect_true(is.na(inf$identified.approx))
  expect_match(inf$identification$note, "C is signed")
  ## Q = 2 with tau2 > 0: the isotropic random effects allow only a permutation
  expect_identical(inf$identification$by.assumption, isTRUE(f$tau2 > 0))
  expect_output(print(summary(inf)), "Identification:")
  ## C >= 0: the data diagnostic, exactly the shared check on (X, C)
  f2 <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = "kmeans", verbose = FALSE, C.signed = FALSE)
  inf2 <- suppressWarnings(nmfre.inference(f2, Y, A, wild.B = 5))
  expect_false(is.na(inf2$identified.approx))
  expect_identical(inf2$identified.approx, nmfkc:::.ident.verdict(f2$X, f2$C, A)$identified.approx)
})

test_that("nmfre(): a sign-free fit starts on a Y with negative entries", {
  skip_unless_full()
  ## Before 1.0.3 the nmfkc() initialisation stopped on min(Y) < 0 although
  ## the sign-free model allows it.  Now the basis starts on the shifted Y and C
  ## by sign-free least squares on the original Y.
  set.seed(7); P <- 6; N <- 60; Q <- 2
  A  <- rbind(intercept = 1, x = runif(N, 0, 2))
  Xt <- matrix(rgamma(P * Q, 1) + 0.2, P, Q); Xt <- sweep(Xt, 2, colSums(Xt), "/")
  Th <- cbind(intercept = c(2, 1), x = c(-1, 0.5))
  Y  <- Xt %*% (Th %*% A + matrix(rnorm(Q * N, 0, 0.5), Q, N)) +
          matrix(rnorm(P * N, 0, 0.1), P, N)
  expect_lt(min(Y), 0)
  f <- nmfre(Y, A = A, rank = Q, seed = 1, X.init = "kmeans", verbose = FALSE)
  expect_s3_class(f, "nmfre")
  expect_true(all(is.finite(f$C)))
  expect_true(all(f$X >= 0))
  ## C.signed = FALSE still needs Y >= 0
  expect_error(nmfre(Y, A = A, rank = Q, seed = 1, verbose = FALSE, C.signed = FALSE))
})

## data with an exact anchor row for each basis listed in `anchors` (row -> basis)
## block > 0 adds `block` to one basis per row (rows 1, 4, 7 -> basis 1, ...): the
## bases stay inside the orthant (a region of positive width) but are distinct.
## Without it the columns of Xt are alike (cosine 0.91 for seed 11), so a fit
## near the truth is already flagged degenerate and X.init = "mix2" restarts.
anchor_data <- function(seed = 11, P = 9, N = 120, Q = 3, anchors = c(`7` = 3), shift = 0,
                        block = 0) {
  set.seed(seed)
  A  <- rbind(intercept = 1, x = runif(N, 0, 2), z = rbinom(N, 1, 0.5))
  Xt <- matrix(runif(P * Q, 0.2, 1), P, Q)
  if (block > 0) Xt[cbind(seq_len(P), (seq_len(P) - 1) %% Q + 1)] <- Xt[cbind(seq_len(P), (seq_len(P) - 1) %% Q + 1)] + block
  for (j in names(anchors)) { Xt[as.integer(j), ] <- 0; Xt[as.integer(j), anchors[[j]]] <- 1 }
  Xt <- sweep(Xt, 2, colSums(Xt), "/")
  Th <- matrix(c(2, 1, 3, -1, 0.5, 0.8, 0.3, -0.4, 0.2)[seq_len(Q * 3)], Q, 3)
  Y  <- Xt %*% (Th %*% A + matrix(rnorm(Q * N, 0, 0.5), Q, N)) +
          matrix(rnorm(P * N, 0, 0.05), P, N) - shift
  rownames(Y) <- paste0("y", seq_len(P))
  list(Y = Y, A = A, Xt = Xt)
}

test_that("nmfre(): X.anchor = 'spa' works on a Y with negative entries", {
  skip_unless_full()
  ## Before 1.0.3 the anchor SPA ran on Y itself and stopped ("SPA needs a
  ## non-negative Y"), while X.init = "spa" ran on the shifted Y.
  d <- anchor_data(Q = 2, anchors = c(`3` = 1, `7` = 2), shift = 0.3)
  expect_lt(min(d$Y), 0)
  f <- suppressWarnings(nmfre(d$Y, d$A, rank = 2, X.anchor = "spa", seed = 1, verbose = FALSE))
  J <- f$X.anchor
  expect_length(J, 2)
  for (q in 1:2) {
    expect_true(all(f$X[J[q], -q] == 0))
    expect_gt(f$X[J[q], q], 0)
  }
})

test_that("nmfre(): fewer anchor rows than bases anchor the first bases only", {
  skip_unless_full()
  d <- anchor_data(Q = 3, anchors = c(`7` = 3))
  f <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.anchor = 7, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE))
  expect_identical(unname(f$X.anchor), 7L)
  expect_identical(names(f$X.anchor), colnames(f$X)[1])
  expect_true(all(f$X[7, -1] == 0))
  expect_gt(f$X[7, 1], 0)
  expect_true(all(f$X[-7, ] > 0))                  # no other zeros imposed
  ## the basis of the anchor row is chosen by likelihood over all 3 assignments
  s <- f$X.anchor.search
  expect_identical(nrow(s), 3L)
  expect_identical(which(s$chosen), which.min(s$nll))
  expect_equal(min(s$nll),
               nmfkc:::.nmfre.marginal.nll(d$Y, f$X, f$C, d$A, f$sigma2, f$tau2), tolerance = 1e-8)
  expect_true(is.numeric(f$nll.free))
  expect_identical(f$call[[1]], quote(nmfre))
  ## by row name, on a Y with negative entries, with an X.init method (used, no message)
  d2 <- anchor_data(Q = 3, anchors = c(`7` = 3), shift = 0.5)
  expect_lt(min(d2$Y), 0)
  msgs <- character(0)
  f2 <- withCallingHandlers(suppressWarnings(nmfre(d2$Y, d2$A, rank = 3, epsilon = 1e-4, X.anchor = "y7",
                                                   X.init = "kmeans++", seed = 2, verbose = FALSE)),
          message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_false(any(grepl("X.init was not used", msgs)))
  expect_true(all(f2$X[7, -1] == 0))
  ## two of three
  f3 <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.anchor = c("y7", "y1"), seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE))
  expect_identical(unname(f3$X.anchor), c(7L, 1L))
  expect_true(all(f3$X[7, -1] == 0) && all(f3$X[1, -2] == 0))
  expect_identical(nrow(f3$X.anchor.search), 6L)   # 3 x 2 ordered choices of bases
  ## refused as for full anchors
  expect_error(nmfre(d$Y, d$A, rank = 3, X.anchor = c(7, 7), verbose = FALSE), "repeats a row")
  expect_error(nmfre(d$Y, d$A, rank = 3, X.anchor = 12, verbose = FALSE), "between 1 and nrow")
  expect_error(nmfre(d$Y, d$A, rank = 3, X.anchor = "y99", verbose = FALSE), "not in Y")
})

test_that("identification region (Q = 3): exact range, arc ends, a point under an anchor", {
  skip_unless_full()
  d1 <- anchor_data()                                   # x7 loads on basis 3 only
  Th <- matrix(c(2, 1, 3, -1, 0.5, 0.8, 0.3, -0.4, 0.2), 3, 3)
  M  <- 0.9 * diag(3) + 0.05 * (matrix(1, 3, 3) - diag(3))   # mix the bases: no zero left
  Xm <- d1$Xt %*% M
  r <- nmfkc:::.nmfre.ident.region(Xm, Th)
  expect_true(r$computed && !r$full)
  expect_lt(r$arc[1], 0); expect_gt(r$arc[2], 0)
  V <- nmfkc:::.nmfre.O3.basis()
  ## O(phi) is orthogonal, fixes 1, and the ends of the arc are where X O first turns negative
  O <- nmfkc:::.nmfre.O3(0.3, V)
  expect_equal(crossprod(O), diag(3), tolerance = 1e-12)
  expect_equal(as.vector(O %*% rep(1, 3)), rep(1, 3), tolerance = 1e-12)
  expect_gte(min(Xm %*% nmfkc:::.nmfre.O3(r$arc[2], V)), -1e-9)
  expect_lt(min(Xm %*% nmfkc:::.nmfre.O3(r$arc[2] + 1e-3, V)), 0)
  expect_lt(min(Xm %*% nmfkc:::.nmfre.O3(r$arc[1] - 1e-3, V)), 0)
  ## the exact range of each entry of O'C over the arc matches a fine scan
  ph <- seq(r$arc[1], r$arc[2], length.out = 20001)
  path <- sapply(ph, function(p) as.vector(crossprod(nmfkc:::.nmfre.O3(p, V), Th)))
  expect_equal(as.vector(r$lower), apply(path, 1, min), tolerance = 1e-6)
  expect_equal(as.vector(r$upper), apply(path, 1, max), tolerance = 1e-6)
  pathX <- sapply(ph, function(p) as.vector(Xm %*% nmfkc:::.nmfre.O3(p, V)))
  expect_equal(as.vector(r$X.lower), apply(pathX, 1, min), tolerance = 1e-6)
  expect_equal(as.vector(r$X.upper), apply(pathX, 1, max), tolerance = 1e-6)
  expect_gte(min(r$X.lower), -1e-9)                     # X stays non-negative on the arc
  ## an anchor row (x7 on basis 3) pins the rotation: the region is a point
  r1 <- nmfkc:::.nmfre.ident.region(d1$Xt, Th)
  expect_lt(max(r1$width), 1e-6)
  ## rows far inside the orthant: the whole circle (cyclic relabellings included)
  d <- anchor_data(anchors = c())
  expect_true(nmfkc:::.nmfre.ident.region(d$Xt, Th)$full)
  ## Q = 2: a point up to relabelling; Q = 4 without anchors and without the search: not computed
  expect_equal(nmfkc:::.nmfre.ident.region(d$Xt[, 1:2] / rep(colSums(d$Xt[, 1:2]), each = 9), Th[1:2, ])$width,
               Th[1:2, ] * 0)
  expect_false(nmfkc:::.nmfre.ident.region(matrix(1 / 9, 9, 4), matrix(1, 4, 2),
                                           ident = list(search = FALSE))$computed)
})

test_that("identification region (Q >= 4): bases pinned by anchors one at a time, then exact or searched", {
  skip_unless_full()
  set.seed(3)
  nrm <- function(X) sweep(X, 2, colSums(X), "/")
  mk <- function(P, Q) nrm(matrix(runif(P * Q, 0.2, 1), P, Q))
  C4 <- matrix(rnorm(12), 4, 3)
  ## anchors on bases 1 and 2 pin them; the two left are a point up to relabelling
  Xa <- mk(12, 4); Xa[1, -1] <- 0; Xa[2, -2] <- 0; Xa <- nrm(Xa)
  ra <- nmfkc:::.nmfre.ident.region(Xa, C4)
  expect_identical(ra$method, "anchors"); expect_identical(ra$pinned, 1:2)
  expect_equal(max(ra$width), 0)
  ## a relative anchor: row 3 is non-zero on bases 1 and 3 only, so once basis 1 is pinned it pins basis 3
  Xr <- mk(12, 4); Xr[1, -1] <- 0; Xr[3, c(2, 4)] <- 0; Xr <- nrm(Xr)
  rr <- nmfkc:::.nmfre.ident.region(Xr, C4)
  expect_identical(rr$pinned, c(1L, 3L)); expect_identical(rr$method, "anchors")
  ## the non-negative variant: a column of C on one basis pins it as well
  Cd <- rbind(c(2, 0, 0), c(0, 3, 0), c(0, 0, 1.5), c(0, 0, 0.7))
  rd <- nmfkc:::.nmfre.ident.region(mk(12, 4), Cd, nonneg = TRUE)
  expect_identical(rd$pinned, 1:2); expect_identical(rd$method, "anchors")
  ## one anchor only: the three bases left form an arc, computed exactly.  A search over all four
  ## (no pinning) stays inside it; it finds only part of it, since the arc is a one-dimensional
  ## subset of a three-dimensional family -- which is why the anchors are removed first
  Xc <- matrix(0.08, 10, 4); Xc[1, ] <- c(1, 0, 0, 0)
  Xc[2:4, 2] <- 1; Xc[5:7, 3] <- 1; Xc[8:10, 4] <- 1; Xc <- nrm(Xc)
  rc <- nmfkc:::.nmfre.ident.region(Xc, C4)
  expect_identical(rc$method, "anchors+exact"); expect_identical(rc$pinned, 1L)
  expect_gt(max(rc$width), 0); expect_equal(max(rc$width[1, ]), 0)
  sc <- nmfkc:::.nmfre.ident.search(Xc, C4, n.dir = 300L, n.global = 3000L)
  expect_true(all(sc$lower >= rc$lower - 1e-8 & sc$upper <= rc$upper + 1e-8))
  ## no anchor: searched along feasible paths from the fit (approximate ranges)
  rb <- nmfkc:::.nmfre.ident.region(mk(12, 4), C4, ident = list(search = TRUE, n.dir = 200L, n.global = 2000L))
  expect_identical(rb$method, "search"); expect_gt(max(rb$width), 0); expect_gt(rb$search$local.feasible, 0)
  ## (counterexample, comment of 2026-11-02) the local region is a point but another component of the region
  ## has feasible transformations labelled near the identity: they must not enter the local ranges
  eta <- 0.01; rows <- list()
  for (q in 1:4) for (i in setdiff(1:4, q)) for (s in c(-1, 1)) {
    x <- rep(1, 4); x[q] <- 0; x[i] <- x[i] + s * eta; rows[[length(rows) + 1]] <- x
  }
  Xe <- do.call(rbind, rows) / 18
  re <- nmfkc:::.nmfre.ident.region(Xe, diag(4), ident = list(search = TRUE, n.dir = 100L, n.global = 5000L, seed = 1))
  expect_identical(re$method, "search"); expect_equal(max(re$width), 0)
  expect_gt(re$search$global.nonperm, 0)                 # not unique globally, reported apart
  ## (counterexample) tiny positive entries are not zeros: no pinning unless a tolerance is asked for
  Xd <- (1 - 1e-12) * rbind(diag(4), 0) + 1e-12 / 5
  expect_length(nmfkc:::.nmfre.ident.pin(Xd, diag(4)), 0)
  expect_identical(nmfkc:::.nmfre.ident.pin(Xd, diag(4), zero.tol = 1e-11), 1:4)
  ## through the fit: the identification list
  fo <- nmfkc:::.nmfre.orient.fit(list(X = Xa, C = C4, U = matrix(0, 4, 5)))
  expect_identical(fo$identification$method, "anchors"); expect_equal(fo$identification$width, 0)
  expect_match(fo$identification$orient$note, "pinned by an anchor")
})

test_that("nmfre(): geomin orientation within the identification region, likelihood and Xi unchanged", {
  skip_unless_full()
  d <- anchor_data()
  Th <- matrix(c(2, 1, 3, -1, 0.5, 0.8, 0.3, -0.4, 0.2), 3, 3,
               dimnames = list(paste0("Basis", 1:3), rownames(d$A)))
  N <- ncol(d$Y)
  M  <- 0.9 * diag(3) + 0.05 * (matrix(1, 3, 3) - diag(3))
  Xm <- d$Xt %*% M
  ## a representation inside the region, away from the geomin member
  V <- nmfkc:::.nmfre.O3.basis(); r <- nmfkc:::.nmfre.ident.region(Xm, Th)
  expect_gt(max(r$width), 0)
  O <- nmfkc:::.nmfre.O3(r$arc[2] * 0.6, V)
  fit <- list(X = Xm %*% O, C = crossprod(O, Th), U = matrix(0, 3, N), sigma2 = 0.05^2, tau2 = 0.25)
  fo <- nmfkc:::.nmfre.orient.fit(fit, orient = "geomin")   # nmfre(orient = "geomin") at its end
  nll <- function(o) nmfkc:::.nmfre.marginal.nll(d$Y, o$X, o$C, d$A, o$sigma2, o$tau2)
  expect_identical(fo$identification$orient$rule, "geomin")
  expect_gt(fo$identification$orient$moved, 0)
  expect_equal(nll(fo), nll(fit), tolerance = 1e-8)
  expect_equal(fo$X %*% fo$C, fit$X %*% fit$C, tolerance = 1e-10)
  expect_equal(unname(colSums(fo$X)), rep(1, 3), tolerance = 1e-8)
  expect_lte(fo$identification$orient$criterion, fo$identification$orient$criterion.start + 1e-12)
  ## the reported ranges contain the representation; the width is that of the start
  expect_true(all(fo$C.ident.lower - 1e-8 <= fo$C & fo$C <= fo$C.ident.upper + 1e-8))
  expect_true(all(fo$X.ident.lower - 1e-8 <= fo$X & fo$X <= fo$X.ident.upper + 1e-8))
  expect_equal(fo$C.ident.width, fo$C.ident.upper - fo$C.ident.lower)
  expect_equal(fo$identification$width, max(r$width), tolerance = 1e-6)
  expect_identical(dimnames(fo$C.ident.lower), dimnames(fo$C))
  ## orient = "none" keeps the representation
  fn <- nmfkc:::.nmfre.orient.fit(fit, orient = "none")
  expect_equal(fn$C, fit$C)
  expect_identical(fn$identification$orient$rule, "none")
  expect_equal(fn$identification$width, fo$identification$width, tolerance = 1e-6)
  ## summary() reports the width, the orientation and the coefficients of undetermined sign
  fs <- structure(c(fo, list(B.blup = matrix(0, 3, N), sigma2 = 0.05^2, tau2 = 0.25, lambda = 0.01,
                             ICC = 0.5, dfU = 10, converged = TRUE)), class = c("nmfre", "nmf"))
  out <- capture.output(nmfkc:::.nmfre.print.ident(fs, 3))
  expect_match(out[1], "width of Theta up to")
  expect_match(paste(out, collapse = " "), "geomin member")
})

test_that("nmfre(): orient defaults to \"geomin\"; orient = \"none\" keeps the fit as it stopped", {
  skip_unless_full()
  d <- anchor_data(anchors = c(), block = 1)            # rows inside the orthant: a region
  f  <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE))
  fn <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE, orient = "none"))
  fg <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE, orient = "geomin"))
  expect_identical(f$identification$orient$rule, "geomin")
  expect_identical(f$X, fg$X); expect_identical(f$C, fg$C); expect_identical(f$U, fg$U)
  expect_identical(fn$identification$orient$rule, "none")
  expect_equal(fn$X %*% fn$C, f$X %*% f$C, tolerance = 1e-8)
  expect_equal(fn$identification$width, f$identification$width, tolerance = 1e-6)
})

test_that("nmfre(): anchored fit has width 0; nmfre.inference() keeps it; optional profile intervals", {
  skip_unless_full()
  d <- anchor_data(Q = 3, anchors = c(`7` = 3))
  f <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.anchor = 7, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE))
  expect_lt(f$identification$width, 1e-6)
  expect_equal(f$identification$orient$angle, 0)
  expect_true(all(f$X[7, -1] == 0))                  # the anchor zeros stay exact
  expect_output(print(summary(f)), "Identification region: a point")
  fn <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.anchor = 7, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE, orient = "none"))
  expect_equal(fn$C, f$C)
  fi <- nmfre.inference(f, d$Y, d$A, wild.bootstrap = FALSE)
  expect_equal(fi$C, f$C)
  expect_equal(fi$identification$width, f$identification$width)   # kept
  expect_false(is.null(fi$identification$by.assumption))            # diagnostic added
  expect_equal(fi$C.ident.width, f$C.ident.width)
  expect_null(fi$C.profile.lower)
  fp <- nmfre.inference(f, d$Y, d$A, wild.bootstrap = FALSE, profile = TRUE, profile.which = 4)
  expect_identical(sum(!is.na(fp$C.profile.lower)), 1L)
  expect_true(fp$C.profile.reached[4])
  expect_lt(fp$C.profile.lower[4], fp$C[4]); expect_gt(fp$C.profile.upper[4], fp$C[4])
  ## close to the Wald interval when the region is a point
  se <- fp$C.se[4]
  expect_lt(abs((fp$C.profile.upper[4] - fp$C.profile.lower[4]) / (2 * qnorm(0.975) * se) - 1), 0.5)
})

test_that("nmfre.inference(): the identification note follows the region computed by nmfre()", {
  skip_unless_full()
  ## an anchor row on one basis only, yet the region is a point: locally unique
  d <- anchor_data(Q = 3, anchors = c(`7` = 3))
  f <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.anchor = 7, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE))
  fi <- nmfre.inference(f, d$Y, d$A, wild.bootstrap = FALSE)
  expect_false(all(fi$identification$anchor))
  expect_lt(fi$identification$width, 1e-6)
  expect_identical(fi$identification$by.assumption, isTRUE(f$tau2 > 0))
  expect_match(fi$identification$note, "locally unique")
  ## no anchor row: a region of positive width, and the note gives it
  d0 <- anchor_data(Q = 3, anchors = NULL, block = 1)
  f0 <- suppressWarnings(nmfre(d0$Y, d0$A, rank = 3, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE))
  expect_gt(f0$identification$width, 1e-6)
  fi0 <- nmfre.inference(f0, d0$Y, d0$A, wild.bootstrap = FALSE)
  expect_false(fi0$identification$by.assumption)
  expect_match(fi0$identification$note, "identification region of width")
  expect_false(grepl("not even", fi0$identification$note))
  ## C >= 0: when non-negativity alone does not pin the factorization but the region is a
  ## point, the summary says so instead of "another (X T, T^-1 C) may fit equally well"
  fn <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.anchor = 7, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE, C.signed = FALSE))
  fin <- nmfre.inference(fn, d$Y, d$A, wild.bootstrap = FALSE)
  out <- paste(capture.output(print(summary(fin))), collapse = "\n")
  if (identical(fin$identified.approx, FALSE) && isTRUE(fin$identification$by.assumption)) {
    expect_match(fin$identification$note, "not by non-negativity alone")
    expect_false(grepl("may fit equally well", out))
  }
  ## the printer: a FALSE verdict with a note prints the note
  fs <- fin; fs$identified.approx <- FALSE; fs$identification$note <- "not by non-negativity alone; test"
  expect_output(print(summary(fs)), "Identification: not by non-negativity alone; test")
})

test_that(".nmfre.nnls() solves the non-negative quadratic program exactly (KKT)", {
  set.seed(3)
  for (r in 1:25) {
    n <- sample(2:9, 1); M <- matrix(rnorm(3 * n * n), 3 * n, n); H <- crossprod(M); g <- 3 * rnorm(n)
    x <- nmfkc:::.nmfre.nnls(H, g); gr <- as.vector(H %*% x - g); s <- max(1, abs(g))
    expect_true(all(x >= 0))
    expect_lt(max(c(0, abs(gr[x > 0]))), 1e-8 * s)                 # stationary on the free set
    expect_true(all(gr[x == 0] >= -1e-8 * s))                       # no descent into the bound
  }
  ## an interior solution is the unconstrained one
  H <- diag(2); expect_equal(nmfkc:::.nmfre.nnls(H, c(1, 2)), c(1, 2))
  expect_equal(nmfkc:::.nmfre.nnls(H, c(-1, 2)), c(0, 2))
})

test_that(".nmfre.rows.nnls() and .nmfre.simplex.X() solve their quadratic programs exactly", {
  set.seed(5)
  for (r in 1:10) {
    Q <- sample(1:4, 1); P <- sample(3:12, 1)
    B <- matrix(rnorm(Q * 30), Q, 30); G <- B %*% t(B) + diag(0.1, Q); Cm <- matrix(3 * rnorm(P * Q), P, Q)
    X1 <- nmfkc:::.nmfre.rows.nnls(G, Cm)
    X2 <- t(apply(Cm, 1, function(cc) nmfkc:::.nmfre.nnls(G, cc)))
    if (Q == 1) X2 <- t(X2)
    expect_equal(X1, X2, tolerance = 1e-8)
    sx <- nmfkc:::.nmfre.simplex.X(G, Cm)
    expect_lt(max(abs(colSums(sx$X) - 1)), 1e-10); expect_true(all(sx$X >= 0))
    ## KKT on the simplices: within a column, the gradient 2(XG - C) equals -2 mu_q on X > 0 and is >= it on X = 0
    gr <- sx$X %*% G - Cm
    for (q in seq_len(Q)) {
      pos <- sx$X[, q] > 1e-12
      if (any(pos)) {
        expect_lt(diff(range(gr[pos, q])), 1e-7 * max(1, abs(Cm)))
        expect_true(all(gr[!pos, q] >= min(gr[pos, q]) - 1e-7 * max(1, abs(Cm))))
      }
    }
  }
})

test_that("nmfre(X.update = \"exact\"): no worse than \"semi\", columns sum to one exactly, recorded", {
  skip_unless_full()
  d <- anchor_data(Q = 3, anchors = c(`7` = 3))
  ## both at the default tolerance (the gradient check below has a wide margin: 2.5e-4
  ## against 0.08; 1e-8 for both took 90 s)
  fs <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.update = "semi", seed = 1, X.init = "kmeans", verbose = FALSE))
  fe <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.update = "exact", seed = 1, X.init = "kmeans", verbose = FALSE))
  expect_identical(fs$X.update, "semi"); expect_identical(fe$X.update, "exact")
  ## the stationarity check of a "semi" fit (anchor rows skip it, so use a free fit)
  d0 <- anchor_data(Q = 3, anchors = NULL, block = 1)
  f0 <- suppressMessages(suppressWarnings(nmfre(d0$Y, d0$A, rank = 3, X.update = "semi", seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE)))
  expect_true(is.numeric(f0$kkt$improvement) && f0$kkt$improvement >= -1e-6)
  expect_null(suppressWarnings(nmfre(d0$Y, d0$A, rank = 3, X.update = "semi", seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE,
                                     kkt.check = FALSE))$kkt)
  expect_null(suppressWarnings(nmfre(d0$Y, d0$A, rank = 3, seed = 1, X.init = "kmeans", epsilon = 1e-4, verbose = FALSE))$kkt)   # "exact" (default): no check
  expect_null(fe$kkt)
  expect_output(print(summary(f0)), "Stationarity check")
  expect_lt(max(abs(colSums(fe$X) - 1)), 1e-9)
  nll <- function(o) nmfkc:::.nmfre.marginal.nll(d$Y, o$X, o$C, d$A, o$sigma2, o$tau2)
  expect_lte(nll(fe), nll(fs) + 1e-6)
  ## sign-free C: dl/dC = 0 at the exact fit
  S <- fe$tau2 * tcrossprod(fe$X) + fe$sigma2 * diag(nrow(d$Y))
  g <- -crossprod(fe$X, solve(S, d$Y - fe$X %*% fe$C %*% d$A)) %*% t(d$A)
  expect_lt(max(abs(g)), 1e-4 * max(1, abs(crossprod(fe$X, solve(S, d$Y)) %*% t(d$A))))
  expect_error(nmfre(d$Y, d$A, rank = 3, X.update = "exact", X.L2.ortho = 1, verbose = FALSE), "does not support")
  ## the default falls back to "semi" with the X penalties
  expect_message(fp <- suppressWarnings(nmfre(d$Y, d$A, rank = 3, X.init = "kmeans", X.L2.ortho = 1, verbose = FALSE, maxit = 200)),
                 "need X.update")
  expect_identical(fp$X.update, "semi")
})

test_that("nmfre(C.signed = FALSE): exact updates by default stop at a KKT point; 'mu' is the legacy path", {
  skip_unless_full()
  set.seed(4); P <- 12; N <- 120; Q <- 2
  A  <- rbind(a0 = 1, a1 = rbinom(N, 1, 0.5), a2 = runif(N))
  Xt <- matrix(runif(P * Q, 0.1, 1), P, Q); Xt <- sweep(Xt, 2, colSums(Xt), "/")
  Th <- matrix(c(3, 3, 6, -1, 2, 5), Q, 3)                            # one negative effect: C >= 0 puts it at 0
  Y  <- pmax(Xt %*% (Th %*% A + matrix(rnorm(Q * N, 0, 0.5), Q, N)) +
               matrix(rnorm(P * N, 0, 0.05), P, N), 0)                # nmfkc() starts need Y >= 0
  ## with a true zero instead, the bound became active only at 1e-8 (12 s; 1e-10: 39 s)
  f  <- suppressWarnings(nmfre(Y, A, rank = Q, C.signed = FALSE, seed = 1, X.init = "kmeans", verbose = FALSE,
                               epsilon = 1e-6, epsilon.outer = 1e-6, maxit = 5e5))
  expect_identical(f$C.nonneg.update, "exact")
  expect_gt(sum(f$C <= 1e-8), 0)                                      # the bound is active somewhere
  expect_true(all(f$C >= 0))
  ## KKT of the marginal likelihood in C: dl/dC = -X' S^-1 (Y - X C A) A'
  S <- f$tau2 * tcrossprod(f$X) + f$sigma2 * diag(P)
  g <- -crossprod(f$X, solve(S, Y - f$X %*% f$C %*% A)) %*% t(A)
  g0 <- -crossprod(f$X, solve(S, Y)) %*% t(A)                        # scale of the gradient at C = 0
  expect_lt(max(c(0, abs(g[f$C > 1e-8]))), 1e-3 * max(abs(g0)))
  expect_true(all(g[f$C <= 1e-8] >= -1e-3 * max(abs(g0))))
  ## the legacy multiplicative updates are still available and do not do better
  fm <- suppressWarnings(nmfre(Y, A, rank = Q, C.signed = FALSE, C.nonneg.update = "mu", seed = 1, X.init = "kmeans", verbose = FALSE))
  expect_identical(fm$C.nonneg.update, "mu")
  nll <- function(o) nmfkc:::.nmfre.marginal.nll(Y, o$X, o$C, A, o$sigma2, o$tau2)
  expect_lte(nll(f), nll(fm) + 1e-6)
  expect_null(suppressWarnings(nmfre(Y, A, rank = Q, seed = 1, X.init = "kmeans", verbose = FALSE))$C.nonneg.update)
})

test_that("nmfre(X.init = \"mix2\", the default): the better of the k-means and SPA starts", {
  skip_unless_full()
  set.seed(5); P <- 6; N <- 60; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  Y  <- pmax(Xt %*% (matrix(abs(rnorm(Q * 2)), Q, 2) %*% A + matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)
  nll <- function(o) nmfkc:::.nmfre.marginal.nll(Y, o$X, o$C, A, o$sigma2, o$tau2)
  f  <- nmfre(Y, A, rank = Q, seed = 1, verbose = FALSE)
  fk <- nmfre(Y, A, rank = Q, seed = 1, X.init = "kmeans", verbose = FALSE)
  fs <- nmfre(Y, A, rank = Q, seed = 1, X.init = "spa", verbose = FALSE)
  st <- f$X.init.search
  expect_identical(st$start, c("kmeans", "spa"))
  expect_equal(st$seed, c(1, 1))
  expect_equal(st$nll, c(nll(fk), nll(fs)), tolerance = 1e-10)
  expect_identical(which(st$chosen), if (nll(fs) < nll(fk)) 2L else 1L)
  expect_equal(f$C, if (nll(fs) < nll(fk)) fs$C else fk$C, tolerance = 1e-10)
  expect_null(fk$X.init.search)                                      # a named start fits once
  expect_equal(nmfre(Y, A, rank = Q, seed = 1, X.init = "mix2", verbose = FALSE)$C, f$C)
  expect_output(print(summary(f)), "Starts (mix2)", fixed = TRUE)
  ## only the kept fit's warnings are passed on (both starts stop at maxit here)
  wn <- function(expr) { w <- 0L; withCallingHandlers(expr, warning = function(c) {
    w <<- w + 1L; invokeRestart("muffleWarning") }); w }
  fm <- NULL
  n_mix <- wn(fm <- nmfre(Y, A, rank = Q, seed = 1, maxit = 5, verbose = FALSE))
  n_one <- wn(nmfre(Y, A, rank = Q, seed = 1, maxit = 5, verbose = FALSE,
                    X.init = fm$X.init.search$start[fm$X.init.search$chosen]))
  expect_gt(n_one, 0L)
  expect_identical(n_mix, n_one)
})

test_that("nmfre(X.init = \"mix2\"): more k-means starts while the kept fit is degenerate", {
  skip_unless_full()
  set.seed(5); P <- 6; N <- 60; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  Y  <- pmax(Xt %*% (matrix(abs(rnorm(Q * 2)), Q, 2) %*% A + matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)
  ## every fit flagged degenerate: the restarts run to the cap (no early stop)
  local_mocked_bindings(.nmfre.degeneracy = function(X, C, cos.tol = 0.9, canc.tol = 5)
    list(max.cosine = 0.95, pair = c(1L, 2L), cancellation = 10, kappa.X = 1, flag = TRUE,
         cos.tol = cos.tol, canc.tol = canc.tol))
  expect_message(f <- nmfre(Y, A, rank = Q, seed = 3, mix.restarts = 3, verbose = FALSE),
                 "still degenerate after 3")
  st <- f$X.init.search
  expect_identical(st$start, c("kmeans", "spa", "kmeans", "kmeans", "kmeans"))
  expect_equal(st$seed, c(3, 3, 4, 5, 6))
  expect_identical(which(st$chosen), which.min(st$nll))
  f0 <- nmfre(Y, A, rank = Q, seed = 3, mix.restarts = 0, verbose = FALSE)
  expect_identical(nrow(f0$X.init.search), 2L)
})

test_that("nmfre(): the variances start from the residual covariance unless given", {
  skip_unless_full()
  set.seed(5); P <- 6; N <- 60; Q <- 2
  A  <- rbind(intercept = 1, male = rbinom(N, 1, 0.5))
  Xt <- matrix(abs(rnorm(P * Q)), P, Q)
  Y  <- pmax(Xt %*% (matrix(abs(rnorm(Q * 2)), Q, 2) %*% A + matrix(rnorm(Q * N, 0, 0.4), Q, N)) +
               matrix(rnorm(P * N, 0, 0.3), P, N), 0)
  f <- nmfre(Y, A, rank = Q, seed = 1, X.init = "kmeans", verbose = FALSE)
  R <- Y - Y %*% t(A) %*% solve(tcrossprod(A), A)
  l <- eigen(tcrossprod(R) / (N - nrow(A)), symmetric = TRUE, only.values = TRUE)$values
  expect_identical(f$variance.start$rule, "residual eigenvalues")
  expect_equal(f$variance.start$sigma2, mean(l[(Q + 1):P]), tolerance = 1e-8)
  expect_true(is.finite(f$variance.start$tau2) && f$variance.start$tau2 > 0)
  ## given: used as given (the start of nmfkc 1.0.3 and earlier)
  g <- nmfre(Y, A, rank = Q, seed = 1, X.init = "kmeans", sigma2 = 1, tau2 = 1, verbose = FALSE)
  expect_identical(g$variance.start$rule, "given")
  expect_equal(c(g$variance.start$sigma2, g$variance.start$tau2), c(1, 1))
  ## a variance not updated keeps the start 1 unless given
  h <- nmfre(Y, A, rank = Q, seed = 1, X.init = "kmeans", sigma2.update = FALSE, verbose = FALSE)
  expect_equal(h$sigma2, 1)
  expect_identical(h$variance.start$rule, "residual eigenvalues")
})

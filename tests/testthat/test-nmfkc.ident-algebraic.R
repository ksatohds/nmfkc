## nmfkc.ident(): the algebraic check on the fitted X and C (any A, also none).
## Exact matrices with a known answer, the exact rectangle for Q = 2, fits with
## and without covariates.

mk_ident <- function() {
  set.seed(11)
  P <- 8; N <- 20
  rdir <- function(n, k) { g <- matrix(stats::rgamma(n * k, 1), n, k); g / rowSums(g) }
  Xa <- rbind(diag(3), rdir(P - 3, 3)); Xa <- sweep(Xa, 2, colSums(Xa), "/")
  Ba <- cbind(diag(5, 3), matrix(stats::rgamma(3 * (N - 3), 2), 3, N - 3))
  Xd <- matrix(stats::rgamma(P * 3, 2), P, 3) + 0.2; Xd <- sweep(Xd, 2, colSums(Xd), "/")
  Bd <- matrix(stats::rgamma(3 * N, 2), 3, N) + 0.2
  ## one possible remixing only: X has an anchor row for every basis, C pure
  ## columns for bases 1 and 2, and row 3 of C inside the support of row 1 only
  Br <- cbind(c(5, 0, 0), c(0, 5, 0), c(3, 0, 2), matrix(stats::rgamma(3 * (N - 3), 2), 3, N - 3))
  Br[1, -2] <- pmax(Br[1, -2], 0.5)
  list(Xa = Xa, Ba = Ba, Xd = Xd, Bd = Bd, Br = Br)
}

test_that("two-sided anchors: identified, cone {0}, nothing moves", {
  m <- mk_ident()
  a <- nmfkc:::.nid.algebraic(m$Xa, m$Ba, n.dir = 300)
  expect_true(a$sufficient)
  expect_true(a$cone.trivial)
  expect_equal(unname(a$moved), c(0, 0, 0))
  expect_match(a$verdict, "^identified")
  expect_equal(a$width.X, 0); expect_equal(a$width.C, 0)
})

test_that("strictly positive X and C: not identified, every kind of direction moves", {
  m <- mk_ident()
  a <- nmfkc:::.nid.algebraic(m$Xd, m$Bd, n.dir = 300)
  expect_false(a$sufficient)
  expect_false(a$cone.trivial)
  expect_match(a$verdict, "^not identified")
  expect_true(all(a$moved > 0))
  ## the column sums of C do not change, and the bounds contain the fit
  expect_true(all(a$C.lower <= m$Bd + 1e-12) && all(a$C.upper >= m$Bd - 1e-12))
})

test_that("a single possible remixing: found by the pairs and the cone, missed by random directions", {
  m <- mk_ident()
  a <- nmfkc:::.nid.algebraic(m$Xa, m$Br, n.dir = 2000)
  expect_match(a$verdict, "^not identified")
  expect_false(a$cone.trivial)
  expect_equal(a$moved[["random"]], 0)
  expect_equal(a$moved[["pair"]], 1)
  expect_gte(a$moved[["cone"]], 1)
  ## the one pair: basis 3 can take a share of basis 1, only upwards
  p <- a$pairs[a$pairs$a.max - a$pairs$a.min > 0, ]
  expect_equal(nrow(p), 1)
  expect_equal(c(p$basis, p$from), c("Basis3", "Basis1"))
  expect_equal(p$a.min, 0)
  expect_gt(p$a.max, 0)
})

test_that("the cone is decided exactly (Farkas / NNLS) on small zero patterns", {
  ## no zeros: every off-diagonal direction is possible
  cn <- nmfkc:::.nid.alg.cone(matrix(0, 0, 6), 6)
  expect_true(all(cn$feasible)); expect_false(cn$trivial)
  ## anchors in X for every basis and in C for every basis: nothing is possible
  m <- mk_ident()
  M <- nmfkc:::.nid.alg.M(m$Xa, m$Ba)
  cn <- nmfkc:::.nid.alg.cone(M, 6)
  expect_true(cn$trivial)
  ## every certificate is a direction of the cone
  M2 <- nmfkc:::.nid.alg.M(m$Xa, m$Br)
  cn2 <- nmfkc:::.nid.alg.cone(M2, 6)
  for (r in cn2$rays) expect_true(all(M2 %*% r >= -1e-8))
  expect_equal(sum(cn2$feasible), 1)
})

test_that("Q = 2: the pairwise intervals are the exact rectangle", {
  fit <- nmfkc(t(as.matrix(iris[, 1:4])), rank = 2, epsilon = 1e-8, maxit = 1e5, verbose = FALSE)
  a <- nmfkc:::.nid.algebraic(fit$X, fit$C, n.dir = 200)
  sn <- nmfkc:::.nid.alg.snap(fit$X, fit$C, 1e-3); X <- sn$X; B <- sn$C
  ## T = [[1-a, b], [a, 1-b]]: column 1 -> (1-a) x1 + a x2, column 2 -> (1-b) x2 + b x1
  x1 <- X[, 1]; x2 <- X[, 2]; s <- colSums(B); b1 <- B[1, ] / s; b2 <- B[2, ] / s
  up <- function(num, den) { k <- den > 0; if (any(k)) min(num[k] / den[k]) else Inf }
  ex <- c(aL = -up(x1, x2 - x1), aU = min(up(x1, x1 - x2), min(b2)),
          bL = -up(x2, x1 - x2), bU = min(up(x2, x2 - x1), min(b1)))
  pa <- a$pairs[a$pairs$basis == "Basis1", ]; pb <- a$pairs[a$pairs$basis == "Basis2", ]
  expect_equal(unname(c(pa$a.min, pa$a.max, pb$a.min, pb$a.max)), unname(ex), tolerance = 1e-6)
})

test_that("nmfkc.ident(fit) without covariates gives the algebraic check only, and prints", {
  fit <- nmfkc(t(as.matrix(iris[, 1:4])), rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(fit, ident.n = 200)
  expect_s3_class(id, "nmfkc.ident")
  expect_null(id$test); expect_null(id$Xi.hat)
  expect_match(id$algebraic$verdict, "^not identified")
  out <- capture.output(print(id))
  expect_true(any(grepl("Algebraic check on the fitted X and C", out)))
  expect_true(any(grepl("no A given", out)))
  ## Y without A is accepted (the test needs A); A without Y is not
  expect_null(nmfkc.ident(fit, t(as.matrix(iris[, 1:4])), ident.n = 50)$test)
  expect_error(nmfkc.ident(fit, A = diag(2)), "Y is needed")
})

test_that("with covariates the earlier components are unchanged and the check is added", {
  set.seed(1)
  X <- cbind(c(0.4, 0.3, 0.2, 0.1, 0), c(0, 0.1, 0.2, 0.3, 0.4))
  A <- rbind(g1 = rep(c(1, 0), each = 20), g2 = rep(c(0, 1), each = 20))
  Y <- (X %*% diag(c(10, 8)) %*% A) * exp(matrix(rnorm(200, sd = 0.1), 5, 40))
  fit <- nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(fit, Y, A, region.n = 2000)
  expect_false(is.null(id$test)); expect_false(is.null(id$Xi.hat))
  expect_true(is.list(id$algebraic))
  ## the check on the fitted X and C needs neither Y nor A, and agrees
  id0 <- nmfkc.ident(fit)
  expect_equal(id0$algebraic, id$algebraic)
  ## the caller's random stream is left alone
  set.seed(5); u1 <- runif(1); set.seed(5); invisible(nmfkc.ident(fit)); u2 <- runif(1)
  expect_equal(u1, u2)
})

test_that("plot(): X, B and C from the algebraic check, bases and at", {
  fit <- nmfkc(t(as.matrix(iris[, 1:4])), rank = 3, epsilon = 1e-8, maxit = 1e5, verbose = FALSE)
  id <- nmfkc.ident(fit, ident.n = 200)
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  expect_invisible(plot(id, type = "X"))
  expect_invisible(plot(id, type = "X", bases = c("Basis1", "Basis3")))
  expect_invisible(plot(id))                         # 150 columns: type "B" by default
  hit <- integer(0)
  expect_invisible(plot(id, type = "B", at = seq(0, 1, length.out = 150), panel = function(q) hit <<- c(hit, q)))
  ## by default "B" draws the bases that move
  cs <- colSums(id$C)
  mv <- which(apply(sweep(id$algebraic$C.upper - id$algebraic$C.lower, 2, cs, "/"), 1, max) > 1e-3)
  expect_equal(hit, unname(mv))
  expect_error(plot(id, type = "B", at = 1:3), "one position per column")
  ## shares of the column sums: B.prob (membership probabilities) and C.prob
  expect_invisible(plot(id, type = "B.prob"))
  sh <- sweep(id$algebraic$C.upper, 2, colSums(id$C), "/")
  expect_true(all(sh <= 1 + 1e-9) && all(sweep(id$algebraic$C.lower, 2, colSums(id$C), "/") >= -1e-9))
  expect_equal(colSums(id$C), colSums(fit$B))
  expect_error(plot(id, bases = "nope"), "bases must")
  ## a fit with covariates and only the fitted matrices: the C grid
  set.seed(1)
  X <- cbind(c(0.4, 0.3, 0.2, 0.1, 0), c(0, 0.1, 0.2, 0.3, 0.4))
  A <- rbind(g1 = rep(c(1, 0), each = 20), g2 = rep(c(0, 1), each = 20))
  Y <- (X %*% diag(c(10, 8)) %*% A) * exp(matrix(rnorm(200, sd = 0.1), 5, 40))
  fa <- nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  expect_invisible(plot(nmfkc.ident(fa), type = "C"))
  expect_invisible(plot(nmfkc.ident(fa, Y, A, region.n = 2000), type = "C", bases = 2))
  expect_invisible(plot(nmfkc.ident(fa, Y, A, region.n = 2000, ci = TRUE), type = "C.prob"))
  grDevices::dev.off(); unlink(f)
})

test_that("plot(col = ): the colours reach the drawing, NA leaves a part out, bad input stops", {
  set.seed(1)
  X <- cbind(c(0.4, 0.3, 0.2, 0.1, 0), c(0, 0.1, 0.2, 0.3, 0.4))
  A <- rbind(g1 = rep(c(1, 0), each = 20), g2 = rep(c(0, 1), each = 20))
  Y <- (X %*% diag(c(10, 8)) %*% A) * exp(matrix(rnorm(200, sd = 0.1), 5, 40))
  fa <- nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE)
  id <- nmfkc.ident(fa, Y, A, region.n = 2000, ci = TRUE)
  ## an uncompressed PDF sets each colour as "r g b SCN" (stroke) or "r g b scn" (fill;
  ## "RG" / "rg" with colormodel = "rgb"): match the three numbers
  ops <- function(...) {
    f <- tempfile(fileext = ".pdf"); grDevices::pdf(f, compress = FALSE)
    plot(id, ...); grDevices::dev.off()
    on.exit(unlink(f)); txt <- readLines(f, warn = FALSE)
    txt[!grepl("CreationDate|ModDate", txt, useBytes = TRUE)]
  }
  has <- function(txt, s) any(grepl(s, txt, fixed = TRUE, useBytes = TRUE))   # the PDF has binary lines
  blue <- "0.000 0.000 1.000 "; red3 <- "0.804 0.000 0.000 "
  gray55 <- "0.549 0.549 0.549 "; gray90 <- "0.898 0.898 0.898 "
  d <- ops(type = "C")
  expect_true(has(d, red3) && has(d, gray55) && has(d, gray90) && !has(d, blue))
  b <- ops(type = "C", col = c(estimate = "blue", shade = NA))
  expect_true(has(b, blue) && has(b, gray55) && !has(b, red3) && !has(b, gray90))
  expect_true(has(ops(type = "X", col = c(range = "blue")), blue))
  expect_true(has(ops(type = "B", col = c(estimate = "blue")), blue))
  expect_false(has(ops(type = "B", col = c(estimate = NA)), red3))
  ## the legend names a colour given by name, otherwise only the mark
  w <- .nid.plot.col(c(estimate = "blue", range = "#1f77b4", interval = "grey40"))$word
  expect_identical(unname(w[c("range", "estimate", "interval")]), c("", "blue", "grey"))
  expect_identical(.nid.plot.col(NULL), .nid.plot.col(c(estimate = "red3")))
  expect_error(plot(id, col = "blue"), "named vector")
  expect_error(plot(id, col = c(colour = "blue")), "named vector")
  expect_error(plot(id, col = c(estimate = "nocolour")), "not a colour")
  ## cex and lwd: the defaults given explicitly change nothing, other values do;
  ## a line width lwd is written as "0.75 * lwd w"
  expect_identical(ops(type = "C", cex = c(estimate = 1.6, range = 0.8), lwd = c(range = 6, interval = 2.5, estimate = 2)), d)
  expect_identical(ops(type = "X", cex = c(estimate = 0.9, range = 0.45), lwd = c(range = 4, estimate = 1.5)),
                   ops(type = "X"))
  expect_identical(ops(type = "B", cex = c(estimate = 0.5)), ops(type = "B"))
  for (tp in c("C", "X", "B"))
    expect_false(identical(ops(type = tp, cex = c(estimate = 3)), ops(type = tp)))
  expect_false(identical(ops(type = "X", lwd = c(estimate = 4)), ops(type = "X")))
  expect_false(has(d, "7.50 w"))
  expect_true(has(ops(type = "C", lwd = c(range = 10)), "7.50 w"))
  expect_true(has(ops(type = "C", lwd = c(interval = 10)), "7.50 w"))
  expect_true(has(ops(type = "X", lwd = c(estimate = 10)), "7.50 w"))   # (the X ranges of id are points)
  expect_true(has(ops(type = "B", lwd = c(range = 10)), "7.50 w"))
  expect_identical(ops(type = "X", lwd = c(interval = 9)), ops(type = "X"))   # not drawn in "X": ignored
  expect_error(plot(id, cex = 0.5), "named numeric vector")
  expect_error(plot(id, cex = c(interval = 1)), "named numeric vector")
  expect_error(plot(id, lwd = c(range = "2")), "named numeric vector")
  expect_error(plot(id, lwd = c(range = 0)), "positive")
})

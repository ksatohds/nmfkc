## Gram-object input to nmfkc(): the large-N route for Nystroem kernel
## covariates.  nmfkc.kernel.gram() chooses landmarks (k-means++ on a
## subsample by default) and accumulates S = A A^T and G0 = Y A^T over
## column blocks so the M x N matrix never exists; nmfkc() accepts the
## object in place of A.  These tests pin down (i) that S, G0 equal the
## one-shot products, (ii) that the fit equals the explicit-matrix fit,
## (iii) the documented restrictions, and (iv) that a matrix A still takes
## the old path exactly.

make_nys_case <- function(seed = 5, p = 4, N = 120, P = 3) {
  set.seed(seed)
  centers <- matrix(stats::rnorm(p * P, sd = 3), p, P)
  lab <- sample.int(P, N, replace = TRUE)
  U <- centers[, lab] + matrix(stats::rnorm(p * N), p, N)
  Y <- matrix(0, P, N); Y[cbind(lab, seq_len(N))] <- 1
  rownames(Y) <- paste0("c", seq_len(P))
  list(U = U, Y = Y, lab = lab)
}

test_that("block-wise S and G0 equal the one-shot Nystroem products", {
  skip_unless_full()
  d <- make_nys_case()
  ## block.size = 25 with N = 120 -> 5 blocks, the last one ragged (20 columns)
  expect_message(
    g <- nmfkc.kernel.gram(d$Y, d$U, V = 9, block.size = 25, seed = 2),
    "nearest-landmark median")
  expect_s3_class(g, "nmfkc.gram")
  expect_identical(g$type, "nystrom"); expect_false(g$signed)
  expect_equal(dim(g$landmarks), c(4L, 9L))
  expect_equal(g$n.blocks, 5L); expect_equal(g$D, 9L); expect_equal(g$N, 120L)
  A <- nmfkc.kernel(g$landmarks, d$U, beta = g$beta)          # 9 x 120
  expect_equal(dim(A), c(9L, 120L))
  expect_true(all(A >= 0))
  expect_equal(g$S,  unname(tcrossprod(unclass(A))),   tolerance = 1e-12)
  expect_equal(unname(g$G0), unname(tcrossprod(d$Y, unclass(A))), tolerance = 1e-12)
  expect_equal(g$A.block(c(1L, 60L, 120L)), unname(unclass(A))[, c(1L, 60L, 120L)],
               tolerance = 1e-12)
  expect_output(print(g), "Nystroem")
})

test_that("landmark selection is reproducible and accepts a supplied matrix", {
  skip_unless_full()
  d <- make_nys_case()
  g1 <- suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 6, seed = 9))
  g2 <- suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 6, seed = 9))
  expect_equal(g1$landmarks, g2$landmarks)
  expect_identical(g1$landmarks.method, "kmeans++")
  ## Other selection methods run and give p x M landmarks
  gk <- suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 6, seed = 9, landmarks = "kmeans"))
  gr <- suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 6, seed = 9, landmarks = "random"))
  expect_equal(dim(gk$landmarks), c(4L, 6L)); expect_equal(dim(gr$landmarks), c(4L, 6L))
  ## Supplied landmarks + explicit beta: no message, same S as the one-shot
  g3 <- nmfkc.kernel.gram(d$Y, d$U, V = g1$landmarks, beta = 0.2, block.size = 40)
  expect_identical(g3$landmarks.method, "supplied")
  A3 <- nmfkc.kernel(g1$landmarks, d$U, beta = 0.2)
  expect_equal(g3$S, unname(tcrossprod(unclass(A3))), tolerance = 1e-12)
  ## Seeding does not disturb the caller's stream
  set.seed(77); r_before <- stats::runif(1)
  set.seed(77); invisible(suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 4)))
  expect_identical(stats::runif(1), r_before)
  ## Argument checks
  expect_error(nmfkc.kernel.gram(d$Y, d$U, V = matrix(1, 3, 2)), "nrow\\(U\\)")
  expect_error(nmfkc.kernel.gram(d$Y, d$U, V = 500, sample.size = 50), "raise sample.size")
})

test_that("Gram-object fit equals the explicit-matrix fit in nmfkc()", {
  skip_unless_full()
  d <- make_nys_case()
  g <- suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 9, block.size = 30, seed = 2))
  A <- nmfkc.kernel(g$landmarks, d$U, beta = g$beta)
  fg <- nmfkc(d$Y, A = g, rank = 3, verbose = FALSE)
  fm <- nmfkc(d$Y, A = A, rank = 3, verbose = FALSE)
  expect_s3_class(fg, "nmfkc")
  expect_equal(fg$iter, fm$iter)
  expect_equal(fg$X,  fm$X,  tolerance = 1e-8)
  expect_equal(fg$C,  fm$C,  tolerance = 1e-8)
  expect_equal(fg$B,  fm$B,  tolerance = 1e-8)
  expect_equal(fg$XB, fm$XB, tolerance = 1e-8)
  expect_equal(fg$objfunc,      fm$objfunc,      tolerance = 1e-8)
  expect_equal(fg$objfunc.iter, fm$objfunc.iter, tolerance = 1e-8)
  expect_equal(fg$r.squared, fm$r.squared, tolerance = 1e-8)
  expect_equal(fg$sigma, fm$sigma, tolerance = 1e-8)
  expect_identical(fg$A.attr$function.name, "nmfkc.gram")
  expect_equal(fg$A.attr$dim, c(9L, 120L))
  expect_match(fg$dims, "Gram input")
  ## S3 verbs and prediction on new data
  expect_output(print(summary(fg)))
  A.new <- nmfkc.kernel(g$landmarks, d$U[, 1:15], beta = g$beta)
  expect_equal(predict(fg, newA = A.new), fm$XB[, 1:15], tolerance = 1e-8)
  cls <- predict(fg, newA = A.new, type = "class")
  expect_length(cls, 15L)
  ## The well-separated clusters are classified correctly
  expect_gt(mean(predict(fg, newA = A, type = "class") == rownames(d$Y)[d$lab]), 0.95)
})

test_that("Gram input with X penalties and C.L1 matches the matrix path", {
  skip_unless_full()
  d <- make_nys_case()
  g <- suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 9, seed = 2))
  A <- nmfkc.kernel(g$landmarks, d$U, beta = g$beta)
  for (args in list(list(C.L1 = 0.05), list(X.L2.ortho = 0.1), list(X.L2.smooth = 0.1),
                    list(X.restriction = "colSqSums"), list(X.restriction = "fixed"))) {
    fg <- do.call(nmfkc, c(list(d$Y, A = g, rank = 3, verbose = FALSE, maxit = 200), args))
    fm <- do.call(nmfkc, c(list(d$Y, A = A, rank = 3, verbose = FALSE, maxit = 200), args))
    expect_equal(fg$C, fm$C, tolerance = 1e-6, info = paste(names(args), collapse = ","))
    expect_equal(fg$objfunc.iter, fm$objfunc.iter, tolerance = 1e-6,
                 info = paste(names(args), collapse = ","))
  }
})

test_that("nmfkc() enforces the Gram-input restrictions", {
  skip_unless_full()
  d <- make_nys_case()
  g <- suppressMessages(nmfkc.kernel.gram(d$Y, d$U, V = 9, seed = 2))
  expect_error(nmfkc(d$Y, A = g, rank = 3, method = "KL", verbose = FALSE), "EU")
  W <- matrix(1, nrow(d$Y), ncol(d$Y)); W[1, 1] <- 0
  expect_error(nmfkc(d$Y, A = g, rank = 3, Y.weights = W, verbose = FALSE), "Y.weights")
  Yna <- d$Y; Yna[1, 1] <- NA
  expect_error(nmfkc(Yna, A = g, rank = 3, verbose = FALSE), "NA")
  expect_error(nmfkc(d$Y[, 1:100], A = g, rank = 3, verbose = FALSE), "N = 120")
  expect_error(nmfkc(d$Y[1:2, ], A = g, rank = 2, verbose = FALSE), "different Y")
  ## A signed (RFF) Gram object belongs to nmfkc.signed()
  gs <- nmfkc.signed.rff.gram(d$Y, d$U, beta = 0.3, D = 8, seed = 1)
  expect_error(nmfkc(d$Y, A = gs, rank = 3, verbose = FALSE), "nmfkc.signed")
  ## ... while nmfkc.signed() accepts the non-negative one too
  fs <- suppressMessages(nmfkc.signed(d$Y, A = g, rank = 3, maxit = 100, verbose = FALSE))
  expect_equal(dim(fs$C), c(3L, 9L))
  ## Fold-based helpers refuse a Gram object with a pointer to the alternative
  expect_error(nmfkc.cv(d$Y, A = g, rank = 3), "cannot be split into folds")
  expect_error(nmfkc.ecv(d$Y, A = g, rank = 2:3), "cannot be split into folds")
  expect_error(nmfkc.rank(d$Y, A = g, rank = 2:3, plot = FALSE), "cannot be split into folds")
})

test_that("a matrix A still takes the original nmfkc() path (regression guard)", {
  skip_unless_full()
  d <- make_nys_case()
  A <- nmfkc.kernel(d$U[, 1:9], d$U, beta = 0.2)
  fit <- nmfkc(d$Y, A = A, rank = 3, verbose = FALSE)
  expect_identical(fit$A.attr$`function.name`, "nmfkc.kernel")
  expect_equal(dim(fit$B), c(3L, 120L))
  expect_equal(fit$objfunc, sum((d$Y - fit$XB)^2), tolerance = 1e-10)
  f0 <- nmfkc(d$Y, rank = 3, verbose = FALSE)             # no covariates
  expect_null(f0$A.attr)
  expect_equal(dim(f0$B), c(3L, 120L))
})

## nmfkc.net(diag.exclude) and nmfkc.net.ecv(diag.exclude, folds, seeds, pred,
## r.squared.cv), 2026-10-09; and the two ECV inputs that used to go wrong
## (NA in Y, a Y.weights replaced by the fold mask).

net_data <- function(N = 20, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::runif(N * 2), N, 2)
  Y <- X %*% matrix(c(2, .3, .3, 1), 2) %*% t(X) + matrix(stats::rexp(N * N, 50), N)
  Y <- (Y + t(Y)) / 2
  dimnames(Y) <- list(paste0("n", 1:N), paste0("n", 1:N))
  Y
}
q <- function(expr) suppressWarnings(suppressMessages(expr))
drop_fields <- function(f, nm) f[setdiff(names(f), nm)]

test_that("nmfkc.net(diag.exclude = TRUE) is a zero-diagonal Y.weights", {
  skip_unless_full()
  Y <- net_data(); N <- nrow(Y)
  W0 <- matrix(1, N, N); diag(W0) <- 0
  for (tp in c("tri", "bi", "signed")) {
    a <- q(nmfkc.net(Y, rank = 2, type = tp, diag.exclude = TRUE))
    b <- q(nmfkc.net(Y, rank = 2, type = tp, Y.weights = W0))
    expect_true(isTRUE(a$diag.exclude))
    expect_identical(drop_fields(a, c("call", "runtime", "diag.exclude")),
                     drop_fields(b, c("call", "runtime")), label = tp)
  }
  ## combined with an NA mask: both the NA pair and the diagonal are masked
  Yna <- Y; Yna[2, 5] <- Yna[5, 2] <- NA
  W1 <- W0; W1[2, 5] <- W1[5, 2] <- 0
  a <- q(nmfkc.net(Yna, rank = 2, diag.exclude = TRUE))
  b <- q(nmfkc.net(Y, rank = 2, Y.weights = W1))
  expect_identical(a$X, b$X)
  ## the default leaves the fit object as it was
  expect_null(q(nmfkc.net(Y, rank = 2))$diag.exclude)
})

test_that("ECV with diag.exclude never holds out, fits or scores the diagonal", {
  skip_unless_full()
  Y <- net_data(); N <- nrow(Y)
  e <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 4, diag.exclude = TRUE, pred = TRUE))
  idx <- sort(unlist(e$folds))
  expect_identical(idx, which(upper.tri(matrix(0, N, N))))
  expect_true(isTRUE(e$diag.exclude))
  expect_true(all(is.na(diag(e$pred[, , 1]))))
  ## the same split and fits as a zero-diagonal Y.weights
  W0 <- matrix(1, N, N); diag(W0) <- 0
  w <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 4, Y.weights = W0))
  expect_identical(e$folds, w$folds)
  expect_equal(e$objfunc, w$objfunc, tolerance = 1e-12)
})

test_that("pred and r.squared.cv agree with the per-fold MSE and the usual R^2", {
  skip_unless_full()
  Y <- net_data(); N <- nrow(Y)
  e <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 4, diag.exclude = TRUE, pred = TRUE))
  ut <- upper.tri(matrix(0, N, N))
  for (h in 1:2) {
    P <- e$pred[, , h]
    M <- P; M[is.na(M)] <- 0
    expect_true(isSymmetric(unname(M)))
    r2 <- 1 - sum((Y[ut] - P[ut])^2) / sum((Y[ut] - mean(Y[ut]))^2)
    expect_equal(unname(e$r.squared.cv[h]), r2, tolerance = 1e-12)
    for (k in seq_along(e$folds)) {
      f <- e$folds[[k]]
      expect_equal(e$objfunc.fold[[h]][k], mean((Y[f] - P[f])^2), tolerance = 1e-12)
    }
  }
  ## without pred = TRUE the object stays small
  expect_null(q(nmfkc.net.ecv(Y, rank = 1, nfolds = 3))$pred)
})

test_that("folds = reuses a split; bad folds are refused", {
  skip_unless_full()
  Y <- net_data(); N <- nrow(Y)
  e <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 3))
  r <- q(nmfkc.net.ecv(Y, rank = 1:2, folds = e$folds))
  expect_identical(r$objfunc, e$objfunc)
  expect_identical(r$folds, e$folds)
  ## a lower-triangle index stands for its upper-triangle mirror
  fl <- e$folds
  i <- fl[[1]][1]; r1 <- (i - 1) %% N + 1; c1 <- (i - 1) %/% N + 1
  if (r1 != c1) {
    fl[[1]][1] <- (r1 - 1) * N + c1
    expect_identical(q(nmfkc.net.ecv(Y, rank = 1:2, folds = fl))$objfunc, e$objfunc)
  }
  ## another implementation's split (upper triangle, no diagonal)
  ut <- which(upper.tri(matrix(0, N, N)))
  set.seed(5); ext <- split(ut[sample.int(length(ut))], rep_len(1:5, length(ut)))
  x <- q(nmfkc.net.ecv(Y, rank = 1:2, folds = ext, diag.exclude = TRUE))
  expect_equal(x$nfolds, 5L)
  expect_true(all(is.finite(x$objfunc)))
  ## refusals
  dup <- e$folds; dup[[2]] <- c(dup[[2]], dup[[1]][1])
  expect_error(nmfkc.net.ecv(Y, rank = 1, folds = dup), "more than one fold")
  dg <- e$folds; dg[[1]] <- c(dg[[1]], 1L)        # (1, 1), on the diagonal
  dg[[2]] <- setdiff(dg[[2]], 1L); dg[[3]] <- setdiff(dg[[3]], 1L)
  expect_error(nmfkc.net.ecv(Y, rank = 1, folds = dg, diag.exclude = TRUE),
               "diag.exclude")
  expect_error(nmfkc.net.ecv(Y, rank = 1, folds = e$folds, nfolds = 5),
               "differs from length")
  expect_error(nmfkc.net.ecv(Y, rank = 1, folds = e$folds, seeds = 1:2),
               "cannot be combined with folds")
})

test_that("seeds repeats the split; each seed's results are those of seed = s", {
  skip_unless_full()
  Y <- net_data()
  e0 <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 3))
  e1 <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 3, seeds = 123))
  expect_identical(e1$objfunc, e0$objfunc)
  expect_null(e1$objfunc.rep)
  s <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 3, seeds = c(5, 9)))
  s9 <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 3, seed = 9))
  expect_identical(unname(s$objfunc.rep[, 2]), unname(s9$objfunc))
  expect_identical(unname(s$r.squared.cv.rep[, 2]), unname(s9$r.squared.cv))
  expect_equal(unname(s$objfunc), unname(rowMeans(s$objfunc.rep)))
  expect_equal(unname(s$objfunc.sd), unname(apply(s$objfunc.rep, 1, stats::sd)))
  expect_identical(s$folds, s$folds.rep[[1]])
  expect_false(identical(s$folds.rep[[1]], s$folds.rep[[2]]))
  expect_output(print(s), "2 seeds")
})

test_that("ECV: NA in Y is never held out and stays masked; Y.weights is used", {
  skip_unless_full()
  Y <- net_data(); N <- nrow(Y)
  Yna <- Y; Yna[2, 5] <- Yna[5, 2] <- NA
  e <- q(nmfkc.net.ecv(Yna, rank = 1:2, nfolds = 3))
  expect_true(all(is.finite(e$objfunc)))
  expect_false(((5 - 1) * N + 2) %in% unlist(e$folds))
  W0 <- matrix(1, N, N); diag(W0) <- 0
  a <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 3))
  b <- q(nmfkc.net.ecv(Y, rank = 1:2, nfolds = 3, Y.weights = W0))
  expect_false(identical(a$objfunc, b$objfunc))
  expect_false(any(diag(matrix(seq_len(N * N), N)) %in% unlist(b$folds)))
})

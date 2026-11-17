## Tests for nmfkc.net.consensus() (consensus clustering of the nodes of an
## nmfkc.net() fit) and for the shared criteria of nmfkc.consensus().

make_net_blocks <- function(seed = 1, n = 10) {
  set.seed(seed)
  g <- rep(1:3, each = n); N <- length(g)
  P <- ifelse(outer(g, g, "=="), 0.8, 0.05)
  Y <- matrix(stats::rbinom(N * N, 1, P), N, N)
  Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]; diag(Y) <- 0
  dimnames(Y) <- list(paste0("n", 1:N), paste0("n", 1:N))
  Y
}

test_that("nmfkc.net.consensus equals a direct loop over nmfkc.net fits", {
  skip_unless_full()
  Y <- make_net_blocks()
  cs <- suppressMessages(nmfkc.net.consensus(Y, rank = 2:3, nrun = 6, seed = 7,
                                             keep.consensus = TRUE))
  for (ki in 1:2) {
    k <- (2:3)[ki]
    M <- Reduce("+", lapply(1:6, function(r) {
      f <- nmfkc.net(Y, rank = k, type = "tri", X.init = "kmeans++", seed = 7 + 1000L * ki + r)
      outer(f$X.cluster, f$X.cluster, "==") + 0
    })) / 6
    dimnames(M) <- dimnames(Y)
    expect_identical(cs$consensus[[ki]], M)
    lt <- M[lower.tri(M)]
    expect_identical(unname(cs$dispersion[ki]), sum(4 * (M - 0.5)^2) / nrow(M)^2)
    expect_identical(unname(cs$pac[ki]), mean(lt > 0.1 & lt < 0.9))
    d <- as.dist(1 - M)
    if (any(d > 0))
      expect_identical(unname(cs$cophenetic[ki]),
                       cor(as.vector(d), as.vector(cophenetic(hclust(d, "average")))))
  }
})

test_that("three planted blocks: rank 3 is the most stable of ranks 2 to 6, for every type", {
  skip_unless_full()
  Y <- make_net_blocks()
  for (tp in c("tri", "bi", "signed")) {
    cs <- suppressMessages(nmfkc.net.consensus(Y, rank = 2:6, type = tp, nrun = 20))
    expect_equal(cs$rank[which.max(cs$dispersion)], 3)
    expect_true(all(cs$dispersion[-2] < cs$dispersion[2]))   # no tie with rank 2
    expect_lt(cs$pac[2], 0.05)
  }
})

test_that("all three types run; class, ranges, print and plot", {
  skip_unless_full()
  Y <- make_net_blocks(n = 6)
  for (tp in c("tri", "bi", "signed")) {
    cs <- suppressMessages(nmfkc.net.consensus(Y, rank = 2:3, type = tp, nrun = 4,
                                               keep.consensus = TRUE))
    expect_s3_class(cs, "nmfkc.consensus")
    expect_identical(cs$type, tp)
    expect_true(all(cs$dispersion >= 0 & cs$dispersion <= 1))
    expect_true(all(cs$pac >= 0 & cs$pac <= 1))
    expect_true(all(cs$consensus[[2]] >= 0 & cs$consensus[[2]] <= 1))
    expect_identical(rownames(cs$consensus[[1]]), rownames(Y))
    expect_output(print(cs), sprintf("nmfkc.net\\(type = \"%s\"\\)", tp))
  }
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(cs))
  expect_no_error(plot(cs, type = "heatmap", rank = 3))
})

test_that("options reach nmfkc.net; input checks; the caller's stream is untouched", {
  skip_unless_full()
  Y <- make_net_blocks(n = 6)
  a <- suppressMessages(nmfkc.net.consensus(Y, rank = 2, nrun = 3, diag.exclude = TRUE, X.init = "nndsvd",
                                            keep.consensus = TRUE))
  M <- Reduce("+", lapply(1:3, function(r) {
    f <- nmfkc.net(Y, rank = 2, X.init = "nndsvd", diag.exclude = TRUE, seed = 123 + 1000L + r)
    outer(f$X.cluster, f$X.cluster, "==") + 0
  })) / 3
  dimnames(M) <- dimnames(Y)
  expect_identical(a$consensus[[1]], M)
  expect_error(nmfkc.net.consensus(Y[, -1], rank = 2, nrun = 2), "square")
  set.seed(5); u1 <- runif(2)
  set.seed(5); invisible(suppressMessages(nmfkc.net.consensus(Y, rank = 2, nrun = 2))); u2 <- runif(2)
  expect_identical(u1, u2)
})

test_that("cores = 2 gives the same result as cores = 1", {
  skip_unless_full()
  skip_on_os("windows")     # PSOCK workers would load the installed nmfkc
  Y <- make_net_blocks(n = 6)
  a <- suppressMessages(nmfkc.net.consensus(Y, rank = 2:3, nrun = 4, cores = 1, keep.consensus = TRUE))
  b <- suppressMessages(nmfkc.net.consensus(Y, rank = 2:3, nrun = 4, cores = 2, keep.consensus = TRUE))
  expect_identical(a, b)
})

# Display options of the *.DOT writers, taken through `...`:
# nmfkc.DOT(Y.order, X.order, edge.fontsize) and edge.fontsize in every writer.

dot_node_labels <- function(d, prefix) {
  hits <- regmatches(d, gregexpr(paste0(prefix, '_[0-9]+ \\[label="[^"]+"\\]'), d))[[1]]
  sub('.*label="([^"]+)".*', "\\1", hits)
}
dot_edges <- function(d) {
  sort(regmatches(d, gregexpr("[A-Z]+_[0-9]+ -> [A-Z]+_[0-9]+ \\[[^]]+\\]", d))[[1]])
}

test_that("nmfkc.DOT: Y.order and X.order change the written order only", {
  skip_unless_full()
  Yi  <- t(as.matrix(iris[, 1:4]))
  res <- nmfkc(Yi, rank = 2, verbose = FALSE)
  d0  <- nmfkc.DOT(res, type = "YX")
  d1  <- nmfkc.DOT(res, type = "YX", Y.order = 4:1, X.order = 2:1)

  ## the nodes are written in the requested order ...
  expect_identical(dot_node_labels(d0, "Y"), rownames(Yi))
  expect_identical(dot_node_labels(d1, "Y"), rev(rownames(Yi)))
  expect_identical(dot_node_labels(d1, "X"), rev(colnames(res$X)))
  ## ... and the edges (ids, labels, widths) are the same set
  expect_identical(dot_edges(d1), dot_edges(d0))

  ## by name gives the same script as by index
  d2 <- nmfkc.DOT(res, type = "YX", Y.order = rev(rownames(Yi)),
                  X.order = rev(colnames(res$X)))
  expect_identical(as.character(d2), as.character(d1))

  ## hidden (isolated) rows are skipped, the rest keep the requested order
  d3 <- nmfkc.DOT(res, type = "YX", threshold = 0.3, Y.order = 4:1)
  shown <- rownames(Yi)[apply(res$X, 1, max) >= 0.3]
  expect_identical(dot_node_labels(d3, "Y"), rev(shown))

  ## invalid orders
  expect_error(nmfkc.DOT(res, Y.order = c(1, 1, 2, 3)), "Y.order")
  expect_error(nmfkc.DOT(res, Y.order = 1:3), "Y.order")
  expect_error(nmfkc.DOT(res, X.order = c("Basis1", "nope")), "X.order")
})

test_that("nmfkc.DOT: Y.order and X.order with covariates (type YXA)", {
  skip_unless_full()
  Y <- matrix(cars$dist, nrow = 1)
  set.seed(1)
  Y2 <- rbind(Y, Y + abs(rnorm(ncol(Y))), 2 * Y)
  A  <- rbind(1, cars$speed)
  res <- nmfkc(Y2, A, rank = 2, verbose = FALSE)
  d0 <- nmfkc.DOT(res, type = "YXA")
  d1 <- nmfkc.DOT(res, type = "YXA", Y.order = 3:1, X.order = 2:1)
  expect_identical(dot_edges(d1), dot_edges(d0))
  expect_identical(dot_node_labels(d1, "X"), rev(dot_node_labels(d0, "X")))
})

test_that("edge.fontsize: default script unchanged, value reaches every writer", {
  skip_unless_full()
  Yi  <- t(as.matrix(iris[, 1:4]))
  res <- nmfkc(Yi, rank = 2, verbose = FALSE)
  d0  <- nmfkc.DOT(res)
  expect_match(d0, "fontsize=8, arrowhead", fixed = TRUE)
  expect_identical(as.character(nmfkc.DOT(res, edge.fontsize = 8)), as.character(d0))
  expect_match(nmfkc.DOT(res, edge.fontsize = 14), "fontsize=14, arrowhead", fixed = TRUE)
  expect_match(nmfkc.DOT(res, edge.fontsize = 10.5), "fontsize=10.5, arrowhead", fixed = TRUE)
  expect_error(nmfkc.DOT(res, edge.fontsize = -1), "edge.fontsize")
  expect_error(nmfkc.DOT(res, edge.fontsize = c(8, 9)), "edge.fontsize")

  set.seed(1)
  N  <- 20
  Y1 <- matrix(abs(rnorm(4 * N)) + 0.5, 4, N)
  Y2 <- matrix(abs(rnorm(3 * N)) + 0.5, 3, N)
  SM <- 0.5 * (function(M) M + t(M))(matrix(abs(rnorm(36)), 6, 6)); diag(SM) <- 0
  fits <- list(
    net = function(...) nmfkc.net.DOT(suppressWarnings(nmfkc.net(SM, rank = 2, verbose = FALSE)), ...),
    ar  = function(...) {
      a <- nmfkc.ar(Y1, degree = 1, intercept = TRUE)
      nmfkc.ar.DOT(suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE)), ...)
    },
    rrr = function(...) nmf.rrr.DOT(suppressWarnings(nmf.rrr(Y1, Y2, rank1 = 2, maxit = 200, verbose = FALSE)), ...),
    ffb = function(...) nmf.ffb.DOT(suppressWarnings(nmf.ffb(Y1, Y2, rank = 2, verbose = FALSE)), ...)
  )
  for (nm in names(fits)) {
    expect_match(fits[[nm]](), "fontsize=8, arrowhead", fixed = TRUE, info = nm)
    expect_match(fits[[nm]](edge.fontsize = 14), "fontsize=14, arrowhead", fixed = TRUE, info = nm)
  }
})

test_that("unknown options in ... are ignored with a warning", {
  skip_unless_full()
  Yi  <- t(as.matrix(iris[, 1:4]))
  res <- nmfkc(Yi, rank = 2, verbose = FALSE)
  expect_warning(d <- nmfkc.DOT(res, edge.font = 14), "edge.font")
  expect_identical(as.character(d), as.character(nmfkc.DOT(res)))
})

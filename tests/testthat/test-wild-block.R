# ---- block wild bootstrap (wild.unit = "block") -------------------------------
# Added 2026-10-18.  Columns are time points; a block wild bootstrap shares one
# multiplier over max(1, round(n^(1/3))) consecutive columns, so residual serial
# dependence within a block is kept.  A block length of 1 must reproduce
# wild.unit = "column"; the tests force other lengths by mocking the internal
# .wild.block.default().

test_that(".wild.block.multipliers is constant within blocks and +-1", {
  set.seed(11)
  w <- nmfkc:::.wild.block.multipliers(23, 4, "rademacher")
  expect_length(w, 23)
  expect_true(all(w %in% c(-1, 1)))
  r <- rle(w)                      # runs of equal values are unions of blocks
  expect_true(all(r$lengths >= 1))
  ## every block of 4 (after the random offset) is constant: the number of sign
  ## changes is at most the number of block edges
  expect_lte(sum(diff(w) != 0), ceiling(23 / 4))
})

test_that(".wild.block.multipliers with L = 1 draws exactly .wild.multipliers", {
  for (d in c("rademacher", "exp", "mammen")) {
    set.seed(5); a <- nmfkc:::.wild.multipliers(17, d)
    set.seed(5); b <- nmfkc:::.wild.block.multipliers(17, 1, d)
    expect_identical(a, b)
  }
  expect_error(nmfkc:::.wild.block.multipliers(10, 0), "positive integer")
  expect_identical(nmfkc:::.wild.block.default(320), 7L)
  expect_identical(nmfkc:::.wild.block.default(1), 1L)
})

ar_fixture <- function(n) {
  set.seed(3)
  Y <- matrix(abs(rnorm(3 * n)) + 1, 3, n)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  list(a = a, lat = nmfkc.ar.latent(f))
}

test_that("nmfkc.ar.latent.inference: block of length 1 equals column", {
  skip_unless_full()
  x <- ar_fixture(25)
  i_col <- suppressWarnings(nmfkc.ar.latent.inference(x$lat, x$a$Y, x$a$A, wild.B = 8, wild.seed = 1))
  local_mocked_bindings(.wild.block.default = function(n) 1L)
  i_b1  <- suppressWarnings(nmfkc.ar.latent.inference(x$lat, x$a$Y, x$a$A, wild.B = 8, wild.seed = 1,
                                                       wild.unit = "block"))
  expect_identical(i_col$coefficients, i_b1$coefficients)
  expect_identical(i_col$spectral.radius.ci.lower, i_b1$spectral.radius.ci.lower)
  expect_identical(i_b1$wild.unit, "block")
  expect_null(i_b1$wild.block)     # not stored: the length is fixed by n
})

test_that("nmfkc.ar.latent.inference: block bootstrap runs; its length changes the draws", {
  skip_unless_full()
  x <- ar_fixture(30)
  ## n = 29 time points -> length round(29^(1/3)) = 3
  expect_identical(nmfkc:::.wild.block.default(ncol(x$a$Y)), 3L)
  i_b <- suppressWarnings(nmfkc.ar.latent.inference(x$lat, x$a$Y, x$a$A, wild.B = 8, wild.seed = 1,
                                                     wild.unit = "block"))
  expect_true(all(is.finite(i_b$coefficients$BSE)))
  expect_output(print(i_b), "unit = \"block\"")
  ## a longer block draws other multipliers
  local_mocked_bindings(.wild.block.default = function(n) 6L)
  i_b6 <- suppressWarnings(nmfkc.ar.latent.inference(x$lat, x$a$Y, x$a$A, wild.B = 8, wild.seed = 1,
                                                      wild.unit = "block"))
  expect_false(identical(i_b6$coefficients$BSE, i_b$coefficients$BSE))
})

test_that("nmfkc.inference(method = 'refit'): block of length 1 equals column", {
  skip_unless_full()
  set.seed(2)
  Y <- matrix(abs(rnorm(4 * 20)) + 1, 4, 20)
  A <- rbind(1, runif(20))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, verbose = FALSE))
  i_col <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.unit = "column",
                                            wild.B = 10, wild.seed = 4))
  i_b3  <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.unit = "block",
                                            wild.B = 10, wild.seed = 4))   # length round(20^(1/3)) = 3
  expect_true(all(is.finite(i_b3$coefficients$BSE)))
  expect_false(identical(i_col$coefficients, i_b3$coefficients))
  ## the result records the unit (the block length is fixed by N, not stored)
  expect_identical(i_b3$wild.unit, "block")
  expect_null(i_b3$wild.block)
  expect_identical(i_col$wild.unit, "column")
  i_def <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.B = 5))
  expect_identical(i_def$wild.unit, "column")
  i_el <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.unit = "element",
                                           wild.B = 5))
  expect_identical(i_el$wild.unit, "element")
  i_os <- suppressWarnings(nmfkc.inference(f, Y, A))          # onestep: no unit
  expect_identical(i_os$wild.unit, NA_character_)
  local_mocked_bindings(.wild.block.default = function(n) 1L)
  i_b1  <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.unit = "block",
                                            wild.B = 10, wild.seed = 4))
  expect_identical(i_col$coefficients, i_b1$coefficients)
})

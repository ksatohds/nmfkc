# ---- block wild bootstrap (wild.unit = "block") -------------------------------
# Added 2026-10-18.  Columns are time points; a block wild bootstrap shares one
# multiplier over wild.block consecutive columns, so residual serial dependence
# within a block is kept.  wild.block = 1 must reproduce wild.unit = "column".

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

test_that("nmfkc.ar.latent.inference: block with wild.block = 1 equals column", {
  skip_unless_full()
  set.seed(3)
  Y <- matrix(abs(rnorm(3 * 25)) + 1, 3, 25)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  lat <- nmfkc.ar.latent(f)
  i_col <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 8, wild.seed = 1))
  i_b1  <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 8, wild.seed = 1,
                                                       wild.unit = "block", wild.block = 1))
  expect_identical(i_col$coefficients, i_b1$coefficients)
  expect_identical(i_col$spectral.radius.ci.lower, i_b1$spectral.radius.ci.lower)
  expect_null(i_col$wild.block)
  expect_identical(i_b1$wild.unit, "block")
})

test_that("nmfkc.ar.latent.inference: block bootstrap runs and records its length", {
  skip_unless_full()
  set.seed(3)
  Y <- matrix(abs(rnorm(3 * 30)) + 1, 3, 30)
  a <- nmfkc.ar(Y, degree = 1, intercept = TRUE)
  f <- suppressWarnings(nmfkc(a$Y, a$A, rank = 2, verbose = FALSE))
  lat <- nmfkc.ar.latent(f)
  i_b <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 8, wild.seed = 1,
                                                     wild.unit = "block"))
  expect_identical(i_b$wild.block, nmfkc:::.wild.block.default(ncol(a$Y)))
  expect_true(all(is.finite(i_b$coefficients$BSE)))
  expect_output(print(i_b), "unit = \"block\" \\(length")
  ## default length here is round(29^(1/3)) = 3; a longer block draws other multipliers
  i_b6 <- suppressWarnings(nmfkc.ar.latent.inference(lat, a$Y, a$A, wild.B = 8, wild.seed = 1,
                                                      wild.unit = "block", wild.block = 6))
  expect_identical(i_b6$wild.block, 6)
  expect_false(identical(i_b6$coefficients$BSE, i_b$coefficients$BSE))
})

test_that("nmfkc.inference(method = 'refit'): block with wild.block = 1 equals column", {
  skip_unless_full()
  set.seed(2)
  Y <- matrix(abs(rnorm(4 * 20)) + 1, 4, 20)
  A <- rbind(1, runif(20))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, verbose = FALSE))
  i_col <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.unit = "column",
                                            wild.B = 10, wild.seed = 4))
  i_b1  <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.unit = "block",
                                            wild.block = 1, wild.B = 10, wild.seed = 4))
  expect_identical(i_col$coefficients, i_b1$coefficients)
  i_b4  <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.unit = "block",
                                            wild.block = 4, wild.B = 10, wild.seed = 4))
  expect_true(all(is.finite(i_b4$coefficients$BSE)))
})

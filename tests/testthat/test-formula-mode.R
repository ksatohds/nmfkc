# test-formula-mode.R
# Tests for Formula Mode support across nmfkc functions

# --- Helper data ---
# Create a small non-negative dataset for testing
make_test_data <- function() {
  set.seed(42)
  data.frame(
    Y1 = rpois(20, 5),
    Y2 = rpois(20, 10),
    Y3 = rpois(20, 8),
    A1 = abs(rnorm(20, 5, 1)),
    A2 = abs(rnorm(20, 3, 1))
  )
}

# ============================================================
# 1. .nmfkc_parse_formula (internal function)
# ============================================================

test_that(".nmfkc_parse_formula: MODE 1 (data provided) basic", {
  df <- make_test_data()
  res <- nmfkc:::.nmfkc_parse_formula(Y1 + Y2 + Y3 ~ A1 + A2, data = df)

  expect_true(is.matrix(res$Y))
  expect_true(is.matrix(res$A))
  # Y is P x N (3 variables x 20 samples)
  expect_equal(nrow(res$Y), 3)
  expect_equal(ncol(res$Y), 20)
  # A is R x N (2 covariates x 20 samples)
  expect_equal(nrow(res$A), 2)
  expect_equal(ncol(res$A), 20)
  # Column names returned
  expect_equal(res$Y_cols, c("Y1", "Y2", "Y3"))
  expect_equal(res$A_cols, c("A1", "A2"))
})

test_that(".nmfkc_parse_formula: MODE 1 dot notation on RHS", {
  df <- make_test_data()
  res <- nmfkc:::.nmfkc_parse_formula(Y1 + Y2 + Y3 ~ ., data = df)

  expect_equal(nrow(res$A), 2)  # A1 and A2 remain
  expect_equal(res$A_cols, c("A1", "A2"))
})

test_that(".nmfkc_parse_formula: MODE 1 dot notation on LHS", {
  df <- make_test_data()
  res <- nmfkc:::.nmfkc_parse_formula(. ~ A1 + A2, data = df)

  expect_equal(nrow(res$Y), 3)  # Y1, Y2, Y3 remain
})

test_that(".nmfkc_parse_formula: MODE 1 A omitted (standard NMF)", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc:::.nmfkc_parse_formula(Y1 + Y2 + Y3 ~ 0, data = df))

  expect_null(res$A)
  expect_equal(nrow(res$Y), 3)
})

test_that(".nmfkc_parse_formula: MODE 2 (direct matrix eval)", {
  # Both matrices: rows=samples, cols=variables. After t(): variables x samples.
  Y_mat <- matrix(1:6, nrow = 3, ncol = 2)   # 3 samples, 2 variables
  A_mat <- matrix(1:6, nrow = 3, ncol = 2)   # 3 samples, 2 covariates
  res <- nmfkc:::.nmfkc_parse_formula(Y_mat ~ A_mat)

  expect_true(is.matrix(res$Y))
  expect_true(is.matrix(res$A))
  # Mode 2: Y_cols and A_cols are NULL
  expect_null(res$Y_cols)
  expect_null(res$A_cols)
})

# ============================================================
# 2. nmfkc() with formula
# ============================================================

test_that("nmfkc: formula mode matches matrix mode", {
  df <- make_test_data()
  Y_mat <- t(as.matrix(df[, c("Y1", "Y2", "Y3")]))
  A_mat <- t(as.matrix(df[, c("A1", "A2")]))

  set.seed(1)
  res_mat <- suppressMessages(nmfkc::nmfkc(Y_mat, A_mat, rank = 2,
                                           epsilon = 1e-4, maxit = 100,
                                           print.dims = FALSE, seed = 1))
  set.seed(1)
  res_frm <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                                           rank = 2, epsilon = 1e-4, maxit = 100,
                                           print.dims = FALSE, seed = 1))

  # Same optimization result
  expect_equal(res_mat$r.squared, res_frm$r.squared, tolerance = 1e-3)
  expect_equal(res_mat$rank, res_frm$rank)
})

test_that("nmfkc: formula mode stores formula.meta", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                                       rank = 2, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  expect_false(is.null(res$formula.meta))
  expect_true(inherits(res$formula.meta$formula, "formula"))
  expect_equal(res$formula.meta$Y_cols, c("Y1", "Y2", "Y3"))
  expect_equal(res$formula.meta$A_cols, c("A1", "A2"))
})

test_that("nmfkc: formula mode without A (standard NMF)", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ 0, data = df,
                                       rank = 2, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  expect_false(is.null(res$formula.meta))
  expect_null(res$formula.meta$A_cols)
})

test_that("nmfkc: formula mode rejects negative values", {
  df <- make_test_data()
  df$Y1[1] <- -1
  expect_error(
    suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                                  rank = 2, print.dims = FALSE)),
    "non-negative"
  )
})

test_that("nmfkc: matrix mode formula.meta is NULL", {
  Y <- matrix(rpois(30, 5), nrow = 3)
  res <- suppressMessages(nmfkc::nmfkc(Y, rank = 2, epsilon = 1e-4,
                                       maxit = 100, print.dims = FALSE))
  expect_null(res$formula.meta)
})

# ============================================================
# 3. predict.nmfkc() with newdata
# ============================================================

test_that("predict: newdata works for formula mode models", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                                       rank = 2, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  new_df <- data.frame(A1 = c(5, 6), A2 = c(3, 4))
  pred <- predict(res, newdata = new_df)

  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), 3)   # P = 3 variables
  expect_equal(ncol(pred), 2)   # 2 new samples
})

test_that("predict: newdata errors for matrix-mode models", {
  Y <- matrix(rpois(30, 5), nrow = 3)
  res <- suppressMessages(nmfkc::nmfkc(Y, rank = 2, epsilon = 1e-4,
                                       maxit = 100, print.dims = FALSE))
  expect_error(
    predict(res, newdata = data.frame(A1 = 1)),
    "formula"
  )
})

test_that("predict: newdata errors when columns are missing", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                                       rank = 2, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  expect_error(
    predict(res, newdata = data.frame(A1 = c(5, 6))),
    "missing required columns"
  )
})

test_that("predict: newdata + newA gives warning and uses newdata", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                                       rank = 2, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  new_df <- data.frame(A1 = c(5, 6), A2 = c(3, 4))
  dummy_A <- matrix(1, nrow = 2, ncol = 2)

  expect_warning(
    pred <- predict(res, newA = dummy_A, newdata = new_df),
    "newdata"
  )
  expect_true(is.matrix(pred))
})

test_that("predict: type='class' works with newdata", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                                       rank = 2, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  new_df <- data.frame(A1 = c(5, 6), A2 = c(3, 4))
  pred <- predict(res, newdata = new_df, type = "class")

  expect_true(is.character(pred))
  expect_equal(length(pred), 2)
})

test_that("predict: backward compatible (newA still works)", {
  df <- make_test_data()
  Y_mat <- t(as.matrix(df[, c("Y1", "Y2", "Y3")]))
  A_mat <- t(as.matrix(df[, c("A1", "A2")]))

  res <- suppressMessages(nmfkc::nmfkc(Y_mat, A_mat, rank = 2,
                                       epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  newA_mat <- matrix(abs(rnorm(4)), nrow = 2, ncol = 2)
  pred <- predict(res, newA = newA_mat)

  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 2)
})

# ============================================================
# 4. CV functions with formula
# ============================================================

test_that("nmfkc.cv: formula mode works", {
  df <- make_test_data()
  res <- suppressWarnings(suppressMessages(
    nmfkc::nmfkc.cv(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                     Q = 2, div = 3, maxit = 100,
                     print.dims = FALSE, print.trace = FALSE)
  ))
  expect_true(is.list(res))
  expect_true(!is.null(res$objfunc))
  expect_true(is.finite(res$objfunc))
})

test_that("nmfkc.ecv: formula mode works", {
  df <- make_test_data()
  res <- suppressWarnings(suppressMessages(
    nmfkc::nmfkc.ecv(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                      Q = 1:2, div = 3, maxit = 100,
                      print.dims = FALSE, print.trace = FALSE,
                      save.time = TRUE)
  ))
  expect_true(is.list(res))
  expect_equal(length(res$objfunc), 2)
})

test_that("nmfkc.rank: formula mode works", {
  df <- make_test_data()
  res <- suppressWarnings(suppressMessages(
    nmfkc::nmfkc.rank(Y1 + Y2 + Y3 ~ A1 + A2, data = df,
                       rank = 1:2, save.time = TRUE,
                       plot = FALSE, maxit = 100,
                       print.dims = FALSE, print.trace = FALSE)
  ))
  expect_true(is.list(res))
  expect_true(!is.null(res$rank.best))
})

# ============================================================
# 5. Edge cases
# ============================================================

test_that("nmfkc: single Y variable formula", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 ~ A1 + A2, data = df,
                                       rank = 1, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  expect_equal(nrow(res$X), 1)  # P = 1
})

test_that("predict: newdata with formula mode (no A = standard NMF) errors", {
  df <- make_test_data()
  res <- suppressMessages(nmfkc::nmfkc(Y1 + Y2 + Y3 ~ 0, data = df,
                                       rank = 2, epsilon = 1e-4, maxit = 100,
                                       print.dims = FALSE))
  expect_error(
    predict(res, newdata = data.frame(A1 = 5)),
    "no A columns"
  )
})


test_that("SE is the sandwich SE and BSE the bootstrap one, not the same number", {
  set.seed(4)
  Y <- matrix(abs(rnorm(5 * 40)) + 1, 5, 40)
  A <- rbind(1, abs(rnorm(40)))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE))
  inf <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit",
                                          wild.B = 40, wild.seed = 1))
  ct <- inf$coefficients
  expect_true("on.boundary" %in% names(ct))
  expect_false(isTRUE(all.equal(ct$SE, ct$BSE)))
  expect_equal(ct$SE, as.vector(inf$C.se.sandwich), tolerance = 1e-12)
  expect_equal(ct$BSE, as.vector(inf$C.se.boot), tolerance = 1e-12)
  ## z still uses the primary SE, which is the bootstrap one under "refit"
  expect_equal(ct$z_value, ifelse(ct$BSE > 0, ct$Estimate / ct$BSE, NA_real_),
               tolerance = 1e-12)
  ## and the boundary flag matches the estimates
  expect_identical(ct$on.boundary, ct$Estimate <= 1e-3)
})

test_that("boundary.tol controls the on.boundary flag", {
  set.seed(4)
  Y <- matrix(abs(rnorm(5 * 40)) + 1, 5, 40)
  A <- rbind(1, abs(rnorm(40)))
  f <- suppressWarnings(nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE))
  inf <- suppressWarnings(nmfkc.inference(f, Y, A, method = "refit", wild.B = 20,
                                          wild.seed = 1, boundary.tol = 1))
  expect_true(all(inf$coefficients$on.boundary[inf$coefficients$Estimate <= 1]))
})

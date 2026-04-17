## Tests for nmfkc.kernel.beta.nearest.med() candidates option

test_that("candidates argument controls the beta grid", {
  set.seed(1)
  U <- matrix(stats::rnorm(5 * 50), 5, 50)

  ## default = "7points"
  info7 <- nmfkc.kernel.beta.nearest.med(U)
  expect_length(info7$beta_candidates, 7)
  expect_equal(info7$beta_candidates[4], info7$beta, tolerance = 1e-12)

  ## "4points"
  info4 <- nmfkc.kernel.beta.nearest.med(U, candidates = "4points")
  expect_length(info4$beta_candidates, 4)

  ## custom numeric t-grid
  t_custom <- c(-0.5, 0, 0.5)
  info_c <- nmfkc.kernel.beta.nearest.med(U, candidates = t_custom)
  expect_length(info_c$beta_candidates, 3)
  expect_equal(info_c$beta_candidates,
               info_c$beta * 10^(-2 * t_custom), tolerance = 1e-12)

  ## invalid character -> error
  expect_error(nmfkc.kernel.beta.nearest.med(U, candidates = "bogus"),
               "7points|4points|numeric")

  ## works with Uk supplied as well (7 candidates by default)
  Uk <- U[, 1:5, drop = FALSE]
  infoUk <- nmfkc.kernel.beta.nearest.med(U, Uk = Uk)
  expect_length(infoUk$beta_candidates, 7)
})

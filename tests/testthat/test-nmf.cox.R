test_that("nmf.cox fits, follows house-style returns, and S3 methods work", {
  skip_if_not_installed("survival")
  d <- survival::veteran
  d$y <- survival::Surv(d$time, d$status)

  fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
                 rank = 2, X.L2.smooth = 1000, verbose = FALSE, maxit = 20)
  expect_s3_class(fit, "nmf.cox")
  ## house-style fields
  expect_true(all(c("C", "X", "beta.t", "objfunc", "objfunc.iter", "iter",
                    "converged", "stop.reason", "runtime", "dims",
                    "X.cluster", "gamma") %in% names(fit)))
  expect_null(fit$Theta)          # parameter matrix is C, not Theta
  expect_null(fit$alpha)          # duplicate alias removed
  ## shapes and constraints
  expect_true(all(fit$X >= 0))
  expect_equal(dim(fit$C), c(2L, ncol(fit$beta.t)))
  expect_equal(unname(fit$beta.t), unname(fit$X %*% fit$C))
  expect_equal(colnames(fit$X), paste0("Basis", 1:2))
  ## S3
  expect_named(coef(fit))
  expect_equal(fitted(fit), fit$beta.t)
  expect_s3_class(summary(fit), "summary.nmf.cox")
  expect_silent(plot(fit))
})

test_that("nmf.cox.inference adds the C coefficients table and Wald tests", {
  skip_if_not_installed("survival")
  d <- survival::veteran
  d$y <- survival::Surv(d$time, d$status)

  fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
                 rank = 2, X.L2.smooth = 1000, verbose = FALSE, maxit = 20)
  inf <- nmf.cox.inference(fit, y ~ trt + age, data = d, A = ~ karno + celltype)
  expect_s3_class(inf, "nmf.cox.inference")
  expect_s3_class(inf, "nmf.cox")
  ## house-style coefficients table
  expect_true(all(c("Basis", "Covariate", "Estimate", "SE", "z_value", "p_value")
                  %in% names(inf$coefficients)))
  expect_equal(nrow(inf$coefficients), prod(dim(fit$C)))
  expect_true(any(grepl("^Basis", inf$coefficients$Basis)))
  ## Wald + pointwise SE
  expect_equal(nrow(inf$wald), ncol(fit$beta.t))
  expect_equal(dim(inf$se.beta.t), dim(fit$beta.t))
  expect_true(all(is.finite(inf$se.beta.t)))
})

test_that("nmf.cox.cv sweeps rank/X.L2.smooth vectors (house style)", {
  skip_if_not_installed("survival")
  d <- survival::veteran
  d$y <- survival::Surv(d$time, d$status)

  cv <- nmf.cox.cv(y ~ trt + age, data = d, A = ~ karno + celltype,
                   rank = 2, X.L2.smooth = c(300, 1000), criterion = "aic",
                   verbose = FALSE, maxit = 15)
  expect_s3_class(cv, "nmf.cox.cv")
  expect_true(all(c("rank.best", "X.L2.smooth.best") %in% names(cv)))
  expect_true(cv$X.L2.smooth.best %in% c(300, 1000))
  expect_equal(dim(cv$aic), c(1L, 2L))
})

test_that("nmf.cox.cf and nmf.cox.phtest run with reproducible folds", {
  skip_if_not_installed("survival")
  d <- survival::veteran
  d$y <- survival::Surv(d$time, d$status)

  cf <- nmf.cox.cf(y ~ trt + age, data = d, A = ~ karno + celltype,
                   rank = 2, X.L2.smooth = 1000, nfolds = 3, seed = 123, maxit = 15)
  expect_s3_class(cf, "nmf.cox.cf")
  expect_named(cf$gamma)
  expect_true(all(is.finite(cf$se)))
  expect_null(cf$alpha)

  ph <- nmf.cox.phtest(y ~ trt + age, data = d, A = ~ karno + celltype,
                       rank = 2, X.L2.smooth = 1000, nfolds = 3, seed = 123, maxit = 15)
  expect_s3_class(ph, "nmf.cox.phtest")
  expect_true(all(is.finite(ph$wald$chisq)))
  expect_equal(dim(ph$C), c(2L, ncol(ph$beta.t)))  # parameter matrix is C (Q x R)
  expect_null(ph$Theta)
})

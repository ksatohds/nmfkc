# Random Fourier Features for nmfkc.signed()

Generates RFF random parameters \\\omega_d \sim \mathcal{N}(0, 2\beta
I_p)\\, \\b_d \sim \mathrm{Uniform}(0, 2\pi)\\ (Rahimi & Recht, 2007)
and applies the RFF transform \$\$z_d(u) = \sqrt{2/D}\\
\cos(\omega_d^\top u + b_d)\$\$ to each column of `U`, yielding a
sign-unrestricted \\D \times N\\ feature matrix \\Z\\ such that \\Z^\top
Z \approx K\\, the Gaussian kernel matrix with bandwidth `beta`.

The return value is a list with the feature matrix `Z` and the
generating parameters `pars = list(omega, b, D, beta)` so that the same
random map can be re-applied to new data (by passing `pars` back) and
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
can record the parameters for downstream summary.

## Usage

``` r
nmfkc.signed.rff(U, beta = NULL, D = ceiling(ncol(U)/2), seed = NULL, ...)
```

## Arguments

- U:

  A \\p \times N\\ numeric matrix; columns are data points.

- beta:

  Positive scalar. Gaussian kernel bandwidth parameter. Can be obtained
  via
  [`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md).
  May be `NULL` only when `pars` is supplied through `...`.

- D:

  Integer. Number of random features. Defaults to
  `ceiling(ncol(U) / 2)`. This default is intended for the
  **training-time fresh generation** only; for test data, always supply
  `pars` via `...` to inherit the training \\D\\ together with \\\omega,
  b\\. For very large \\N\\ the default may be excessive (direct MU cost
  is \\O(QD^2)\\); choose a smaller \\D\\ manually. For very small
  \\N\\, RFF is not recommended; use a full kernel matrix with
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
  instead.

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before generating
  \\\omega, b\\, for reproducibility. Ignored when `pars` is supplied.

- ...:

  Hidden option `pars`: a list `list(omega, b, D, beta)` obtained from a
  previous call (`Ztrain$pars`). When supplied, \\\omega, b\\ are reused
  and the `beta`, `D`, `seed` arguments are ignored. Use this to apply
  the same random map to test data.

## Value

A list with two elements:

- `Z`:

  A \\D \times N\\ sign-unrestricted numeric matrix. Pass this to
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  as its `A` argument.

- `pars`:

  A list `list(omega, b, D, beta)`. Pass this to
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  via its `pars` argument (for summary display) and to subsequent
  `nmfkc.signed.rff()` calls (to reuse the same random map on new data).

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md),
[`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md)

## Examples

``` r
# \donttest{
## Iris 3-class classification with RFF + direct MU (Ding-Li-Jordan)
data(iris)
set.seed(1)
idx <- sample(nrow(iris), 100)    # 100 training, 50 test

## Scale features using TRAINING mean/sd; transpose to p x N layout
mn      <- colMeans(iris[idx, 1:4])
sc      <- apply(iris[idx, 1:4], 2, sd)
U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))   # 4 x 100
U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))   # 4 x  50

## One-hot encode training labels as a Q_obs x N target matrix
levs      <- levels(iris$Species)
Y.train   <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
rownames(Y.train) <- levs                                       # 3 x 100
lab.train <- iris$Species[idx]
lab.test  <- iris$Species[-idx]

## Beta candidates from nearest-neighbour median heuristic
beta_info <- nmfkc.kernel.beta.nearest.med(U.train)
betas     <- beta_info$beta_candidates

## CV over beta candidates: for each beta, generate RFF, fit, and
## evaluate column-wise CV-MSE on training data
cv_mse <- numeric(length(betas))
for (i in seq_along(betas)) {
  rff_i   <- nmfkc.signed.rff(U.train, beta = betas[i], D = 50, seed = 1)
  cv_i    <- nmfkc.signed.cv(Y.train, A = rff_i$Z, rank = 3, seed = 123)
  cv_mse[i] <- cv_i$objfunc
}
beta_best <- betas[which.min(cv_mse)]

## Generate signed RFF features with the best beta
rff.train <- nmfkc.signed.rff(U.train, beta = beta_best, D = 50, seed = 1)
rff.test  <- nmfkc.signed.rff(U.test, pars = rff.train$pars)

## Fit on training data only
res <- nmfkc.signed(Y.train, A = rff.train$Z, rank = 3,
                    pars = rff.train$pars, verbose = FALSE)

## Predict on training and test separately
pred.train <- predict(res, newA = rff.train$Z, type = "class")
pred.test  <- predict(res, newA = rff.test$Z,  type = "class")
mean(pred.train == as.character(lab.train))
#> [1] 0.98
mean(pred.test  == as.character(lab.test))
#> [1] 1
# }
```

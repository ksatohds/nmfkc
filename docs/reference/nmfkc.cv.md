# Perform k-fold cross-validation for NMF with kernel covariates

`nmfkc.cv` performs k-fold cross-validation for the tri-factorization
model \\Y \approx X C A = X B\\, where

- \\Y(P,N)\\ is the observation matrix,

- \\A(R,N)\\ is the covariate (or kernel) matrix,

- \\X(P,Q)\\ is the basis matrix,

- \\C(Q,R)\\ is the parameter matrix, and

- \\B(Q,N)\\ is the coefficient matrix (\\B = C A\\).

Given \\Y\\ (and optionally \\A\\), \\X\\ and \\C\\ are fitted on each
training split and predictive performance is evaluated on the held-out
split.

## Usage

``` r
nmfkc.cv(Y, A = NULL, rank = 2, data, ...)
```

## Arguments

- Y:

  Observation matrix, or a formula (see
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) for
  Formula Mode).

- A:

  Covariate matrix. If `NULL`, the identity matrix is used. Ignored when
  `Y` is a formula.

- rank:

  Rank of the basis matrix \\X\\. Default is 2.

- data:

  A data frame (required when `Y` is a formula with column names).

- ...:

  Additional arguments controlling CV and the internal
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) call:

  `Y.weights`

  :   Non-negative weight matrix or vector
      ([`lm()`](https://rdrr.io/r/stats/lm.html)-style: loss \\\sum W \\
      r^2\\). Binary `{0,1}` masks (`TRUE` / `FALSE` also accepted) are
      the typical ECV usage – 0/`FALSE` excludes an element. See
      [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) for
      full details.

  `div`

  :   Number of folds (\\k\\); default: `5`.

  `seed`

  :   Integer seed for reproducible partitioning; default: `123`.

  `shuffle`

  :   Logical. If `TRUE` (default), randomly shuffles samples (standard
      CV); if `FALSE`, splits sequentially (block CV; recommended for
      time series).

  `Q`

  :   (Deprecated) Alias for `rank`.

  *Arguments passed to* [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)

  :   e.g., `C.L1`, `epsilon`, `maxit`, `method` (`"EU"` or `"KL"`),
      `X.restriction`, `X.init`, etc.

## Value

A list with components:

- `objfunc`:

  Mean loss per valid entry over all folds (MSE for `method="EU"`).

- `sigma`:

  Residual standard error (RMSE). Available only if `method="EU"`; on
  the same scale as `Y`.

- `objfunc.block`:

  Loss for each fold.

- `block`:

  Vector of fold indices (1, …, `div`) assigned to each column of \\Y\\.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.cv.md),
[`nmfkc.ar.degree.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.degree.cv.md)

## Examples

``` r
# Example 1 (with explicit covariates):
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
res <- nmfkc.cv(Y, A, rank = 1)
res$objfunc
#> [1] 286.3378

# Example 2 (kernel A and beta sweep):
Y <- matrix(cars$dist, nrow = 1)
U <- matrix(c(5, 10, 15, 20, 25), nrow = 1)
V <- matrix(cars$speed, nrow = 1)
betas <- 25:35/1000
obj <- numeric(length(betas))
for (i in seq_along(betas)) {
  A <- nmfkc.kernel(U, V, beta = betas[i])
  obj[i] <- nmfkc.cv(Y, A, rank = 1, nfolds = 10)$objfunc
}
betas[which.min(obj)]
#> [1] 0.031
```

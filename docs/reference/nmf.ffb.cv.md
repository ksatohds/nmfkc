# Cross-Validation for NMF-FFB

Performs K-fold cross-validation to evaluate the equilibrium mapping of
the NMF-FFB model.

For each fold, `nmf.ffb` is fitted on the training samples, yielding an
equilibrium mapping \\\hat Y_1 = M\_{\mathrm{model}} Y_2\\. The held-out
endogenous variables \\Y_1\\ are then predicted from \\Y_2\\ using this
mapping, and the mean absolute error (MAE) over all entries in the test
block is computed. The returned value is the average MAE across folds.

This implements the hyperparameter selection strategy described in the
paper: hyperparameters are chosen by predictive cross-validation rather
than direct inspection of the internal structural matrices.

**`method = "fiml"`.** When `method = "fiml"` is passed (through `...`),
the column-wise scheme above is not run: the equilibrium prediction
\\M\_{\mathrm{model}} Y_2\\ is a reduced-form quantity, and the reduced
form is the same with and without feedback, so predicting held-out
columns cannot select \\\Theta_1\\ (that is what the BIC path and the
bootstrap-calibrated LR test in
[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md) /
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)
are for). The only tunable quantity is the stage-1 rank, and the call is
delegated to the element-wise CV of the feed-forward fit,
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)`(Y = Y1, A = Y2, rank = rank, ...)`,
whose object is returned (`$objfunc` holds the CV error per rank; `rank`
may be a vector, default `1:min(3, P1)`). `C1.L1` and `C2.L1` are
ignored in that case.

## Usage

``` r
nmf.ffb.cv(
  Y1,
  Y2,
  rank = NULL,
  X.init = "nndsvd",
  X.L2.ortho = 100,
  C1.L1 = 1,
  C2.L1 = 0.1,
  epsilon = 1e-06,
  maxit = 5000,
  ...
)
```

## Arguments

- Y1:

  A non-negative numeric matrix of endogenous variables with **rows =
  variables (P1), columns = samples (N)**.

- Y2:

  A non-negative numeric matrix of exogenous variables with **rows =
  variables (P2), columns = samples (N)**. Must satisfy
  `ncol(Y1) == ncol(Y2)`.

- rank:

  Integer; rank (number of latent factors) passed to `nmf.ffb`. If
  `NULL`, `nmf.ffb` decides the effective rank (via `...` or
  `nrow(Y2)`).

- X.init:

  Initialization strategy for `X`, forwarded to
  [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md).
  One of `"nndsvd"` (default), `"kmeans"`, `"kmeansar"`, `"runif"`, a
  numeric \\P_1 \times Q\\ matrix, or `NULL` (alias for `"nndsvd"`). See
  [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md) for
  details.

- X.L2.ortho:

  L2 orthogonality penalty for `X`.

- C1.L1:

  L1 sparsity penalty for `C1` (\\\Theta_1\\).

- C2.L1:

  L1 sparsity penalty for `C2` (\\\Theta_2\\).

- epsilon:

  Convergence threshold for `nmf.ffb`.

- maxit:

  Maximum number of iterations for `nmf.ffb`.

- ...:

  Additional arguments passed to `nmf.ffb` (except for `rank`, `seed`,
  `div`, `shuffle`, which are handled here). Also accepts: `nfolds`
  (number of folds, default 5; `div` also accepted), `seed` (master
  random seed, default `NULL`), `shuffle` (logical, default `TRUE`),
  `cores` (integer; number of parallel workers for the fold loop,
  default `getOption("mc.cores", 1L)`). Fold partitions and per-fold
  seeds are drawn before the loop and each `nmf.ffb` fit self-seeds, so
  results are identical to the sequential run for any `cores`.

## Value

A numeric scalar: mean MAE across CV folds (the fits use `method = "mu"`
unless `method` is given). With `method = "fiml"`, the
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)
object of the stage-1 rank sweep.

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md)

## Examples

``` r
Y <- t(iris[, -5])
Y1 <- Y[1:2, ]
Y2 <- Y[3:4, ]
mae <- nmf.ffb.cv(Y1, Y2, rank = 2, maxit = 500, nfolds = 3)
#> Warning: maximum iterations (500) reached...
#> Warning: maximum iterations (500) reached...
#> Warning: maximum iterations (500) reached...
mae
#> [1] 1.656582
# \donttest{
# rank selection for the likelihood-based estimator
ecv <- nmf.ffb.cv(Y1, Y2, rank = 1:2, method = "fiml", nfolds = 3)
#> Performing Element-wise CV for Q = 1,2 (3-fold)...
ecv$objfunc
#>      Q=1      Q=2 
#> 2.915778 2.915701 
# }
```

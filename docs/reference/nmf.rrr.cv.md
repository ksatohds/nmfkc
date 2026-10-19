# Sample-wise k-fold Cross-Validation for nmfae

`nmfae.cv` performs k-fold cross-validation by splitting columns
(samples) of \\Y_1\\ and \\Y_2\\ into `div` folds. For each fold, the
model \\Y_1 \approx X_1 \Theta X_2 Y_2\\ is fitted on the training
samples and predictive performance is evaluated on the held-out samples.

When `Y2` is a kernel matrix created by
[`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md)
(detected via attributes), the symmetric kernel splitting convention is
used: `Y2[train, train]` for training and `Y2[train, test]` for
prediction.

## Usage

``` r
nmf.rrr.cv(Y1, Y2 = Y1, rank1 = 2, rank2 = NULL, ...)
```

## Arguments

- Y1:

  Output matrix \\Y_1\\ (P1 x N). Non-negative.

- Y2:

  Input matrix \\Y_2\\ (P2 x N), or a kernel matrix (N x N). Default is
  `Y1` (autoencoder).

- rank1:

  Integer. Rank of the response basis. Default is 2.

- rank2:

  Integer. Rank of the covariate basis. Default (`NULL`) = `rank1`.

- ...:

  Additional arguments passed to
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md)
  (e.g., `epsilon`, `maxit`, `Y1.weights`). Also accepts: `nfolds`
  (number of folds, default 5; `div` also accepted), `seed` (integer
  seed, default 123), `shuffle` (logical, default `TRUE`), and `cores`
  (fit the folds in parallel; default `getOption("mc.cores", 1L)`, PSOCK
  on Windows / forking elsewhere). Each fold is an independent
  self-seeded fit and the per-fold errors are summed in order, so the
  result is identical for any `cores`. For backward compatibility, `Q`,
  `R` are accepted as aliases for `rank`, `rank.encoder`.

  Rank aliases accepted here for backward compatibility: `Q` for
  `rank1`, `R` for `rank2`.

## Value

A list with components:

- objfunc:

  Mean squared error per valid element over all folds.

- sigma:

  Residual standard error (RMSE), same scale as \\Y_1\\.

- objfunc.block:

  Per-fold squared error totals.

- block:

  Integer vector of fold assignments (1, ..., `div`) for each column.

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`nmf.rrr.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md),
[`nmf.rrr.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.kernel.beta.cv.md),
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md)

## Examples

``` r
Y <- t(iris[1:30, 1:4])
res <- nmf.rrr.cv(Y, rank1 = 2, rank2 = 2, nfolds = 5, maxit = 500)
res$sigma
#> [1] 0.1429115
```

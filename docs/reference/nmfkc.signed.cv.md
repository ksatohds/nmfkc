# Column-wise k-fold cross-validation for nmfkc.signed

Column-wise k-fold CV by held-out samples: for each fold, the model is
fit on the training columns and evaluated on the held-out columns by
solving for new-sample coefficients via a weighted refit with \\X\\
fixed.

## Usage

``` r
nmfkc.signed.cv(Y, A, rank = 2, ...)
```

## Arguments

- Y:

  Real-valued \\Q\_{\mathrm{obs}} \times N\\ response matrix (signed
  entries allowed).

- A:

  Real-valued \\D \times N\\ covariate matrix (signed).

- rank:

  Integer \\Q\\.

- ...:

  Passed to
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md);
  also accepts `nfolds` (default 5; `div` alias), `seed` (default 123),
  `shuffle` (default `TRUE`).

## Value

A list with `objfunc` (mean squared prediction error), `sigma` (RMSE),
`objfunc.block` (per-fold MSE vector), `block` (integer fold assignment
of length \\N\\). Field names match
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md).

## Lifecycle

This function is **experimental**.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md),
[`nmfkc.signed.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.ecv.md)

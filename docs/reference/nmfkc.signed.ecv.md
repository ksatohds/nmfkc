# Element-wise cross-validation for nmfkc.signed

Element-wise k-fold CV (Wold's CV): held-out elements are masked via
`Y.weights = 0` during fitting, and the RMSE on those elements is
reported. Loops over candidate `rank` values.

## Usage

``` r
nmfkc.signed.ecv(Y, A, rank = 1:3, ...)
```

## Arguments

- Y:

  Real-valued \\Q\_{\mathrm{obs}} \times N\\ response matrix (signed
  entries allowed).

- A:

  Real-valued \\D \times N\\ covariate matrix (signed).

- rank:

  Integer vector of candidate ranks (default `1:3`).

- ...:

  Passed to
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md);
  also accepts `nfolds` (default 5; `div` alias), `seed` (default 123).
  `Y.weights` (a matrix, one weight per column, or a scalar) is kept in
  every fit; cells of weight 0 are never held out, and the held-out loss
  is the weighted mean. (Before 1.0.0 a `Y.weights` given here was
  dropped.)

## Value

A list with `objfunc` (MSE per rank), `sigma` (RMSE), `objfunc.fold`
(per-fold per-rank), `folds`, `Q.grid`.

## Lifecycle

This function is **experimental**.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md),
[`nmfkc.signed.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.cv.md)

# Element-wise Cross-Validation for Signed-Bottleneck NMF-AE

Element-wise k-fold cross-validation for
[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)
to select the decoder / encoder ranks \\(Q, R)\\. Mirrors
[`nmf.rrr.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md)
but uses the **weighted** Signed-Bottleneck NMF-AE fit path
(`Y1.weights`): test-fold elements are zero-weighted during fitting, and
held-out MSE is computed on those elements.

## Usage

``` r
nmf.rrr.signed.ecv(Y1, Y2 = Y1, rank1 = 1:2, rank2 = NULL, ...)
```

## Arguments

- Y1:

  Output matrix (P1 x N).

- Y2:

  Input matrix (P2 x N). Default `Y1`.

- rank1:

  Integer vector of candidate response-basis ranks. Default `1:2`.

- rank2:

  Integer vector of candidate covariate-basis ranks, or `NULL` (default:
  pair rank2 = rank1, diagonal grid).

- ...:

  Additional arguments:

  `nfolds` / `div`

  :   Number of folds. Default 5.

  `seed`

  :   RNG seed for fold assignment. Default 123.

  `Y1.weights`

  :   A matrix, one weight per column, or a scalar; kept in every fit.
      Cells of weight 0 are never held out, and the held-out loss is the
      weighted mean. (Before 1.0.0 it was silently replaced by the fold
      mask.)

  `cores`

  :   Evaluate the \\(Q,R)\\-pair \\\times\\ fold grid in parallel.
      Default `getOption("mc.cores", 1L)` (PSOCK on Windows, forking
      elsewhere). Each task is an independent self-seeded fit and
      results are aggregated in order, so the returned object is
      identical for any `cores`.

  `nstart`

  :   Number of random restarts per fit. Default 1. Signed models have
      more local minima (the bottleneck can carry both signs), so
      `nstart >= 10` is recommended for reproducible rank selection.

  Other args

  :   `epsilon`, `maxit`, `warm.start`, etc.\\ are passed to
      [`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md).

  Rank aliases accepted here for backward compatibility: `Q` for
  `rank1`, `R` for `rank2`.

## Value

An object of class `c("nmfae.signed.ecv", "nmfae.ecv")` with `objfunc`
(MSE per pair), `sigma` (RMSE), `objfunc.fold` (per-fold MSE), `folds`,
`QR`, `paired`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md),
[`nmf.rrr.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md)

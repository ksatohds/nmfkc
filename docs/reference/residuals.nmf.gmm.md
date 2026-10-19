# Residuals from an NMF-GMM fit

\\Y - \hat Y\\, matching
[`residuals.nmf`](https://ksatohds.github.io/nmfkc/reference/residuals.nmf.md).
Without this method
[`residuals()`](https://rdrr.io/r/stats/residuals.html) fell through to
[`stats::residuals.default`](https://rdrr.io/r/stats/residuals.html) and
returned `NULL`, even though
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) worked.

## Usage

``` r
# S3 method for class 'nmf.gmm'
residuals(object, Y, ...)
```

## Arguments

- object:

  An object of class `"nmf.gmm"`.

- Y:

  The observation matrix the model was fitted on.

- ...:

  Ignored.

## Value

The residual matrix.

## Development status

Part of the **experimental** NMF-GMM family: see
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md). The
interface may change in future versions — argument names, defaults and
the contents of the returned object are not yet stable.

## See also

[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md),
[`fitted.nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.gmm.md)

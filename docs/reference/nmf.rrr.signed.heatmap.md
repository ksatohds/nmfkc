# Heatmap visualization of nmfae.signed factor matrices

Displays the factor blocks of a
[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)
fit as side-by-side heatmaps. Non-negative blocks (\\X_1, C\_{+},
C\_{-}, X_2\\) use the white-orange-red palette; the signed combined
bottleneck \\C = C\_{+} - C\_{-}\\ is rendered with a diverging
blue-white-red palette so positive and negative weights are visually
distinguishable.

## Usage

``` r
nmf.rrr.signed.heatmap(
  x,
  Y1.label = NULL,
  X1.label = NULL,
  X2.label = NULL,
  Y2.label = NULL,
  palette.pos = NULL,
  palette.signed = NULL,
  show.C = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `"nmfae.signed"`.

- Y1.label:

  Character vector for rows of \\X_1\\.

- X1.label:

  Decoder basis labels.

- X2.label:

  Encoder basis labels.

- Y2.label:

  Input variable labels.

- palette.pos:

  Palette for non-negative blocks. Default white-orange-red.

- palette.signed:

  Palette for signed \\C\\. Default blue-white-red.

- show.C:

  Logical. If `TRUE` (default), shows the combined signed \\C = C\_{+} -
  C\_{-}\\ as a separate panel.

- ...:

  Not used.

## Value

Invisible `NULL`. Called for its side effect (plot).

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md),
[`nmf.rrr.heatmap`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.heatmap.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y1 <- matrix(abs(rnorm(12)), 3, 4)
Y2 <- matrix(abs(rnorm(20)), 5, 4)
res <- nmf.rrr.signed(Y1, Y2, rank1 = 2, rank2 = 2, maxit = 200)
#> Warning: maximum iterations (200) reached...
nmf.rrr.signed.heatmap(res)

# }
```

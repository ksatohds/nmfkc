# Plot method for nmfae.ecv objects

Visualizes element-wise cross-validation results. When `rank.encoder`
was `NULL` (paired), a line plot of sigma vs rank is drawn. When
`rank.encoder` was explicitly specified (grid), a heatmap of sigma over
the (rank, rank.encoder) grid is drawn.

## Usage

``` r
# S3 method for class 'nmfae.ecv'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"nmfae.ecv"` returned by
  [`nmf.rrr.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md).

- ...:

  Additional graphical parameters (currently unused).

## Value

Invisible `NULL`. Called for its side effect of producing a plot.

## See also

[`nmf.rrr.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md)

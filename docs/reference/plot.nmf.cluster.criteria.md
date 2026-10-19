# Plot clustering-quality criteria across a sequence of fits

Line plot of `silhouette`, `CPCC`, and `dist.cor` against the result
index, for an object from
[`nmf.cluster.criteria`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.criteria.md).
X-axis ticks default to each result's `$rank` (overridable via the
`names` argument of
[`nmf.cluster.criteria`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.criteria.md)).

## Usage

``` r
# S3 method for class 'nmf.cluster.criteria'
plot(
  x,
  main = "Clustering quality across rank",
  xlab = "rank (Q)",
  ylab = "criterion",
  lwd = 2,
  ...
)
```

## Arguments

- x:

  An object of class `"nmf.cluster.criteria"`.

- main:

  Plot title.

- xlab, ylab:

  Axis labels.

- lwd:

  Line width.

- ...:

  Further arguments passed to the initial
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## See also

[`nmf.cluster.criteria`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.criteria.md)

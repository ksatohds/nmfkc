# Plot a cluster-flow (alluvial) diagram

Draws the alluvial / Sankey-style cluster-flow diagram for an object
created by
[`nmf.cluster.flow`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.flow.md).

## Usage

``` r
# S3 method for class 'nmf.cluster.flow'
plot(
  x,
  col = NULL,
  lwd = 1,
  xlab = "rank (Q)",
  ylab = "individuals",
  main = "Cluster flow across rank",
  ...
)
```

## Arguments

- x:

  An object of class `"nmf.cluster.flow"`.

- col:

  Optional colour vector indexed by **reference** cluster id (`col[k]`
  colours every individual whose reference cluster is `k`); recycled if
  shorter than the number of reference clusters. Defaults to the
  object's palette.

- lwd:

  Line width of the flow segments.

- xlab, ylab:

  Axis labels.

- main:

  Plot title.

- ...:

  Further arguments passed to the initial
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) call.

## Value

`x`, invisibly.

## See also

[`nmf.cluster.flow`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.flow.md)

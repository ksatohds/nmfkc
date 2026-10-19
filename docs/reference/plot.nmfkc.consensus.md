# Plot a consensus rank-selection (nmfkc.consensus) object

Two views of a
[`nmfkc.consensus`](https://ksatohds.github.io/nmfkc/reference/nmfkc.consensus.md)
result:

- `type = "criteria"` (default): the stability curves `cophenetic`
  (blue) and `dispersion` (red) against rank. No "best" marker is drawn
  – a stable rank tends to be a coarse one, so the curves are shown for
  inspection rather than as an optimum.

- `type = "heatmap"`: the consensus matrix of each requested `rank`,
  reordered by average-linkage hierarchical clustering of \\1 - \bar C\\
  (blue = 0 / different cluster, red = 1 / same cluster). Requires
  `keep.consensus = TRUE` at compute time.

## Usage

``` r
# S3 method for class 'nmfkc.consensus'
plot(
  x,
  type = c("criteria", "heatmap"),
  rank = NULL,
  mfrow = NULL,
  mar = NULL,
  col = grDevices::hcl.colors(50, "Blue-Red"),
  main = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `"nmfkc.consensus"`.

- type:

  `"criteria"` or `"heatmap"`.

- rank:

  For `type = "heatmap"`, the rank(s) to display. `NULL` (default) shows
  every rank stored in `x`.

- mfrow:

  Panel layout for multiple heatmaps, as `c(nrow, ncol)`. `NULL`
  (default) uses [`n2mfrow`](https://rdrr.io/r/grDevices/n2mfrow.html)
  for a near-square grid.

- mar:

  Per-panel margins `c(b, l, t, r)` for the heatmap view. `NULL`
  (default) uses tight margins for a grid and wider margins (for sample
  labels) when a single rank is shown.

- col:

  Heatmap colour palette (length-50 blue-to-red by default).

- main:

  Plot title for the `"criteria"` view (heatmap panels are titled per
  rank).

- ...:

  Further arguments passed to the underlying
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) /
  [`image`](https://rdrr.io/r/graphics/image.html).

## Value

`x`, invisibly.

## See also

[`nmfkc.consensus`](https://ksatohds.github.io/nmfkc/reference/nmfkc.consensus.md)

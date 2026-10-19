# Plot a rank-selection (nmf.rank) object

Draws the concise three-criterion rank-selection plot for an object
returned by
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)
(or `nmfkc.net.rank`, `nmfkc.signed.rank`, `nmfae.rank`,
`nmfae.signed.rank`): `r.squared` (red) and the broken-stick-corrected
effective-rank index `effective.rank.index` (green, legend label
`eff.rank`; **not** the raw `effective.rank.ratio`) on the left \\\[0,
1\]\\ axis, and the cross-validation error `sigma.ecv` (blue) on the
right axis, each with points, rank-number labels and a highlighted best
marker.

## Usage

``` r
# S3 method for class 'nmf.rank'
plot(
  x,
  main = NULL,
  xlab = "Rank (Q)",
  ylab = "R-squared / eff.rank (0-1)",
  lwd = 3,
  ...
)
```

## Arguments

- x:

  An object of class `"nmf.rank"`.

- main:

  Plot title (defaults to the title stored in `x`).

- xlab, ylab:

  Axis labels.

- lwd:

  Line width of the criterion curves.

- ...:

  Further arguments passed to the initial
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## See also

[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)

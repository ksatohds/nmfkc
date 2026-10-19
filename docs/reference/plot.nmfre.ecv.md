# Plot method for nmfre.ecv objects

Draws the element-wise CV score (held-out RMSE) against the rank,
marking the minimizing rank as “Best (Min)”.

## Usage

``` r
# S3 method for class 'nmfre.ecv'
plot(
  x,
  main = "NMF-RE element-wise CV",
  xlab = "Rank (Q)",
  ylab = "sigma.ecv (held-out RMSE)",
  ...
)
```

## Arguments

- x:

  An object of class `"nmfre.ecv"`.

- main:

  Plot title.

- xlab, ylab:

  Axis labels.

- ...:

  Passed to [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The input object, invisibly.

## See also

[`nmfre.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfre.ecv.md)

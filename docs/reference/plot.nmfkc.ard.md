# Plot method for nmfkc.ard objects

Bar plot of the per-component relevance (descending), with a line at the
pruning threshold; bars above it (the estimated rank) are highlighted.

## Usage

``` r
# S3 method for class 'nmfkc.ard'
plot(x, main = NULL, ...)
```

## Arguments

- x:

  An object of class `"nmfkc.ard"`.

- main:

  Plot title.

- ...:

  Passed to [`barplot`](https://rdrr.io/r/graphics/barplot.html).

## Value

`x`, invisibly.

# `plot.nmfae` displays the convergence trajectory of the objective function across iterations. The title shows the achieved \\R^2\\.

`plot.nmfae` displays the convergence trajectory of the objective
function across iterations. The title shows the achieved \\R^2\\.

## Usage

``` r
# S3 method for class 'nmfae'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"nmfae"` returned by
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md).

- ...:

  Additional graphical parameters passed to `plot`.

## Value

Invisible `NULL`. Called for its side effect (plot).

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`nmf.rrr.heatmap`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.heatmap.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(runif(20), nrow = 4)
res <- nmf.rrr(Y, rank1 = 2)
plot(res)

# }
```

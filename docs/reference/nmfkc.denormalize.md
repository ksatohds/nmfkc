# Denormalize a matrix from \\\[0,1\]\\ back to its original scale

`nmfkc.denormalize` rescales a matrix with values in \\\[0,1\]\\ back to
its original scale using the column-wise minima and maxima of a
reference matrix.

## Usage

``` r
nmfkc.denormalize(x, ref = x)
```

## Arguments

- x:

  A numeric matrix (or vector) with values in \\\[0,1\]\\ to be
  denormalized.

- ref:

  A reference matrix used to obtain the original column-wise minima and
  maxima. Must have the same number of columns as `x`.

## Value

A numeric matrix with values transformed back to the original scale.

## See also

[`nmfkc.normalize`](https://ksatohds.github.io/nmfkc/reference/nmfkc.normalize.md)

## Examples

``` r
x <- nmfkc.normalize(iris[, -5])
x_recovered <- nmfkc.denormalize(x, iris[, -5])
apply(x_recovered - iris[, -5], 2, max)
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#> 0.000000e+00 4.440892e-16 0.000000e+00 2.220446e-16 
```

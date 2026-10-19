# Normalize a matrix to the range \\\[0,1\]\\

`nmfkc.normalize` rescales the values of a matrix to lie between 0 and 1
using the column-wise minimum and maximum values of a reference matrix.

## Usage

``` r
nmfkc.normalize(x, ref = x)
```

## Arguments

- x:

  A numeric matrix (or vector) to be normalized.

- ref:

  A reference matrix from which the column-wise minima and maxima are
  taken. Default is `x`.

## Value

A matrix of the same dimensions as `x`, with each column rescaled to the
\\\[0,1\]\\ range.

## See also

[`nmfkc.denormalize`](https://ksatohds.github.io/nmfkc/reference/nmfkc.denormalize.md)

## Examples

``` r
# Example.
x <- nmfkc.normalize(iris[,-5])
apply(x,2,range)
#>      Sepal.Length Sepal.Width Petal.Length Petal.Width
#> [1,]            0           0            0           0
#> [2,]            1           1            1           1
```

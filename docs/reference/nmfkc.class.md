# Create a class (one-hot) matrix from a categorical vector

`nmfkc.class` converts a categorical or factor vector into a class
matrix (one-hot encoded representation), where each row corresponds to a
category and each column corresponds to an observation.

## Usage

``` r
nmfkc.class(x)
```

## Arguments

- x:

  A categorical vector or a factor.

## Value

A binary matrix with one row per unique category and one column per
observation. Each column has exactly one entry equal to 1, indicating
the category of the observation.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)

## Examples

``` r
# Example.
Y <- nmfkc.class(iris$Species)
Y[,1:6]
#>            1 2 3 4 5 6
#> setosa     1 1 1 1 1 1
#> versicolor 0 0 0 0 0 0
#> virginica  0 0 0 0 0 0
```

# Plot method for objects of class `nmfkc`

`plot.nmfkc` produces a diagnostic plot for the return value of `nmfkc`,
showing the objective function across iterations.

## Usage

``` r
# S3 method for class 'nmfkc'
plot(x, ...)
```

## Arguments

- x:

  An object of class `nmfkc`, i.e., the return value of `nmfkc`.

- ...:

  Additional arguments passed to the base
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function.

## Value

Called for its side effect (a plot). Returns `NULL` invisibly.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`summary.nmfkc`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
plot(result)

```

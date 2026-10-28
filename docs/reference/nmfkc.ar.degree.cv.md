# Optimize lag order for the autoregressive model

`nmfkc.ar.degree.cv` selects the optimal lag order for an autoregressive
model by applying cross-validation over candidate degrees.

This function accepts both standard matrices (Variables x Time) and `ts`
objects (Time x Variables). `ts` objects are automatically transposed
internally.

## Usage

``` r
nmfkc.ar.degree.cv(
  Y,
  rank = 1,
  degree = 1:2,
  intercept = TRUE,
  plot = TRUE,
  ...
)
```

## Arguments

- Y:

  Observation matrix \\Y(P,N)\\ or a `ts` object.

- rank:

  Rank of the basis matrix. For backward compatibility, `Q` is accepted
  via `...`.

- degree:

  A vector of candidate lag orders to be evaluated.

- intercept:

  Logical. If TRUE (default), an intercept is added to the covariate
  matrix.

- plot:

  Logical. If TRUE (default), a plot of the objective function values is
  drawn.

- ...:

  Additional arguments passed to `nmfkc.cv`. Also accepts `cores`
  (`getOption("mc.cores", 1L)`) to evaluate the candidate degrees in
  parallel; results are identical to the sequential loop for any `cores`
  (PSOCK cluster on Windows, forking elsewhere).

## Value

A list with components:

- degree:

  The lag order that minimizes the cross-validation objective function.

- degree.max:

  Maximum recommended lag order, computed as \\10 \log\_{10}(N)\\
  following the `ar` function in the stats package.

- objfunc:

  Objective function values for each candidate lag order.

## See also

[`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md),
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md)

## Examples

``` r
# Example using ts object directly
d <- AirPassengers

# Selection of degree (using ts object)
# Note: Y is automatically transposed if it is a ts object
nmfkc.ar.degree.cv(Y=d, rank=1, degree=11:14)
#> degree=11...
#> 0sec
#> degree=12...
#> 0.1sec
#> degree=13...
#> 0.1sec
#> degree=14...
#> 0.1sec

#> $degree
#> [1] 12
#> 
#> $degree.max
#> [1] 21
#> 
#> $objfunc
#>       11       12       13       14 
#> 507.7053 271.4061 276.0141 278.1065 
#> 
```

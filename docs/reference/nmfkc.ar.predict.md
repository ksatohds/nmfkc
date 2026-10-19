# Forecast future values for NMF-VAR model

`nmfkc.ar.predict` computes multi-step-ahead forecasts for a fitted
NMF-VAR model using recursive forecasting.

If the fitted model contains time series property information (from
`nmfkc.ar`), the forecasted values will have appropriate time-based
column names.

## Usage

``` r
nmfkc.ar.predict(x, Y, degree = NULL, n.ahead = 1)
```

## Arguments

- x:

  An object of class `nmfkc` (the fitted model).

- Y:

  The historical observation matrix used for fitting (or at least the
  last `degree` columns).

- degree:

  Optional integer. Lag order (D). If `NULL` (default), it is inferred
  from `x$A.attr` (when available) or from the dimensions of `x$C`.

- n.ahead:

  Integer (\>=1). Number of steps ahead to forecast.

## Value

A list with components:

- pred:

  A \\P \times n.ahead\\ matrix of predicted values. Column names are
  future time points if time information is available.

- time:

  A numeric vector of future time points corresponding to the columns of
  `pred`.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md)

## Examples

``` r
# Forecast AirPassengers
d <- AirPassengers
ar_data <- nmfkc.ar(d, degree = 2)
result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#> Y(1,142)~X(1,1)C(1,3)A(3,142)=XB(1,142)...
#> 0sec
pred <- nmfkc.ar.predict(result, Y = matrix(d, nrow = 1), degree = 2, n.ahead = 3)
pred$pred
#>          1961 1961.08333 1961.16667
#> [1,] 427.6732   424.8679   422.1085
```

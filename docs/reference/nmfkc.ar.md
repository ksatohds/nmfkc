# Construct observation and covariate matrices for a vector autoregressive model

`nmfkc.ar` generates the observation matrix and covariate matrix
corresponding to a specified autoregressive lag order.

If the input `Y` is a `ts` object, its time properties are preserved in
the `"tsp_info"` attribute, adjusted for the lag. Additionally, the
column names of `Y` and `A` are set to the corresponding time points.

## Usage

``` r
nmfkc.ar(Y, degree = 1, intercept = TRUE)
```

## Arguments

- Y:

  An observation matrix (P x N) or a `ts` object. If `Y` is a `ts`
  object (typically N x P), it is automatically transposed to match the
  (P x N) format.

- degree:

  The lag order of the autoregressive model. The default is 1.

- intercept:

  Logical. If TRUE (default), an intercept term is added to the
  covariate matrix.

## Value

A list containing:

- Y:

  Observation matrix (P x N_A) used for NMF. Includes adjusted
  `"tsp_info"` attribute and time-based column names.

- A:

  Covariate matrix (R x N_A) constructed according to the specified lag
  order. Includes adjusted `"tsp_info"` attribute and time-based column
  names.

- A.columns:

  Index matrix used to generate `A`.

- degree.max:

  Maximum lag order.

## References

Satoh, K. (2025). Applying non-negative matrix factorization with
covariates to multivariate time series data as a vector autoregression
model. *Japanese Journal of Statistics and Data Science*.
arXiv:2501.17446.
[doi:10.1007/s42081-025-00314-0](https://doi.org/10.1007/s42081-025-00314-0)

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.ar.degree.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.degree.cv.md),
[`nmfkc.ar.stationarity`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.stationarity.md),
[`nmfkc.ar.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.DOT.md)

## Examples

``` r
# Example using AirPassengers (ts object)
d <- AirPassengers
ar_data <- nmfkc.ar(d, degree = 2)
dim(ar_data$Y)
#> [1]   1 142
dim(ar_data$A)
#> [1]   3 142

# Example using matrix input
Y <- matrix(1:20, nrow = 2)
ar_data <- nmfkc.ar(Y, degree = 1)
ar_data$degree.max
#> [1] 9
```

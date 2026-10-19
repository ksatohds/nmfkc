# Stationarity of an NMF-VAR model

Reports whether a fitted NMF-VAR is stationary, and how cheaply that
could be decided. Three numbers are returned, in increasing order of
sharpness: `colsum.max` \\=\max_j c_j\\, the largest column sum of
\\\sum_d\Xi_d = X\sum_d\Theta_d\\, certifies stationarity when below 1
with no eigenvalue computation at all (a sufficient condition only, so a
value above 1 is merely inconclusive); `spectral.radius.sum`
\\=\rho(\sum_d G_d)\\ is below 1 exactly when the model is stationary
(and is bracketed by \\\min_j c_j\\ and \\\max_j c_j\\); and
`spectral.radius` is the spectral radius of the companion matrix itself.

Note that the \\c_j\\ are weighted by \\\bm 1'X\\, i.e.
`colSums(X) %*% Lambda`. Using the raw column sums of \\\sum_d\Theta_d\\
is correct only when \\\bm 1'X=\bm 1'\\ and can certify stationarity for
a non-stationary fit otherwise.

This function answers the yes/no question only. For the dynamics behind
the number — the transition matrices, the roots, the impulse responses,
the long-run mean — use
[`nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md),
and for the standard error of \\\rho\\ use
[`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md).

## Usage

``` r
nmfkc.ar.stationarity(x, method = c("latent", "companion"))
```

## Arguments

- x:

  A fitted `"nmfkc"` object obtained from an
  [`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md)
  design, or an object from
  [`nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md).

- method:

  Route used for `spectral.radius`. `"latent"` (default) uses the
  \\DQ\times DQ\\ companion matrix of the \\G_d\\; `"companion"` uses
  the \\PD\times PD\\ companion matrix of the \\\Xi_d\\. The two
  companion matrices share the same non-zero eigenvalues, so the results
  agree up to floating-point rounding, but `"latent"` costs
  \\O((DQ)^3)\\ instead of \\O((PD)^3)\\ — for a 47-variable, 7-lag,
  rank-4 model that is a \\28\times28\\ eigenproblem instead of a
  \\329\times329\\ one. Use `"companion"` to reproduce the value
  bit-for-bit as computed by earlier versions of the package.

## Value

An object of class `"nmfkc.ar.stationarity"`: a list with

- spectral.radius:

  Spectral radius of the companion matrix (route selected by `method`).
  A value below 1 indicates stationarity.

- stationary:

  Logical; `TRUE` when `spectral.radius < 1`.

- method:

  The route actually used.

- spectral.radius.sum:

  \\\rho(\sum_d G_d)\\; below 1 iff stationary.

- colsum:

  Column sums \\c_j\\ of \\\sum_d\Xi_d\\, i.e. `colSums(X) %*% Lambda`
  (length \\P\\).

- colsum.max:

  \\\max_j c_j\\; below 1 it certifies stationarity.

- dims:

  Named vector of `P`, `Q`, `D`.

## References

Satoh, K. (2025). Applying non-negative matrix factorization with
covariates to multivariate time series data as a vector autoregression
model. *Japanese Journal of Statistics and Data Science*.
arXiv:2501.17446.
[doi:10.1007/s42081-025-00314-0](https://doi.org/10.1007/s42081-025-00314-0)

## See also

[`nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md),
[`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md),
[`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md)

## Examples

``` r
# Check stationarity of fitted AR model
d <- AirPassengers
ar_data <- nmfkc.ar(d, degree = 2)
result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#> Y(1,142)~X(1,1)C(1,3)A(3,142)=XB(1,142)...
#> 0sec
st <- nmfkc.ar.stationarity(result)
st$spectral.radius
#> [1] 0.968468
st$stationary
#> [1] TRUE
st                 # full report
#> rho(companion, 2x2, method="latent") = 0.9685  ->  stationary
#> rho(sum_d G_d)        = 0.9675   (< 1 iff stationary)
#> bracket: min_j c_j = 0.9675  <=  rho(sum_d G_d)  <=  max_j c_j = 0.9675
#>    max_j c_j < 1: stationarity certified without an eigenvalue
#> 
#> For the dynamics behind this number see nmfkc.ar.latent(), for a test of H0: rho >= 1 see nmfkc.ar.latent.inference().
```

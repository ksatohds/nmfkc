# Plot the stationarity bracket of an NMF-VAR model

Draws the Perron-Frobenius bracket \\\min_j c_j\le\rho(\sum_d
G_d)\le\max_j c_j\\ against the unit boundary, where \\c_j\\ are the
column sums of \\\sum_d\Theta_d\\. It shows at a glance whether the
cheap sufficient condition settles stationarity or an eigenvalue
computation is needed. For the model's dynamics see
[`plot.nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.ar.latent.md).

## Usage

``` r
# S3 method for class 'nmfkc.ar.stationarity'
plot(x, ...)
```

## Arguments

- x:

  An object from
  [`nmfkc.ar.stationarity`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.stationarity.md).

- ...:

  Passed to the underlying plot call.

## Value

`invisible(NULL)`; called for the plot.

## See also

[`nmfkc.ar.stationarity`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.stationarity.md),
[`plot.nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.ar.latent.md)

## Examples

``` r
set.seed(1)
Y   <- matrix(abs(rnorm(4 * 60)) + 1, 4, 60)
ar  <- nmfkc.ar(Y, degree = 2)
fit <- nmfkc(ar$Y, ar$A, rank = 2, verbose = FALSE)
plot(nmfkc.ar.stationarity(fit))
```

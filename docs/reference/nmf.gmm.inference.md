# Statistical inference for an NMF-GMM fit (given basis)

This function is **experimental and still under development**. The
interface may change in future versions: argument names, defaults and
the contents of the returned object are not yet stable. Code written
against it today may need adjusting after an update. The rest of the
package does not carry this caveat.

`nmf.gmm.inference` adds Wald inference on the covariate-coefficient
matrix \\C\\ (\\=\Theta\\) of a fitted
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md)
object, conditional on the estimated basis \\\hat X\\ and mixture. It
reports the analytic outer-product (mixture-information) standard error
— the coverage-calibrated choice at \\K\>1\\ — together with a
**wild-bootstrap** standard error and percentile confidence interval,
mirroring
[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md)
/
[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md).
Supported for the tied-covariance variant.

## Usage

``` r
nmf.gmm.inference(object, Y, A = object$A, ...)
```

## Arguments

- object:

  A fitted `"nmf.gmm"` object (with `cov = "tied"`).

- Y:

  The data matrix used in the fit.

- A:

  The covariate matrix used in the fit. Default `object$A`.

- ...:

  Additional arguments, named as in the other wild-bootstrap inference
  functions of the package: `wild.B` (bootstrap replicates, default
  500), `wild.level` (confidence level, default 0.95) and `seed`
  (default 123).

## Value

`object` with class `c("nmf.gmm.inference", "nmf.gmm")` and added
components: `coefficients` (data frame with `Basis`, `Covariate`,
`Estimate`, `SE`, `BSE`, `z_value`, `p_value`, `CI_low`, `CI_high`),
`C.se.outer`, `C.se.sandwich`, `C.se.boot` (Q x R matrices) and `Xi`
(\\= X C\\).

## See also

[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md),
[`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(abs(rnorm(20 * 60)) + 1, 20, 60); A <- rbind(1, rnorm(60))
fit <- nmf.gmm(Y, A, rank = 2, K = 2)
fit <- nmf.gmm.inference(fit, Y, A, wild.B = 300)
head(fit$coefficients)
#>    Basis Covariate    Estimate        SE       BSE     z_value   p_value
#> 1 Basis1      Cov1 36.27851926 0.7824778 0.3225164 46.36364255 0.0000000
#> 2 Basis2      Cov1  0.14127734 0.5034434 0.3956258  0.28062206 0.7790003
#> 3 Basis1      Cov2 -0.06061948 0.7665946 0.2878565 -0.07907633 0.9369719
#> 4 Basis2      Cov2  0.27859573 0.5495305 0.2928593  0.50697046 0.6121755
#>       CI_low    CI_high
#> 1 35.6346689 36.9039131
#> 2 -0.6222280  0.8555773
#> 3 -0.5826109  0.4764388
#> 4 -0.3218638  0.8484674
# }
```

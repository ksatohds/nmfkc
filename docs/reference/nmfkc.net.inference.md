# Statistical Inference for Symmetric NMF Parameters

Performs statistical inference on the parameter matrix \\C\\ of a
symmetric NMF model (\\Y \approx XCX^\top\\).

This is a wrapper around
[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
that automatically sets the covariate matrix to \\A = X^\top\\, which is
the defining property of symmetric NMF.

## Usage

``` r
nmfkc.net.inference(object, Y, wild.bootstrap = TRUE, ...)
```

## Arguments

- object:

  A fitted model returned by
  [`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) with
  `Y.symmetric = "bi"` or `"tri"`, or the newer
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  with `type = "bi"` or `"tri"`.

- Y:

  The original symmetric matrix (P x P).

- wild.bootstrap:

  Logical; if TRUE (default), perform wild bootstrap inference for the C
  matrix.

- ...:

  Additional arguments passed to
  [`nmfkc.inference()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
  (e.g., `wild.B`, `wild.seed`, `C.p.side`).

## Value

The object augmented with inference fields (same as `nmfkc.inference`),
with class
`c("nmfkc.net.inference", "nmfkc.inference", "nmf.inference", "nmfkc", "nmf")`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md),
[`nmfkc.net.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.DOT.md),
[`summary.nmfkc.net.inference`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.net.inference.md)

## Examples

``` r
# \donttest{
library(nmfkc)
Y <- matrix(c(0,1,1,0,0,0,
              1,0,1,0,0,0,
              1,1,0,1,0,0,
              0,0,1,0,1,1,
              0,0,0,1,0,1,
              0,0,0,1,1,0), 6, 6)
res <- nmfkc.net(Y, rank = 2, type = "tri", nstart = 20)
res_inf <- nmfkc.net.inference(res, Y)
summary(res_inf)
#> 
#> Call:
#> nmfkc.net(Y = Y, rank = 2, type = "tri", nstart = 20)
#> 
#> Dimensions: 6 
#> Rank (Q):    2 
#> Runtime:     0.0sec 
#> Iterations: 25 / 5000 (converged)  epsilon = 0.0001
#> 
#> Statistics:
#>   Objective function:      5.193
#>   R-squared (cor^2):       0.3957
#>   R-squared (uncentered):  0.6291
#>   R-squared (centered):    0.3768
#>   Residual Std Error:      0.3798
#>   Mean Absolute Error:     0.3031
#> 
#> Structure Diagnostics:
#>   Basis (X) Sparsity:      33.3% (< 1e-4)
#>   Separability:            1 (1 = every basis has a row of its own)
#> 
#> Symmetric NMF type: tri
#> Inference (conditional on X):
#>   sigma^2:   0.1623 
#>   C matrix:  2 x 2 (max asymmetry: 2.1e-17)
#> 
#> C coefficients: 1 total, 0 significant
#>       Row:Col  Estimate Std. Error (Boot) z value   Pr(>z) 
#> Basis1:Basis2     0.053      0.222  0.104    0.24   0.4062  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
# }
```

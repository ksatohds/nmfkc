# Statistical Inference for NMF-AE Parameter Matrix

Performs post-estimation inference for \\\Theta\\ in the three-layer NMF
model \\Y_1 \approx X_1 \Theta X_2 Y_2\\, conditional on \\(\hat{X}\_1,
\hat{X}\_2)\\. Uses sandwich covariance estimation and one-step wild
bootstrap with non-negative projection.

## Usage

``` r
nmf.rrr.inference(object, Y1, Y2 = Y1, wild.bootstrap = TRUE, ...)
```

## Arguments

- object:

  An object of class `"nmfae"` returned by
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md).

- Y1:

  Output matrix \\Y_1\\ (P1 x N). Must match the data used in
  [`nmf.rrr()`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md).

- Y2:

  Input matrix \\Y_2\\ (P2 x N). Default is `Y1` (autoencoder).

- wild.bootstrap:

  Logical. If `TRUE` (default), performs wild bootstrap for bootstrap SE
  and confidence intervals. If `FALSE`, only sandwich SE and z-test
  p-values are computed (faster).

- ...:

  Additional arguments:

  `wild.B`

  :   Number of bootstrap replicates. Default is 500.

  `wild.seed`

  :   Seed for bootstrap. Default is 123.

  `wild.level`

  :   Confidence level for bootstrap CI. Default is 0.95.

  `sandwich`

  :   Logical. Use sandwich covariance. Default is `TRUE`.

  `C.p.side`

  :   P-value type: `"one.sided"` (default) or `"two.sided"`.

  `cov.ridge`

  :   Ridge stabilization for information matrix inversion. Default is
      1e-8.

  `print.trace`

  :   Logical. If `TRUE`, prints progress. Default is `FALSE`.

## Value

An object of class `c("nmfae.inference", "nmfae")`, inheriting all
components from the input `object`, with additional inference
components:

- sigma2.used:

  Estimated \\\sigma^2\\ used for inference.

- C.se:

  Sandwich standard errors for \\\Theta\\ (Q x R matrix).

- C.se.boot:

  Bootstrap standard errors for \\\Theta\\ (Q x R matrix).

- C.ci.lower:

  Lower CI bounds for \\\Theta\\ (Q x R matrix).

- C.ci.upper:

  Upper CI bounds for \\\Theta\\ (Q x R matrix).

- coefficients:

  Data frame with Estimate, SE, BSE, z, p-value for each element of
  \\\Theta\\.

- C.p.side:

  P-value type used.

- identified.approx:

  Whether the factorization is unique up to permutations of both latent
  spaces, i.e. under \\X_1\to X_1T_1\\, \\\Theta\to T_1^{-1}\Theta
  T_2^{-1}\\, \\X_2\to T_2X_2\\. A sufficient condition: the diagnostic
  of
  [`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
  passes on \\(X_1, \Theta X_2)\\ (pinning \\T_1\\) and on \\(X_1\Theta,
  X_2)\\ (pinning \\T_2\\), and \\Y_2\\ has full row rank. Zeros are
  judged relative only (`pin.tol`, `sep.tol` through `...`). An
  algebraic check on the point estimate: no re-fits.

- identification:

  The two checks, `left` and `right`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`summary.nmfae.inference`](https://ksatohds.github.io/nmfkc/reference/summary.nmfae.inference.md)

## Examples

``` r
Y <- matrix(c(1,0,1,0, 0,1,0,1, 1,1,0,0), nrow=3, byrow=TRUE)
res <- nmf.rrr(Y, rank1=2, rank2=2)
res2 <- nmf.rrr.inference(res, Y)
summary(res2)
#> 
#> Call:
#> nmf.rrr(Y1 = Y, rank1 = 2, rank2 = 2)
#> 
#> Dimensions:
#>   Y1:              3 x 4 
#>   Y2:              3 x 4   (autoencoder) 
#>   Response rank Q: 2 
#>   Covariate rank R: 2 
#>   Parameters:      X1(3x2) + C(2x2) + X2(2x3) = 16 
#> 
#> Convergence:
#>   Iterations:      23 / 5000 (converged)  epsilon = 0.0001
#>   Runtime:         0.0 secs 
#> 
#> Goodness of fit:
#>   Objective function:      0.7558
#>   R-squared (cor^2):       0.7843
#>   R-squared (uncentered):  0.874
#>   R-squared (centered):    0.7481
#>   Residual Std Error:      0.251
#>   Mean Absolute Error:     0.2031
#>   Effective Rank:          2.00 / 2  (100.0%)
#> 
#> Structure Diagnostics:
#>   Response (X1) Sparsity:  33.3% (< 1e-4)
#>   Bottleneck (C) Sparsity: 0.0% (< 1e-4)
#>   Covariate (X2) Sparsity: 33.3% (< 1e-4)
#> 
#> Coefficients (conditional on X1, X2): 4 total, 2 significant
#>   Cov:Resp  Estimate Std. Error (Boot) z value   Pr(>z) 
#> Cov1:Resp1     1.657      0.102  0.083   16.28   <2e-16 ***
#> Cov1:Resp2     0.000      0.189  0.067    0.00   0.4993  
#> Cov2:Resp1     0.000      0.189  0.066    0.00   0.4993  
#> Cov2:Resp2     1.657      0.102  0.090   16.28   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#>   Identification: unique up to a permutation of the bases; the zeros it rests on hold to 0.0002, so an equivalent fit mixes the bases by at most about 0.02% (both latent spaces)
```

# Statistical inference for the coefficient matrix C from NMF-RE

`nmfre.inference` performs statistical inference on the coefficient
matrix \\C\\ (\\\Theta\\) from a fitted `nmfre` model, conditional on
the estimated basis matrix \\\hat{X}\\ and random effects \\\hat{U}\\.

Under the working model \\Y^\* = Y - X\hat{U} \approx X C A +
\varepsilon\\, inference is conducted via sandwich covariance estimation
and one-step wild bootstrap with non-negative projection.

The result is compatible with
[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md)
for visualization (pass the result directly as `x` with `type = "YXA"`).

## Usage

``` r
nmfre.inference(object, Y, A = NULL, wild.bootstrap = TRUE, ...)
```

## Arguments

- object:

  An object of class `"nmfre"` returned by
  [`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md).

- Y:

  Observation matrix (P x N). Must match the data used in
  [`nmfre()`](https://ksatohds.github.io/nmfkc/reference/nmfre.md).

- A:

  Covariate matrix (K x N). Default is `NULL` (intercept only).

- wild.bootstrap:

  Logical. If `TRUE` (default), performs wild bootstrap for confidence
  intervals and bootstrap standard errors.

- ...:

  Additional arguments:

  `wild.B`

  :   Number of bootstrap replicates. Default is 500.

  `wild.seed`

  :   Seed for bootstrap. Default is 123.

  `wild.level`

  :   Confidence level for bootstrap CI. Default is 0.95.

  `C.p.side`

  :   P-value type: `"one.sided"` (default) or `"two.sided"`.

  `cov.ridge`

  :   Ridge stabilization. Default is 1e-8.

  `print.trace`

  :   Logical. Default is `FALSE`.

## Value

The input `object` with additional inference components:

- sigma2.used:

  Estimated \\\sigma^2\\ used for inference.

- C.vec.cov:

  Full covariance matrix for \\vec(C)\\.

- C.se:

  Sandwich standard errors for \\C\\.

- C.se.boot:

  Bootstrap standard errors for \\C\\.

- C.ci.lower:

  Lower CI bounds for \\C\\.

- C.ci.upper:

  Upper CI bounds for \\C\\.

- coefficients:

  Data frame with Basis, Covariate, Estimate, SE, BSE, z_value, p_value,
  CI_low, CI_high.

- C.p.side:

  P-value type used.

- identified.approx:

  Whether \\(X, C)\\ is unique up to a permutation of the bases, by the
  data. With `C.signed = FALSE` it is the diagnostic of
  [`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
  on \\(X, C)\\. With a signed \\C\\ (the default) non-negativity cannot
  pin the rotation, so it is `NA`: only the working model's isotropic
  random effects \\U\sim N(0,\tau^2 I)\\ do, by allowing orthogonal
  \\T\\ alone, which with the column sums of \\X\\ fixed is a
  permutation when \\Q=2\\ or every basis has an anchor row. That is
  identification by assumption. An algebraic check on the point
  estimate: no re-fits.

- identification:

  The parts of the diagnostic, plus `by.assumption` (the isotropic
  argument above holds, with \\\tau^2\>0\\) and, for a signed \\C\\,
  `note`.

## References

Satoh, K. (2026). Wild Bootstrap Inference for Non-Negative Matrix
Factorization with Random Effects. arXiv:2603.01468.
<https://arxiv.org/abs/2603.01468>

## See also

[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md),
[`summary.nmfre`](https://ksatohds.github.io/nmfkc/reference/summary.nmfre.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(intercept = 1, speed = cars$speed)
res <- nmfre(Y, A, rank = 1, wild.bootstrap = FALSE)
res2 <- nmfre.inference(res, Y, A)
res2$coefficients
#>    Basis Covariate   Estimate        SE       BSE   z_value     p_value
#> 1 Basis1 intercept -17.579095 5.5981348 5.5578503 -3.140170 0.001688498
#> 2 Basis1     speed   3.932409 0.4027284 0.3937048  9.764418 0.000000000
#>       CI_low   CI_high
#> 1 -28.482843 -6.692752
#> 2   3.171203  4.772687
```

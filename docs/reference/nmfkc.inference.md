# Statistical inference for the parameter matrix C (Theta)

`nmfkc.inference` performs statistical inference on the parameter matrix
\\C\\ (\\\Theta\\) from a fitted `nmfkc` model, conditional on the
estimated basis matrix \\\hat{X}\\.

Under the working model \\Y = X C A + \varepsilon\\ where
\\\varepsilon\_{pn} \stackrel{iid}{\sim} N(0, \sigma^2)\\, inference is
conducted via sandwich covariance estimation and one-step wild bootstrap
with non-negative projection.

## Usage

``` r
nmfkc.inference(object, Y, A = NULL, wild.bootstrap = TRUE, ...)
```

## Arguments

- object:

  An object of class `"nmfkc"` returned by
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md).

- Y:

  Observation matrix (P x N). Must match the data used in
  [`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md).

- A:

  Covariate matrix (K x N). Default is `NULL` (same as identity; in this
  case \\B = C\\ and inference is on \\B\\ directly).

- wild.bootstrap:

  Logical. If `TRUE` (default), performs wild bootstrap for confidence
  intervals and bootstrap standard errors. Set to `FALSE` to skip
  bootstrap (faster, only sandwich SE is computed).

- ...:

  Additional arguments:

  `method`

  :   Bootstrap engine. `"onestep"` (default): a one-step Newton
      linearisation around the fit using the inverse information `Hinv`
      (fast; the sandwich SE is primary). `"refit"`: a residual wild
      bootstrap that re-estimates \\C\\ to convergence with the basis
      \\X\\ held FIXED. The "refit" engine uses no information matrix,
      so it stays valid when the Fisher information \\AA'\\ is singular
      (over-parameterised covariates) or \\C\\ sits on the \\\ge 0\\
      boundary; the bootstrap SE/CI become primary and the p-value is a
      two-sided bootstrap p-value.

  `wild.dist`

  :   Multiplier distribution (orthogonal to `method`): `"rademacher"`
      (\\\pm 1\\), `"mammen"` (Mammen 1993 two-point), or `"exp"`
      (\\\mathrm{Exp}(1)-1\\). Default `"exp"` for `method="onestep"`
      (backward compatible), `"rademacher"` for `"refit"`.

  `wild.unit`

  :   For `method="refit"`: `"element"` (default, i.i.d. multiplier per
      matrix cell), `"column"` (one multiplier per sample column, shared
      over rows) or `"block"` (one multiplier per block of
      `L = max(1, round(N^(1/3)))` consecutive columns, shared over
      rows; `N = ncol(Y)`). Use `"block"` when the columns are time
      points with serially dependent residuals, e.g. the design of
      [`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md):
      residuals at nearby time points then share a multiplier. `L` is 3
      for `N = 27`, 5 for 100, 8 for 500 and 1 (the same draws as
      `"column"`) for `N <= 3`. The block grid starts at a random offset
      in `0, ..., L-1`, drawn anew in each replicate, so the first and
      last blocks may be shorter. `L` is fixed by `N`; it cannot be set
      and is not stored in the result.

  `refit.epsilon`, `refit.maxit`

  :   For `method="refit"`: convergence tolerance and iteration cap of
      the re-fits, defaulting to `1e-8` and `100000`. These are
      deliberately far tighter than
      [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)'s
      own `1e-4`: under multiplicative updates a coefficient heading for
      the non-negativity boundary approaches 0 slowly, so at a loose
      tolerance every replicate stops while that coefficient is still
      spuriously positive. The replicates then pile up above the point
      estimate, the percentile interval can exclude it entirely, and the
      bootstrap SE is understated. Loosen only if the re-fits are too
      slow, and check that the intervals still contain the estimates.

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

  `sep.tol`, `pin.tol`

  :   Tolerances of the identification diagnostic (defaults `0.999` and
      `0.01`); see `identified.approx` below and the Details of
      [`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md).
      Coefficients at or below `boundary.tol` also count as zero there.

## Value

An object of class `c("nmfkc.inference", "nmfkc")`, inheriting all
components from the input `object`, with additional inference
components:

- sigma2.used:

  Estimated \\\sigma^2\\ used for inference.

- C.se:

  Sandwich standard errors for \\C\\ (Q x K matrix).

- C.se.boot:

  Bootstrap standard errors for \\C\\ (Q x K matrix).

- C.ci.lower:

  Lower CI bounds for \\C\\ (Q x K matrix).

- C.ci.upper:

  Upper CI bounds for \\C\\ (Q x K matrix).

- coefficients:

  Data frame with one row per element of \\C\\. `SE` is always the
  sandwich standard error and `BSE` always the bootstrap one, so the two
  columns differ; `z_value` and `p_value` use whichever is primary for
  the chosen `method` (bootstrap under `"refit"`, sandwich under
  `"onestep"`). `on.boundary` flags coefficients at or below
  `boundary.tol` (default `1e-3`): their interval and p-value behave
  correctly in the "do not reject" direction, but the value is not a
  calibrated p-value — a one-sided bootstrap p piles up near \\1/2\\ at
  the boundary.

- C.p.side:

  P-value type used.

- wild.unit:

  Multiplier unit of the `method = "refit"` bootstrap (`"element"`,
  `"column"` or `"block"`); `NA` otherwise.

- identified.approx:

  Logical, an approximate diagnostic of whether the factorization is
  unique up to a permutation of the bases. The intervals are conditional
  on \\\hat X\\; when this is `TRUE` they concern the same \\C\\
  whichever admissible representative \\(XT,\\T^{-1}C)\\ is taken, and
  when `FALSE` the entries of \\C\\ describe this representative. `TRUE`
  needs \\A\\ of full row rank and either \\Q=1\\,
  `X.restriction = "fixed"`, or a fixed per-column scale (`"colSums"` /
  `"colSqSums"`) with one of two zero patterns: an anchor row of \\X\\
  for every basis and the rows of \\C\\ pinned one at a time, or the
  mirror condition. The argument is in the Details of
  [`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md).
  An algebraic check on the point estimate: no re-fits.

- identification:

  Its parts: `anchor`, `pinned`, `pin.order`, `pure`, `X.pinned`,
  `X.pin.order`, `ok`, `route`, `A.full.rank` and the margins
  `anchor.margin`, `pin.margin`, `pure.margin`, `X.pin.margin` (0 =
  exact).

## References

Satoh, K. (2026). Wild Bootstrap Inference for Non-Negative Matrix
Factorization with Random Effects. arXiv:2603.01468.
<https://arxiv.org/abs/2603.01468>

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`summary.nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.inference.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(intercept = 1, speed = cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
result2 <- nmfkc.inference(result, Y, A)
summary(result2)
#> 
#> Call:
#> nmfkc(Y = Y, A = A, rank = 1)
#> 
#> Dimensions: Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50) 
#> Rank (Q):    1 
#> Runtime:     0.0sec 
#> Method:      EU 
#> Iterations: 67 / 5000 (converged)  epsilon = 0.0001  last change = 9.8e-05
#> Missing:     0 (0.0%) 
#> 
#> Statistics:
#>   Objective function:      12982
#>   R-squared (cor^2):       0.6511
#>   R-squared (uncentered):  0.8961
#>   R-squared (centered):    0.601
#>   Residual Std Error:      16.17
#>   Mean Absolute Error:     12.71
#> 
#> Structure Diagnostics:
#>   Basis (X) Sparsity:      0.0% (< 1e-4)
#>   Coef (B) Sparsity:       0.0% (< 1e-4)
#> 
#> Inference (conditional on X):
#>   sigma^2:   270.5 
#> 
#> Coefficients (conditional on X): 2 total, 1 significant
#>        Cov:Basis  Estimate Std. Error (Boot) z value   Pr(>z) 
#> intercept:Basis1     0.156      6.512  3.351    0.02   0.4904  
#>     speed:Basis1     2.900      0.464  0.439    6.25 2.102e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
```

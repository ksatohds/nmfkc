# Summary method for objects of class `nmfre`

Displays a concise summary of an NMF-RE model fit, including dimensions,
convergence, variance components, and a coefficient table following
standard R regression output conventions.

## Usage

``` r
# S3 method for class 'nmfre'
summary(
  object,
  ci.show = FALSE,
  max.coef = 20,
  by = c("covariate", "basis"),
  ...
)

# S3 method for class 'summary.nmfre'
print(x, ...)
```

## Arguments

- object:

  An object of class `nmfre`, returned by
  [`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md).

- ci.show:

  Logical. If `TRUE`, show confidence interval columns (default
  `FALSE`). Named object-first to match the package style (`C.signed`,
  `X.init`, ...). Legacy `show_ci` is accepted via `...`.

- max.coef:

  Largest number of coefficient rows to print (default 20). The table
  has one row per (basis, covariate) pair, so it grows as \\Q\times K\\;
  beyond the cap the significant rows are preferred and the omission is
  reported.

- by:

  Grouping order of the coefficient rows: `"covariate"` (default; list
  all bases for each covariate) or `"basis"` (list all covariates for
  each basis).

- ...:

  Additional arguments (currently unused).

- x:

  An object of class `"summary.nmfre"`.

## Value

The input object, invisibly.

## See also

[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(intercept = 1, speed = cars$speed)
res <- nmfre(Y, A, rank = 1, maxit = 5000)
summary(res)
#> NMF-RE: Y(1,50) = X(1,1) [C(1,2) A + U(1,50)]
#> Iterations: 43 / 5000 (converged)  epsilon = 1e-05  outer 18 / 500
#> R-squared (cor^2):    0.9328 (XB+blup), 0.6511 (XB)
#> R-squared (uncentered):     0.9799 (XB+blup), 0.9091 (XB)
#> R-squared (centered): 0.9227 (XB+blup), 0.6511 (XB)
#> 
#> Variance components:
#>   sigma2 = 106.9  (residual)
#>   tau2   = 120.2  (random effect)
#>   lambda = 0.8891  (sigma2 / tau2)
#>   ICC    = 0.5294  (tau2*tr(X'X) / (tau2*tr(X'X) + sigma2*P))
#>   dfU    = 26.47
#> 
#> Coefficients (Theta): run nmfre.inference(fit, Y, A) for SE / p-values.
```

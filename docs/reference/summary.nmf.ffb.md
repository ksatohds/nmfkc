# Summary method for nmf.ffb objects

Produces a formatted summary of a fitted NMF-FFB model, including matrix
dimensions, convergence, stability diagnostics, fit statistics, and
inference results (if available). For `nmf.ffb(method = "fiml")` objects
the fit statistics are the log-likelihoods, parameter counts and BIC of
the feed-forward null, the unpenalized feedback fit and the BIC-selected
model, the two likelihood-ratio statistics, and – after
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)
– their bootstrap p-values and the false-selection rate under the null.

## Usage

``` r
# S3 method for class 'nmf.ffb'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"nmf.ffb"` returned by
  [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md).

- ...:

  Not used.

## Value

An object of class `"summary.nmf.ffb"` (the fitted model tagged for
printing); printed by
[`print.summary.nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmf.ffb.md).

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)

## Examples

``` r
Y <- t(iris[, -5])
Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
result <- nmf.ffb(Y1, Y2, rank = 2, maxit = 500)
summary(result)
#> NMF-FFB (FIML, X fixed from stage 1): Y1(2,N) = X(2,2) [C1(2,2) Y1 + C2(2,2) Y2]
#> Stage 1 NMF | stage 2 FIML: 391 / 500 (converged)  epsilon = 1e-06  fiml 87 / 3000
#> 
#> Stability diagnostics:
#>   Spectral radius(XC1): 0.0000 (stable)
#>   ||XC1||_1:            0.0000
#> 
#> Likelihood (FF null | unpenalized feedback | BIC-selected):
#>   loglik:     -346.165 |   -346.165 |   -346.165
#>   npar:              9 |         11 |          9
#>   BIC:          737.43 |     747.45 |     737.43
#>   LR vs null: full = 0.000 (free entries = 2), selected = 0.000 (nnz = 0)
#>   selected: C1.L1 = Inf, nnz = 0, rho(XC1) = 0.0000, restriction = union
#> 
#>   Support candidates (best 2 of 2 by BIC; one parameter costs log(N) = 5.01):
#>     nnz=0  rho=0.000 BIC=      737.43 (selected)  feed-forward (no feedback)
#>     nnz=2  rho=0.000 BIC=      747.45 (+10.02)    Factor2<-Sepal.Length, Factor1<-Sepal.Width
#> 
#> Fit statistics:
#>   MAE (mean absolute error):       1.3106
#>   Effective Rank:                  1.58 / 2  (78.9%)
```

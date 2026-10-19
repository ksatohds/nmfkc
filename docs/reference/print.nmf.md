# Print method for fitted NMF models

Prints the call, the model shape and whether the fit converged, in the
spirit of `print.lm`: enough to see what was fitted and whether to trust
it, with everything else left to
[`summary()`](https://rdrr.io/r/base/summary.html). Registered on the
shared `"nmf"` class, so it covers every optimizer that returns one
([`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md) /
`nmf.rrr`,
[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md) and
the signed variants); `nmf.gmm` has its own methods.

Without it, printing a fit dumps the whole list — over 350 lines for a
\\4\times82\\ problem, and unusable for a large one.

Unlike `lm`, the coefficient matrix is **not** shown: it is \\Q\times
K\\ with \\K\\ the number of covariates — the number of samples when
there is no covariate matrix — so it is routinely far too large for a
print method. Use [`coef`](https://rdrr.io/r/stats/coef.html) for it.

## Usage

``` r
# S3 method for class 'nmf'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  A fitted model object inheriting class `"nmf"`.

- digits:

  Number of significant digits.

- ...:

  Additional arguments (currently unused).

## Value

`x`, invisibly.

## See also

[`summary.nmfkc`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.md),
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
result
#> 
#> Call:
#> nmfkc(Y = Y, A = A, rank = 1)
#> 
#> Model:       Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)
#> Convergence: 67 / 5000 (converged)  epsilon = 0.0001  last change = 9.8e-05
#> R-squared:   0.6511
#> 
#> Use coef() for the coefficient matrix, summary() for diagnostics.
#> 
coef(result)   # the coefficient matrix C
#>             Cov1    Cov2
#> Basis1 0.1563808 2.89964
```

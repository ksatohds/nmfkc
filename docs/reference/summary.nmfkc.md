# Summary method for objects of class `nmfkc`

Produces a summary of an `nmfkc` object, including matrix dimensions,
runtime, fit statistics, and diagnostics.

## Usage

``` r
# S3 method for class 'nmfkc'
summary(object, ...)
```

## Arguments

- object:

  An object of class `nmfkc`, i.e., the return value of `nmfkc`.

- ...:

  Additional arguments (currently unused).

## Value

An object of class `summary.nmfkc`, containing summary statistics. Among
them, `separability`, `anchor.purity` and `anchor.row`; see Details.

## Details

**Separability** measures how close the basis is to having an anchor row
for every basis – a row of \\X\\ that carries that basis alone. For each
basis \\q\\ take the largest share it has in any row, \\\max_j \tilde
x\_{jq} / \sum_r \tilde x\_{jr}\\, with \\\tilde X\\ the basis rescaled
so that its columns sum to one (so the value does not depend on
`X.restriction`); these are `anchor.purity`, attained at rows
`anchor.row`. The separability is their minimum: 1 when every basis has
a row of its own, \\1/Q\\ when every row mixes all bases equally, and
`NA` for `rank = 1`. The same three values are stored on the fit, in
`fit$criterion`.

Anchors for every basis, together with \\\mathrm{rank}(B) = Q\\, make
the factorization unique up to scale, so a value near 1 says the
solution is nearly unique. A fit with `X.anchor` scores exactly 1. The
measure looks at \\X\\ only, however: zeros in \\B\\ or in \\\Theta\\
can pin the solution as well, so a low value does not show that it is
not unique.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md),
[`plot.nmfkc`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
summary(result)
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
```

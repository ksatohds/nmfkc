# Rename decoder and encoder bases

Assigns user-specified names to the decoder (X1 columns) and encoder (X2
rows) bases of an `nmfae` object. The names propagate to \\\Theta\\, the
coefficients table, and all downstream displays such as `summary`,
`nmfae.DOT`, and `nmfae.heatmap`.

## Usage

``` r
nmf.rrr.rename(x, X1.colnames = NULL, X2.rownames = NULL)
```

## Arguments

- x:

  An object of class `"nmfae"` returned by
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md).

- X1.colnames:

  Character vector of length \\Q\\ for decoder bases (columns of \\X_1\\
  / rows of \\\Theta\\). If `NULL` (default), the decoder names are left
  unchanged.

- X2.rownames:

  Character vector of length \\R\\ for encoder bases (rows of \\X_2\\ /
  columns of \\\Theta\\). If `NULL` (default), the encoder names are
  left unchanged.

## Value

A modified copy of `x` with updated names.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(runif(15), nrow = 3)
res <- nmf.rrr(Y, rank1 = 2, rank2 = 2)
res <- nmf.rrr.rename(res,
  X1.colnames = c("Resp1", "Resp2"),
  X2.rownames = c("Cov1", "Cov2"))
summary(res)
#> 
#> Call:
#> nmf.rrr(Y1 = Y, rank1 = 2, rank2 = 2)
#> 
#> Dimensions:
#>   Y1:              3 x 5 
#>   Y2:              3 x 5   (autoencoder) 
#>   Response rank Q: 2 
#>   Covariate rank R: 2 
#>   Parameters:      X1(3x2) + C(2x2) + X2(2x3) = 16 
#> 
#> Convergence:
#>   Iterations:      229 / 5000 (converged)  epsilon = 0.0001
#>   Runtime:         0.0 secs 
#> 
#> Goodness of fit:
#>   Objective function:      0.1071
#>   R-squared (cor^2):       0.9132
#>   R-squared (uncentered):  0.9795
#>   R-squared (centered):    0.8988
#>   Residual Std Error:      0.08448
#>   Mean Absolute Error:     0.05914
#>   Effective Rank:          1.34 / 2  (67.0%)
#> 
#> Structure Diagnostics:
#>   Response (X1) Sparsity:  16.7% (< 1e-4)
#>   Bottleneck (C) Sparsity: 0.0% (< 1e-4)
#>   Covariate (X2) Sparsity: 0.0% (< 1e-4)
#> 
# }
```

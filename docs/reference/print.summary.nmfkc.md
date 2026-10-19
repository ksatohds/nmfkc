# Print method for `summary.nmfkc` objects

Prints a formatted summary of an `nmfkc` model fit.

## Usage

``` r
# S3 method for class 'summary.nmfkc'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `summary.nmfkc`.

- digits:

  Minimum number of significant digits to be used.

- ...:

  Additional arguments (currently unused).

## Value

Called for its side effect (printing). Returns `x` invisibly.

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
print(summary(result))
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

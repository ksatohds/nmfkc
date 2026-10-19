# Extract fitted values from NMF models

Returns the reconstructed matrix \\\hat{Y} = X B\\ from a fitted NMF
model.

For `nmf.ffb` objects, returns the equilibrium prediction \\\hat{Y}\_1 =
M\_{model} Y_2\\ if available. Supply `Y1` and `Y2` to get the direct
reconstruction \\X (C_1 Y_1 + C_2 Y_2)\\ instead.

## Usage

``` r
# S3 method for class 'nmf'
fitted(object, ...)

# S3 method for class 'nmfae'
fitted(object, ...)

# S3 method for class 'nmfre'
fitted(object, type = c("blup", "fixed"), ...)

# S3 method for class 'nmf.ffb'
fitted(object, ...)
```

## Arguments

- object:

  A fitted model object of class `"nmf"`, `"nmfkc"`, `"nmfae"`,
  `"nmfre"`, or `"nmf.ffb"`.

- ...:

  For `nmf.ffb`: optionally `Y1` and `Y2`.

- type:

  For `nmfre` objects: `"blup"` (default, the fit \\X(\Theta A + U)\\
  including random effects) or `"fixed"` (\\X\Theta A\\, fixed effects
  only). Chosen so that `Y - fitted(object, type)` equals
  `residuals(object, Y, type)`.

## Value

The fitted matrix \\X B\\.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`residuals.nmf`](https://ksatohds.github.io/nmfkc/reference/residuals.nmf.md)

## Examples

``` r
result <- nmfkc(matrix(runif(50), 5, 10), rank = 2)
#> Y(5,10)~X(5,2)B(2,10)...
#> 0sec
fitted(result)
#>             [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 0.283545147 0.5544245 0.3665403 0.4393818 0.5465527 0.5492524 0.3401122
#> [2,] 0.573324096 0.3797386 0.1954475 0.1570642 0.5411677 0.3675262 0.2439152
#> [3,] 0.697082792 0.3857293 0.1817058 0.1160069 0.6001821 0.3707009 0.2510807
#> [4,] 0.256756109 0.6385265 0.4323792 0.5325227 0.5987466 0.6341660 0.3896859
#> [5,] 0.008325715 0.7966175 0.5851916 0.7827758 0.6097011 0.7983124 0.4771441
#>            [,8]      [,9]     [,10]
#> [1,] 0.40027636 0.4492140 0.6933691
#> [2,] 0.72684369 0.6125685 0.8818754
#> [3,] 0.87528451 0.7144869 1.0189297
#> [4,] 0.37764976 0.4612217 0.7236190
#> [5,] 0.09860759 0.3244998 0.5678558
```

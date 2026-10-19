# Extract residuals from NMF models

Returns the residual matrix \\Y - \hat{Y}\\ from a fitted NMF model.
Requires the original observation matrix `Y` to be supplied.

For `nmfre` objects, residuals are computed from the BLUP reconstruction
(\\Y - X(B\_{blup})\\) by default. Set `type = "fixed"` to use
fixed-effects only.

## Usage

``` r
# S3 method for class 'nmf'
residuals(object, Y, ...)

# S3 method for class 'nmfae'
residuals(object, Y, ...)

# S3 method for class 'nmfre'
residuals(object, Y, type = c("blup", "fixed"), ...)

# S3 method for class 'nmf.ffb'
residuals(object, Y, ...)
```

## Arguments

- object:

  A fitted model object.

- Y:

  The original observation matrix used for fitting.

- ...:

  Not used.

- type:

  For `nmfre` objects: `"blup"` (default) or `"fixed"`.

## Value

The residual matrix.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`fitted.nmf`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.md)

## Examples

``` r
Y <- matrix(runif(50), 5, 10)
result <- nmfkc(Y, rank = 2)
#> Y(5,10)~X(5,2)B(2,10)...
#> 0sec
residuals(result, Y)
#>            [,1]         [,2]         [,3]        [,4]        [,5]        [,6]
#> [1,]  0.3853578 -0.163275964 -0.097967791  0.01995718  0.02734230  0.06501486
#> [2,] -0.2339840 -0.026251849 -0.031518900  0.42130815  0.03967501 -0.20229696
#> [3,]  0.1593010  0.103888577 -0.110554162 -0.39241640  0.09115974  0.01881988
#> [4,] -0.2164343  0.078131349  0.002078884  0.14999504  0.01795844 -0.11960495
#> [5,] -0.2263132 -0.001433012  0.246505947 -0.11394744 -0.16940061  0.17737861
#>             [,7]        [,8]        [,9]       [,10]
#> [1,]  0.02926708 -0.37985832  0.24888837 -0.27154213
#> [2,]  0.33737676  0.09290804 -0.16548991 -0.19477019
#> [3,] -0.04780274  0.17631770 -0.09389127  0.16367036
#> [4,]  0.14231008  0.18015257 -0.16957629  0.02358155
#> [5,] -0.39176047  0.03814927  0.08426470  0.30873492
```

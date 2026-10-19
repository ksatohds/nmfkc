# Extract coefficients from NMF models

Returns the `coefficients` data frame from a fitted NMF model that has
been passed through an inference function
([`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md),
[`nmf.rrr.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.inference.md),
[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md)).

If inference has not been run, returns the parameter matrix \\C\\
(\\\Theta\\) directly.

For `nmf.ffb` objects, returns \\C_2\\ (exogenous block) as fallback.

## Usage

``` r
# S3 method for class 'nmf'
coef(object, ...)

# S3 method for class 'nmf.ffb'
coef(object, ...)
```

## Arguments

- object:

  A fitted model object of class `"nmf"`, `"nmfkc"`, `"nmfae"`,
  `"nmfre"`, or `"nmf.ffb"`.

- ...:

  Not used.

## Value

A data frame of coefficients (if inference was performed), or the
parameter matrix \\C\\.

## See also

[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md),
[`nmf.rrr.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.inference.md),
[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md),
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
coef(result)  # returns C matrix
#>             Cov1    Cov2
#> Basis1 0.1563808 2.89964

result2 <- nmfkc.inference(result, Y, A)
coef(result2)  # returns coefficients data frame
#>    Basis Covariate  Estimate        SE       BSE    z_value      p_value
#> 1 Basis1      Cov1 0.1563808 6.5123713 3.3507470 0.02401288 4.904212e-01
#> 2 Basis1      Cov2 2.8996401 0.4642225 0.4391579 6.24622937 2.102398e-10
#>     CI_low   CI_high on.boundary
#> 1 0.000000 10.740750       FALSE
#> 2 2.129304  3.942492       FALSE
```

# Predict method for nmfre objects

Predicts the response from a fitted NMF-RE model. With `newA` supplied,
returns the **fixed-effect** prediction \\X\\\Theta\\A\_{new}\\ (new
units carry no estimated random effect \\U\\, so only the population
mean is predicted). Without `newA`, returns the in-sample BLUP fit
\\X(\Theta A + U)\\ (`object$XB.blup`), matching
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

## Usage

``` r
# S3 method for class 'nmfre'
predict(object, newA = NULL, ...)
```

## Arguments

- object:

  An object of class `"nmfre"` from
  [`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md).

- newA:

  Optional covariate matrix (K x M) for new units. If `NULL` (default),
  the in-sample BLUP fit is returned.

- ...:

  Not used.

## Value

A matrix of predicted (P x M) values.

## See also

[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`fitted.nmf`](https://ksatohds.github.io/nmfkc/reference/fitted.nmf.md)

# Prediction method for objects of class `nmfkc`

`predict.nmfkc` generates predictions from an object of class `nmfkc`,
either using the fitted covariates or a new covariate matrix.

When the model was fitted using a formula (Formula Mode), a `newdata`
data frame can be supplied instead of `newA`; the covariate matrix is
then constructed automatically from the stored formula metadata.

## Usage

``` r
# S3 method for class 'nmfkc'
predict(
  object,
  newA = NULL,
  newdata = NULL,
  type = c("response", "prob", "class"),
  ...
)
```

## Arguments

- object:

  An object of class `nmfkc`, i.e., the return value of `nmfkc`.

- newA:

  Optional. A new covariate matrix to be used for prediction.

- newdata:

  Optional data frame. Only available when the model was fitted using a
  formula. Covariate columns are extracted automatically using the
  stored formula metadata. If both `newdata` and `newA` are supplied,
  `newdata` takes precedence (with a warning).

- type:

  Type of prediction to return. Options are "response" (fitted values
  matrix), "prob" (soft-clustering probabilities), or "class"
  (hard-clustering labels based on row names of X).

- ...:

  Further arguments passed to or from other methods.

## Value

Depending on `type`: a numeric matrix (`"response"` or `"prob"`) or a
character vector of class labels (`"class"`).

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md)

## Examples

``` r
# Prediction with newA
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
newA <- rbind(1, c(10, 20, 30))
predict(result, newA = newA)
#>          [,1]     [,2]     [,3]
#> [1,] 29.15278 58.14918 87.14558
```

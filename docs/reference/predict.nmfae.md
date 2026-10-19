# Predict method for nmfae objects

`predict.nmfae` computes fitted or predicted values from a three-layer
NMF model. Without `newY2`, returns the in-sample fitted values \\X_1
\Theta X_2 Y_2\\. With `newY2`, computes out-of-sample predictions \\X_1
\Theta X_2 \cdot \mathrm{newY2}\\.

When `type = "class"`, each column is classified to the row with the
maximum predicted value (useful when \\Y_1\\ is a one-hot class matrix
from
[`nmfkc.class`](https://ksatohds.github.io/nmfkc/reference/nmfkc.class.md)).

If `Y1` (actual values) is provided, it is stored as an attribute so
that `plot.predict.nmfae` can produce an observed-vs-predicted scatter
plot (for `type = "response"`) or a confusion matrix heatmap (for
`type = "class"`).

## Usage

``` r
# S3 method for class 'nmfae'
predict(object, newY2 = NULL, Y1 = NULL, type = c("response", "class"), ...)
```

## Arguments

- object:

  An object of class `"nmfae"` returned by
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md).

- newY2:

  Optional new input matrix (P2 x M) for prediction. If `NULL`, returns
  in-sample fitted values.

- Y1:

  Optional actual output matrix for comparison plotting.

- type:

  Character. `"response"` (default) returns the predicted matrix.
  `"class"` returns a factor of predicted class labels (row with max
  value).

- ...:

  Not used.

## Value

For `type = "response"`: a matrix of class `"predict.nmfae"`. For
`type = "class"`: a factor of class `"predict.nmfae"` with predicted
class labels. If `Y1` was provided, actual classes are stored in
`attr(result, "actual")`.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`plot.predict.nmfae`](https://ksatohds.github.io/nmfkc/reference/plot.predict.nmfae.md),
[`nmfkc.class`](https://ksatohds.github.io/nmfkc/reference/nmfkc.class.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(runif(20), nrow = 4)
res <- nmf.rrr(Y, rank1 = 2)
pred <- predict(res)
# }
```

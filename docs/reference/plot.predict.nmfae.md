# Plot method for predict.nmfae objects

For `type = "response"`: if actual values \\Y_1\\ were stored, displays
an observed-vs-predicted scatter plot with \\R^2\\ in the title.
Otherwise, displays the predicted matrix as a heatmap.

For `type = "class"`: if actual classes were stored, displays a
confusion matrix heatmap with accuracy (ACC) in the title.

## Usage

``` r
# S3 method for class 'predict.nmfae'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"predict.nmfae"` returned by
  [`predict.nmfae`](https://ksatohds.github.io/nmfkc/reference/predict.nmfae.md).

- ...:

  Additional graphical parameters passed to `plot` or `image`.

## Value

Invisible `NULL`. Called for its side effect (plot).

## See also

[`predict.nmfae`](https://ksatohds.github.io/nmfkc/reference/predict.nmfae.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(runif(20), nrow = 4)
res <- nmf.rrr(Y, rank1 = 2)
pred <- predict(res)
plot(pred)

# }
```

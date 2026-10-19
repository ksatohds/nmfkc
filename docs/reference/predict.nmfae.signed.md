# Predict method for nmfae.signed

Computes \\\hat Y_1 = X_1 (C\_{+} - C\_{-}) X_2 Y_2^{\mathrm{new}}\\.
Since \\\Theta = C\_{+} - C\_{-}\\ is signed, predictions may contain
negative entries even when \\Y_1 \ge 0\\ in training.

## Usage

``` r
# S3 method for class 'nmfae.signed'
predict(object, newY2 = NULL, Y1 = NULL, type = c("response", "class"), ...)
```

## Arguments

- object:

  A fitted `"nmfae.signed"` object.

- newY2:

  New input matrix (P2 x N_new). If `NULL`, returns the training fitted
  values.

- Y1:

  Optional reference Y1 for scatter / confusion plot.

- type:

  Output: `"response"` (raw signed) or `"class"`.

- ...:

  Unused.

## Value

A numeric matrix (`"response"`) or factor (`"class"`).

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)

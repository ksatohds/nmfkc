# Print/plot methods for nmfkc cross-validation objects

Compact print and (for rank sweeps) a score-vs-rank plot with a “Best
(Min)” marker, for the objects returned by
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md),
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md) and
[`nmfkc.bicv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.bicv.md).

## Usage

``` r
# S3 method for class 'nmfkc.ecv'
print(x, ...)

# S3 method for class 'nmfkc.ecv'
plot(x, main = "nmfkc element-wise CV", ylab = NULL, ...)

# S3 method for class 'nmfkc.bicv'
print(x, ...)

# S3 method for class 'nmfkc.bicv'
plot(
  x,
  main = "nmfkc bi-cross-validation",
  ylab = "sigma (held-out RMSE)",
  ...
)

# S3 method for class 'nmfkc.cv'
print(x, ...)

# S3 method for class 'nmfkc.cv'
plot(
  x,
  main = "nmfkc CV: per-fold held-out loss",
  ylab = "objfunc (block)",
  ...
)
```

## Arguments

- x:

  A `"nmfkc.ecv"` / `"nmfkc.cv"` / `"nmfkc.bicv"` object.

- ...:

  Passed to [`plot`](https://rdrr.io/r/graphics/plot.default.html).

- main, ylab:

  Plot labels.

## Value

`x`, invisibly.

# Plot method for nmfae.cv objects

Displays a bar chart of per-fold cross-validation errors from
[`nmf.rrr.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.cv.md).
The overall RMSE (sigma) is shown in the title.

## Usage

``` r
# S3 method for class 'nmfae.cv'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"nmfae.cv"` returned by
  [`nmf.rrr.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.cv.md).

- ...:

  Additional graphical parameters passed to `barplot`.

## Value

Invisible `NULL`. Called for its side effect (plot).

## See also

[`nmf.rrr.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.cv.md)

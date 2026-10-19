# Print method for nmf.rrr element-wise CV results

A two-line summary, matching `print.nmfkc.ecv`. Without it the object
printed as a raw list (33 lines of `$objfunc` / `$folds`).

## Usage

``` r
# S3 method for class 'nmfae.ecv'
print(x, ...)
```

## Arguments

- x:

  An object of class `"nmfae.ecv"` (or `"nmfae.signed.ecv"`).

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`nmf.rrr.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md),
[`plot.nmfae.ecv`](https://ksatohds.github.io/nmfkc/reference/plot.nmfae.ecv.md)

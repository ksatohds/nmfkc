# Print method for summary.nmfae objects

Prints a formatted summary of an NMF-AE model fit.

## Usage

``` r
# S3 method for class 'summary.nmfae'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  max.coef = 20,
  by = c("covariate", "basis"),
  ...
)
```

## Arguments

- x:

  An object of class `"summary.nmfae"`.

- digits:

  Minimum number of significant digits to be used.

- max.coef:

  Maximum number of coefficient rows to display. If the table has more
  rows, only significant rows (p \< 0.05) are shown. Default is 20.

- by:

  Character; grouping order of the coefficients table. `"covariate"`
  (default) lists all bases within each covariate (1-1, 1-2, ...);
  `"basis"` lists all covariates within each basis (1-1, 2-1, ...).

- ...:

  Additional arguments (currently unused).

## Value

Called for its side effect (printing). Returns `x` invisibly.

## See also

[`summary.nmfae`](https://ksatohds.github.io/nmfkc/reference/summary.nmfae.md)

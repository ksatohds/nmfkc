# Print method for summary.nmfkc.inference objects

Prints a formatted summary including the coefficients table.

## Usage

``` r
# S3 method for class 'summary.nmfkc.inference'
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

  An object of class `"summary.nmfkc.inference"`.

- digits:

  Minimum number of significant digits.

- max.coef:

  Maximum coefficient rows to display. Default is 20.

- by:

  Character; grouping order of the coefficients table. `"covariate"`
  (default) lists all bases within each covariate (1-1, 1-2, ...);
  `"basis"` lists all covariates within each basis (1-1, 2-1, ...).

- ...:

  Additional arguments (currently unused).

## Value

Called for its side effect (printing). Returns `x` invisibly.

## See also

[`summary.nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.inference.md)

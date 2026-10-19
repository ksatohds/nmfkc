# Print method for summary.nmfkc.net.inference objects

Prints a formatted summary of a symmetric NMF model with inference
results. The coefficients table uses `Basis.row` and `Basis.col` labels
reflecting the symmetric structure \\C\_{qr}\\.

## Usage

``` r
# S3 method for class 'summary.nmfkc.net.inference'
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

  An object of class `"summary.nmfkc.net.inference"`.

- digits:

  Minimum number of significant digits.

- max.coef:

  Largest number of coefficient rows to print (default 20). The table
  has one row per basis pair, so it grows as \\Q^2\\; beyond the cap the
  significant rows are preferred and the omission is reported.

- by:

  Character; grouping order of the coefficients table. Because the
  symmetric model sets \\A = X^\top\\, the column factor `Basis.col`
  plays the covariate role: `"covariate"` (default) lists all
  `Basis.row` within each `Basis.col` (1-1, 2-1, ...), while `"basis"`
  lists all `Basis.col` within each `Basis.row` (1-1, 1-2, ...).

- ...:

  Additional arguments (currently unused).

## Value

Called for its side effect (printing). Returns `x` invisibly.

## See also

[`summary.nmfkc.net.inference`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.net.inference.md)

# Print method for summary.nmfae.signed.inference objects

Prints the Signed-Bottleneck NMF-AE summary followed by the coefficients
table of Theta.

## Usage

``` r
# S3 method for class 'summary.nmfae.signed.inference'
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

  An object of class `"summary.nmfae.signed.inference"`.

- digits:

  Minimum number of significant digits.

- by:

  Character; grouping order of the coefficients table. `"covariate"`
  (default) lists all bases within each covariate (1-1, 1-2, ...);
  `"basis"` lists all covariates within each basis (1-1, 2-1, ...).

- ...:

  Additional arguments (currently unused).

## Value

Called for its side effect (printing). Returns `x` invisibly.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`summary.nmfae.signed.inference`](https://ksatohds.github.io/nmfkc/reference/summary.nmfae.signed.inference.md)

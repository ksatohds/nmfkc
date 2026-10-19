# Summary method for nmfae.signed

Produces a summary with dimensions, convergence, fit statistics, and
structure diagnostics (sparsity and negative-mass ratio).

## Usage

``` r
# S3 method for class 'nmfae.signed'
summary(object, ...)
```

## Arguments

- object:

  An `nmfae.signed` object.

- ...:

  Unused.

## Value

An object of class `"summary.nmfae.signed"`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)

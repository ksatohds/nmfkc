# Summary method for nmfae.signed.inference objects

Produces a summary of a fitted Signed-Bottleneck NMF-AE model with
inference results. Extends
[`summary.nmfae.signed`](https://ksatohds.github.io/nmfkc/reference/summary.nmfae.signed.md)
by attaching the `coefficients` table and p-value side from
[`nmf.rrr.signed.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.inference.md).

## Usage

``` r
# S3 method for class 'nmfae.signed.inference'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"nmfae.signed.inference"`.

- ...:

  Additional arguments (currently unused).

## Value

An object of class `"summary.nmfae.signed.inference"`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr.signed.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.inference.md),
[`summary.nmfae.signed`](https://ksatohds.github.io/nmfkc/reference/summary.nmfae.signed.md)

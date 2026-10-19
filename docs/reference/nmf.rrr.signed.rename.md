# Rename Resp/Cov labels on nmfae.signed objects

Replaces the default `"Resp1", "Resp2", ...` (response basis / X1
columns and Cp/Cn/C rows) and `"Cov1", "Cov2", ...` (covariate basis /
X2 rows and Cp/Cn/C columns) with user-supplied labels. Propagates to
coefficients tables if present (e.g., from `nmfae.inference`).

## Usage

``` r
nmf.rrr.signed.rename(x, X1.colnames = NULL, X2.rownames = NULL)
```

## Arguments

- x:

  An `"nmfae.signed"` object.

- X1.colnames:

  Character vector of length Q for decoder labels.

- X2.rownames:

  Character vector of length R for encoder labels.

## Value

The renamed object (same class).

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)

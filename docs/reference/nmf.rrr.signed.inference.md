# Statistical Inference for Signed-Bottleneck NMF-AE Signed Bottleneck

Post-estimation inference for the **signed** bottleneck \\\Theta =
C\_{+} - C\_{-}\\ in the Signed-Bottleneck NMF-AE model \\Y_1 \approx
X_1 \Theta X_2 Y_2\\, conditional on \\(\hat X_1, \hat X_2)\\. Uses
sandwich covariance and wild bootstrap **without** the non-negativity
projection that
[`nmf.rrr.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.inference.md)
applies (because \\\Theta\\ is unconstrained in sign here).

## Usage

``` r
nmf.rrr.signed.inference(object, Y1, Y2 = Y1, wild.bootstrap = TRUE, ...)
```

## Arguments

- object:

  A fitted `"nmfae.signed"` object.

- Y1:

  Output matrix used during fitting.

- Y2:

  Input matrix used during fitting. Default `Y1`.

- wild.bootstrap:

  Logical. Default `TRUE`.

- ...:

  Additional arguments:

  `wild.B`

  :   Bootstrap replicates. Default 500.

  `wild.seed`

  :   RNG seed. Default 123.

  `wild.level`

  :   CI confidence level. Default 0.95.

  `sandwich`

  :   Use sandwich covariance. Default `TRUE`.

  `C.p.side`

  :   P-value type: `"two.sided"` (default for Signed-Bottleneck NMF-AE)
      or `"one.sided"`.

  `cov.ridge`

  :   Ridge stabilization. Default `1e-8`.

  `print.trace`

  :   Logical. Default `FALSE`.

## Value

An object of class
`c("nmfae.signed.inference", "nmfae.inference", "nmfae.signed", "nmfae", "nmf")`
with added fields:

- sigma2.used:

  Estimated \\\sigma^2\\.

- C.se, C.se.boot:

  Sandwich / bootstrap SEs for \\\Theta\\ (Q x R).

- C.ci.lower, C.ci.upper:

  Bootstrap CIs.

- coefficients:

  Data frame with Estimate, SE, BSE, z, p-value, CI.

- C.p.side:

  P-value side used.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md),
[`nmf.rrr.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.inference.md)

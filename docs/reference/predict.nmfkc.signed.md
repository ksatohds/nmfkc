# Predict method for nmfkc.signed

Computes \\\widehat Y = X \\ C \\ A\_{\mathrm{new}}\\ (\\= X (C\_{+} -
C\_{-})(A\_{+}^{\mathrm{new}} - A\_{-}^{\mathrm{new}})\\). For
`type = "response"` the raw prediction is returned (possibly signed).

With signed covariates the scores \\\bm b_n = C\bm a_n\\ can be
negative, so the NMF-LAB recipe of normalizing \\\bm b_n\\ to a
membership vector does not apply. Since \\X\\ is column-stochastic,
\\\widehat{\bm y}\_n = X\bm b_n\\ still sums to \\\sum_q b\_{qn}\\
(close to one for one-hot targets) and is the least-squares estimate of
the class indicator; for `type = "prob"` it is mapped to the probability
simplex by the **Euclidean projection** \$\$\widehat{\bm p}\_n =
\mathop{\mathrm{arg\\min}}\_{\bm p \ge 0,\\ \bm 1^\top \bm p = 1} \lVert
\bm p - \widehat{\bm y}\_n \rVert^2 = \max(\widehat{\bm y}\_n - \tau_n
\bm 1, 0),\$\$ the closest probability vector in the same Frobenius
geometry the fit minimizes (Duchi et al. 2008; Wang & Carreira-Perpinan
2013; \\O(P \log P)\\ per column by sorting). The former rule – clip
negatives to zero and renormalize – is kept as `prob.method = "clip"`.
`type = "class"` is \\\arg\max_p \widehat y\_{pn}\\, which both rules
preserve.

## Usage

``` r
# S3 method for class 'nmfkc.signed'
predict(object, newA = NULL, type = c("response", "prob", "class"), ...)
```

## Arguments

- object:

  A fitted `"nmfkc.signed"` object.

- newA:

  Real-valued \\D \times N\_{\mathrm{new}}\\ covariate matrix.

- type:

  Output: `"response"` (raw signed), `"prob"`, or `"class"`.

- ...:

  Hidden option `prob.method`: how `type = "prob"` maps \\\widehat{\bm
  y}\_n\\ to the simplex, `"simplex"` (default; Euclidean projection) or
  `"clip"` (clip negatives, renormalize).

## Value

A numeric matrix (`"response"` or `"prob"`) or a character vector
(`"class"`).

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

Duchi, J., Shalev-Shwartz, S., Singer, Y., & Chandra, T. (2008).
Efficient projections onto the l1-ball for learning in high dimensions.
*Proceedings of the 25th International Conference on Machine Learning
(ICML)*, 272–279.

Wang, W., & Carreira-Perpinan, M. A. (2013). Projection onto the
probability simplex: an efficient algorithm with a simple proof, and an
application. *arXiv:1309.1541*.

## See also

[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md),
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.signed.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.cv.md)

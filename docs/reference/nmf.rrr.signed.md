# Signed-Bottleneck NMF-AE: Three-Layer NMF-AE with Signed Bottleneck

`nmfae.signed` fits a three-layer non-negative matrix factorization
autoencoder with a **signed** bottleneck, solving \$\$Y_1 \approx X_1
(C\_{+} - C\_{-}) X_2 Y_2, \quad X_1 \ge 0,\\ C\_{+} \ge 0,\\ C\_{-} \ge
0,\\ X_2 \ge 0,\$\$ where \\\Theta = C\_{+} - C\_{-}\\ is the signed
bottleneck. The basis matrices \\X_1\\ (columns sum to 1) and \\X_2\\
(rows sum to 1) retain their non-negative "parts-based"
interpretability, while \\\Theta\\ can express anti-correlations (e.g.,
refractive index up vs. Abbe number down).

The algorithm uses **Direct Multiplicative Updates** derived from Ding
et al. (2010) sign-splitting technique, applied block-wise to the four
non-negative blocks \\(C\_{+}, C\_{-}, X_1, X_2)\\. Each block update
monotonically decreases the true objective \\\\Y_1 - X_1(C\_{+} -
C\_{-})X_2 Y_2\\\_F^2\\ (Lee-Seung auxiliary function method).

**Relation to
[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md):**
When \\\Theta \ge 0\\ suffices (the `nmfae` case), `nmfae.signed`
reduces to `nmfae` up to the \\C\_{+} - C\_{-}\\ parameterization. Use
`nmfae.signed` when the data exhibit negative cross-property
correlations that tri-NMF-AE cannot express (e.g., high refractive index
\<-\> low Abbe number trade-off).

## Usage

``` r
nmf.rrr.signed(
  Y1,
  Y2 = Y1,
  rank1 = 2,
  rank2 = NULL,
  epsilon = 1e-04,
  maxit = 5000,
  verbose = FALSE,
  ...
)
```

## Arguments

- Y1:

  Output matrix \\Y_1\\ (P1 x N). Negative entries allowed
  (`Y.signed = TRUE` is auto-detected).

- Y2:

  Input matrix \\Y_2\\ (P2 x N). Must be non-negative. Default `Y1`
  (autoencoder).

- rank1:

  Integer. Response-basis rank Q. Default 2.

- rank2:

  Integer. Covariate-basis rank R. Default (`NULL`) = `rank1`.

- epsilon:

  Relative convergence tolerance on the objective. Default `1e-4`.

- maxit:

  Maximum iterations. Default 5000.

- verbose:

  Logical. Print progress. Default `FALSE`.

- ...:

  Additional arguments:

  `warm.start`

  :   One of `TRUE` (default, **hybrid**: warm-start \\X_1, X_2\\ from
      [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md)
      but initialize \\C\_{+}, C\_{-}\\ randomly), `"full"` (warm-start
      everything including \\C\_{+} = C\_{\mathrm{tri}}\\, \\C\_{-} =
      \delta\\), or `FALSE` (random for all blocks). The hybrid default
      avoids the \\C\_{-} = 0\\ local-minimum trap inherited from
      tri-NMF-AE while still benefiting from good \\X_1, X_2\\
      initialization. Ignored when \\Y_1\\ has negative entries.

  `nstart`

  :   Integer, default 1 (cf.
      [`kmeans`](https://rdrr.io/r/stats/kmeans.html)). Number of random
      restarts for \\C\_{+}, C\_{-}\\. Each restart uses seed
      `seed + 7919 * (s-1)`. Returns the best run by final objective.
      **Signed models have more local minima than non-negative ones**
      because the bottleneck \\\Theta = C\_{+} - C\_{-}\\ can take both
      positive and negative values; during exploration a larger `nstart`
      (e.g., 10-50) reduces the chance of being trapped at an inferior
      stationary point (particularly the \\C\_{-} = 0\\ trap from
      warm-start from non-negative tri-NMF-AE). Use the default 1 for
      fast development and raise for publication-grade runs.

  `cores`

  :   Run the `nstart` restarts in parallel (only when `nstart > 1`).
      Default `getOption("mc.cores", 1L)` (PSOCK cluster on Windows,
      forking elsewhere). Each restart is fully self-seeded and the best
      run is selected by first-minimum objective, so the returned fit is
      identical for any `cores`.

  `X1.L2.ortho`, `X2.L2.ortho`

  :   Non-negative L2 orthogonality penalties (default 0) on the
      **columns** of \\X_1\\ and the **rows** of \\X_2\\, penalizing
      \\(\lambda/2)\lVert\mathrm{offdiag}(X_1^\top X_1)\rVert^2\\ and
      \\(\lambda/2)\lVert\mathrm{offdiag}(X_2 X_2^\top)\rVert^2\\
      respectively. Same convention as
      [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md);
      encourage more distinct (less overlapping) response / covariate
      bases.

  `C.L2`

  :   Non-negative ridge penalty (default 0) on the signed bottleneck
      \\C = C\_{+} - C\_{-}\\, adding \\\lambda\lVert C\_{+} -
      C\_{-}\rVert^2\\. Shrinks \\\Theta\\ toward zero (with zero
      gradient on the unidentified common mode \\C\_{+} + C\_{-}\\),
      injected into both the unweighted and weighted
      \\C\_{+}\\/\\C\_{-}\\ updates.

  `Y1.weights`

  :   Optional non-negative weight matrix (P1 x N) or vector (length N)
      for \\Y_1\\, analogous to the `weights` argument of
      [`lm`](https://rdrr.io/r/stats/lm.html). Loss becomes \\\sum
      W\_{ij} \\ (Y\_{1,ij} - \hat Y\_{1,ij})^2\\
      ([`lm()`](https://rdrr.io/r/stats/lm.html)-style, **linear** in
      \\W\\). Logical matrices (`TRUE` / `FALSE`) are also accepted.
      Used by
      [`nmf.rrr.signed.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.ecv.md)
      to hold out test elements via a binary mask \\W \in \\0,1\\\\;
      real-valued weights for importance weighting are also supported.
      Default: if `Y1` has `NA`, a binary mask is auto-generated (0 for
      `NA`, 1 elsewhere).

  `Cp.init`, `Cn.init`

  :   Explicit \\Q \times R\\ non-negative matrices for initialization.
      Overrides warm.start.

  `C.init`

  :   Explicit signed \\Q \times R\\ matrix, internally split into
      \\(C\_{+}, C\_{-})\\.

  `X1.init`, `X2.init`

  :   Explicit basis matrices.

  `seed`

  :   RNG seed. Default 123.

  `print.trace`

  :   Logical. Print iteration trace. Default `FALSE`.

  `prefix.dec`, `prefix.enc`

  :   Label prefixes for the response/covariate bases. Default `"Resp"`,
      `"Cov"`.

  Rank aliases accepted here for backward compatibility: `Q` for
  `rank1`, `R` for `rank2`.

## Value

An object of class `c("nmfae.signed", "nmfae", "nmf")` with:

- X1:

  Decoder basis (P1 x Q), column sum 1.

- Cp, Cn:

  Non-negative parts of \\\Theta\\ (each Q x R).

- C:

  Signed bottleneck \\\Theta = C\_{+} - C\_{-}\\ (Q x R).

- X2:

  Encoder basis (R x P2), row sum 1.

- Y1hat:

  Fitted values \\X_1 (C\_{+} - C\_{-}) X_2 Y_2\\.

- B1:

  Decoder-side scores \\B_1 = (C\_{+} - C\_{-}) X_2 Y_2\\ (Q x N), so
  that \\\widehat Y_1 = X_1 B_1\\. **Signed**, since \\C\\ is.

- B2:

  Encoder-side scores \\B_2 = X_2 Y_2\\ (R x N), i.e. \\B_1\\ before the
  \\C\\ map. Non-negative.

- H:

  **Deprecated**; identical to `B1`. Kept so that code written against
  earlier releases keeps working. Use `B1`.

- B.prob, B.cluster:

  Sample-level soft/hard clustering from a non-negative clip of the
  (possibly signed) \\B_1\\; use with care since \\C\\ may be signed.
  There are deliberately no `B1.prob` / `B2.prob` counterparts: a
  column-normalized membership is not a distribution when the scores can
  be negative.

- X1.prob:

  Row-normalized \\X_1\\ (P1 x Q): each RESPONSE VARIABLE's soft
  membership over the Q response groups. Well defined (\\X_1 \ge 0\\).

- X1.cluster:

  Hard response-group label per response variable.

- X2.prob:

  Column-normalized \\X_2\\ (R x P2): each COVARIATE VARIABLE's soft
  membership over the R covariate groups. Well defined (\\X_2 \ge 0\\).

- X2.cluster:

  Hard covariate-group label per covariate variable.

- rank:

  `c(Q = Q, R = R)`.

- dims:

  `c(P1, P2, N)`.

- objfunc, objfunc.iter:

  Final and per-iteration objective values.

- r.squared:

  \\\mathrm{cor}(Y_1, \widehat Y_1)^2\\ (Pearson; in \\\[0,1\]\\).

- r.squared.uncentered:

  Uncentered \\R^2 = 1 - \\Y_1 - \widehat Y_1\\\_F^2 / \\Y_1\\\_F^2\\
  (baseline = zero matrix).

- r.squared.centered:

  Row-mean centered \\1 - \\Y_1 - \widehat Y_1\\\_F^2 / \\Y_1 - \bar
  Y\_{p\cdot}\\\_F^2\\.

- sigma, mae:

  Residual SE and mean absolute error.

- niter, runtime:

  Iterations and elapsed seconds.

- Y.signed:

  Logical; whether \\Y_1\\ contained negative entries.

- call:

  Matched call.

## Lifecycle

This function is **experimental**; interface may change.

## References

Satoh, K. and Tokuda, Y. (2026). Co-clustering of Response and Covariate
Variables by Tri-Factorizing Their Non-negative Regression Coefficient
Matrix. *arXiv preprint* arXiv:2607.27474.
[doi:10.48550/arXiv.2607.27474](https://doi.org/10.48550/arXiv.2607.27474)

Ding, C.H.Q., Li, T., and Jordan, M.I. (2010). Convex and
Semi-Nonnegative Matrix Factorizations. *IEEE TPAMI*, 32(1), 45-55.

Satoh, K. (2026). Signed-Bottleneck NMF-AE: Signed-Bottleneck 3-Layer
NMF (research memo, 2026-04-18).

Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
semi-nonnegative matrix factorizations. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 32(1), 45–55.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`predict.nmfae.signed`](https://ksatohds.github.io/nmfkc/reference/predict.nmfae.signed.md),
[`summary.nmfae.signed`](https://ksatohds.github.io/nmfkc/reference/summary.nmfae.signed.md),
[`nmf.rrr.signed.rename`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.rename.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y1 <- matrix(abs(rnorm(12)), 3, 4)
Y2 <- matrix(abs(rnorm(20)), 5, 4)
res <- nmf.rrr.signed(Y1, Y2, rank1 = 2, rank2 = 2, maxit = 500)
summary(res)
#> 
#> Call:
#> nmf.rrr.signed(Y1 = Y1, Y2 = Y2, rank1 = 2, rank2 = 2, maxit = 500)
#> 
#> Signed-Bottleneck NMF-AE (Direct MU; signed bottleneck Theta = Cp - Cn)
#> 
#> Dimensions:
#>   P1 (output rows): 3
#>   P2 (input rows):  5
#>   N  (samples):     4
#>   Q (decoder rank): 2
#>   R (encoder rank): 2
#>   Parameters:       24
#> 
#> Convergence:
#>   Iterations:       16 / 500 (converged)  epsilon = 0.0001
#>   Runtime (secs):   0.03
#>   Final objfunc:    2.128
#> 
#> Goodness of fit:
#>   R-squared (cor^2):    0.04789
#>   R-squared (uncentered):     0.7371
#>   R-squared (centered): 0.03128
#>   Sigma (RMSE):         0.4211
#>   MAE:                  0.344
#> 
#> Structure (range / sparsity / negative mass):
#>   X1 (P1 x Q):     range [6.232e-139, 0.9987]   sparsity  16.7%
#>   Cp (Q x R):      range [0, 2.459]   sparsity  75.0%
#>   Cn (Q x R):      range [0, 0.01288]   sparsity  50.0%
#>   C=Cp-Cn (Q x R): range [-0.01288, 2.459]   sparsity  25.0%   neg-mass   0.6%
#>   X2 (R x P2):     range [1.027e-30, 0.8354]   sparsity  30.0%
#>   Encoding H (Q x N):                                   neg-mass   0.1%
#> 
# }
```

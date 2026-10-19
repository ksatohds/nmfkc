# Three-Layer Non-negative Matrix Factorization (NMF-AE)

`nmfae` fits a three-layer nonnegative matrix factorization model \\Y_1
\approx X_1 \Theta X_2 Y_2\\, where \\X_1\\ is a decoder basis (column
sum 1), \\\Theta\\ is a bottleneck parameter matrix, \\X_2\\ is an
encoder basis (row sum 1), and \\Y_2\\ is the input matrix.

When `Y2 = Y1`, the model acts as a non-negative autoencoder. When
`Y1 != Y2`, it acts as a heteroencoder.

Initialization uses a three-step NMF procedure via
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md): (1)
`nmfkc(Y1, rank=Q)` to obtain \\X_1\\, (2) `nmfkc(Y1, A=Y2, rank=Q)`
with fixed \\X_1\\ to obtain \\C = \Theta X_2\\, (3) `nmfkc(Y2, rank=R)`
to factor \\C\\ into \\\Theta\\ and \\X_2\\.

## Usage

``` r
nmf.rrr(
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

## Source

Satoh, K. (2025). Applying Non-negative Matrix Factorization with
Covariates to Multivariate Time Series. *Japanese Journal of Statistics
and Data Science*.

## Arguments

- Y1:

  Output matrix \\Y_1\\ (P1 x N). Non-negative. May contain `NA`s
  (handled via `Y1.weights`).

- Y2:

  Input matrix \\Y_2\\ (P2 x N). Non-negative. Default is `Y1`
  (autoencoder).

- rank1:

  Integer. Rank of the response basis \\X_1\\ (P1 x Q). Default is 2.

- rank2:

  Integer. Rank of the covariate basis \\X_2\\ (R x P2). Default
  (`NULL`) sets `rank2 = rank1`.

- epsilon:

  Positive convergence tolerance. Default is `1e-4`.

- maxit:

  Maximum number of multiplicative update iterations. Default is 5000.

- verbose:

  Logical. If `TRUE`, prints progress messages during fitting. Default
  is `FALSE`.

- ...:

  Additional arguments:

  `method`

  :   Objective function: Euclidean distance `"EU"` (default) or
      Kullback-Leibler divergence `"KL"`. Both use Lee-Seung
      multiplicative updates for the three factors \\X_1, \Theta, X_2\\
      of \\Y_1 \approx X_1 \Theta X_2 Y_2\\; for `"KL"` the residual SE
      `sigma` is `NA` (not on the data scale).

  `Y1.weights`

  :   Optional non-negative weight matrix (P1 x N) or vector for
      \\Y_1\\, analogous to the `weights` argument of
      [`lm`](https://rdrr.io/r/stats/lm.html). Loss becomes \\\sum
      W\_{ij} \\ (Y\_{1,ij} - \hat Y\_{1,ij})^2\\
      ([`lm()`](https://rdrr.io/r/stats/lm.html)-style, **linear** in
      \\W\\). Logical matrices (`TRUE` / `FALSE`) are also accepted.
      Typical ECV / CV usage passes a binary mask \\W \in \\0,1\\\\ for
      held-out elements; real-valued weights for importance weighting
      are also supported. Default: if `Y1` has `NA`, a binary mask is
      auto-generated (0 for `NA`, 1 elsewhere).

  `C.L1`

  :   L1 regularization parameter for \\C\\. Default is 0.

  `X1.L2.ortho`

  :   L2 orthogonality regularization for \\X_1\\ columns. Default is 0.

  `X2.L2.ortho`

  :   L2 orthogonality regularization for \\X_2\\ rows. Default is 0.

  `seed`

  :   Integer seed for reproducibility. Default is 123.

  `nstart`

  :   Number of random restarts for the
      [`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
      initialisation steps (passed to the k-means initialiser of the
      \\X_1\\ and \\X_2\\ factorisations). Default `1` (single start;
      the historical behaviour). A larger value (e.g.\\ 10-20) gives a
      more stable initialisation and is recommended before inference.

  `print.trace`

  :   Logical. If `TRUE`, prints progress. Default is `FALSE`.

  `X.init`

  :   Initialization method for the two bases, passed to the internal
      [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
      steps: `"kmeans"` (default), `"nndsvd"`, `"spa"`, ...

  `X.anchor`

  :   Anchor rows of the response basis \\X_1\\, imposed as in
      [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md):
      `NULL` (default), `"spa"`, or one row of `Y1` per basis. The zeros
      are kept for the whole fit. The scores of \\X_1\\ are \\\Theta X_2
      Y_2\\, so `rank1 <= min(rank2, nrow(Y2))` is required. The rows
      used are returned as `X.anchor`.

  Rank aliases accepted here for backward compatibility: `Q` for
  `rank1`, `R` for `rank2`.

## Value

An object of class `"nmfae"`, a list with components:

- X1:

  Decoder basis matrix (P1 x Q), column sum 1.

- C:

  Parameter matrix (Q x R).

- X2:

  Encoder basis matrix (R x P2), row sum 1.

- Y1hat:

  Fitted values \\X_1 \Theta X_2 Y_2\\ (P1 x N).

- B1:

  Decoder-side scores \\B_1 = C X_2 Y_2\\ (Q x N), so that \\\widehat
  Y_1 = X_1 B_1\\. Analogue of `B` in
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md).

- H:

  **Deprecated**; identical to `B1`. Kept so that code written against
  earlier releases keeps working. Use `B1`.

- B2:

  Encoder-side scores \\B_2 = X_2 Y_2\\ (R x N), i.e. \\B_1\\ before the
  \\C\\ map.

- B1.prob, B2.prob:

  Column-normalized \\B_1\\ / \\B_2\\: each SAMPLE's soft membership
  over the Q response groups and over the R covariate groups
  respectively.

- B1.cluster, B2.cluster:

  Hard label per sample (argmax over the columns of the corresponding
  `.prob`).

- X1.prob:

  Row-normalized \\X_1\\ (P1 x Q): each RESPONSE VARIABLE's soft
  membership over the Q response groups.

- X1.cluster:

  Hard response-group label per response variable (argmax over `X1.prob`
  rows).

- X2.prob:

  Column-normalized \\X_2\\ (R x P2): each COVARIATE VARIABLE's soft
  membership over the R covariate groups.

- X2.cluster:

  Hard covariate-group label per covariate variable (argmax over
  `X2.prob` columns).

- rank:

  Named integer vector `c(Q, R)`.

- method:

  Objective used (`"EU"` or `"KL"`).

- objfunc:

  Final objective value.

- objfunc.iter:

  Objective values by iteration.

- r.squared:

  \\\mathrm{cor}(Y, \widehat Y)^2\\ (Pearson; in \\\[0,1\]\\).

- r.squared.uncentered:

  Uncentered \\R^2 = 1 - \\Y - \widehat Y\\\_F^2 / \\Y\\\_F^2\\
  (baseline = zero matrix).

- r.squared.centered:

  Row-mean centered \\1 - \\Y - \widehat Y\\\_F^2 / \\Y - \bar
  Y\_{p\cdot}\\\_F^2\\.

- niter:

  Number of iterations performed.

- runtime:

  Elapsed time as a `difftime` object.

- n.missing:

  Number of missing (or zero-weighted) elements in \\Y_1\\.

- n.total:

  Total number of elements in \\Y_1\\ (P1 x N).

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Satoh, K. and Tokuda, Y. (2026). Co-clustering of Response and Covariate
Variables by Tri-Factorizing Their Non-negative Regression Coefficient
Matrix. *arXiv preprint* arXiv:2607.27474.
[doi:10.48550/arXiv.2607.27474](https://doi.org/10.48550/arXiv.2607.27474)

Lee, D. D. and Seung, H. S. (2001). Algorithms for Non-negative Matrix
Factorization. *Advances in Neural Information Processing Systems*, 13.

Saha, S. et al. (2022). Hierarchical Deep Learning Neural Network
(HiDeNN): An Artificial Intelligence (AI) Framework for Computational
Science and Engineering. *Computer Methods in Applied Mechanics and
Engineering*, 399.

## See also

[`nmf.rrr.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.inference.md),
[`predict.nmfae`](https://ksatohds.github.io/nmfkc/reference/predict.nmfae.md),
[`nmf.rrr.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.ecv.md),
[`nmf.rrr.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.DOT.md),
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)

## Examples

``` r
# Autoencoder example
Y <- matrix(c(1,0,1,0, 0,1,0,1, 1,1,0,0), nrow=3, byrow=TRUE)
res <- nmf.rrr(Y, rank1=2, rank2=2)
res$r.squared
#> [1] 0.7842934

# Heteroencoder example
Y1 <- matrix(c(1,0,0,1), nrow=2)
Y2 <- matrix(runif(8), nrow=4)
res2 <- nmf.rrr(Y1, Y2, rank1=2, rank2=2)
```

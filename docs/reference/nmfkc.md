# Optimize NMF with kernel covariates (Full Support for Missing Values)

`nmfkc` fits a nonnegative matrix factorization with kernel covariates
under the tri-factorization model \\Y \approx X C A = X B\\.

This function supports two major input modes:

1.  **Matrix Mode (Existing)**: `nmfkc(Y=matrix, A=matrix, ...)`

2.  **Formula Mode (New)**:
    `nmfkc(formula=Y_vars ~ A_vars, data=df, rank=Q, ...)`

The rank of the basis matrix can be specified using either the `rank`
argument (preferred for formula mode) or the hidden `Q` argument (for
backward compatibility).

## Usage

``` r
nmfkc(
  Y,
  A = NULL,
  rank = NULL,
  data,
  epsilon = 1e-04,
  maxit = 5000,
  verbose = TRUE,
  ...
)
```

## Arguments

- Y:

  Observation matrix (P x N), OR a formula object for Formula Mode. In
  Formula Mode, use `Y1 + Y2 ~ A1 + A2` with `data`, or
  `Y_matrix ~ A_matrix` for direct matrix evaluation. Supports dot
  notation (`. ~ A1 + A2`) when `data` is supplied.

- A:

  Covariate matrix. Default is `NULL` (no covariates). For large \\N\\,
  `A` may instead be a Gram object of class `"nmfkc.gram"` from
  [`nmfkc.kernel.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.gram.md)
  (block-wise Nystr\\om kernel covariates), which holds only \\S =
  AA^\top\\ (\\D \times D\\) and \\G_0 = YA^\top\\ (\\P \times D\\) so
  that the \\D \times N\\ matrix never has to exist in memory. The
  Euclidean multiplicative updates are unchanged (they only use those
  two products) and the fit equals the explicit-matrix fit up to
  floating-point summation order. Gram input is restricted to
  `method = "EU"` without `Y.weights` or `NA` in `Y`; `B`, `XB` and the
  fit statistics are rebuilt block by block after the fit; and the
  resampling helpers
  ([`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md),
  [`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md),
  [`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md))
  do not accept it. Gram objects built from signed features
  ([`nmfkc.signed.rff.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.gram.md))
  are refused here; they belong to
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md).
  Ignored when `Y` is a formula.

- rank:

  Integer. The rank of the basis matrix \\X\\ (Q). Preferred over `Q`.

- data:

  Optional. A data frame from which variables in the formula should be
  taken.

- epsilon:

  Positive convergence tolerance. From the 10th iteration the fit stops
  when \\\|f_i - f\_{i-1}\| / \max(\|f_i\|, 1) \le\\ `epsilon`, where
  \\f_i\\ is the objective at iteration \\i\\. The rule is relative
  while the objective exceeds 1 and *absolute* below it, so on data of
  small magnitude (an objective well below 1) the default `1e-4` can
  stop early, far from the minimum. Rescale `Y` or lower `epsilon`
  (e.g.\\ `1e-8`) in that case. Coefficients that sit on the
  non-negativity boundary also need a tight `epsilon`; see `C.init`
  under `...`.

- maxit:

  Maximum number of iterations.

- verbose:

  Logical. If `TRUE` (default), prints matrix dimensions and elapsed
  time.

- ...:

  Additional arguments passed for fine-tuning regularization,
  initialization, constraints, and output control. This includes the
  backward-compatible arguments `Q` and `method`.

  - `Y.weights`: Optional weight matrix (P x N) or vector (length N)
    with non-negative entries, analogous to the `weights` argument of
    [`lm`](https://rdrr.io/r/stats/lm.html). When supplied, the
    objective becomes \\\sum W\_{ij} \\ (Y\_{ij} - (XB)\_{ij})^2\\
    (i.e.\\ **linear** in \\W\\;
    [`lm()`](https://rdrr.io/r/stats/lm.html)-style weighted least
    squares). Logical matrices (`TRUE` / `FALSE`) are also accepted and
    coerced to 1 / 0. The primary use case is missing-value masking for
    ECV / CV, where \\W\_{ij} \in \\0, 1\\\\ (`FALSE` / `TRUE`)
    indicates held-out vs.\\ used elements; real-valued weights for
    observation-level importance weighting are also supported. Default
    `NULL`: if `Y` contains `NA` a binary mask is auto-constructed (0
    for `NA`, 1 elsewhere); otherwise no weighting.

  - `X.L2.ortho`: Nonnegative penalty parameter for the orthogonality of
    \\X\\ (default: 0). It minimizes the off-diagonal elements of the
    Gram matrix \\X^\top X\\, reducing the correlation between basis
    vectors (conceptually minimizing \\\\ X^\top X -
    \mathrm{diag}(X^\top X) \\\_F^2\\). (Formerly `lambda.ortho`). The
    penalty \\(\lambda/2)\sum\_{q \ne r} (x_q^\top x_r)^2\\ is added to
    the objective as it stands, so how strongly a given \\\lambda\\ acts
    depends on the scale of `Y`; with column-sum-normalized \\X\\ each
    term is also of order \\1/P^2\\, so small values such as \\\lambda
    \le 0.1\\ usually do nothing. To state \\\lambda\\ relative to the
    fit, pass `X.L2.ortho = lambda0 * sum(Y^2)`.

  - `X.L2.smooth`: Nonnegative penalty parameter for row-smoothness of
    \\X\\ (default: 0). Adds \\\lambda\\\mathrm{tr}(X^\top L X)\\ with
    \\L\\ the path-graph Laplacian over the \\P\\ rows, i.e.\\ it
    penalizes squared differences between adjacent rows
    \\\sum_q\sum\_{j\ge 2}(x\_{jq}-x\_{j-1,q})^2\\, yielding
    gently-varying (smooth) bases. Integrates with the multiplicative
    updates (non-negativity and monotone descent preserved). Useful when
    the rows of \\Y\\ have a natural order (e.g.\\ time points).

  - `C.L1`: Nonnegative penalty parameter for L1 regularization on \\C\\
    (default: 0). Promotes **sparsity in the parameter matrix**
    \\\Theta\\ (variable selection over basis-covariate links).
    (Formerly `lambda`).

  - `Q`: Backward-compatible name for the rank of the basis matrix (Q).

  - `method`: Objective function: Euclidean distance `"EU"` (default) or
    Kullback–Leibler divergence `"KL"`.

  - `X.restriction`: Constraint for columns of \\X\\. Options:
    `"colSums"` (default), `"colSqSums"`, `"totalSum"`, `"none"`, or
    `"fixed"`. `"none"` applies no normalization to \\X\\ after each
    update, allowing it to absorb the scale freely.

  - `X.init`: Method for initializing the basis matrix \\X\\. Options:
    `"kmeans"` (default), `"kmeansar"`, `"kmeans++"`, `"runif"`,
    `"nndsvd"`, or a user-specified matrix. `"kmeansar"` applies
    \\k\\-means initialization and then fills zero entries with
    `Uniform(0, mean(Y)/100)`, analogous to NNDSVDar. `"kmeans++"` seeds
    the \\k\\-means centres by \\D^2\\ weighting (Arthur &
    Vassilvitskii, 2007) before Lloyd refinement, giving a more careful,
    stable initialization (`nstart` is not used in this case). If the
    initialization produces a zero basis column – typically `"kmeans"`
    or `"kmeans++"` on data with two or more all-zero observation
    columns, which k-means may gather into a cluster whose centre is the
    zero vector – and `X.restriction` cannot normalize it, the
    initialization is rerun on the non-zero columns of `Y` and a message
    says so. The fit itself still uses every column. A user-supplied
    matrix with a zero column stops with an error instead. `"spa"`
    starts from the anchor rows that the successive projection algorithm
    (SPA; Gillis & Vavasis, 2014) finds in `Y`: \\B\\ is started at
    those rows of `Y` and \\X\\ at the non-negative least-squares fit to
    them. The anchor rows' zeros are filled with small positive values,
    so this is only a starting value and does not make the factorization
    unique; `X.anchor` does. Choosing the rows and the least-squares fit
    involve no randomness; the filled values are drawn as for
    `"kmeansar"`, under `seed`.

  - `X.anchor`: Anchor rows of \\X\\, imposed as constraints (default
    `NULL`: none). An anchor row carries one basis alone, \\X\_{j_q r} =
    0\\ for \\r \ne q\\. The zeros are kept for the whole fit, and they
    make the factorization unique up to the column scale that
    `X.restriction` fixes: another solution \\XR\\ with the same zeros
    needs a diagonal \\R\\. With covariates this identifies \\C\\
    (\\\Theta\\) too, provided `rank <= nrow(A)` and `A` has full row
    rank; a larger `rank` is refused, because \\B = CA\\ then has rank
    below `rank`. `"spa"` chooses the rows by SPA, which recovers them
    when `Y` is (nearly) separable – each basis has a row of `Y` that it
    alone loads – and \\B\\ has full rank. SPA returns `rank` rows even
    when `Y` has no such rows: compare the objective with a fit from
    `X.init = "spa"`, which starts from the same rows without the
    constraint, to see what the assumption costs. Alternatively give the
    rows, one per basis, as indices or row names; basis \\q\\ is then
    the one anchored at the \\q\\-th row given, and the bases are not
    reordered after the fit. The start is built from the anchors, so
    `X.init` is not used (a message says so if it was given). The rows
    used are returned as `X.anchor`.

  - `nstart`: Number of random starts for initialization of \\X\\
    (default: 1). Used by `kmeans` (when `X.init = "kmeans"` or
    `"kmeansar"`) and by the multi-start evaluation (when
    `X.init = "runif"`).

  - `seed`: Integer seed for reproducibility (default: 123).

  - `C.init`: Optional numeric matrix giving the initial value of the
    parameter matrix \\C\\ (i.e., \\\Theta\\). If `A` is `NULL`, `C` has
    dimension \\Q \times N\\ (equivalently \\B\\); otherwise, `C` has
    dimension \\Q \times K\\ where \\K = nrow(A)\\. Default initializes
    all entries to 1. The multiplicative updates keep an exact zero at
    zero, so a zero in `C.init` (or in a user-supplied `X.init`) is a
    structural zero that holds throughout the fit. A coefficient whose
    optimum lies on the boundary, by contrast, is approached only slowly
    and may remain visibly positive at the default `epsilon`. To impose
    a zero, put it in the initial value rather than relying on the
    updates to reach it.

  - `Y.symmetric`: **Removed.** Symmetric NMF (\\Y \approx X X^\top\\ or
    \\X C X^\top\\) has moved to the dedicated
    [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
    function (types `"tri"`, `"bi"`, `"signed"`), which uses the correct
    Frobenius bilateral-gradient updates. Passing `Y.symmetric` to
    `nmfkc()` now stops with a message pointing to
    [`nmfkc.net()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md).

  - `prefix`: Prefix for column names of \\X\\ and row names of \\B\\
    (default: "Basis").

  - `print.trace`: Logical. If `TRUE`, prints progress every 10
    iterations (default: `FALSE`).

  - `print.dims`: Deprecated. Use `verbose` instead.

  - `detail`: Level of post-fit criterion computation. `"fast"`
    (**default**) computes everything except the \\O(N^2)\\
    sample-clustering criteria; `"full"` adds `silhouette`, `CPCC` and
    `dist.cor`, which need two \\N\times N\\ distance matrices and a
    cophenetic correlation — at \\N=2000\\ that was 26 times the cost of
    the rest of the call; `"minimal"` returns only information criteria
    and skips \\XB\\. The default changed from `"full"` in 0.8.9:
    nothing consumed those three (see the `criterion` entry under
    Value), so every call was paying for them. For backward
    compatibility, `save.time = TRUE` maps to `"fast"` and
    `save.memory = TRUE` maps to `"minimal"`.

## Value

A list with components:

- call:

  The matched call, as captured by
  [`match.call()`](https://rdrr.io/r/base/match.call.html).

- dims:

  A character string summarizing the matrix dimensions of the model.

- runtime:

  A character string summarizing the computation time.

- X:

  Basis matrix. Column normalization depends on `X.restriction`.

- B:

  Coefficient matrix \\B = C A\\.

- XB:

  Fitted values for \\Y\\.

- C:

  Parameter matrix.

- B.prob:

  Soft-clustering probabilities derived from columns of \\B\\.

- B.cluster:

  Hard-clustering labels (argmax over \\B.prob\\ for each column).

- X.prob:

  Row-wise soft-clustering probabilities derived from \\X\\.

- X.cluster:

  Hard-clustering labels (argmax over \\X.prob\\ for each row).

- X.anchor:

  Present only when `X.anchor` was used: the anchor row of each basis,
  as a named integer vector.

- X.restriction:

  The constraint that was applied to the columns of \\X\\. It decides
  how much of the \\(X,C)\to(XT,T^{-1}C)\\ freedom is pinned down, which
  inference on the latent parameters needs to know.

- A.attr:

  List of attributes of the input covariate matrix `A`, containing
  metadata like lag order and intercept status if created by `nmfkc.ar`
  or `nmfkc.kernel`.

- formula.meta:

  If fitted via Formula Mode, a list with `formula`, `Y_cols`, and
  `A_cols`; otherwise `NULL`.

- objfunc:

  Final objective value.

- objfunc.iter:

  Objective values by iteration. The trace is trimmed to `[10:end]` so
  that the plot is not dominated by the first few steps, so its length
  is 9 short of `iter` whenever the fit ran at least ten iterations; use
  `iter`, not `length(objfunc.iter)`, for the iteration count.

- iter:

  Number of iterations performed.

- maxit:

  The iteration cap that was in force.

- epsilon:

  The convergence tolerance that was in force.

- epsilon.iter:

  Relative change of the objective at the last step — the quantity the
  stopping rule compares with `epsilon`. Measured as \\\|f_i -
  f\_{i-1}\| / \max(\|f_i\|, 1)\\; the floor of 1 means that for an
  objective below 1 the test is an absolute one, which keeps a
  near-exact fit from iterating forever but makes the number depend on
  the scale of \\Y\\.
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  divides by \\\|f\_{i-1}\|\\ instead, so the two are comparable only
  while the objective exceeds 1.

- objfunc.increases:

  Number of steps at which the objective rose, counted on the recorded
  (trimmed) trace. The multiplicative update is monotone by itself, so a
  positive count means something outside it — a restriction imposed by
  projection, or a numerical problem — is pushing the iterate back, in
  which case the fit can oscillate and exhaust `maxit` instead of
  converging.

- converged:

  Logical; whether `epsilon.iter` met `epsilon` before `maxit`. A run
  that merely exhausted `maxit` is `FALSE`.

- r.squared:

  \\R^2 = \mathrm{cor}(Y, XB)^2\\ (Pearson; scale-invariant;
  \\\[0,1\]\\).

- r.squared.uncentered:

  Uncentered \\R^2 = 1 - \\Y - XB\\\_F^2 / \\Y\\\_F^2\\ (baseline = zero
  matrix; natural for non-negative factorizations without an intercept).

- r.squared.centered:

  Row-mean centered \\R^2 = 1 - \\Y - XB\\\_F^2 / \\Y - \bar
  Y\_{p\cdot}\\\_F^2\\, the multivariate regression \\R^2\\.

- method:

  Character string indicating the optimization method used (`"EU"` or
  `"KL"`).

- n.missing:

  Number of missing (or zero-weighted) elements in \\Y\\.

- n.total:

  Total number of elements in \\Y\\.

- rank:

  The rank \\Q\\ used in the factorization.

- sigma:

  The residual standard error, representing the typical deviation of the
  observed values \\Y\\ from the fitted values \\X B\\.

- mae:

  Mean Absolute Error between \\Y\\ and \\X B\\.

- criterion:

  A list with `effective.rank`, `effective.rank.index`, and the
  separability of the basis: `separability` (one value; 1 = every basis
  has a row of its own), `anchor.purity` (per basis) and `anchor.row`
  (where each basis is purest). See the Details of
  [`summary.nmfkc`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.md)
  for their definition and reading.

  `effective.rank.index` is `effective.rank` rescaled onto \\\[0,1\]\\
  by the broken-stick correction \\(\hat r_Q-E_Q)/(Q-E_Q)\\ with
  \\E_Q=\exp(H_Q-1)\\, the effective rank a random split of the variance
  would produce; it is the quantity
  [`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)
  plots, and `summary` prints it in place of the crispness. Read it as
  **how evenly the across-sample coefficient variance is shared**, which
  is not the same as how useful the factors are: a factor whose
  coefficient is large but constant carries no variance and pulls the
  index down (correctly — the raw value can go negative and is clamped
  at 0, meaning the variance is spread over fewer factors than chance
  would give), while two duplicated factors split the variance evenly
  and score near 1. `NA` at \\Q=1\\, where \\E_1=Q\\.

  With `detail = "full"` the list additionally carries the
  **sample-clustering** criteria `silhouette` (mean silhouette width of
  the hard clustering, computed in the original data space `dist(t(Y))`
  with the per-sample labels), `CPCC` (cophenetic correlation of a
  hierarchical clustering of the coefficient distances `dist(t(B))`) and
  `dist.cor` (correlation between original-data and coefficient
  distances). They are **absent by default** because they cost
  \\O(N^2)\\ and nothing consumes them: they are no longer part of rank
  selection,
  [`summary.nmfkc`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.md)
  does not print them, and
  [`nmf.cluster.criteria`](https://ksatohds.github.io/nmfkc/reference/nmf.cluster.criteria.md)
  — the right entry point when several ranks are to be compared —
  recomputes them from the fits it is given.

  `effective.rank` is the **effective rank**: \\\exp\\ of the Shannon
  entropy of the explained-variance distribution \\p_k =
  \mathrm{var}(B\_{k\cdot}) / \sum_j \mathrm{var}(B\_{j\cdot})\\. By the
  trace identity \\\sum_k \mathrm{var}(B\_{k\cdot}) =
  \mathrm{tr}(\mathrm{Cov}(B))\\, \\p_k\\ is the exact fraction of the
  total coefficient variance carried by factor \\k\\, so the entropy
  measures how that variance is spread across factors. It ranges in
  \\\[1, Q\]\\ (1 when one factor carries all the variance, \\Q\\ when
  all contribute equally) and counts the number of latent factors that
  actively shape across-sample variation. This is the PCA-style
  explained-variance / effective-dimensionality measure and reuses the
  \\\exp(\mathrm{entropy})\\ functional form of Roy & Vetterli (2007).

## References

Satoh, K. (2026). Applying Non-negative Matrix Factorization with
Covariates to the Longitudinal Data as Growth Curve Model. *American
Journal of Mathematical and Management Sciences*. In press.
arXiv:2403.05359. <https://arxiv.org/abs/2403.05359>

Satoh, K. (2025). Applying non-negative matrix factorization with
covariates to multivariate time series data as a vector autoregression
model. *Japanese Journal of Statistics and Data Science*.
arXiv:2501.17446.
[doi:10.1007/s42081-025-00314-0](https://doi.org/10.1007/s42081-025-00314-0)

Satoh, K. (2026). Applying non-negative matrix factorization with
covariates to label matrix for classification. *Japanese Journal of
Statistics and Data Science*.
[doi:10.1007/s42081-026-00349-x](https://doi.org/10.1007/s42081-026-00349-x)

Ding, C., Li, T., Peng, W., & Park, H. (2006). Orthogonal Nonnegative
Matrix Tri-Factorizations for Clustering. In *Proc. 12th ACM SIGKDD*
(pp. 126–135).
[doi:10.1145/1150402.1150420](https://doi.org/10.1145/1150402.1150420)

Roy, O., & Vetterli, M. (2007). The effective rank: A measure of
effective dimensionality. In *15th European Signal Processing Conference
(EUSIPCO)* (pp. 606–610).

Gillis, N., & Vavasis, S. A. (2014). Fast and robust recursive
algorithms for separable nonnegative matrix factorization. *IEEE
Transactions on Pattern Analysis and Machine Intelligence*, 36(4),
698–714.
[doi:10.1109/TPAMI.2013.226](https://doi.org/10.1109/TPAMI.2013.226)

## See also

[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
(standard errors and p-values for \\C\\),
[`summary.nmfkc`](https://ksatohds.github.io/nmfkc/reference/summary.nmfkc.md),
[`plot.nmfkc`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.md),
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md),
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md),
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md),
[`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md),
[`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md),
[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md),
[`predict.nmfkc`](https://ksatohds.github.io/nmfkc/reference/predict.nmfkc.md)

## Examples

``` r
# Example 1. Matrix Mode (Existing)
X <- cbind(c(1,0,1),c(0,1,0))
B <- cbind(c(1,0),c(0,1),c(1,1))
Y <- X %*% B
rownames(Y) <- paste0("P",1:nrow(Y))
colnames(Y) <- paste0("N",1:ncol(Y))
print(X); print(B); print(Y)
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> [3,]    1    0
#>      [,1] [,2] [,3]
#> [1,]    1    0    1
#> [2,]    0    1    1
#>    N1 N2 N3
#> P1  1  0  1
#> P2  0  1  1
#> P3  1  0  1
res <- nmfkc(Y,rank=2,epsilon=1e-6)
#> Y(3,3)~X(3,2)B(2,3)...
#> 0sec
res$X
#>    Basis1      Basis2
#> P1      0 0.498047869
#> P2      1 0.003904261
#> P3      0 0.498047869
res$B
#>              N1           N2        N3
#> Basis1 0.000000 0.9999995176 0.9920988
#> Basis2 2.007838 0.0001206861 2.0079012

# Example 2. Formula Mode
set.seed(1)
dummy_data <- data.frame(Y1=rpois(10,5), Y2=rpois(10,10),
                         A1=abs(rnorm(10,5)), A2=abs(rnorm(10,3)))
res_f <- nmfkc(Y1 + Y2 ~ A1 + A2, data=dummy_data, rank=2)
#> A created using columns: A1, A2
#> Y(2,10)~X(2,2)C(2,2)A(2,10)=XB(2,10)...
#> 0sec

# Example 3. Anchors: in Example 1, rows P1 and P2 of X carry one basis
# each.  X.anchor imposes such rows (here found by SPA), which makes the
# factorization unique; X.init = "spa" only starts from them.
res_a <- nmfkc(Y, rank = 2, X.anchor = "spa", epsilon = 1e-6)
#> Y(3,3)~X(3,2)B(2,3)...
#> 
#> nmfkc(): anchors chosen by SPA at rows P1, P2 ($X.anchor records the basis each one anchors); X keeps zeros elsewhere in these rows.
#> 0sec
res_a$X.anchor
#> Basis1 Basis2 
#>      1      2 

# For symmetric NMF (Y approximated by X X^T or X C X^T),
# use \code{\link{nmfkc.net}()} instead.
```

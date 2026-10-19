# NMF-FFB: non-negative matrix factorization with latent feedback

Fits the NMF-FFB model \$\$ Y_1 = X B + E, \qquad B = \Theta_1 Y_1 +
\Theta_2 Y_2 + U, \$\$ in which a non-negative basis \\X\\ (\\P_1 \times
Q\\) generates the endogenous block \\Y_1\\ from latent scores \\B\\
that are driven by the exogenous block \\Y_2\\ (feed-forward,
\\\Theta_2\\) and fed back by \\Y_1\\ itself (latent feedback,
\\\Theta_1\\). The function returns the estimated basis, the structural
coefficient matrices and the implied equilibrium (input-output) mapping.

At equilibrium, the model can be written as \$\$ Y_1 \approx (I - X
\Theta_1)^{-1} X \Theta_2 Y_2 \equiv M\_{\mathrm{model}} Y_2, \$\$ where
\\M\_{\mathrm{model}} = (I - X \Theta_1)^{-1} X \Theta_2\\ is a
Leontief-type cumulative-effect operator in latent space.

Internally, the latent feedback and exogenous loading matrices are
stored as `C1` and `C2`, corresponding to \\\Theta_1\\ and \\\Theta_2\\,
respectively.

Two estimators are available through `method`:

- `"fiml"` (default):

  A two-stage likelihood-based estimator. **Stage 1** estimates the
  basis with the feed-forward fit
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)`(Y1, A = Y2)`
  (`X.init`, `X.L2.ortho`, `epsilon`, `maxit`, `seed` are forwarded),
  unless a basis is supplied through `X`. Columns are normalized to unit
  sum. **Stage 2** holds \\\hat X\\ fixed and fits the Gaussian working
  model \\U \sim N(0, \Phi)\\, \\E \sim N(0, \mathrm{diag}(\psi))\\ by
  full-information maximum likelihood (L-BFGS-B with an analytic
  gradient) under the non-negativity of \\\Theta_1, \Theta_2\\ and an
  exclusion restriction on \\\Theta_1\\. It fits (a) the feed-forward
  null \\\Theta_1 = 0\\ (a non-negative MIMIC factor model with
  correlated factors), (b) the unpenalized feedback model on the
  admitted entries, and (c) an L1 path over `C1.L1.path`. The penalized
  problem is non-convex, so every point of the path is fitted from
  several starting points (`starts`); every distinct support proposed by
  any (penalty, start) pair is re-estimated without penalty (itself from
  three starts, keeping the best log-likelihood) and BIC is minimized
  over all distinct candidates together with the null and the
  unpenalized model. The support with the smallest BIC is the *selected*
  model, reported in `C1`, `C2`, `Phi`, `psi`. `C1.L1` and `C2.L1` are
  not used by this estimator. The likelihood-ratio statistics against
  the null are returned without p-values: \\\Theta_1 \ge 0\\ puts the
  null on the boundary of the parameter space and the BIC refit is a
  post-selection statistic, so a chi-square reference is not valid.
  Their calibration by parametric bootstrap is done by
  [`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md).

- `"mu"`:

  The legacy estimator: joint multiplicative updates of \\X, \Theta_1,
  \Theta_2\\ for the structural-form squared error \\\lVert Y_1 -
  X(\Theta_1 Y_1 + \Theta_2 Y_2)\rVert_F^2\\ with an orthogonality
  penalty on \\X\\ and L1 penalties `C1.L1`, `C2.L1`. Kept,
  bit-identical, so that analyses published with it reproduce. Because
  the structural and reduced forms have the same fit once \\X\\ is free,
  this estimator cannot separate \\\Theta_1\\ from \\\Theta_2\\; prefer
  `"fiml"`.

## Usage

``` r
nmf.ffb(
  Y1,
  Y2,
  rank = NULL,
  X.init = "nndsvd",
  X.L2.ortho = 100,
  C1.L1 = 1,
  C2.L1 = 0.1,
  epsilon = 1e-06,
  maxit = 5000,
  seed = 123,
  ...,
  method = c("fiml", "mu"),
  X = NULL,
  C1.restriction = c("union", "none"),
  C1.restriction.threshold = 0.05,
  Phi.restriction = c("full", "diag"),
  C1.L1.path = NULL,
  select = c("BIC", "none")
)
```

## Arguments

- Y1:

  A non-negative numeric matrix of endogenous variables with **rows =
  variables (P1), columns = samples (N)**.

- Y2:

  A non-negative numeric matrix of exogenous variables with **rows =
  variables (P2), columns = samples (N)**. Must satisfy
  `ncol(Y1) == ncol(Y2)`.

- rank:

  Integer; number of latent factors \\Q\\. If `NULL`, \\Q\\ is taken
  from a hidden argument in `...` or defaults to `nrow(Y2)`.

- X.init:

  Initialization strategy for the basis matrix `X` (\\P_1 \times Q\\).
  One of:

  - `"nndsvd"` (default): Non-negative Double SVD with additive
    randomness (NNDSVDar; Boutsidis & Gallopoulos 2008), computed
    internally via `.nndsvdar(Y1, Q)`. Requires \\Q \le \min(P_1, N)\\
    (over-rank case falls back to `"runif"`). Uses a full SVD of
    \\Y_1\\, so for very large \\Y_1\\ consider switching to `"kmeans"`
    to avoid SVD memory / compute cost.

  - `"kmeans"`: k-means on the columns of \\Y_1\\ (samples clustered
    into \\Q\\ groups); the transposed cluster centers become \\X\\.
    Scales well for large \\Y_1\\; this is the default of
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md).

  - `"kmeansar"`: `"kmeans"` followed by filling zero entries of \\X\\
    with \\\mathrm{Uniform}(0, \bar Y_1 / 100)\\ (NNDSVDar-style
    additive randomness to escape trivial stationary points).

  - `"runif"`: Uniform random entries in \\\[0, 1\]\\.

  - `"spa"`: the anchor rows the successive projection algorithm finds
    in \\Y_1\\ (see
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)); a
    starting value only. `method = "fiml"` only.

  - A numeric \\P_1 \times Q\\ matrix supplied by the user; negative
    entries are projected to 0.

  - `NULL`: backward-compatible alias for `"nndsvd"`.

  To impose anchor rows on the stage-1 basis instead, pass `X.anchor`
  (`"spa"` or one row of \\Y_1\\ per basis) through `...`, as in
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md): stage
  2 holds the basis fixed, so the zeros stay. It needs
  `rank <= nrow(Y2)`, is not available with `method = "mu"` or a
  supplied `X`, and the rows used are returned as `X.anchor`. In all
  cases the result is column-normalized to `colSums(X) = 1` before
  iteration. The menu mirrors
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)'s
  `X.init` option for consistency across the package.

- X.L2.ortho:

  L2 orthogonality penalty for `X`. This controls the penalty term
  \\\lambda_X \lVert X^\top X - \mathrm{diag}(X^\top X) \rVert_F^2\\.
  Default: `100`.

- C1.L1:

  L1 sparsity penalty for `C1` (i.e., \\\Theta_1\\). Default: `1.0`.
  **Scale note**: this function adds `C1.L1` to the multiplicative
  denominator, whereas
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md) and
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  add `C.L1 / 2`. The same nominal value is therefore twice as strong
  here. The difference is retained so that published `nmf.ffb` fits
  reproduce; halve the value to match the other models.

- C2.L1:

  L1 sparsity penalty for `C2` (i.e., \\\Theta_2\\). Default: `0.1`.
  Same scale note as `C1.L1`.

- epsilon:

  Relative convergence threshold for the objective function. Iterations
  stop when the relative change in **reconstruction loss** falls below
  this value. Default: `1e-6`. **Note**: the test is on the unpenalized
  loss (`objfunc`), not on the penalized objective the updates actually
  minimize (`objfunc.penalized`). Every other optimizer in the package
  tests the penalized value, so with a large `X.L2.ortho` or `C*.L1`
  this function can stop while the quantity being optimized is still
  moving. Both traces are returned; compare them if the penalties are
  strong.

- maxit:

  Maximum number of iterations for the multiplicative updates. Default:
  `5000` (matches
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) and
  other MU functions in the package).

- seed:

  Random seed used to initialize `X`, `C1`, and `C2`. Default: `123`.
  For `method = "fiml"` the seed only reaches the stage-1
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) fit;
  the FIML stage is deterministic.

- ...:

  Additional hidden arguments. For `method = "fiml"`: `fiml.maxit`
  (L-BFGS-B iteration cap, default `3000`), `factr` (`optim` tolerance,
  default `1e3`), and `Q` (alias of `rank`). For `method = "mu"` the
  following control the optional feedforward baseline (used both as an
  \\X\\ warm-start and as the reference for `SC.map`, the input-output
  structural fidelity defined in Satoh (2025) §4.SC.map):

  `nmfkc.baseline`

  :   Controls whether a feedforward
      [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)(Y1,
      A = Y2) fit is used as baseline. Possible values:

      - Default (not given) — `nmf.ffb` runs
        [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
        **internally** when `X.init` is a string method (`"nndsvd"`,
        `"kmeans"`, ...) or `NULL`, forwarding `X.init`, `X.L2.ortho`,
        `epsilon`, `maxit`, `seed`. The fitted \\X\\ of the baseline is
        then used as warm-start for the nmf.ffb MU iterations, and
        `SC.map` is computed. This means `nmf.ffb(Y1, Y2, rank = Q)`
        runs end-to-end without a prior `nmfkc` call.

      - `TRUE` — same as above, but force the internal
        [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
        call even when `X.init` is a user-supplied matrix (the matrix is
        overridden).

      - `FALSE` — opt out; no internal call, `SC.map = NA` (pre-v0.6.8
        behavior).

      - An
        [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
        result (list with `$X` and `$C`) — use as the baseline directly
        (no internal call); also adopted as `X.init` when the latter is
        a string / NULL.

  `M.simple`

  :   Optional \\P_1 \times P_2\\ pre-computed baseline mapping. Takes
      precedence over `nmfkc.baseline` for the SC.map calculation but
      does not affect warm-start.

  `Q`

  :   Backward-compat alias for `rank`.

- method:

  `"fiml"` (default) or `"mu"`; see Description.

- X:

  Optional basis for `method = "fiml"`: a \\P_1 \times Q\\ non-negative
  matrix, or an
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) /
  `nmf.ffb` object whose `$X` is used. When supplied, stage 1 is skipped
  and `rank` is taken from `ncol(X)`.

- C1.restriction:

  Exclusion restriction on \\\Theta_1\\ for `method = "fiml"`: `"union"`
  (default), `"none"`, or a \\Q \times P_1\\ 0/1 matrix (1 = free). See
  the section *Exclusion restriction*.

- C1.restriction.threshold:

  Loading threshold for `C1.restriction = "union"`. Default `0.05`.

- Phi.restriction:

  Covariance of the latent disturbance \\U\\ for `method = "fiml"`:
  `"full"` (default; positive definite via Cholesky) or `"diag"`.

- C1.L1.path:

  Numeric vector of L1 penalties on \\\Theta_1\\ defining the path
  (`method = "fiml"`). Default `NULL`, meaning
  `N * c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)`.

- select:

  `"BIC"` (default): the support with the smallest BIC among all
  candidates proposed along the path is the selected model; `"none"`:
  the unpenalized feedback fit is returned as the selected model and the
  path is skipped.

## Value

An object of class `c("nmf.ffb", "nmf")`, a list with components:

- X:

  Estimated basis matrix (\\P_1 \times Q\\).

- C1:

  Estimated latent feedback matrix (\\\Theta_1\\, \\Q \times P_1\\); for
  `"fiml"` the BIC-selected estimate.

- C2:

  Estimated exogenous loading matrix (\\\Theta_2\\, \\Q \times P_2\\).

- XC1:

  Feedback matrix \\X \Theta_1\\.

- XC2:

  Direct-effect matrix \\X \Theta_2\\.

- XC1.radius:

  Spectral radius \\\rho(X \Theta_1)\\.

- XC1.norm1:

  Induced 1-norm \\\lVert X \Theta_1 \rVert\_{1,\mathrm{op}}\\.

- Leontief.inv:

  Leontief-type inverse \\(I - X \Theta_1)^{-1}.\\

- M.model:

  Equilibrium mapping \\M\_{\mathrm{model}} = (I - X \Theta_1)^{-1} X
  \Theta_2\\.

- amplification:

  Latent amplification factor \\\lVert M\_{\mathrm{model}}
  \rVert\_{1,\mathrm{op}} / \bigl\lVert X \Theta_2
  \bigr\rVert\_{1,\mathrm{op}}\\; meaningful only when `XC1.radius > 0`
  (it equals 1 when no feedback is selected).

- amplification.bound:

  Geometric-series upper bound \\1 / (1 - \lVert X \Theta_1
  \rVert\_{1,\mathrm{op}})\\ if \\\lVert X \Theta_1
  \rVert\_{1,\mathrm{op}} \< 1\\, otherwise `Inf`.

- rank:

  Effective latent dimension used in the fit.

- SC.cov:

  Correlation between sample and model-implied covariance (flattened) of
  \\Y_1\\. See *second-moment fidelity* in Satoh (2025).

- SC.map:

  Correlation between the equilibrium operator \\M\_{\mathrm{model}}\\
  and a feedforward baseline mapping \\M\_{\mathrm{simple}} = X_0
  \Theta_0\\, computed only when the baseline is supplied via `M.simple`
  or `nmfkc.baseline` in `...`; otherwise `NA`. See *input-output
  structural fidelity* in Satoh (2025).

- mae:

  Mean absolute error between \\Y_1\\ and its equilibrium prediction
  \\\hat Y_1 = M\_{\mathrm{model}} Y_2\\.

- objfunc:

  Vector of reconstruction losses per iteration (`"mu"`); `NULL` for
  `"fiml"`.

- objfunc.penalized:

  Vector of penalized objective values per iteration (`"mu"`); `NULL`
  for `"fiml"`.

- iter, maxit, epsilon, converged:

  Convergence bookkeeping. For `"fiml"`: the number of objective
  evaluations of the selected fit, the L-BFGS-B cap (`fiml.maxit`), the
  stage-1 tolerance, and `optim()$convergence == 0`.

- method:

  `"fiml"` or `"mu"`.

The following are present for `method = "fiml"` only (`SC.cov` and
`SC.map` are then `NULL`):

- Phi, psi, loglik, npar:

  Latent-disturbance covariance (\\Q \times Q\\), unique variances
  (length \\P_1\\), log-likelihood and parameter count of the selected
  model.

- null:

  The feed-forward null: list `C2, Phi, psi, loglik, npar, M.model`.

- full:

  The unpenalized feedback fit: list
  `C1, C2, Phi, psi, loglik, npar, XC1.radius`.

- path:

  Data frame with one row per (`C1.L1.path`, start) of the L1 path
  (`C1.L1.path = 0` is the unpenalized fit re-estimated on its non-zero
  entries, `Inf` the null):
  `C1.L1.path, start, support_id, nnz, rho, loglik, BIC, MAE` (each
  after re-estimation on the proposed support), `pen.value` (the
  penalized objective reached by that start, smaller is better) and
  `duplicate` (`TRUE` when the same support was already proposed by an
  earlier row).

- candidates, supports, support.selected:

  One row per distinct support
  (`support_id, nnz, rho, loglik, BIC, MAE, lambda1.first, start.first, selected`;
  the null and the model with every admitted entry free are always
  candidates, so `candidates` can hold a support that appears in no row
  of `path`); the supports themselves (list of logical \\Q \times P_1\\
  matrices indexed by `support_id`); and the id of the selected one.

- C1.free, C1.L1.path, C1.L1.selected, support:

  The free-entry matrix used (\\Q \times P_1\\ 0/1), the path, the
  smallest penalty at which the selected support was proposed (`0` for
  the unpenalized model, `Inf` for the null), the selected support
  (logical \\Q \times P_1\\) and the starts used.

- LR, LR.df:

  Likelihood-ratio statistics
  `c(full = 2(l_full - l_null), selected = 2(l_sel - l_null))` and the
  naive degrees of freedom
  `c(full = sum(C1.restriction), selected = nnz)` (also stored as
  `attr(LR, "df")`). No p-value is attached; see
  [`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md).

- BIC, AIC:

  Named vectors `c(null, full, selected)`.

- call:

  The matched call, from which
  [`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)
  inherits the design.

## Exclusion restriction

Feedback is identified only through exclusion restrictions: an outcome
may not feed back into a factor on which it has more than a negligible
loading, because such an entry is nearly equivalent to a change of the
outcome's loading and uniqueness. With `C1.restriction = "union"`
(default) entry \\(q, i)\\ of \\\Theta_1\\ is excluded if \\q =
\arg\max\_{q'} X\_{i q'}\\ (the dominant factor) or \\X\_{iq} \ge\\
`C1.restriction.threshold`. Both halves are needed: blocking the
dominant factor alone would leave an outcome with a substantial second
loading free to feed that factor, and blocking only the factors above
the threshold would leave an outcome whose largest loading is below the
threshold free to feed its own. `"none"` frees every entry (not
recommended: the model is then identified only through the
non-negativity and the covariance structure). A user-supplied \\Q \times
P_1\\ 0/1 matrix is used as given.

The restriction is derived from the estimated basis and therefore from
the same \\Y_1\\ that is subsequently tested; see the *Calibration*
section of
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)
for what this implies.

## Lifecycle

`method = "fiml"` became the default in version 0.9.8, as did
`C1.restriction = "union"` (earlier fiml fits blocked the dominant
factor only), and
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)
gained the `calibration` argument. `method = "mu"` is the legacy
estimator, kept for the reproducibility of published analyses; it will
be deprecated in a later release.

## References

Satoh, K. (2025). Applying non-negative matrix factorization with
covariates to structural equation modeling for blind input-output
analysis. arXiv:2512.18250. <https://arxiv.org/abs/2512.18250>

## See also

[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md),
[`nmf.ffb.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.cv.md),
[`nmf.ffb.split`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.split.md),
[`nmf.ffb.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.DOT.md),
[`summary.nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/summary.nmf.ffb.md)

## Examples

``` r
# Simple NMF-FFB with iris data (non-negative)
Y <- t(iris[, -5])
Y1 <- Y[1:2, ]  # Sepal
Y2 <- Y[3:4, ]  # Petal
result <- nmf.ffb(Y1, Y2, rank = 2)
result$LR          # feedback vs feed-forward null (no p-value here)
#>         full     selected 
#> 4.837375e-10 0.000000e+00 
#> attr(,"df")
#>     full selected 
#>        2        0 
result$BIC
#>     null     full selected 
#> 737.4257 747.4470 737.4257 
result$path
#>    C1.L1 start support_id nnz rho   loglik      BIC      MAE pen.value
#> 1   0.00  full          1   0   0 -346.165 737.4257 1.310584        NA
#> 2   0.30  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 3   0.75  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 4   1.50  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 5   3.00  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 6   7.50  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 7  15.00  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 8  30.00  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 9  75.00  full          1   0   0 -346.165 737.4257 1.310584   346.165
#> 10   Inf  null          1   0   0 -346.165 737.4257 1.310584        NA
#>    duplicate
#> 1      FALSE
#> 2       TRUE
#> 3       TRUE
#> 4       TRUE
#> 5       TRUE
#> 6       TRUE
#> 7       TRUE
#> 8       TRUE
#> 9       TRUE
#> 10      TRUE

# Legacy multiplicative-update estimator
result.mu <- nmf.ffb(Y1, Y2, rank = 2, maxit = 500, method = "mu")
#> Warning: maximum iterations (500) reached...
result.mu$mae
#> [1] 1.692159
```

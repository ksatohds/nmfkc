# Non-negative Matrix Factorization with Random Effects

Estimates the NMF-RE model \$\$Y = X(\Theta A + U) + \mathcal{E}\$\$
where \\Y\\ (\\P \times N\\) is a non-negative observation matrix, \\X\\
(\\P \times Q\\) is a non-negative basis matrix learned from the data,
\\\Theta\\ (\\Q \times K\\) is a non-negative coefficient matrix
capturing systematic covariate effects on latent scores, \\A\\ (\\K
\times N\\) is a covariate matrix, and \\U\\ (\\Q \times N\\) is a
random effects matrix capturing unit-specific deviations in the latent
score space.

NMF-RE can be viewed as a mixed-effects latent-variable model defined on
a reconstruction (mean) structure. The non-negativity constraint on
\\X\\ induces sparse, parts-based loadings, achieving measurement-side
variable selection without an explicit sparsity penalty. Inference on
\\\Theta\\ provides covariate-side variable selection by identifying
which covariates significantly affect which components.

Estimation is an outer–inner ECM: an inner loop (fixed \\\lambda\\) runs
a ridge/BLUP update for \\U\\, a complete-EM semi-NMF update for \\X\\,
and a fixed-effect update for \\\Theta\\; an outer loop runs the EM
M-steps for \\(\sigma^2, \tau^2)\\. The variance components are
estimated from the data; the effective degrees of freedom \\df_U\\ are
reported only as a diagnostic.

`nmfre()` performs **optimization only**. Hypothesis tests and standard
errors for \\\Theta\\ are obtained separately with
[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md)
(sandwich SE + wild bootstrap), mirroring the
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) /
`nmfkc.inference` split.

## Usage

``` r
nmfre(
  Y,
  A = NULL,
  rank = 2,
  C.signed = TRUE,
  epsilon = 1e-05,
  maxit = 5000,
  ...
)
```

## Arguments

- Y:

  Observation matrix (P x N), non-negative.

- A:

  Covariate matrix (K x N). Default is a row of ones (intercept only).

- rank:

  Integer. Rank of the basis matrix \\X\\. Default is 2. For backward
  compatibility, `Q` is accepted via `...`.

- C.signed:

  Logical. Whether the fixed-effect coefficients \\C\\ (\\= \Theta\\ in
  the paper) are sign-free; the single switch that selects the whole
  estimation scheme. `TRUE` (default, recommended, matches the paper)
  treats \\C\\ as real-valued and updates it by exact least squares,
  with the basis \\X\\ estimated by the complete-EM semi-NMF step and a
  two-sided test (interior null). `FALSE` constrains \\C \ge 0\\ via a
  multiplicative update, with \\X\\ estimated by the positive-part
  multiplicative update and a one-sided test (boundary null). The basis
  \\X\\ is always non-negative. A character value (`"signed"` /
  `"nonneg"`) is also accepted for backward compatibility.

- epsilon:

  Convergence tolerance for relative change in objective (default 1e-5).

- maxit:

  Maximum number of iterations. Default `5000` (matches
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) and the
  other MU functions in the package). When the cap is hit without
  meeting the relative- tolerance criterion, a
  `"maximum iterations (...) reached..."` warning is emitted so users
  notice unconverged fits.

- ...:

  Additional arguments for initialization, variance control, dfU
  control, optimization, and inference settings.

  - `X.init`: Initial basis matrix (P x Q), or `NULL`. When `NULL`,
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) is
    called internally to generate initial values. A method name
    (`"kmeans"`, `"nndsvd"`, `"spa"`, ...) is forwarded to that
    [`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
    initialization.

  - `X.anchor`: anchor rows of \\X\\, imposed as in
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md):
    `NULL` (default), `"spa"`, or one row of `Y` per basis. The zeros
    are kept for the whole fit, including through the floor this
    function applies to \\X\\. There is no `rank <= nrow(A)` limit here:
    the random effects give the scores full rank. The rows used are
    returned as `X.anchor`.

  - `C.init`: Initial coefficient matrix (Q x K), or `NULL`. When
    `NULL`,
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) is
    called internally to generate initial values.

  - `U.init`: Initial random effects matrix (Q x N), or `NULL` (all
    zeros).

  - `prefix`: Prefix for basis names (default `"Basis"`).

  - `sigma2`: Initial residual variance (default 1).

  - `sigma2.update`: Logical. Update \\\sigma^2\\ during iterations
    (default `TRUE`).

  - `tau2`: Initial random effect variance (default 1).

  - `tau2.update`: Logical. Update \\\tau^2\\ by moment matching
    (default `TRUE`).

  - `x.postvar`: Logical. Include the posterior-variance term in the
    semi-NMF \\X\\-step (default `TRUE`; advanced). Applies only to the
    sign-free / semi-NMF path (`C.signed = TRUE`).

  - `X.L2.smooth`: Non-negative row-smoothness penalty on the basis
    \\X\\ (default 0), adding \\(\lambda/2)\\\mathrm{tr}(X^\top L X)\\
    with \\L\\ the path-graph Laplacian over rows (squared adjacent-row
    differences). Well suited to ordered rows (e.g.\\ longitudinal /
    spatial bases). Injected into the \\X\\-step for both
    `C.signed = TRUE/FALSE`; leaves the random-effect variance updates
    unchanged.

  - `X.L2.ortho`: Non-negative column-orthogonality penalty on \\X\\
    (default 0), penalizing \\(\lambda/2)\lVert\mathrm{offdiag}(X^\top
    X)\rVert^2\\.

  - `C.L2`: Non-negative ridge penalty on \\\Theta = C\\ (default 0),
    adding \\\lambda\lVert C\rVert^2\\. A Gaussian prior on the fixed
    effects; for `C.signed = TRUE` the \\C\\-step stays a closed-form
    Sylvester ridge-LS solve, for `C.signed = FALSE` it is added to the
    multiplicative-update denominator.

  - `dfU.control`: Deprecated and inert. The algorithm imposes no cap on
    \\df_U\\ (`"off"`, the only behaviour); \\df_U\\ is reported as a
    diagnostic only.

  - `print.trace`: Logical. If `TRUE`, print progress every 100
    iterations (default `FALSE`).

  - `seed`: Integer seed for reproducibility (default 1).

  - `nstart`: Number of random restarts for the
    [`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
    initialisation step (passed to the k-means initialiser). Default `1`
    (single start; historical behaviour). A larger value (e.g.\\ 10-20)
    gives a more stable initialisation.

  - `inner.maxit`, `outer.maxit`: Maximum inner (fixed-\\\lambda\\
    block-coordinate) and outer (EM variance) iterations (defaults
    `10000` and `500`).

  - `epsilon.outer`: Convergence tolerance for the outer EM loop on
    \\\lambda\\ (default `1e-6`).

## Value

A list of class `"nmfre"` with components. The model is \\Y = X(\Theta
A + U) + \mathcal{E}\\.

**Core matrices**

- `X`:

  Basis matrix \\X\\ (\\P \times Q\\), columns normalized to sum to 1.

- `X.prob`:

  Row-wise soft-clustering probabilities from the non-negative \\X\\
  (each row normalized to sum to 1), as in
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md).

- `X.cluster`:

  Hard-clustering label for each row of \\X\\ (argmax over `X.prob`).

- `C`:

  Coefficient matrix \\\Theta\\ (\\Q \times K\\).

- `U`:

  Random effects matrix \\U\\ (\\Q \times N\\).

**Variance components**

- `sigma2`:

  Residual variance \\\hat{\sigma}^2\\.

- `tau2`:

  Random effect variance \\\hat{\tau}^2\\.

- `lambda`:

  Ridge penalty \\\lambda = \sigma^2 / \tau^2\\.

**Convergence diagnostics**

- `converged`:

  Logical. Whether the algorithm converged.

- `stop.reason`:

  Character string describing why iteration stopped.

- `iter`:

  Number of iterations performed.

- `maxit`:

  Maximum iterations setting used.

- `epsilon`:

  Convergence tolerance used.

- `objfunc`:

  Final objective function value \\\\Y - X(\Theta A + U)\\^2 + \lambda
  \\U\\^2\\.

- `rel.change.final`:

  Final relative change in objective.

- `objfunc.iter`:

  Numeric vector of the fixed-\\\lambda\\ penalized objective \\\\Y -
  X(\Theta A + U)\\^2 + \lambda\\U\\^2\\ per iteration. Monotone within
  an inner loop but *not* across outer iterations (the penalty jumps
  when \\\lambda\\ is updated).

- `rss.trace`:

  Numeric vector of \\\\Y - X(\Theta A + U)\\^2\\ per iteration.

- `nll.trace`:

  Numeric vector of the marginal negative log-likelihood
  \\\ell(X,\Theta,\sigma^2,\tau^2)\\ per iteration (random effects
  integrated out). This is the ECM-monotone quantity and is what
  [`plot.nmfre`](https://ksatohds.github.io/nmfkc/reference/plot.nmfre.md)
  displays.

**Effective degrees of freedom (dfU) diagnostics**

- `dfU`:

  Final effective degrees of freedom \\\mathrm{df}\_U = N \sum_q d_q /
  (d_q + \lambda)\\, where \\d_q\\ are eigenvalues of \\X'X\\.

- `dfU.cap`:

  Upper bound imposed on \\\mathrm{df}\_U\\.

- `dfU.cap.rate`:

  Rate used to compute the cap.

- `lambda.enforced`:

  Final \\\lambda\\ enforced to satisfy the cap.

- `dfU.hit.cap`:

  Logical. Whether the cap was binding.

- `dfU.hit.iter`:

  Iteration at which the cap first bound.

- `dfU.frac`:

  \\\mathrm{df}\_U / (NQ)\\, fraction of maximum df.

- `dfU.cap.frac`:

  \\\mathrm{df}\_U^{\mathrm{cap}} / (NQ)\\.

**Fitted matrices**

- `B`:

  Fixed-effect scores \\\Theta A\\ (\\Q \times N\\).

- `B.prob`:

  Column-normalized probabilities from \\\max(\Theta A, 0)\\.

- `B.blup`:

  BLUP scores \\\Theta A + U\\ (\\Q \times N\\).

- `B.blup.pos`:

  Non-negative BLUP scores \\\max(\Theta A + U, 0)\\ (\\Q \times N\\).

- `B.blup.prob`:

  Column-normalized probabilities from \\\max(\Theta A + U, 0)\\.

- `XB`:

  Fitted values from fixed effects \\X \Theta A\\ (\\P \times N\\).

- `XB.blup`:

  Fitted values including random effects \\X(\Theta A + U)\\ (\\P \times
  N\\).

**Fit statistics**

- `r.squared`:

  Pearson \\\mathrm{cor}(Y, X(\Theta A + U))^2\\ (BLUP prediction).

- `r.squared.uncentered`:

  Uncentered \\1 - \\Y - X(\Theta A + U)\\\_F^2 / \\Y\\\_F^2\\ (BLUP;
  baseline = zero matrix).

- `r.squared.centered`:

  Row-mean centered \\1 - \\Y - X(\Theta A + U)\\\_F^2 / \\Y - \bar
  Y\_{p\cdot}\\\_F^2\\ (BLUP; baseline = per-row mean).

- `r.squared.fixed`:

  Pearson \\\mathrm{cor}(Y, X\Theta A)^2\\ (fixed-only prediction).

- `r.squared.fixed.uncentered`, `r.squared.fixed.centered`:

  Uncentered and centered \\R^2\\ for the fixed-only prediction.

- `ICC`:

  Trace-based Intraclass Correlation Coefficient. In the NMF-RE model,
  the conditional covariance of the \\n\\-th observation column is
  \\\mathrm{Var}(Y_n) = \tau^2 X X^\top + \sigma^2 I_P\\, a \\P \times
  P\\ matrix. Unlike a standard random intercept model where the design
  matrix \\Z\\ is a simple indicator (so the ICC reduces to \\\tau^2 /
  (\sigma^2 + \tau^2)\\), the basis matrix \\X\\ plays the role of \\Z\\
  in a random slopes model, making the variance contribution of \\U\\
  depend on \\X\\. To obtain a scalar summary, we take the trace of each
  component: \$\$\mathrm{ICC} = \frac{\tau^2 \\ \mathrm{tr}(X^\top X)}
  {\tau^2 \\ \mathrm{tr}(X^\top X) + \sigma^2 P}.\$\$ This equals the
  average (over \\P\\ dimensions) proportion of per-column variance
  attributable to the random effects.

**Sign convention**

- `C.signed`:

  Logical. Whether \\C\\ was sign-free (`TRUE`) or non-negative
  (`FALSE`).

Standard errors, z-values, p-values, and confidence intervals for
\\\Theta\\ are **not** computed here; obtain them by passing the fit to
[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md).

## References

Satoh, K. (2026). Wild Bootstrap Inference for Non-Negative Matrix
Factorization with Random Effects. arXiv:2603.01468.
<https://arxiv.org/abs/2603.01468>

## See also

[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md),
[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md),
[`summary.nmfre`](https://ksatohds.github.io/nmfkc/reference/summary.nmfre.md)

## Examples

``` r
# Example 1. cars data
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(intercept = 1, speed = cars$speed)
res <- nmfre(Y, A, rank = 1, maxit = 5000)
summary(res)
#> NMF-RE: Y(1,50) = X(1,1) [C(1,2) A + U(1,50)]
#> Iterations: 43 / 5000 (converged)  epsilon = 1e-05  outer 18 / 500
#> R-squared (cor^2):    0.9328 (XB+blup), 0.6511 (XB)
#> R-squared (uncentered):     0.9799 (XB+blup), 0.9091 (XB)
#> R-squared (centered): 0.9227 (XB+blup), 0.6511 (XB)
#> 
#> Variance components:
#>   sigma2 = 106.9  (residual)
#>   tau2   = 120.2  (random effect)
#>   lambda = 0.8891  (sigma2 / tau2)
#>   ICC    = 0.5294  (tau2*tr(X'X) / (tau2*tr(X'X) + sigma2*P))
#>   dfU    = 26.47
#> 
#> Coefficients (Theta): run nmfre.inference(fit, Y, A) for SE / p-values.

# \donttest{
# Example 2. Orthodont data (nlme)
if (requireNamespace("nlme", quietly = TRUE)) {
  Y <- matrix(nlme::Orthodont$distance, 4, 27)
  male <- ifelse(nlme::Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
  A <- rbind(intercept = 1, male = male)

  # Fit (sign-free Theta by default; variances estimated by EM/ECM)
  res <- nmfre(Y, A, rank = 1)
  summary(res)
}# }
#> NMF-RE: Y(4,27) = X(4,1) [C(1,2) A + U(1,27)]
#> Iterations: 41 / 5000 (converged)  epsilon = 1e-05  outer 20 / 500
#> R-squared (cor^2):    0.8246 (XB+blup), 0.4167 (XB)
#> R-squared (uncentered):     0.9974 (XB+blup), 0.9915 (XB)
#> R-squared (centered): 0.7594 (XB+blup), 0.2134 (XB)
#> 
#> Variance components:
#>   sigma2 = 1.933  (residual)
#>   tau2   = 48.19  (random effect)
#>   lambda = 0.04011  (sigma2 / tau2)
#>   ICC    = 0.6100  (tau2*tr(X'X) / (tau2*tr(X'X) + sigma2*P))
#>   dfU    = 23.28
#> 
#> Coefficients (Theta): run nmfre.inference(fit, Y, A) for SE / p-values.
```

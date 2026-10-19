# Fit NMF-GMM: a Gaussian-mixture latent-class extension of NMF with covariates

This function is **experimental and still under development**. The
interface may change in future versions: argument names, defaults and
the contents of the returned object are not yet stable. Code written
against it today may need adjusting after an update. The rest of the
package does not carry this caveat.

`nmf.gmm` fits the model \$\$\bm b_n\mid(z_n=k)\sim N_Q(C\bm
a_n+\bm\mu_k,\Sigma_k),\qquad \bm y_n = X\bm b_n+\bm\varepsilon_n,\quad
X\ge 0,\$\$ by a generalized EM algorithm: a \\K\\-component Gaussian
mixture is placed on the latent NMF scores \\\bm b_n\\, whose
class-conditional mean combines a covariate regression \\C\bm a_n\\
(\\C=\Theta\\ in the paper, the \\Q\times R\\ coefficient matrix) with a
class shift \\\bm\mu_k\\. The shared basis \\X\\ stays non-negative and
column-normalized. Clustering is model based, through the posterior
responsibilities. \\K=1\\ (tied covariance) reduces exactly to NMF with
random effects
([`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md)).

Like the rest of the package, `nmf.gmm` performs **optimization only**;
use
[`nmf.gmm.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md)
for standard errors and tests of \\C\\, and
[`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md)
to choose \\K\\.

## Usage

``` r
nmf.gmm(Y, A = NULL, rank, K = 1, ...)
```

## Arguments

- Y:

  Data matrix \\Y\\ (P x N).

- A:

  Covariate matrix \\A\\ (R x N). Default `NULL` uses an intercept-only
  \\1\times N\\ matrix (a plain Gaussian mixture on the scores, with no
  covariate adjustment). A row of ones should be present for the
  identifiability constraint \\\sum_k \xi_k\bm\mu_k = 0\\ (see
  `intercept`).

- rank:

  Integer rank \\Q\\ of the basis.

- K:

  Integer number of mixture components. Default 1 (= NMF-RE).

- ...:

  Additional arguments:

  - `cov`: score-covariance structure. `"tied"` (default, shared
    diagonal \\\Sigma_k=\mathrm{diag}(\tau^2_1,\dots,\tau^2_Q)\\, \\Q\\
    variances); `"free"` (per-class diagonal, \\Q\times K\\); or
    `"scalar"` (isotropic \\\Sigma_k=\tau^2 I_Q\\, a single variance —
    the most parsimonious variant; at \\K=1\\ it is the
    [`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md)
    model, and `tau2` is returned as one number).

  - `X.init`: initial basis. A \\P\times Q\\ non-negative matrix, or an
    initialization method name: `"nndsvd"` (default), `"kmeans++"`,
    `"kmeans"`, ... (passed to the shared initializer), or `"nmf"`, the
    basis of an ordinary NMF fitted by
    `nmfkc(Y, rank = rank, seed = seed)` without covariates (requires a
    non-negative `Y`;
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)'s own
    default seed is 123). The EM restarts (`nstart`) all start from this
    one basis and differ only in the seeding of the mixture, so the fit
    can depend on the choice of `X.init`; the basis actually used is
    returned as `X0`. `"spa"` starts from the anchor rows the successive
    projection algorithm finds in `Y` (see
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md);
    requires a non-negative `Y`).

  - `X.anchor`: anchor rows of \\X\\, imposed as in
    [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md):
    `NULL` (default), `"spa"`, or one row of `Y` per basis (requires a
    non-negative `Y`). The EM's multiplicative basis update keeps the
    zeros. There is no `rank <= nrow(A)` limit here: the mixture gives
    the scores full rank.
    [`nmf.gmm.twostage()`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.twostage.md)
    accepts it too, finding the anchors on `Y`. The rows used are
    returned as `X.anchor`.

  - `intercept`: index of the intercept row of `A` (default 1); the
    class-mean average is absorbed into that column of \\C\\.

  - `nstart`: EM restarts (default 1 for `K=1`, else 8).

  - `cores`: run the `nstart` restarts in parallel (default
    `getOption("mc.cores", 1L)`; PSOCK cluster on Windows, forking
    elsewhere). Each restart is an independent self-seeded EM and
    results are combined in order, so the returned fit is identical for
    any `cores`.
    ([`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md)
    parallelizes over `K` instead and runs each inner fit with
    `cores = 1`.)

  - `maxit`, `tol`: outer EM cap / tolerance (2e5, 1e-7). `tol` is the
    relative change of the marginal log-likelihood, the usual EM rule.

  - `ptol`: a second stopping rule, on the relative change of the
    *reported* parameters \\X\\, \\C\\ and \\\bm\mu\\ (default 1e-9;
    `NULL` disables it). The log-likelihood rule alone can fail to fire
    when a nuisance variance drifts to its boundary – \\\tau^2\to 0\\ on
    a sample that carries no random effect – so the objective keeps
    creeping up after the partition, the basis and the covariate effect
    have settled. The returned `stop_by` says which rule fired
    (`"loglik"`, `"parameters"` or `"maxit"`).

  - `vfloor`: lower bound on the variance parameters. The default
    (`NULL`) makes it scale-relative, \\10^{-12}\\\mathrm{mean}(Y^2)\\,
    so the fit is equivariant under \\Y\to cY\\; a fixed bound would
    bind on small-scale data.

  - `fixX`: hold the basis at its initial value instead of updating it
    (default `FALSE`). Used to compare the joint fit with
    [`nmf.gmm.twostage`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.twostage.md)
    on the *same* \\X\\, so that the two routes share the basis. They
    still fit different data: the two-stage route fits the
    reconstituted, shifted residuals.

  - `class.effects`: covariates whose coefficients differ by class
    (default `NULL`: every covariate has one coefficient shared by all
    classes). Row indices of `A`, row names, or — when `A` is a formula
    — term names, so `class.effects = "diet"` frees every indicator row
    of a factor `diet`. The class-\\k\\ mean becomes \\C\bm
    a_n+\Delta_k\bm a\_{n,S}+\bm\mu_k\\ with \\\sum_k\xi_k\Delta_k=0\\,
    so the columns \\S\\ of `C` hold the \\\xi\\-weighted average of the
    class-specific coefficients, which are returned in `C.class`. The
    intercept row cannot be freed: \\\bm\mu_k\\ already makes it
    class-specific. With a continuous covariate, the class-specific
    coefficients are identified when the covariate varies within every
    class. With a categorical covariate, freeing all of its indicator
    rows gives every level its own class means, so the likelihood no
    longer ties the classes of one level to those of another: which
    class of one level corresponds to which class of another is then
    carried by the starting values rather than by the data. Not
    supported by
    [`nmf.gmm.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md)
    or
    [`nmf.gmm.twostage`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.twostage.md).

  - `seed`: RNG seed (default 1). `prefix`: basis-name prefix (default
    `"Basis"`).

  - `data`: a data frame with one row per column of `Y`, required when
    `A` is a formula (see below).

  - `standardize`: when `A` is a formula, center and scale the
    constructed covariate columns (default `TRUE`). Factors are expanded
    to treatment indicators before standardization, so a \\5\\-level
    factor becomes four centered, scaled columns.

  `A` may also be a one-sided **formula** (e.g. `~ size` or `~ diet`),
  evaluated in `data`: the design matrix is built with
  [`model.matrix`](https://rdrr.io/r/stats/model.matrix.html), its
  intercept column is replaced by the package's own intercept row, and
  the remaining columns are standardized by default. The constructed
  numeric \\A\\ is returned in the fit (fields `A`, `A.formula`,
  `A.center`, `A.scale`), so
  [`nmf.gmm.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md)
  works unchanged.

## Value

An object of class `"nmf.gmm"`: a list with `X` (basis), `X0` (the
column-normalized initial basis the EM started from), `C` (\\=\Theta\\,
Q x R), `mu` (class means, Q x K), `tau2`, `sigma2`, `xi` (mixing
proportions), `gamma` (responsibilities, N x K), `cluster` (hard
labels), `loglik`, `BIC`, `ICL`, `n.params`, `entropy`, `Yhat`,
`objfunc(.iter)`, `iter`, `converged`, `stop_by`, `K`, `rank`, `cov`,
`dims` and `runtime`. With `class.effects`, also `class.effects` (the
freed rows of `A`, by name), `class.rows` (their indices) and `C.class`
(a Q x \|S\| x K array of class-specific coefficients; the
class-specific effect in observation space is `X %*% C.class[, , k]`).

## See also

[`nmf.gmm.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.inference.md),
[`nmf.gmm.select`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.select.md),
[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md) (the
\\K=1\\ special case).

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(abs(rnorm(20 * 60)) + 1, 20, 60)
A <- rbind(1, rnorm(60))                 # intercept + one covariate
fit <- nmf.gmm(Y, A, rank = 2, K = 2)
fit$BIC
#> [1] 2398.661
table(fit$cluster)
#> 
#>  1  2 
#> 14 46 
# }
```

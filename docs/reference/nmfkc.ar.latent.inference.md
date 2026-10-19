# Bootstrap inference for the latent VAR of an NMF-VAR model

This function is **experimental**. The interface may change in future
versions: argument names, defaults and the contents of the returned
object are not yet stable.

`nmfkc.ar.latent.inference` attaches a wild bootstrap to the quantities
reported by
[`nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md):
standard errors and intervals for the spectral radius, the long-run
mean, the long-run composition and the entries of the latent transition
matrices \\G_d\\.

## Usage

``` r
nmfkc.ar.latent.inference(object, Y, A, ...)
```

## Arguments

- object:

  An object from
  [`nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md).

- Y, A:

  The observation and covariate matrices the model was fitted on, i.e.
  the `Y` and `A` returned by
  [`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md).

- ...:

  Additional arguments:

  `wild.B`

  :   Number of bootstrap replicates. Default 500.

  `wild.dist`

  :   Multiplier distribution, `"rademacher"` (default) or `"exp"`
      (mean-centred).

  `wild.unit`

  :   `"column"` (default, one multiplier per time point), `"block"`
      (one per block of `max(1, round(n^(1/3)))` consecutive time
      points, `n` the number of time points) or `"element"` (one per
      cell).

  `wild.level`

  :   Confidence level. Default 0.95.

  `wild.seed`

  :   Seed for the bootstrap. Default 123.

  `fit`

  :   The fitted `"nmfkc"` object. Only needed when `object` does not
      carry `X` / `Theta`.

  `truncate`

  :   Logical; clip \\Y^\ast\\ at 0 (default `TRUE`).

  `epsilon`, `maxit`

  :   Convergence control of the re-fits, default `1e-8` and `20000` —
      tight, because a coefficient heading for the non-negativity
      boundary moves slowly and a loose tolerance biases anything that
      is compared with zero.

  `kkt.tol`, `kkt.ratio.min`

  :   Thresholds of the strict complementarity diagnostic (defaults
      `1e-3` and `5`).

  `sep.tol`

  :   Row-normalised loading a column of \\X\\ must reach for that basis
      to count as anchored. Default `0.999`.

  `pin.tol`

  :   A coefficient of \\\Theta\\ at most `pin.tol` times the
      coefficient of the same covariate for the basis being pinned
      counts as zero when the rows of \\\Theta\\ are pinned; so does one
      at most `kkt.tol`. Default `0.01`: a mixing of the bases of up to
      1\\ actually reached is reported –
      [`print()`](https://rdrr.io/r/base/print.html) states it, and
      `identification$pin.margin`, `$pure.margin` and `$X.pin.margin`
      hold the largest leaks that were needed.

  `refit.args`

  :   Named list of extra arguments for the bootstrap re-fits. `method`
      and `X.restriction` are taken from the original fit automatically,
      so the replicates use the same estimator; anything else the
      original call set (`X.init`, penalties, ...) cannot be recovered
      from the fitted object and must be passed here. The function warns
      when it sees such arguments in the recorded call.

  `cores`

  :   Number of workers for the re-fits, default
      `getOption("mc.cores", 1L)` (sequential). The bootstrap
      multipliers are drawn before the loop and each re-fit is
      deterministic given its data, so the result does not depend on
      `cores`.

## Value

An object of class `"nmfkc.ar.latent.inference"`: a list with

- coefficients:

  Data frame of the \\DQ^2\\ entries of the \\G_d\\, with columns `Lag`,
  `Basis` (effect at \\t\\), `Covariate` (cause at \\t-d\\), `Estimate`,
  `BSE`, `CI_low`, `CI_high`, `p_value`. Descriptive unless
  `identified.approx`.

- G, G.se, G.ci.lower, G.ci.upper, G.p.value:

  The same quantities as lists of \\Q\times Q\\ matrices, one per lag.

- spectral.radius:

  The radius of the \\DQ\times DQ\\ latent companion matrix — the same
  quantity
  [`nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md)
  reports, not \\\rho(\sum_d G_d)\\, which differs when \\D\>1\\. With
  `spectral.radius.se`, `spectral.radius.ci.lower` / `.upper` and
  `spectral.radius.boot.mean`.

- prob.nonstationary:

  Bootstrap tail probability \\P^\ast(\rho^\ast\ge1)\\. **Not** a
  p-value for \\H_0:\rho\ge1\\ — see Details.

- complex, complex.frac:

  Whether the point estimate has complex roots, and the fraction of
  replicates that do.

- mu.y:

  With `.ci.lower` / `.ci.upper`; invariant, so inferential.

- p.star:

  With `.ci.lower` / `.ci.upper`; descriptive unless
  `identified.approx`.

- kkt:

  Strict-complementarity diagnostic: `n.boundary`, `all.positive`,
  `ratio`, and the two margins `delta.dual` and `delta.prim`.

- identified.approx:

  Logical: is the factorization unique up to a permutation of the bases,
  to the precision of the fit? `TRUE` for \\Q=1\\, for
  `X.restriction = "fixed"` (where \\X\\ is not estimated, so \\T=I\\),
  or when the per-column scale is fixed (`"colSums"` / `"colSqSums"`)
  and either every column of \\X\\ has an anchor row at `sep.tol` and
  every row of \\\Theta\\ is pinned at `pin.tol` / `kkt.tol`, or the
  mirror condition holds (see Details). The precision is the margin of
  the route that passed: an equivalent fit mixes the bases by at most
  about that much, and [`print()`](https://rdrr.io/r/base/print.html)
  states it. When `FALSE` the entry-level summaries and `p.star` are
  descriptive; \\\rho\\, the eigenvalues and \\\bm\mu_y\\ remain
  inferential either way.

- identification:

  The parts of the diagnostic: for the first route `anchor` (per column
  of \\X\\), `pinned` (per row of \\\Theta\\) and `pin.order`; for the
  mirror route `pure` (a pure covariate per row of \\\Theta\\),
  `X.pinned` (per basis) and `X.pin.order`; `ok` and `route` (the
  route(s) that passed); and the margins `anchor.margin`, `pin.margin`,
  `pure.margin` and `X.pin.margin` saying how far the fit is from an
  exact configuration (0 = exact; the pin margins are the largest
  relative entry the pinning had to treat as zero).

- refit.args:

  The arguments the bootstrap re-fits actually used, so the replicates
  can be checked against the original estimator.

Plus the bookkeeping entries `wild.B`, `wild.B.requested`, `n.fail`,
`truncate.rate`, `wild.unit`, `wild.dist`, `wild.level`, `dims`.

## Details

**What is identified, and why that decides what can be inferred.** The
factorization admits \\(X,\Theta)\to(XT,T^{-1}\Theta)\\ for invertible
\\T\\, so an entry of \\\Theta\\ means nothing across replicates;
[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
deals with this by conditioning on the estimated \\X\\. This function
instead re-estimates the basis in every replicate, propagating its
uncertainty. That splits the output in two:

- **Inferential.** \\\rho\\, the eigenvalues and the observed long-run
  mean \\\bm\mu_y=X\bm\mu_b\\ are invariant under the whole family
  (\\G_d\to T^{-1}G_dT\\ leaves the spectrum alone), so their replicates
  need no alignment and their intervals mean what they say.

- **Conditional on identification.** The entries \\(G_d)\_{qq'}\\ and
  the composition \\\bm p^\ast\\ need \\T\\ to be a permutation, and
  three conditions together deliver that: the column scale fixed
  (`X.restriction = "colSums"`), an **anchor row** for every column of
  \\X\\, and every row of \\\Theta\\ **pinned**. The anchor rows force
  \\T\ge0\\. The rows are then pinned one at a time: row \\k\\ joins the
  pinned set \\S\\ when some covariate (column of \\\Theta\\) is
  non-zero only for bases in \\S\cup\\k\\\\, including \\k\\, and, for
  every \\q\in S\\, such a covariate has a zero coefficient for \\q\\.
  The first row pinned needs a **pure column** (a covariate that loads
  on that basis alone); later rows only need zeros that separate them
  from the rows already pinned, possibly in different columns. Each
  pinned row makes \\e_k\\ a column of \\T\\, so pinning every row makes
  \\T\\ a permutation. (This used to require a pure column for every
  row, which is stronger: every fit that passed before passes now.)
  **Either condition alone is not enough**: with \\X=I_2\\,
  \\\Theta=\left(\begin{smallmatrix}2&3\\1&1\end{smallmatrix}\right)\\
  and
  \\T=\left(\begin{smallmatrix}0.9&0.1\\0.1&0.9\end{smallmatrix}\right)\\
  the basis is perfectly separable and column-normalized, yet
  \\XT\ge0\\, \\T^{-1}\Theta\ge0\\, \\X\Theta\\ is unchanged, \\T\\ is
  no permutation, and \\G\\ differs between the two representatives — no
  row of \\\Theta\\ can be pinned there. Nor are zeros that merely
  separate every pair of rows enough: with anchors and the covariate
  columns \\(0,1,1)\\, \\(1,0,1)\\, \\(1,1,0)\\, \\T=(J-I)/2\\ is
  admissible. A coefficient counts as zero when it is at most `kkt.tol`
  or at most `pin.tol` times the coefficient of the same covariate for
  the basis being pinned (that ratio is the mixing the covariate still
  allows), and a covariate whose coefficients are all numerical noise is
  ignored. **The mirror route** swaps the two factors and is checked as
  well, with nothing to choose: a pure covariate for every basis
  (forcing \\T^{-1}\ge0\\), and the bases pinned one at a time through
  the rows of \\X\\ – the first by an anchor row, the later ones by
  variables whose zero loadings separate them from the bases already
  pinned (an entry of \\X\\ counts as zero when it is at most
  `1 - sep.tol` times the loading of the basis being pinned). Either
  route suffices; `identification$route` names the one that passed.
  `identified.approx` reports the verdict and `identification` its
  parts. **What a pass means, to the precision of the fit**:
  multiplicative updates leave small positive values where the argument
  wants exact zeros, so the check treats values within the tolerances as
  zero and reports how exact they are. The margin of the route that
  passed is the largest such value relative to the entry it was compared
  with, i.e. the rotation slack: every equivalent fit \\(XT,
  T^{-1}\Theta)\\ mixes one basis into another by at most about that
  much, so the entries of \\G_d\\ are determined to about that relative
  precision. [`print()`](https://rdrr.io/r/base/print.html) states it
  (e.g. "the zeros it rests on hold to 0.0024, so an equivalent fit
  mixes the bases by at most about 0.24\\ part that fails. Two things
  lie outside the check: whether the estimated zeros are zeros in the
  population (a question for their intervals), and other fits that are
  only nearly as good (different starting values can reach one).
  Replicates are permutation-aligned to the original basis regardless.
  `X.restriction = "fixed"` is the one clean case for \\Q\>1\\ — the
  basis is not estimated, so \\T=I\\.

**`prob.nonstationary` is not a unit-root test.** It is the bootstrap
tail probability \\P^\ast(\rho^\ast\ge1)\\, with replicates drawn around
the *estimated* model. Nothing imposes \\H_0:\rho=1\\, so it does not
have the interpretation of a p-value for that hypothesis; a proper test
would have to resample under the null.

**The zeros of \\\Theta\\ are only meaningful under misspecification.**
The strict-complementarity diagnostic reports the dual margin
\\\delta\_{\mathrm{dual}}=\min\|\Lambda^\ast\|\\ over the zero
coefficients and the primal margin
\\\delta\_{\mathrm{prim}}=\min\Theta\\ over the positive ones. Both must
be bounded away from zero for the support to be recovered exactly and
for the active coefficients to attain the distribution that knowing the
support would give. If the model is correctly specified,
\\\Xi_0=X\Theta\\ makes \\\Lambda^\ast=0\\ exactly, so strict
complementarity fails at *every* zero: \\P(\hat\theta_k=0)\\ then tends
to a constant in \\(0,1)\\ (one half for an isolated zero) instead of to
one, and no resampling scheme is consistent there. A small `kkt$ratio`
is therefore not a numerical defect — it says the zeros are not
estimating a sign restriction. What makes them estimate something is
misspecification: they recover the coordinates whose *unconstrained*
population coefficient is strictly negative.

**Bootstrap validity at the boundary.** A naive nonparametric bootstrap
is inconsistent for a parameter on the boundary. The scheme used here
resamples the residuals and re-imposes the non-negativity constraint in
every re-fit, which is the practical remedy. `wild.unit = "block"` in
addition keeps the serial dependence of the residuals within a block;
the lagged design itself stays fixed (see below).

**Bootstrap p-values invert the centred distribution.** Every replicate
of a non-negative coefficient is \\\ge 0\\, so comparing the raw
replicates with zero would return \\p = 0\\ for every entry. The
reported `G.p.value` is the basic-bootstrap two-sided \\p\\ obtained
from \\P^\ast(\hat G^\ast \ge 2\hat G)\\.

**Resampling scheme.** With \\\hat Y = X\Theta A\\ and residuals
\\E=Y-\hat Y\\, each replicate forms \\Y^\ast=\hat Y+E\odot W\\ and
re-fits [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
at the same rank. \\W\\ carries one multiplier per time point
(`wild.unit = "column"`, the default, which preserves the
contemporaneous covariance between series), one per block of consecutive
time points (`wild.unit = "block"`, below), or one per cell. With i.i.d.
multipliers over time the intervals ignore residual autocorrelation and
may be too narrow when it is present. Negative entries of \\Y^\ast\\ are
truncated at 0 to keep the response non-negative; the truncation rate is
reported, and a large one is a warning sign for the whole procedure.

**Block wild bootstrap (`wild.unit = "block"`).** One multiplier is
drawn per block of \\L = \max(1, \mathrm{round}(n^{1/3}))\\ consecutive
time points and shared by all series, where \\n\\ is the number of time
points, `ncol(Y)` (after the lags are removed): \\L = 3\\ for \\n =
27\\, 5 for 100, 8 for 500, 10 for 1000, and \\L = 1\\ (the same draws
as `"column"`) for \\n \le 3\\. The block grid starts at a random offset
in \\0, \dots, L-1\\, drawn anew in each replicate, so the first and
last blocks may be shorter and no time point is always at a block edge.
Residuals within a block keep their serial dependence, which matters
when the residuals are autocorrelated (for example when the lag order is
too small); for a correctly specified VAR the default `"column"`
suffices. \\n^{1/3}\\ is the usual rate for the block length of block
resampling. The length is fixed by \\n\\; it cannot be set and is not
stored in the result.

**This is a fixed-design bootstrap**: the covariate matrix \\A\\ is held
at the observed lags rather than rebuilt from \\Y^\ast\\, so the
intervals are conditional on the observed past and do not propagate the
randomness of the lagged design itself. They may therefore be somewhat
optimistic.

## References

Satoh, K. (2025). Applying non-negative matrix factorization with
covariates to multivariate time series data as a vector autoregression
model. *Japanese Journal of Statistics and Data Science*.
arXiv:2501.17446.
[doi:10.1007/s42081-025-00314-0](https://doi.org/10.1007/s42081-025-00314-0)

## See also

[`nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.md),
[`nmfkc.ar.stationarity`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.stationarity.md),
[`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md),
[`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y   <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
ar  <- nmfkc.ar(Y, degree = 1)
fit <- nmfkc(ar$Y, ar$A, rank = 2, verbose = FALSE)
lat <- nmfkc.ar.latent(fit)
nmfkc.ar.latent.inference(lat, ar$Y, ar$A, wild.B = 50)
#> Warning: The factorization is not even approximately unique up to a permutation: basis 1, 2 has no anchor row in X (so T >= 0 is not forced); row 1, 2 of Theta cannot be pinned (no covariate separates it from the other bases); nor does the mirror condition hold (row 1, 2 of Theta has no pure covariate). The entries of G_d and p.star are therefore DESCRIPTIVE; rho, the eigenvalues and mu.y are invariant under the whole family and remain the inferential targets.
#> Bootstrap inference for the latent VAR
#>   B = 50 (0 failed), unit = "column", dist = "rademacher", truncation rate = 0.32%
#> 
#> Latent transition matrices  G_d = Theta_d %*% X   (row = effect at t, column = cause at t-d)
#>   NOTE: the identification diagnostic FAILS, so these
#>        entries (and p.star) are DESCRIPTIVE. rho, the eigenvalues
#>        and mu.y are invariant and inferential.
#> 
#> lag1:
#>        Basis1 Basis2
#> Basis1  0.189  0.140
#> Basis2  0.085  0.083
#> SE:
#>        Basis1 Basis2
#> Basis1  0.152  0.141
#> Basis2  0.172  0.132
#> p-value:
#>        Basis1 Basis2
#> Basis1    0.2   0.28
#> Basis2    1.0   0.64
#> 
#> rho(G) = 0.2573   SE = 0.1484   95% CI = [0.0707, 0.5778]
#>   bootstrap mean = 0.3209  (bias +0.0636)
#>   bootstrap tail P*(rho* >= 1) = 0.000
#>   (a tail probability around the fitted model, NOT a unit-root test: the
#>    replicates are not drawn under an imposed H0: rho = 1)
#> 
#> complex eigenvalues: point estimate no, in 0.0% of replicates
#> 
#> long-run composition p*:
#>        estimate  lower  upper
#> Basis1   0.4749 0.2466 0.7462
#> Basis2   0.5251 0.2538 0.7534
# }
```

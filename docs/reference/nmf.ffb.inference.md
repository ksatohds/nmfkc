# Bootstrap inference for NMF-FFB (X held fixed)

`nmf.ffb.inference` performs statistical inference on the structural
coefficient matrices \\C_1\\ (latent feedback, \\\Theta_1\\) and \\C_2\\
(exogenous loading, \\\Theta_2\\) from a fitted
[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md)
model, *after* the feed-forward null has been rejected by
[`nmf.ffb.test`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md).
The basis \\\hat X\\ is held fixed throughout, which avoids label
switching and gives a clean conditional interpretation: uncertainty of
the structural coefficients given the measurement model. Which
resampling scheme is run depends on `object$method`.

**`method = "fiml"`: a parametric bootstrap from the selected model.**
\\B\\ data sets are drawn from the selected model (with its
\\\Theta_1\\) and \\(\Theta_1, \Theta_2, \Phi, \psi)\\ are re-estimated
with \\X\\ and the selected support held fixed. Confidence intervals are
the centred (basic) percentile intervals \\\[2\hat\theta -
q^\*\_{1-\alpha/2},\\ 2\hat\theta - q^\*\_{\alpha/2}\]\\; `support_rate`
is the share of replicates with \\\|\hat\theta^\*\| \>\\ `threshold` (0
for entries outside the selected support), and `prob.unsupported` is
\\1 - \mathrm{support\\rate}\\.

These are *post-selection* quantities, conditional on the support that
BIC picked, and they can overstate the evidence. They answer "how large
is this entry, given that the procedure kept it", not "is there feedback
at all"; the calibrated statement about the presence of feedback is
[`nmf.ffb.test`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md),
which is the step that belongs before this one. Each replicate seeds
itself (`seed + B + b`), so the result does not depend on `cores`.

**`method = "mu"` (legacy fits, and objects without a `method` field).**
The procedure is a **full pair bootstrap** that holds the basis matrix
\\\hat X\\ from the original fit fixed across all replicates (which
avoids label switching and gives a clean conditional interpretation:
“uncertainty of the structural coefficients given the measurement
model”):

1.  For each replicate \\b = 1, \dots, B\\, resample column indices
    \\(i_1, \dots, i_N)\\ with replacement from \\\\1, \dots, N\\\\ and
    form \\Y_1^{(b)} = Y_1\[, i\]\\, \\Y_2^{(b)} = Y_2\[, i\]\\.

2.  Re-estimate \\(C_1^{(b)}, C_2^{(b)})\\ by running the
    [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md)
    multiplicative updates *with \\X = \hat X\\ held fixed* (no \\X\\
    update; no centroid sort), using the same `C1.L1`, `C2.L1` as the
    original fit.

3.  Discard replicates that violate stationarity (\\\rho(X C_1^{(b)})
    \ge 1\\) or have an amplification ratio exceeding the
    geometric-series bound by more than 1\\

Because \\C_1, C_2 \ge 0\\ are non-negative by construction, exact zeros
are essentially never observed in the bootstrap distribution.
Significance is assessed via a **support rate** at a small display
threshold \\\delta\\ (default `0.01`): \$\$ \mathrm{sup}(c) \\=\\
\frac{1}{\|\mathrm{valid}\|} \sum\_{b \in \mathrm{valid}}
\mathbf{1}\\\left( \hat c^{(b)} \> \delta \right). \$\$ This is a
one-sided counterpart of the classical \\p\\-value: large `support_rate`
indicates strong evidence that the entry is meaningfully positive.
Significance markers follow the lavaan convention with the natural
correspondence \\p = 1 - \mathrm{sup}\\: `*` (sup \> 0.95), `**` (sup \>
0.99), `***` (sup \> 0.999). Cutoffs use strict greater-than so the rule
mirrors the standard R convention for p-values (p \< 0.05 / 0.01 / 0.001
→ */**/***), translated to support_rate via \\\mathrm{sup} = 1 - p\\.

## Usage

``` r
nmf.ffb.inference(
  object,
  Y1,
  Y2,
  B = 1000L,
  threshold = 0.01,
  boot.level = 0.95,
  C1.L1 = 1,
  C2.L1 = 0.1,
  seed = 123L,
  ...
)
```

## Arguments

- object:

  A fitted object returned by
  [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md).
  Must contain `X`, `C1`, `C2`.

- Y1:

  Endogenous variable matrix (P1 x N). Must match the data used in
  [`nmf.ffb()`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md).

- Y2:

  Exogenous variable matrix (P2 x N). Same.

- B:

  Number of bootstrap replicates. Default `1000`, the value the
  published analysis used; the other inference functions in the package
  default their `wild.B` to 500, and this one is deliberately left at
  1000 so the manuscript results reproduce out of the box. Note that
  `***` (sup \> 0.999) is only reachable when *every* replicate clears
  `threshold`, at any `B`: raising `B` does not add a finer grade, it
  makes the same grade stronger evidence (all 1000 rather than all 500).
  Lowering it to 500 roughly halves the running time, since each
  replicate is a re-fit. For `method = "fiml"` a replicate is one FIML
  fit on the selected support, so `B = 1000` is affordable; the
  expensive bootstrap is the null one, which
  [`nmf.ffb.test`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md)
  runs.

- threshold:

  Display threshold \\\delta\\ for the support rate
  \\\Pr\_{\mathrm{boot}}(\hat c^{(b)} \> \delta)\\. Default `0.01`;
  entries below this magnitude are treated as effectively zero in the
  path diagram.

- boot.level:

  Confidence level for the bootstrap CI. Default `0.95`.

- C1.L1, C2.L1:

  L1 sparsity penalties used by the original `method = "mu"` fit. These
  must match the fit's hyperparameters for the bootstrap to estimate the
  correct model. Defaults (`1.0`, `0.1`) match `nmf.ffb`'s defaults but
  you should pass the actual values used. Not used for
  `method = "fiml"`.

- seed:

  Base RNG seed for the bootstrap. For `method = "mu"` each replicate
  uses `seed + b` (resampling) and `seed + 1000 + b` (\\C_1, C_2\\
  initialization); for `method = "fiml"` see Description. Default `123`.

- ...:

  Hidden options. Shared: `cores` (number of parallel workers, default
  `getOption("mc.cores", 1L)` for `"fiml"`, `1` for `"mu"`; `ncores` is
  accepted as an alias) and `print.trace`. For `method = "fiml"`:
  `factr` and `fiml.maxit` (L-BFGS-B tolerance and cap, default the
  values recorded on `object`). The null bootstrap is no longer run here
  at all: see
  [`nmf.ffb.test`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md).
  (Removed: `boot.null`, which used to skip it, and the coefficient
  intervals). For `method = "mu"`:

  `epsilon`

  :   Convergence tolerance for the inner fixed-X MU loop. Default
      `1e-8` – deliberately tighter than a plain fit, because under
      multiplicative updates an entry heading for the non-negativity
      boundary approaches 0 slowly, so a loose tolerance stops every
      replicate while such an entry is still above `threshold`. The
      support rate is then inflated and significance over-declared: on a
      pure-noise design the support of null entries falls from 0.75 to
      0.20 to 0.017 as `epsilon` goes 1e-6, 1e-8, 1e-10, while genuinely
      supported entries stay at 1. Tighten further when near-threshold
      entries matter; loosen only for speed.

  `maxit`

  :   Maximum iterations for the inner MU loop. Default `100000` (the
      tighter tolerance needs the headroom).

## Value

The input `object` with additional bootstrap inference components:

- coefficients:

  Data frame with rows for every entry of \\C_1\\ and \\C_2\\ and
  columns `Type` ("C1" / "C2"), `Basis`, `Covariate`, `Estimate`,
  `CI_low`, `CI_high`, `support_rate`, `prob.unsupported` (\\= 1 -
  \mathrm{support\\rate}\\) and `sig`.

- C1.support.rate, C2.support.rate:

  Per-element support rates (Q x P1 and Q x P2 matrices).

- C1.ci.lower, C1.ci.upper, C2.ci.lower, C2.ci.upper:

  Per-element CI bounds (percentile for `"mu"`, centred percentile for
  `"fiml"`).

- C1.boot.draws, C2.boot.draws:

  Bootstrap distributions: 3D arrays of shape B x Q x P1 (and B x Q x
  P2). Invalid replicates contain `NA`.

- rho.boot.draws:

  Per-replicate spectral radius \\\rho(X C_1^\*)\\.

- AR.boot, iter.boot:

  (`"mu"` only) per-replicate amplification ratio and inner-loop
  iteration count.

\#'

- boot.B, boot.threshold, boot.level:

  Inputs recorded for reproducibility.

- boot.n.valid, boot.n.invalid:

  Validity counts (for `"fiml"`: of the selected-model bootstrap).

- boot.method:

  How the replicates were drawn.

- identified.approx, identification:

  Whether \\(X, \[C_1\\ C_2\])\\ is unique up to a permutation of the
  bases – the diagnostic of
  [`nmfkc.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.inference.md)
  with \\A = \[Y_1; Y_2\]\\ (see the Details of
  [`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md)).
  The inference is conditional on \\\hat X\\; when `TRUE` it concerns
  the same coefficients whichever admissible representative is taken.
  The exclusion restriction's zeros in \\C_1\\ count. An algebraic check
  on the point estimate: no re-fits.

The calibrated test of the feed-forward null – `LR.p.boot`,
`prob.select.null`, `C1.restriction.change.rate` and the null bootstrap
behind them – is **not** returned here. It moved to
[`nmf.ffb.test`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md)
in 0.9.8, so that an object carrying coefficient intervals cannot be
mistaken for a test of feedback (see NEWS).

## Lifecycle

This function's interface changed at v0.6.8: the legacy 1-step Newton
wild bootstrap (with sandwich SE) has been replaced by the full pair
bootstrap described above, following the paper revision. The fields
`sigma2.used`, `C2.se`, `C2.se.boot`, `C2.p.side` that the previous
implementation produced are no longer present.

## References

Satoh, K. (2025). Applying non-negative matrix factorization with
covariates to structural equation modeling for blind input-output
analysis. arXiv:2512.18250. <https://arxiv.org/abs/2512.18250>

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`nmf.ffb.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.DOT.md)

## Examples

``` r
# \donttest{
Y <- t(iris[, -5])
Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
res  <- nmf.ffb(Y1, Y2, rank = 2)
## intervals for the entries the test has already licensed reading.
inf <- nmf.ffb.inference(res, Y1, Y2, B = 20)   # use B = 1000 in practice
head(inf$coefficients)  # estimate, interval and support rate per entry
#>   Type   Basis    Covariate  Estimate    CI_low   CI_high support_rate
#> 1   C1 Factor1 Sepal.Length 0.0000000 0.0000000 0.0000000            0
#> 2   C1 Factor2 Sepal.Length 0.0000000 0.0000000 0.0000000            0
#> 3   C1 Factor1  Sepal.Width 0.0000000 0.0000000 0.0000000            0
#> 4   C1 Factor2  Sepal.Width 0.0000000 0.0000000 0.0000000            0
#> 5   C2 Factor1 Petal.Length 1.4069840 1.3190281 1.8602862            1
#> 6   C2 Factor2 Petal.Length 0.4013679 0.1949612 0.6157461            1
#>   prob.unsupported sig
#> 1                1    
#> 2                1    
#> 3                1    
#> 4                1    
#> 5                0 ***
#> 6                0 ***

res.mu  <- nmf.ffb(Y1, Y2, rank = 2, method = "mu")
#> Warning: maximum iterations (5000) reached...
res.mu2 <- nmf.ffb.inference(res.mu, Y1, Y2, B = 200)
head(res.mu2$coefficients)
#>   Type   Basis    Covariate     Estimate       CI_low      CI_high support_rate
#> 1   C1 Factor1 Sepal.Length 9.624411e-01 9.978009e-01 9.983034e-01            1
#> 2   C1 Factor2 Sepal.Length 4.652216e-03 3.049265e-03 4.505394e-03            0
#> 3   C1 Factor1  Sepal.Width 4.303632e-02 2.426654e-04 4.476438e-04            0
#> 4   C1 Factor2  Sepal.Width 9.903750e-01 9.904116e-01 9.935386e-01            1
#> 5   C2 Factor1 Petal.Length 2.318005e-02 2.237557e-03 2.903178e-03            0
#> 6   C2 Factor2 Petal.Length 3.800995e-09 5.063817e-14 6.615945e-10            0
#>   prob.unsupported sig
#> 1                0 ***
#> 2                1    
#> 3                1    
#> 4                0 ***
#> 5                1    
#> 6                1    
# }
```

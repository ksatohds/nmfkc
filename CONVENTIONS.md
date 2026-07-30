# nmfkc package conventions

Decisions that apply across the whole package. Read this before adding a
function or changing a default; several of the rules below were arrived at by
measurement, and the measurements are recorded so they can be re-checked rather
than re-argued.

This file is excluded from the build (`.Rbuildignore`).

---

## 1. Naming

**Arguments.** Inference functions take their bootstrap controls through `...`
with the `wild.` prefix:

| argument | meaning | default |
|:---|:---|:---|
| `wild.B` | number of bootstrap replicates | 500 |
| `wild.dist` | multiplier distribution | function-specific |
| `wild.unit` | `"element"` or `"column"` | function-specific |
| `wild.level` | confidence level | 0.95 |
| `wild.seed` | bootstrap seed | 123 |

One deliberate exception: `nmf.ffb.inference()` keeps `B = 1000`, the value its
published analysis used, so those results reproduce out of the box. Its help
says so.

**Return values.** Estimate / uncertainty pairs follow
`<quantity>`, `<quantity>.se`, `<quantity>.ci.lower`, `<quantity>.ci.upper`
(e.g. `C.se`, `C.ci.lower`; `spectral.radius.se`, `spectral.radius.ci.lower`).
Coefficient tables use the columns `Basis`, `Covariate`, `Estimate`, `SE`,
`BSE`, `z_value`, `p_value`, `CI_low`, `CI_high` — `nmfkc.DOT()` depends on
`Basis` / `Covariate`.

**Optimization / inference split.** `f()` fits, `f.inference()` does inference,
`f.cv()` / `f.ecv()` / `f.rank()` resample or sweep. Fitters never do inference.

---

## 2. Seeds and the random stream

**Every function that calls `set.seed()` must restore the caller's stream.**
Use the helpers in `R/parallel-utils.R`, right after the seed is resolved:

```r
.rng <- .nmfkc.rng.save(seed)
on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
```

Without this, each call leaves `.Random.seed` at a fixed state, so a user loop
drawing random numbers after a fit silently repeats itself. This was true of 23
entry points before 2026-07-26.

**`seed = NULL` is not expressible through `...`.** Most functions resolve the
seed with

```r
seed <- if (!is.null(extra_args$seed)) extra_args$seed else 123
```

and `list(...)` cannot distinguish "`seed = NULL` was passed" from "`seed` was
not passed". An explicit `NULL` therefore yields the default, not a free-running
stream. Do not write `seed = NULL` in a call and expect it to mean anything.
`NULL` only works where `seed` is a real formal (`nmf.ffb`,
`nmfkc.signed.rff`), and there it means "do not seed" — note that bare
`set.seed(NULL)` re-initializes from the clock, so guard it with
`if (!is.null(seed))`.

**Default seeds are concrete (123), not NULL**, so an analysis is reproducible
without the user having to think about it.

**CV wrappers share one seed across folds, and that is fine.** The concern that
this makes every fold start from an identical initialization does not apply:
the k-means initializer is built from each fold's own training data, so the
initializations differ by fold anyway. Measured on a plain design, they differ
*more* across folds (max abs diff 3.71 column-wise, 0.43 element-wise) than two
different seeds differ on the same data (0.95). Deriving per-fold seeds would
change every CV number for a marginal gain, so we do not.

---

## 3. Changes must not change results

The package is used for published analyses. Any refactor, optimization or
tidy-up is expected to be **bit-identical**, and is verified before it lands:

1. capture the outputs of a representative battery,
2. `git stash` the change, capture the same battery again,
3. compare with `identical()` **and** `all.equal(tolerance = 0)`.

Do the comparison **in one R session**, sourcing the old file into a separate
environment when possible. Comparing `.rds` files written by different sessions
has produced false mismatches here (a stale baseline), which cost real time.

Numerical changes are acceptable only when the old behaviour is a defect, and
then the commit message must say what changed and why.

**Known trap: `tcrossprod(X)` is not `X %*% t(X)`.** The one-argument forms
dispatch to BLAS `dsyrk`, which is not bitwise identical to the `dgemm` of the
two-argument form (observed divergence ~1e-15 at some shapes, e.g. K = 5). Only
the two-argument `tcrossprod(A, B)` / `crossprod(A, B)` are safe swaps under a
bit-identical mandate.

Also not bit-identical, though mathematically equal: reassociating a triple
product, batching per-column `crossprod` into one gemm (gemv vs gemm differ by
~1e-15), and `sum()` (long double) versus a sequential `+` accumulation — use
`Reduce("+", x)` when reproducing a loop's accumulation exactly.

---

## 4. Parallelism is opt-in and result-identical

Sweeps and resampling loops go through `.nmfkc.parlapply()`:

```r
cores <- if (!is.null(extra$cores)) extra$cores else getOption("mc.cores", 1L)
res   <- .nmfkc.parlapply(tasks, run_one, cores = cores)
```

- default is `getOption("mc.cores", 1L)`, i.e. **sequential** — CRAN allows at
  most 2 cores in examples/tests, and a library should not seize every core;
- `cores = 1` takes a plain `lapply`, so it is bitwise the old sequential loop;
- Windows uses a PSOCK cluster, elsewhere forking;
- results are identical for any `cores` **only because each task is
  deterministic given its index** (self-seeded). Do not change per-task seeding
  without revisiting this guarantee;
- nested levels must not both parallelize: the outer level passes `cores = 1`
  to the inner one (`nmf.gmm.select` → `nmf.gmm`).

PSOCK gotchas that cost time here: worker closures need the caller's frame
`clusterExport`ed (lazy promises are forced that way), workers do not inherit
the attached search path (pass every package whose exported symbols a worker
closure or its formulas reference, via `packages =`), and a binding that errors
when forced must be skipped rather than exported.

`nmfkc.cv`'s fold loop is *not* currently parallelized, but not for RNG reasons
— see §2; it simply has not been done.

---

## 5. Convergence tolerance in inference

`nmfkc()` defaults to `epsilon = 1e-4`, which is fine for reconstruction but
**not for inference**. Under multiplicative updates a coefficient heading for
the non-negativity boundary approaches 0 slowly, so a loose tolerance stops
while it is still spuriously positive, and anything compared with zero inherits
the error — in the loose direction, which **over-declares significance**.

The re-fitting inference functions therefore converge tightly:

| function | setting | default |
|:---|:---|:---|
| `nmfkc.inference(method = "refit")` | `refit.epsilon`, `refit.maxit` | `1e-8`, `1e5` |
| `nmf.ffb.inference()` | `epsilon`, `maxit` | `1e-8`, `1e5` |
| `nmfkc.ar.latent.inference()` | `epsilon`, `maxit` | `1e-8`, `2e4` |

Measured consequences of a loose tolerance: on `vars::Canada` the percentile CI
excluded the point estimate for 4 of 10 coefficients at `1e-4` and 0 of 10 at
`1e-8`, and the bootstrap SE was understated about twofold; on a pure-noise
design `nmf.ffb.inference`'s support rate for null entries fell 0.75 → 0.20 →
0.017 as `epsilon` went 1e-6 → 1e-8 → 1e-10 while supported entries stayed at
1.0. Near-threshold entries stay tolerance-sensitive even at `1e-8`; the README
says so.

Fit the model tightly too (`epsilon = 1e-8`) when the fit feeds inference.

---

## 6. Inference under the non-negativity constraint

- **Bootstrap p-values must invert the *centred* replicate distribution.**
  Comparing raw replicates with zero is degenerate when every replicate is
  `>= 0`: it returns p = 0 for every coefficient. Use
  `P*(Chat* >= 2 Chat)` (basic bootstrap) and honour `C.p.side`.
- **The one-sided Wald test is valid at the boundary.** With
  `theta_hat = max(0, theta_tilde)`, `P(theta_hat > t) = P(theta_tilde > t)`
  for every `t > 0`, so its asymptotic size is exactly alpha. The
  `0.5*chisq_0 + 0.5*chisq_1` mixture applies to the *likelihood-ratio*
  statistic, which is a different test — do not invoke it here.
- **A CI that excludes the point estimate is usually a convergence artefact,**
  not an inherent boundary phenomenon; see §5 before reaching for special
  boundary machinery.
- **Prefer invariant targets when the basis is re-estimated.** Entries of
  `Theta` are not identified under `(X, Theta) -> (XT, T^-1 Theta)`;
  `rho(G)`, the eigenvalues of `G`, `Xi_d`, impulse responses and `mu_y` are.
  `nmfkc.ar.latent.inference()` uses this: because `rho` is invariant its
  replicates need no label alignment at all.
- **Do not claim an entry of a rotation-covariant matrix is identified.**
  `X.restriction = "colSums"` removes only the **scale** part of `T` (without it
  `G_qq' -> (t_q'/t_q) G_qq'`, so only the diagonal survives), and separability
  of the **fitted** `X` does not remove the rest. The uniqueness theorem assumes
  the *population* basis is separable and then singles out the separable
  representative; the estimator has no reason to land on it, because the
  objective is flat along the family. Counterexample worth remembering:

      X = I_2,  Theta = [2 3; 1 1],  T = [.9 .1; .1 .9]

  `X` is perfectly separable with unit column sums, yet `XT >= 0`,
  `colSums(XT) = 1`, `T^-1 Theta >= 0`, `X Theta` is unchanged, `T` is no
  permutation, and `G = Theta X` differs between the two representatives. So
  `nmfkc.ar.latent.inference()` sets `identified = TRUE` only for `Q == 1`; for
  `Q >= 2` the entry-level output is labelled **descriptive**.

  Report the invariants (`rho`, the eigenvalues, `mu_y`) as the inferential
  targets. This rule cost two wrong claims before it was written down: first
  "colSums is enough", then "colSums plus separability is enough".
- **A bootstrap tail probability is not a p-value.** Replicates drawn around the
  *estimated* model give `P*(rho* >= 1)`, not a test of `H0: rho = 1` — nothing
  imposes the null. Name such quantities `prob.*`, never `p.*`, and say in the
  help and the print-out what they are not.
- **Perron-Frobenius bounds apply to the matrix whose radius you want.** The
  stationarity bracket uses the column sums of `sum_d Xi_d = X Lambda`, i.e.
  `colSums(X) %*% Lambda`. The unweighted `colSums(Lambda)` is the special case
  `1'X = 1'` and is wrong otherwise: `X = (20, 20)'`, `Lambda = (0.04, 0.04)`
  gives `max_j c_j = 0.04`, which would certify stationarity for a fit with
  `rho = 1.6`.
- **A resampling scheme must re-fit the SAME estimator.** A bootstrap re-fit that
  falls back on `nmfkc()`'s defaults describes a different estimator than the one
  whose uncertainty is being reported. Carry the fit's `method` and
  `X.restriction` (both recorded on the object), inherit a literal `seed` from
  the recorded call, and warn about anything else the call set that cannot be
  recovered — with a `refit.args=` escape hatch.
- **Know which constraint fixes what.** Only `"colSums"` and `"colSqSums"` fix
  the per-column scale of `X`. `"totalSum"` constrains the grand total, so
  columns can still trade scale; `"none"` constrains nothing; `"fixed"` is the
  opposite extreme — `X` is not estimated, so `T = I` and there is no rotational
  freedom at all. Do not lump `"fixed"` in with `"none"`.
- **Do not report a threshold-based diagnostic as a theorem.** The
  anchor-row / pure-column argument needs *exact* anchors and *exact* zeros;
  multiplicative updates leak. Name such flags `*.approx`, return the margins
  saying how far from exact the fit is, and say in the help and print-out that
  this is evidence, not proof.
- **Real roots do not mean monotone decay.** A negative real root alternates in
  sign every period (period 2). Report `alternating` alongside `cycle.period`,
  and note that for `Q = 2, D = 1` with `G >= 0` the discriminant
  `(a-d)^2 + 4bc` is non-negative, so the roots are real *by construction* — the
  absence of a cycle there is a property of the design, not a finding. And
  positive real roots rule out *oscillation*, not *overshoot*: a non-normal `G`
  can send a coordinate up before it comes down. Name the flag after the roots
  (`non.oscillatory`), never after the trajectory (`monotone`).
- **`b_t` follows a VARMA, not a VAR.** With `y_t = X b_t + e_t` the exact
  identity is `b_t = sum_d G_d b_{t-d} + theta + sum_d Theta_d e_{t-d}`. `G_d`
  is the AR operator of the fitted VAR reduced to the latent space; it is not a
  Granger-causality statement between latent coordinates, and `p.star` is the
  fixed point of the *residual-free* composition dynamics, not the limit or the
  mean of `b_t / 1'b_t`.
- **A vanishing KKT dual margin is a statement, not a bug.** Under correct
  specification `Xi_0 = X Theta` makes the population multipliers exactly zero,
  so strict complementarity fails at *every* zero coefficient and
  `P(theta_k = 0)` tends to a constant in (0, 1) — one half for an isolated zero
  — rather than to 1. The estimated zeros only estimate something under
  **misspecification**, where they recover the coordinates whose *unconstrained*
  population coefficient is strictly negative. Report both margins
  (`delta.dual` over the zeros, `delta.prim` over the positives) and say which
  reading applies rather than calling a small ratio "doubtful".
- **A naive nonparametric bootstrap is inconsistent on the boundary.** Resample
  in a way that re-imposes the constraint in every re-fit (as the wild-bootstrap
  re-fits here do). A moving-block scheme would additionally preserve the serial
  dependence that a fixed design discards.

**Function split for NMF-VAR.** `nmfkc.ar.latent()` computes the latent VAR
representation; `plot()` on it draws the dynamics; `nmfkc.ar.latent.inference()`
bootstraps them. `nmfkc.ar.stationarity()` answers the yes/no stationarity
question and nothing else, and `plot()` on it draws the bracket. There is
deliberately no `nmfkc.ar.stationarity.inference()`: `rho` is an eigenvalue of
`sum_d G_d`, so its inference belongs with the other latent quantities, and one
bootstrap (each replicate a full re-fit) serves all of them.

---

## 7. Documentation

- roxygen2 with markdown; `NAMESPACE` and `man/*.Rd` are generated — never
  edit them by hand.
- `\bm` is used throughout the Rd sources (138 files) and ships to CRAN fine;
  keep new pages consistent rather than switching to `\boldsymbol`.
- Vignette `.html` files are git-ignored build products. Editing a `.Rmd`
  does **not** update the rendered page a reader is looking at — rebuild it.
- When a comment explains *why* a value or an argument is what it is, keep it
  true. Three comments in this package described behaviour that did not
  happen, and each cost an investigation.

---

## 8. Known divergences, deliberately not fixed

Two places differ from the rest of the package in ways that are *wrong by the
package's own conventions* but were left alone because correcting them would
move published numbers. Both are stated in the affected help pages; do not
"tidy" them without deciding to accept the change.

- **`nmf.ffb` L1 is twice as strong.** It adds `C1.L1` to the multiplicative
  denominator where `nmfkc` / `nmfae` / `nmfkc.net` add `C.L1 / 2`. The same
  nominal penalty therefore shrinks twice as hard in `nmf.ffb`.
- **`nmf.ffb` and `nmfkc.net` converge on the unpenalized objective.** The
  break test reads the reconstruction loss, not the penalized objective the
  updates minimize; every other optimizer tests the penalized value. With a
  strong penalty they can stop while the optimized quantity is still moving.
  `nmf.ffb` returns both traces (`objfunc`, `objfunc.full`) so the gap is at
  least visible.

By contrast these *were* fixed, because the results they produced were not
merely different but wrong: a KL fit whose relative-change denominator lacked
`abs()` and so "converged" in two iterations at any tolerance
(`nmf.rrr`); `X.restriction = "fixed"` that did not hold `X` fixed and
`X.L2.ortho` that was never forwarded (`nmfkc.net`); and a signed
initialisation that drew every entry of `C0` from `U(-1, 1)`, so with
probability `(1/2)^3` at `Q = 2` the whole matrix came out negative, `Cp` was
identically zero and `X` collapsed to zeros on the first iteration.

---

## 9. Housekeeping

- Commit messages are date-prefixed, e.g. `20260726 <what changed>`.
- Version = total commits / 1000, last digit dropped (883 commits → 0.8.8).
- Run `devtools::check()` before merging `develop` → `main`. The only expected
  NOTE offline is "unable to verify current time".
- Temporary scripts and comparison `.rds` files belong in the session scratch
  directory, never in the repo.

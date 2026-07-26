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
the attached search path (pass `packages = c("nmfkc", "survival")` when a
formula needs `Surv()`), and a binding that errors when forced must be skipped
rather than exported.

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
| `nmfkc.ar.stationarity.inference()` | `epsilon`, `maxit` | `1e-8`, `2e4` |

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
  `nmfkc.ar.stationarity.inference()` uses this: because `rho` is invariant its
  replicates need no label alignment, and only the composition `p*` (invariant
  merely up to basis order) is aligned.

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

## 8. Housekeeping

- Commit messages are date-prefixed, e.g. `20260726 <what changed>`.
- Version = total commits / 1000, last digit dropped (883 commits → 0.8.8).
- Run `devtools::check()` before merging `develop` → `main`. The only expected
  NOTE offline is "unable to verify current time".
- Temporary scripts and comparison `.rds` files belong in the session scratch
  directory, never in the repo.

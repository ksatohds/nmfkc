# NMF-RE element-wise cross-validation for rank selection

Selects the basis rank \\Q\\ for
[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md) by
Wold-style element-wise (entry-holdout) cross-validation that
**exercises the random effects**. For each fold the held-out entries are
hidden and filled by iterative imputation: fit the NMF-RE model on the
current matrix, replace the held-out entries with the BLUP prediction
\\X(\Theta A + U)\\, and repeat. The score is the held-out prediction
RMSE (`sigma`); the selected rank minimizes it.

Because the held-out entries of a column are predicted using the random
effect \\u_n\\ fitted from that column's *retained* entries, this
evaluates the full NMF-RE model (including \\U\\) — unlike
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md),
which masks entries by zero weight and predicts from the fixed-effect
fit \\X\Theta A\\ only. The two scores are therefore *not* directly
comparable.

## Usage

``` r
nmfre.ecv(Y, A = NULL, rank = 1:3, C.signed = TRUE, ...)
```

## Arguments

- Y:

  Observation matrix (P x N), non-negative.

- A:

  Covariate matrix (K x N). Default is a row of ones (intercept only).

- rank:

  Integer vector of ranks \\Q\\ to evaluate (default `1:3`).

- C.signed:

  Logical. Sign convention for \\\Theta\\ passed to each
  [`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md) fit
  (`TRUE` = sign-free, default; `FALSE` = non-negative). The basis
  update rule follows this choice automatically.

- ...:

  Additional arguments:

  - `nfolds`: Number of folds (default 5; legacy `nfold` also accepted).

  - `rounds`: Iterative-imputation rounds per fold (default 4).

  - `seed`: RNG seed for the fold assignment (default 1).

  - `cores`: Integer; number of parallel workers for the rank sweep
    (default `getOption("mc.cores", 1L)`). Because the fold matrix is
    fixed before the sweep and each
    [`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md) fit
    self-seeds, results are identical to the sequential run for any
    `cores`.

  - `print.trace`: Logical; print per-rank scores (default `FALSE`).

  - Convergence controls forwarded to
    [`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md):
    `epsilon` (default `1e-5`), `epsilon.outer` (default `1e-3`),
    `inner.maxit` (default `1500`), `outer.maxit` (default `80`),
    `maxit` (default `40000`). CV does not need the tight tolerances of
    a final fit, so these are loosened by default.

## Value

A list of class `"nmfre.ecv"` with components:

- `rank`:

  The ranks evaluated.

- `sigma`:

  Named numeric vector of held-out RMSE per rank (same field name as
  [`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)).

- `best`:

  The rank minimizing `sigma`.

- `nfolds`, `rounds`, `C.signed`:

  Settings used.

## See also

[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`nmfre.inference`](https://ksatohds.github.io/nmfkc/reference/nmfre.inference.md),
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)

## Examples

``` r
# \donttest{
if (requireNamespace("nlme", quietly = TRUE)) {
  Y <- matrix(nlme::Orthodont$distance, 4, 27)
  male <- ifelse(nlme::Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
  A <- rbind(intercept = 1, male = male)
  cv <- nmfre.ecv(Y, A, rank = 1:3)
  cv$best
  plot(cv)
}

# }
```

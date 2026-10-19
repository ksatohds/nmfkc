# Select the bandwidth of random-feature covariates by cross-validation

Chooses the Gaussian-kernel bandwidth \\\beta\\ (and optionally the
number of random features \\D\\) for Random Fourier Feature covariates
by cross-validation, in the same way that
[`nmfkc.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.cv.md)
does for exact kernel covariates. For every candidate the random
features are regenerated on the columns of `U` and the model is
cross-validated column-wise:

- `type = "signed"`: cosine features from
  [`nmfkc.signed.rff`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md)
  (signed), fitted with
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  through
  [`nmfkc.signed.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.cv.md);

- `type = "positive"`: positive random features from
  [`nmfkc.rff.positive`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.md)
  (non-negative), fitted with
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) through
  [`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md).

The candidate with the smallest cross-validated prediction error is
returned. Since the same random map cannot be shared across folds
without leaking (the map is data-independent, so it can: only the fit is
cross-validated), one draw per candidate is used and its `seed` is fixed
so the comparison across \\\beta\\ is paired.

For large \\N\\ set `sample.size`: the selection then runs on a random
subsample of the columns (default: all columns), after which the chosen
\\\beta\\ is used with
[`nmfkc.signed.rff.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.gram.md)
/
[`nmfkc.rff.positive.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.gram.md)
on the full data.

## Usage

``` r
nmfkc.rff.beta.cv(
  Y,
  rank = 2,
  U,
  beta = NULL,
  D = 500L,
  type = c("signed", "positive"),
  plot = TRUE,
  ...
)
```

## Arguments

- Y:

  Response matrix \\P \times N\\ (non-negative for `type = "positive"`).

- rank:

  Number of bases \\Q\\.

- U:

  Input matrix \\p \times N\\; columns are samples.

- beta:

  Numeric vector of candidate bandwidths. Default `NULL`: the
  seven-point grid of
  [`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md)
  around the median nearest-neighbour heuristic (computed on the
  subsample when `sample.size` is set).

- D:

  Integer vector of candidate numbers of random features. A single value
  (default `500`) selects \\\beta\\ only; several values select over the
  \\(\beta, D)\\ grid.

- type:

  `"signed"` (cosine RFF, default) or `"positive"`.

- plot:

  Logical. If `TRUE` (default), plot the CV objective against `beta`
  (one line per `D`).

- ...:

  Passed to the CV function
  ([`nmfkc.signed.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.cv.md)
  or
  [`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md):
  e.g. `nfolds`, `epsilon`, `maxit`). Also accepts:

  - `sample.size`: number of columns used for the selection (default
    all).

  - `seed`: integer seed for the subsample and the random maps (default
    123; the caller's stream is restored).

  - `cores`: evaluate the candidates in parallel with `.nmfkc.parlapply`
    (default `getOption("mc.cores", 1L)`). Each candidate is
    deterministic given its index, so the result is identical for any
    `cores`.

  - `hyperbolic`: passed to
    [`nmfkc.rff.positive`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.md).

## Value

A list with

- `beta`, `D`:

  The selected bandwidth and feature count.

- `objfunc`:

  Matrix of CV objective values, `length(beta)` rows by `length(D)`
  columns.

- `beta.candidates`, `D.candidates`:

  The grids.

- `type`, `sample.size`, `n.used`:

  Settings used.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## See also

[`nmfkc.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.cv.md)
(exact kernel; also works for Nystr\\om covariates with
`U = landmarks, V = data`),
[`nmfkc.signed.rff.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.gram.md),
[`nmfkc.rff.positive.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.gram.md)

## Examples

``` r
# \donttest{
data(iris)
set.seed(1)
idx <- sample(nrow(iris), 100)
U <- t(scale(iris[idx, 1:4]))
levs <- levels(iris$Species)
Y <- sapply(iris$Species[idx], function(s) as.integer(levs == s)); rownames(Y) <- levs

## signed RFF: choose beta over the median-heuristic grid at D = 50
sel <- nmfkc.rff.beta.cv(Y, rank = 3, U, D = 50, plot = FALSE)
sel$beta
#> [1] 0.2681876
g   <- nmfkc.signed.rff.gram(Y, U, beta = sel$beta, D = 50, seed = 1)
fit <- nmfkc.signed(Y, A = g, rank = 3, verbose = FALSE)

## positive RFF: select over a (beta, D) grid
sel2 <- nmfkc.rff.beta.cv(Y, rank = 3, U, D = c(100, 400), type = "positive",
                          plot = FALSE, verbose = FALSE)
#> Warning: candidate beta = 577.8, D = 100 failed (nmfkc(): the objective is NaN at iteration 10.  A column of X became zero and could not be normalized; this can happen when the rank is larger than the data support.  Try a smaller rank or X.init = "nndsvd".); treated as Inf.
#> Warning: candidate beta = 124.5, D = 100 failed (nmfkc(): the objective is NaN at iteration 10.  A column of X became zero and could not be normalized; this can happen when the rank is larger than the data support.  Try a smaller rank or X.init = "nndsvd".); treated as Inf.
#> Warning: candidate beta = 26.82, D = 100 failed (nmfkc(): the objective is NaN at iteration 10.  A column of X became zero and could not be normalized; this can happen when the rank is larger than the data support.  Try a smaller rank or X.init = "nndsvd".); treated as Inf.
#> Warning: candidate beta = 577.8, D = 400 failed (nmfkc(): the objective is NaN at iteration 10.  A column of X became zero and could not be normalized; this can happen when the rank is larger than the data support.  Try a smaller rank or X.init = "nndsvd".); treated as Inf.
#> Warning: candidate beta = 124.5, D = 400 failed (nmfkc(): the objective is NaN at iteration 10.  A column of X became zero and could not be normalized; this can happen when the rank is larger than the data support.  Try a smaller rank or X.init = "nndsvd".); treated as Inf.
#> Warning: candidate beta = 26.82, D = 400 failed (nmfkc(): the objective is NaN at iteration 10.  A column of X became zero and could not be normalized; this can happen when the rank is larger than the data support.  Try a smaller rank or X.init = "nndsvd".); treated as Inf.
sel2$objfunc
#>                D=100      D=400
#> 577.79268        Inf        Inf
#> 124.48166        Inf        Inf
#>  26.81876        Inf        Inf
#>   5.77793 0.40283493 0.60412928
#>   1.24482 0.24717367 0.55493823
#>   0.26819 0.09430055 0.06380387
#>   0.05778 0.07752821 0.07436292
# }
```

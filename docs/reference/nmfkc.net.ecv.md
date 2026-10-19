# Element-wise cross-validation for nmfkc.net (upper-triangle folds)

k-fold CV with folds taken over the upper triangle of the symmetric
\\Y\\ (mirrored to the lower triangle) to prevent information leakage
through symmetry. A single entry point covers all three symmetric NMF
variants; `type` selects the fitting function:

- `"tri"` (default):
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  with \\C \ge 0\\

- `"bi"`:
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  with \\C = I_Q\\

- `"signed"`:
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  with signed \\C = C\_{+} - C\_{-}\\

## Usage

``` r
nmfkc.net.ecv(Y, rank = 1:3, type = c("tri", "bi", "signed"), ...)
```

## Arguments

- Y:

  Symmetric N x N non-negative matrix.

- rank:

  Integer vector of ranks to evaluate. Default `1:3`.

- type:

  Model type: `"tri"` (default), `"bi"`, or `"signed"`.

- ...:

  Passed to the underlying fitter; also accepts:

  - `nfolds` (default 5; `div` alias) and `seed` (default 123).

  - `diag.exclude` (default `FALSE`): leave the diagonal out of the
    folds and out of every fit (see
    [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)),
    for a network without self-loops. Otherwise its zeros would be
    scored as edges to predict.

  - `folds`: the folds to use, as a list of linear indices into \\Y\\ –
    e.g.\\ the `$folds` of an earlier run, so that two models are
    compared on the same split. A lower-triangle index is read as its
    upper-triangle mirror. An entry in two folds, or one that cannot be
    held out (`NA`, zero weight, or the diagonal under `diag.exclude`),
    is an error. `nfolds` is then `length(folds)`. (Within this function
    the folds depend only on \\N\\, `nfolds`, `seed` and the entries
    that can be held out, so runs that share these already share their
    folds.)

  - `seeds`: a vector of seeds (overriding `seed`), to repeat the CV
    over different splits. Each split refits with its seed, so the
    results for seed \\s\\ are those of `seed = s`. `objfunc`, `sigma`
    and `r.squared.cv` are then averaged over the seeds, with their
    spread in `objfunc.sd` and `r.squared.cv.sd`. Not with `folds`.

  - `pred` (default `FALSE`): also return the held-out predictions (an
    \\N \times N\\ array per rank; kept optional because it is \\N^2\\
    per rank).

  - `Y.weights`: as in
    [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md).
    Entries of weight 0 are never held out, every fit keeps the weights,
    and the held-out loss is the weighted mean. (Before 1.0.0 a
    `Y.weights` given here was silently replaced by the fold mask.)

  - `cores` (`getOption("mc.cores", 1L)`) to evaluate the rank x fold
    grid in parallel; results are identical to the sequential run for
    any `cores` (PSOCK cluster on Windows, forking elsewhere).

  `NA` entries of \\Y\\ are never held out and stay masked in every fit
  (before 1.0.0 they were fitted as zeros and made `objfunc` `NA`).

## Value

A list with `objfunc`, `sigma`, `objfunc.fold`, `folds`, `Q.grid`,
`type`, and

- `r.squared.cv`:

  Held-out \\R^2\\ per rank, pooled over the folds: \\1 - \sum (y - \hat
  y)^2 / \sum (y - \bar y)^2\\ over all held-out upper-triangle entries
  (weighted when `Y.weights` is given). Unlike the `r.squared` of a fit
  (a squared correlation), it can be negative.

- `pred`:

  With `pred = TRUE`: \\N \times N \times\\ (ranks) array of held-out
  predictions, mirrored to the lower triangle; `NA` where an entry was
  never held out.

- `diag.exclude`:

  `TRUE`, present only when used.

- `seeds`, `objfunc.rep`, `objfunc.sd`, `r.squared.cv.rep`,
  `r.squared.cv.sd`, `folds.rep`:

  Present only with two or more `seeds`: the per-seed values (ranks x
  seeds), their standard deviations, and each seed's folds.
  `objfunc.fold`, `folds` and `pred` are those of the first seed.

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)

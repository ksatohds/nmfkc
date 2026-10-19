# Choose the number of mixture components K for NMF-GMM

This function is **experimental and still under development**. The
interface may change in future versions: argument names, defaults and
the contents of the returned object are not yet stable. Code written
against it today may need adjusting after an update. The rest of the
package does not carry this caveat.

`nmf.gmm.select` fits
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md) over
a vector of \\K\\ values and reports the log-likelihood, BIC and ICL for
each, selecting `K.best` by BIC (`K.best.icl` by ICL). If a vector of
known labels is supplied via `truth`, the adjusted Rand index of the
hard clustering is added.

## Usage

``` r
nmf.gmm.select(Y, A = NULL, rank, K = 1:5, ...)
```

## Arguments

- Y, A, rank:

  As in
  [`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md).

- K:

  Integer vector of candidate component counts. Default `1:5`.

- ...:

  Additional arguments passed to
  [`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md)
  (e.g. `cov`, `X.init`, `nstart`, `maxit`, `seed`). Also accepts
  `truth` (a length-N vector of known class labels for the ARI column),
  `verbose` (logical, print the table; default `TRUE`), and `cores`
  (evaluate the `K` candidates in parallel; default
  `getOption("mc.cores", 1L)`). Parallelism uses a PSOCK cluster on
  Windows and forking elsewhere; because each `K` is an independent
  self-seeded fit and results are returned in order, the table and the
  selected `K` are identical for any `cores`.

## Value

An object of class `"nmf.gmm.select"`: a list with the `table` (data
frame over `K`), `K.best` (min BIC), `K.best.icl` (min ICL) and `fits`
(the fitted objects).

## See also

[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(abs(rnorm(20 * 80)) + 1, 20, 80); A <- rbind(1, rnorm(80))
sel <- nmf.gmm.select(Y, A, rank = 2, K = 1:4, verbose = FALSE)
sel$K.best
#> [1] 2
# }
```

# Two-stage (adjust-then-cluster) baseline for NMF-GMM

This function is **experimental and still under development**. The
interface may change in future versions: argument names, defaults and
the contents of the returned object are not yet stable. Code written
against it today may need adjusting after an update. The rest of the
package does not carry this caveat.

`nmf.gmm.twostage` runs the *two-stage* route that
[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md) is
designed to improve on, as a matched baseline: (1) estimate
least-squares scores on the initial basis, (2) regress the scores on the
covariates *blind to the class* and keep the residuals, (3) reconstitute
the residuals in observation space, shift them to non-negativity, and
(4) refit an intercept-only `nmf.gmm` from the *same* basis
initialization. Both routes start from the same basis and use the same
mixture family and settings. They differ in the order of adjustment and
clustering and, unless `fixX = TRUE`, in the basis each refits to its
own data: the two-stage route refits it to the reconstituted, shifted
residuals. When the covariate is associated with the class, the
class-blind regression of stage (2) displaces the class means (the
oracle proposition of Satoh 2026); when the covariate is
(near-)mean-independent of the class, that displacement vanishes in the
population, so both routes target the same class means, although the
fitted partitions need not coincide in a finite sample.

## Usage

``` r
nmf.gmm.twostage(Y, A = NULL, rank, K = 1, ...)
```

## Arguments

- Y:

  Data matrix \\Y\\ (P x N).

- A:

  Covariate matrix \\A\\ (R x N) including an intercept row, or a
  one-sided formula evaluated in `data` (as in
  [`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md)).
  An intercept-only `A` is an error: there is nothing to adjust for.

- rank:

  Integer rank \\Q\\ of the basis.

- K:

  Integer number of mixture components.

- ...:

  Additional arguments as in
  [`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md)
  (`cov`, `X.init`, `nstart`, `maxit`, `seed`, `data`, `standardize`,
  ...); they are applied to both stages.

## Value

An object of class `c("nmf.gmm.twostage", "nmf.gmm")`: the stage-2 fit
(all [`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md)
fields and S3 methods apply), plus a `twostage` list with the
non-negativity `shift`, the covariate matrix `A` that was removed, and
the shared basis initialization `X0`.

## See also

[`nmf.gmm`](https://ksatohds.github.io/nmfkc/reference/nmf.gmm.md) (the
joint route this baselines).

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(abs(rnorm(12 * 40)) + 1, 12, 40)
A <- rbind(1, rnorm(40))
ts <- nmf.gmm.twostage(Y, A, rank = 2, K = 2, nstart = 2, maxit = 100)
table(ts$cluster)
#> 
#>  1  2 
#>  6 34 
# }
```

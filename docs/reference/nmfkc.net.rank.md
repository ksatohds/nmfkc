# Rank selection for nmfkc.net (concise diagnostics)

Fits
[`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
across a range of ranks and reports the three rank-selection criteria –
`r.squared`, the effective rank (utilization), and the element-wise CV
error `sigma.ecv` – with the same concise diagnostics plot as
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md).

## Usage

``` r
nmfkc.net.rank(
  Y,
  rank = 1:5,
  type = c("tri", "bi", "signed"),
  detail = c("full", "fast"),
  plot = TRUE,
  ...
)
```

## Arguments

- Y:

  Symmetric (network) observation matrix.

- rank:

  Integer vector of ranks to evaluate.

- type:

  One of `"tri"`, `"bi"`, `"signed"`.

- detail:

  `"full"` (default) also runs element-wise CV (`sigma.ecv`); `"fast"`
  skips it (plots r.squared and eff.rank only, and recommends the
  R-squared elbow).

- plot:

  Logical; draw the diagnostics plot (default `TRUE`).

- ...:

  Passed on to
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  and
  [`nmfkc.net.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.ecv.md)
  (e.g.\\ `nstart`, `maxit`, `nfolds`, `seed`). Also accepts `cores`
  (`getOption("mc.cores", 1L)`) to fit the per-rank models (and the ECV
  grid) in parallel; results are identical for any `cores` (PSOCK
  cluster on Windows, forking elsewhere).

## Value

A list with `rank.best` (ECV minimum, or the R-squared elbow under
`detail = "fast"`) and `criteria` (data frame: `rank`, `effective.rank`,
`effective.rank.ratio`, `r.squared`, `sigma.ecv`).

## References

Roy, O., & Vetterli, M. (2007). The effective rank: A measure of
effective dimensionality. *Proc. EUSIPCO*, 606–610. (`effective.rank`)
Wold, S. (1978). Cross-validatory estimation of the number of components
in factor and principal components models. *Technometrics*, 20(4),
397–405. (`sigma.ecv`)

## See also

[`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md),
[`nmfkc.net.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.ecv.md),
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)

## Examples

``` r
# \donttest{
Y <- matrix(c(0,1,1,0,0,0, 1,0,1,0,0,0, 1,1,0,1,0,0,
              0,0,1,0,1,1, 0,0,0,1,0,1, 0,0,0,1,1,0), 6, 6)
nmfkc.net.rank(Y, rank = 1:3, type = "tri", nstart = 5, nfolds = 3)
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: did not converge in 100 iterations
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Warning: maximum iterations (5000) reached...
#> Note: sample-clustering quality (silhouette / CPCC / dist.cor) is not part of rank selection; compute it from a list of fits with nmf.cluster.criteria().  See ?nmf.cluster.criteria

# }
```

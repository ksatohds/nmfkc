# Rank selection for nmfkc.signed (concise diagnostics)

Fits
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
across a range of ranks and reports `r.squared`, the effective rank, and
the element-wise CV error `sigma.ecv`, with the same concise plot as
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md).

## Usage

``` r
nmfkc.signed.rank(
  Y,
  A,
  rank = 1:5,
  detail = c("full", "fast"),
  plot = TRUE,
  ...
)
```

## Arguments

- Y:

  Observation matrix (may contain negative entries).

- A:

  Covariate matrix (may be signed).

- rank:

  Integer vector of ranks to evaluate.

- detail:

  `"full"` (default) also runs element-wise CV (`sigma.ecv`); `"fast"`
  skips it (plots r.squared and eff.rank only, and recommends the
  R-squared elbow).

- plot:

  Logical; draw the diagnostics plot (default `TRUE`).

- ...:

  Passed on to
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  and `nmfkc.signed.ecv` (e.g.\\ `maxit`, `nfolds`, `seed`).

## Value

A list with `rank.best` and `criteria` (`rank`, `effective.rank`,
`effective.rank.ratio`, `r.squared`, `sigma.ecv`).

## References

Roy, O., & Vetterli, M. (2007). The effective rank: A measure of
effective dimensionality. *Proc. EUSIPCO*, 606–610. (`effective.rank`)
Wold, S. (1978). Cross-validatory estimation of the number of components
in factor and principal components models. *Technometrics*, 20(4),
397–405. (`sigma.ecv`)

## See also

[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md),
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)

# Rank selection for nmfae.signed (paired rank, concise diagnostics)

Fits
[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)
with a **paired** decoder/encoder rank (\\Q = R\\) across a range of
ranks and reports `r.squared`, the effective rank (of the latent
encoding \\H\\), and the element-wise CV error `sigma.ecv`, with the
same concise plot as
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md).
For a full \\(Q, R)\\ grid use
[`nmf.rrr.signed.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.ecv.md).

## Usage

``` r
nmf.rrr.signed.rank(
  Y1,
  Y2 = Y1,
  rank1 = 1:5,
  detail = c("full", "fast"),
  plot = TRUE,
  ...
)
```

## Arguments

- Y1:

  Endogenous matrix (\\P_1 \times N\\); may be signed.

- Y2:

  Exogenous matrix; defaults to `Y1`.

- rank1:

  Integer vector of (paired) ranks to evaluate (both bases use the same
  value). Legacy `Q` accepted via `...`.

- detail:

  `"full"` (default) also runs element-wise CV (`sigma.ecv`); `"fast"`
  skips it (plots r.squared and eff.rank only, and recommends the
  R-squared elbow).

- plot:

  Logical; draw the diagnostics plot (default `TRUE`).

- ...:

  Passed on to
  [`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md)
  and
  [`nmf.rrr.signed.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.ecv.md).
  Also accepts `cores` to evaluate the rank sweep (and the
  element-wise CV) in parallel; default `getOption("mc.cores", 1L)`.
  Each rank is an independent self-seeded fit and results are gathered
  in order, so the output is identical for any `cores`. Rank aliases
  accepted here for backward compatibility: `Q` for `rank1`.

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

[`nmf.rrr.signed`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.md),
[`nmf.rrr.signed.ecv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.signed.ecv.md),
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)

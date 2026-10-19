# Bi-cross-validation for NMF rank selection

Owen & Perry's (2009) bi-cross-validation (BCV) for choosing the NMF
rank. A lightweight CV engine in the spirit of
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md):
it returns the held-out error per rank and nothing more (no plot, no
table) – pass the result to `which.min(sigma)`, or build your own
diagnostics.

Unlike the element-wise CV of
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)
(which holds out scattered *entries* and refits with weights), BCV holds
out a whole **row-block and column-block** simultaneously: the model is
fit only on the retained block \\D\\, and the held-out block \\A\\ is
predicted by folding the held-out rows/columns onto the fixed
\\D\\-factors via non-negative regression (\\\hat A = L_I R_J\\).
Because the held-out rows and columns never enter the fit, there is no
information leakage. Covariates are ignored (plain NMF). The recommended
setting is to leave out roughly half the rows and half the columns
(`nfolds = 2`).

## Usage

``` r
nmfkc.bicv(Y, rank = 1:3, ...)
```

## Arguments

- Y:

  Observation matrix (\\P \times N\\), non-negative.

- rank:

  Integer vector of ranks to evaluate.

- ...:

  Advanced options, rarely needed (defaults in parentheses): `nfolds`
  (`2`), the number of row and column folds (the grid is
  `nfolds x nfolds`; `2` leaves out half the rows / columns, Owen &
  Perry's recommendation); `seed` (`123`, fold-assignment seed); and
  `nnls.maxit` (`100`, multiplicative-update iterations for the fold-in
  non-negative regressions). Any other arguments are passed to
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) for the
  per-block fits (e.g.\\ `maxit`). `cores` (`getOption("mc.cores", 1L)`)
  evaluates the ranks in parallel; results are identical to the
  sequential loop for any `cores` (PSOCK cluster on Windows, forking
  elsewhere).

## Value

A list (cf.\\
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md))
with:

- objfunc:

  Held-out mean squared error for each rank.

- sigma:

  Its square root (RMSE) for each rank.

- rank:

  The evaluated rank vector.

- nfolds:

  The number of folds used.

## Details

Each fold keeps about \\(1 - 1/\text{nfolds})\\ of the rows and columns,
so the retained block \\D\\ must have more than `rank` rows **and**
columns. The largest testable rank is therefore about \\(1 -
1/\text{nfolds})\min(P, N) - 1\\; with `nfolds = 2` this is roughly
\\\min(P, N)/2 - 1\\. Ranks above this return `NA` and trigger a
[`warning`](https://rdrr.io/r/base/warning.html) that names the limit
and the `nfolds` (or
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md))
that would reach the requested ranks. Raising `nfolds` lifts the limit
at the cost of a smaller hold-out and more compute (\\(\text{nfolds} -
1)^2\\ full fits per rank).

## References

A. B. Owen and P. O. Perry (2009). Bi-cross-validation of the SVD and
the nonnegative matrix factorization. *Ann. Appl. Stat.* 3(2):564–594.
[doi:10.1214/08-AOAS227](https://doi.org/10.1214/08-AOAS227) .

## See also

[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md),
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md),
[`nmfkc.consensus`](https://ksatohds.github.io/nmfkc/reference/nmfkc.consensus.md),
[`nmfkc.ard`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ard.md)

## Examples

``` r
# \donttest{
## rank-3 non-negative data; bi-CV needs enough kept rows/cols per
## fold (> rank), so use a matrix with ample dimensions.
set.seed(1)
X <- matrix(abs(rnorm(30 * 3)), 30, 3)
B <- matrix(abs(rnorm(3 * 40)), 3, 40)
bv <- nmfkc.bicv(X %*% B, rank = 1:6)   # nfolds = 2 (Owen & Perry) by default
#> bi-CV: ranks 1,2,3,4,5,6, 2x2 fold grid (Owen-Perry 2009)...
bv$sigma                  # held-out RMSE per rank
#>     rank=1     rank=2     rank=3     rank=4     rank=5     rank=6 
#> 0.44454653 0.31984467 0.09025809 0.06145124 0.05091540 0.03910144 
bv$rank[which.min(bv$sigma)]
#> [1] 6
# }
```

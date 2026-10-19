# Element-wise Cross-Validation for nmfae (Wold's CV)

`nmfae.ecv` performs k-fold element-wise cross-validation by randomly
holding out individual elements of \\Y_1\\, assigning them a weight of 0
via `Y1.weights`, and evaluating the reconstruction error on those
held-out elements.

This method (also known as Wold's CV) is suitable for determining the
optimal rank pair \\(Q, R)\\ in three-layer NMF. Both `rank1` and
`rank2` accept vector inputs. When `rank2 = NULL` (default), `rank2` is
set equal to `rank1` and pairs are evaluated element-wise (i.e., \\(Q_1,
R_1), (Q_2, R_2), \dots\\). When `rank.encoder` is explicitly specified,
all combinations of `rank` and `rank.encoder` are evaluated via
`expand.grid`.

## Usage

``` r
nmf.rrr.ecv(Y1, Y2 = Y1, rank1 = 1:2, rank2 = NULL, ...)
```

## Arguments

- Y1:

  Output matrix \\Y_1\\ (P1 x N).

- Y2:

  Input matrix \\Y_2\\ (P2 x N). Default is `Y1`.

- rank1:

  Integer vector of response-basis ranks to evaluate. Default is `1:2`.

- rank2:

  Integer vector of covariate-basis ranks to evaluate. Default is
  `NULL`, which sets `rank2 = rank1` and evaluates element-wise pairs.

- ...:

  Additional arguments passed to
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md)
  (e.g., `epsilon`, `maxit`). Also accepts: `nfolds` (number of folds,
  default 5; `div` also accepted), `seed` (integer seed, default 123),
  and `cores` (evaluate the \\(Q,R)\\-pair \\\times\\ fold grid in
  parallel; default `getOption("mc.cores", 1L)`, PSOCK on Windows /
  forking elsewhere). Each task is an independent self-seeded fit and
  results are aggregated in order, so the returned object is identical
  for any `cores`. For backward compatibility, `Q` and `R` are accepted
  as aliases for `rank` and `rank.encoder`. `Y1.weights` (a matrix, one
  weight per column, or a scalar) is kept in every fit; cells of weight
  0 are never held out, and the held-out loss is the weighted mean.
  (Before 1.0.0 a `Y1.weights` given here was silently replaced by the
  fold mask.)

  Rank aliases accepted here for backward compatibility: `Q` for
  `rank1`, `R` for `rank2`.

## Value

A list with components:

- objfunc:

  Named numeric vector of mean MSE for each (Q, R) pair.

- sigma:

  Named numeric vector of RMSE (square root of MSE) for each pair.

- objfunc.fold:

  Named list of per-fold MSE vectors for each pair.

- folds:

  List of length `div` containing the held-out element indices for each
  fold.

- QR:

  Data frame with columns `Q` and `R` listing the evaluated pairs.

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)

## Examples

``` r
Y <- t(iris[1:30, 1:4])
# Default: rank2=NULL -> paired rank1=rank2
res <- nmf.rrr.ecv(Y, rank1 = 1:2, nfolds = 3, maxit = 200)
#> Element-wise CV: 2 (Q,R) pairs, 3-fold, 6 tasks...
#>   Q=1, R=1: MSE=0.021319, sigma=0.1460
#>   Q=2, R=2: MSE=0.021287, sigma=0.1459
res$sigma
#>   Q=1,R=1   Q=2,R=2 
#> 0.1460097 0.1459009 
# Explicit rank2: full grid (kept small so the example stays under CRAN's
# 10s budget -- a real sweep would use wider ranges and a larger maxit)
res2 <- nmf.rrr.ecv(Y, rank1 = 1:2, rank2 = 1:2, nfolds = 3, maxit = 200)
#> Element-wise CV: 4 (Q,R) pairs, 3-fold, 12 tasks...
#>   Q=1, R=1: MSE=0.021319, sigma=0.1460
#>   Q=2, R=1: MSE=0.021319, sigma=0.1460
#>   Q=1, R=2: MSE=0.021321, sigma=0.1460
#>   Q=2, R=2: MSE=0.021287, sigma=0.1459
res2$sigma
#>   Q=1,R=1   Q=2,R=1   Q=1,R=2   Q=2,R=2 
#> 0.1460097 0.1460097 0.1460179 0.1459009 
```

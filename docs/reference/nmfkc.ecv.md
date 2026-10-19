# Perform Element-wise Cross-Validation (Wold's CV)

`nmfkc.ecv` performs k-fold cross-validation by randomly holding out
individual elements of the data matrix (element-wise), assigning them a
weight of 0 via `Y.weights`, and evaluating the reconstruction error on
those held-out elements.

This method (also known as Wold's CV) is theoretically robust for
determining the optimal rank (Q) in NMF. This function supports vector
input for `Q`, allowing simultaneous evaluation of multiple ranks on the
same folds.

For symmetric (network) data use
[`nmfkc.net.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.ecv.md),
which creates upper-triangle folds to prevent information leakage
through the symmetric entries \\Y\_{ij} = Y\_{ji}\\. Passing the old
`Y.symmetric` argument here is no longer supported and stops with a
redirect message.

## Usage

``` r
nmfkc.ecv(Y, A = NULL, rank = 1:3, data, ...)
```

## Arguments

- Y:

  Observation matrix, or a formula (see
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) for
  Formula Mode).

- A:

  Covariate matrix. Ignored when `Y` is a formula.

- rank:

  Vector of ranks to evaluate (e.g., 1:5). For backward compatibility,
  `Q` is accepted via `...`.

- data:

  A data frame (required when `Y` is a formula with column names).

- ...:

  Additional arguments passed to
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) (e.g.,
  `method="EU"`). Also accepts: `nfolds` (number of folds, default 5;
  `div` also accepted), `seed` (integer seed, default 123), and `cores`
  (evaluate the rank x fold task grid in parallel; default
  `getOption("mc.cores", 1L)`). Parallelism uses a PSOCK cluster on
  Windows and forking elsewhere; because each task is a deterministic
  self-seeded fit and results are returned in order, `objfunc`/`sigma`
  are identical for any `cores`. `Y.weights` (as in
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md): a P x
  N matrix, one weight per column, or a scalar) is kept in every fit;
  cells of weight 0 are never held out, and the held-out loss is the
  weighted mean. (Before 1.0.0 a `Y.weights` given here was silently
  replaced by the fold mask.)

## Value

A list with components:

- objfunc:

  Numeric vector containing the Mean Squared Error (MSE) for each Q,
  named `"Q=2"`, `"Q=3"`, .... The names travel with the values:
  [`unname()`](https://rdrr.io/r/base/unname.html) them before combining
  results into a data frame or with
  [`c()`](https://rdrr.io/r/base/c.html) if they are not wanted.

- sigma:

  Numeric vector containing the Residual Standard Error (RMSE) for
  each Q. Only available if method="EU".

- objfunc.fold:

  List of length equal to Q vector. Each element contains the MSE values
  for the k folds.

- folds:

  A list of length `div`, containing the linear indices of held-out
  elements for each fold (shared across all Q).

## References

Wold, S. (1978). Cross-validatory estimation of the number of components
in factor and principal components models. *Technometrics*, 20(4),
397–405.
[doi:10.1080/00401706.1978.10489693](https://doi.org/10.1080/00401706.1978.10489693)
Owen, A. B., & Perry, P. O. (2009). Bi-cross-validation of the SVD and
the nonnegative matrix factorization. *Ann. Appl. Stat.* 3(2), 564–594.
[doi:10.1214/08-AOAS227](https://doi.org/10.1214/08-AOAS227)
(cross-validation of the NMF rank; see also
[`nmfkc.bicv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.bicv.md)).

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md);
other rank-selection criteria:
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md),
[`nmfkc.bicv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.bicv.md),
[`nmfkc.consensus`](https://ksatohds.github.io/nmfkc/reference/nmfkc.consensus.md),
[`nmfkc.ard`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ard.md).

## Examples

``` r
# Element-wise CV to select rank
Y <- t(iris[1:30, 1:4])
res <- nmfkc.ecv(Y, rank = 1:2, nfolds = 3)
#> Performing Element-wise CV for Q = 1,2 (3-fold)...
res$objfunc
#>       Q=1       Q=2 
#> 0.2597611 0.4376291 
```

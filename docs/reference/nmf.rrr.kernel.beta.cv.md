# Optimize kernel beta for nmfae by cross-validation

`nmfae.kernel.beta.cv` selects the optimal `beta` parameter of the
kernel function by evaluating
[`nmf.rrr.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.cv.md)
for each candidate value. The kernel matrix \\A = K(U, V; \beta)\\
replaces \\Y_2\\ in the three-layer NMF model.

When `beta = NULL`, candidate values are automatically generated via
[`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md).

## Usage

``` r
nmf.rrr.kernel.beta.cv(
  Y1,
  rank1 = 2,
  rank2 = NULL,
  U,
  V = NULL,
  beta = NULL,
  plot = TRUE,
  ...
)
```

## Arguments

- Y1:

  Output matrix \\Y_1\\ (P1 x N). Non-negative.

- rank1:

  Integer. Rank of the response basis. Default is 2.

- rank2:

  Integer. Rank of the covariate basis. Default (`NULL`) = `rank1`.

- U:

  Covariate matrix \\U\\ (K x M). Rows are features, columns are samples
  (or knot points for non-symmetric kernels).

- V:

  Covariate matrix \\V\\ (K x N). If `NULL` (default), `V = U` and a
  symmetric kernel is used.

- beta:

  Numeric vector of candidate beta values. If `NULL`, automatically
  determined via
  [`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md).

- plot:

  Logical. If `TRUE` (default), plots the objective function curve.

- ...:

  Additional arguments. Kernel-specific args (`kernel`, `degree`) are
  passed to
  [`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md);
  all others (`div`, `seed`, `shuffle`, `epsilon`, `maxit`, etc.) are
  passed to
  [`nmf.rrr.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.cv.md).
  Also accepts `cores` to evaluate the candidate `beta` values in
  parallel (default `getOption("mc.cores", 1L)`); each inner CV then
  runs sequentially, and because results are gathered in order the
  selected `beta` is identical for any `cores`. For backward
  compatibility, `Q` and `R` are accepted as aliases for `rank` and
  `rank.encoder`.

  Rank aliases accepted here for backward compatibility: `Q` for
  `rank1`, `R` for `rank2`.

## Value

A list with components:

- beta:

  The beta value that minimizes the cross-validation objective.

- objfunc:

  Named numeric vector of objective function values for each candidate
  beta.

## See also

[`nmf.rrr.cv`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.cv.md),
[`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md),
[`nmfkc.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.cv.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
U <- matrix(cars$speed, nrow = 1)
res <- nmf.rrr.kernel.beta.cv(Y, rank1 = 1, rank2 = 1, U = U,
                             beta = c(0.01, 0.02, 0.05), nfolds = 5)
#> beta=0.01...
#> 0.1sec
#> beta=0.02...
#> 0.1sec
#> beta=0.05...
#> 0.1sec

res$beta
#> [1] 0.05
```

# Block-wise Gram accumulation of positive random features for large N

The
[`nmfkc.rff.positive`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.md)
counterpart of
[`nmfkc.signed.rff.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.gram.md):
walks the columns of `U` in blocks, generates the non-negative feature
block \\\Phi_b\\ with the same random map, accumulates \\S =
\Phi\Phi^\top\\ (\\D \times D\\) and \\G_0 = Y\Phi^\top\\ (\\P \times
D\\), and discards the block, so the \\D \times N\\ feature matrix never
exists in memory. Because the features are non-negative, the returned
`"nmfkc.gram"` object is accepted by
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) (as well
as by
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)),
and the NMF-LAB membership interpretation of \\B = \Theta\Phi\\ is
retained.

## Usage

``` r
nmfkc.rff.positive.gram(
  Y,
  U,
  beta = NULL,
  D = NULL,
  seed = NULL,
  block.size = NULL,
  ...
)
```

## Arguments

- Y:

  Non-negative \\P \times N\\ response matrix (e.g. one-hot labels).
  `NA` is not allowed.

- U:

  A \\p \times N\\ numeric matrix of inputs; columns are samples.

- beta:

  Positive scalar. Gaussian kernel bandwidth. Ignored when `pars` is
  supplied.

- D:

  Integer. Number of random features. Required (no \\N/2\\ default
  here). Ignored when `pars` is supplied.

- seed:

  Optional integer seed for \\\omega\\; the caller's stream is restored.
  Ignored when `pars` is supplied.

- block.size:

  Integer. Columns of `U` per block. Default: about 256 MB of features
  per block (`2^25 / D` columns), capped at `N`.

- ...:

  Hidden options `pars` (reuse a random map from
  [`nmfkc.rff.positive`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.md)
  or a previous call) and `hyperbolic` (see
  [`nmfkc.rff.positive`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.md)).

## Value

An object of class `"nmfkc.gram"` with `S`, `G0`, `N`, `D`,
`type = "prf"`, `signed = FALSE`, `pars`, `beta`, `block.size`,
`n.blocks`, `A.block` (block generator) and `rownames = NULL`. Features
for new data: `nmfkc.rff.positive(U.new, pars = g$pars)$Z`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Choromanski, K. et al. (2021). Rethinking attention with Performers.
*ICLR*. arXiv:2009.14794.

## See also

[`nmfkc.rff.positive`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.md),
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.kernel.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.gram.md)
(Nystr\\om),
[`nmfkc.signed.rff.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.gram.md)

## Examples

``` r
# \donttest{
data(iris)
set.seed(1)
idx <- sample(nrow(iris), 100)
mn <- colMeans(iris[idx, 1:4]); sc <- apply(iris[idx, 1:4], 2, sd)
U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))
U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))
levs    <- levels(iris$Species)
Y.train <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
rownames(Y.train) <- levs
beta <- nmfkc.kernel.beta.nearest.med(U.train)$beta

g <- nmfkc.rff.positive.gram(Y.train, U.train, beta = beta, D = 200,
                             seed = 1, block.size = 25)
g
#> Gram summary of positive random features (non-negative) for nmfkc()
#>   N (samples):        100
#>   D (features):       200
#>   Y rows (P):         3
#>   beta (bandwidth):   5.778
#>   blocks:             4 x 25 columns
#>   S  = A A^T:         200 x 200
#>   G0 = Y A^T:         3 x 200
#> Pass as `A` to nmfkc(); use `$pars` with nmfkc.rff.positive() for new data.
res <- nmfkc(Y.train, A = g, rank = 3, verbose = FALSE)   # standard NMF-LAB
Z.test <- nmfkc.rff.positive(U.test, pars = g$pars)$Z
pred   <- predict(res, newA = Z.test, type = "class")
mean(pred == as.character(iris$Species[-idx]))
#> [1] 0.38
# }
```

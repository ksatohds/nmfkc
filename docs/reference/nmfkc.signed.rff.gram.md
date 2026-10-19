# Block-wise Gram accumulation of Random Fourier Features for large N

Builds the two summary matrices that
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
needs from a Random Fourier Feature (RFF) covariate matrix \\Z\\ (\\D
\times N\\) without ever forming \\Z\\ itself: \$\$S = Z Z^\top \quad (D
\times D), \qquad G_0 = Y Z^\top \quad (Q\_{\mathrm{obs}} \times D).\$\$
The input `U` is processed in column blocks; for each block the features
\\Z_b = \sqrt{2/D}\cos(\omega U_b + b)\\ are generated with
[`nmfkc.signed.rff`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md),
added into \\S\\ and \\G_0\\, and discarded. Peak memory is therefore
\\O(D^2 + Q\_{\mathrm{obs}} D + D \cdot \mathtt{block.size})\\ instead
of \\O(DN)\\, which is what makes the RFF route usable when \\N\\ is in
the hundreds of thousands and \\D\\ exceeds the input dimension \\p\\ by
a large factor (the \\D \times N\\ matrix is \\D/p\\ times the size of
the data itself).

The returned object is passed to
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
as its `A` argument. The multiplicative updates are unchanged, since
they only ever use \\S\\ and \\G_0\\; the fit equals the one obtained
from the explicit matrix `nmfkc.signed.rff(U, ...)$Z` with
`warm.start = FALSE` up to floating-point summation order.

## Usage

``` r
nmfkc.signed.rff.gram(
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

  Real-valued \\Q\_{\mathrm{obs}} \times N\\ response matrix (e.g. a
  one-hot label matrix). `NA` is not allowed.

- U:

  A \\p \times N\\ numeric matrix of inputs; columns are samples.

- beta:

  Positive scalar. Gaussian kernel bandwidth \\k(u, u') = \exp(-\beta
  \lVert u - u' \rVert^2)\\. Can be obtained with
  [`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md),
  which itself works on a random subsample. Ignored when `pars` is
  supplied.

- D:

  Integer. Number of random features. Required (there is no `N/2`
  default here: on the large-\\N\\ problems this function is for, that
  default would be enormous, and the per-iteration cost of
  [`nmfkc.signed()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  is \\O(QD^2)\\). Ignored when `pars` is supplied.

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before generating
  \\\omega, b\\; the caller's random stream is restored afterwards.
  Ignored when `pars` is supplied.

- block.size:

  Integer. Number of columns of `U` per block. Default: as many columns
  as keep one block of features near 256 MB (`2^25 / D` columns), capped
  at `N`.

- ...:

  Hidden option `pars`: a list `list(omega, b, D, beta)` from a previous
  [`nmfkc.signed.rff`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md)
  or `nmfkc.signed.rff.gram` call. When supplied, the same random map is
  reused and `beta`, `D`, `seed` are ignored.

## Value

An object of class `"nmfkc.gram"`: a list with

- `S`:

  \\D \times D\\ matrix \\Z Z^\top\\ (signed, positive semi-definite).

- `G0`:

  \\Q\_{\mathrm{obs}} \times D\\ matrix \\Y Z^\top\\ (signed).

- `N`, `D`:

  Number of samples and of features.

- `pars`:

  The RFF parameters `list(omega, b, D, beta)`. Pass them to
  [`nmfkc.signed.rff`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md)
  to generate features for new data, e.g. for
  [`predict.nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/predict.nmfkc.signed.md).

- `block.size`, `n.blocks`:

  Blocking actually used.

- `A.block`:

  A function `A.block(idx)` returning the \\D \times \|\mathtt{idx}\|\\
  feature block for the columns `idx` of `U`.
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  uses it to rebuild \\B = CZ\\ block by block after the fit. It closes
  over `U` and `pars`, so the object keeps `U` alive (no copy is made).

- `rownames`:

  Always `NULL` (RFF features are unnamed); present so that
  [`nmfkc.signed()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)
  can treat the object like a matrix.

## Choosing the bandwidth

With Gram input the fold-based
[`nmfkc.signed.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.cv.md)
is not available (it must split the columns of \\A\\). Select \\\beta\\
as in the NMF-LAB paper: build one Gram object per candidate on the
training columns, fit, and score the validation columns with
[`predict()`](https://rdrr.io/r/stats/predict.html) on features from
`nmfkc.signed.rff(U.val, pars = g$pars)`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
machines. *Advances in NIPS*, 20.

Satoh, K. (2026). Applying non-negative matrix factorization with
covariates to label matrix for classification. *Japanese Journal of
Statistics and Data Science*.
[doi:10.1007/s42081-026-00349-x](https://doi.org/10.1007/s42081-026-00349-x)

## See also

[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md),
[`nmfkc.signed.rff`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md),
[`predict.nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/predict.nmfkc.signed.md),
[`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md)

## Examples

``` r
# \donttest{
## Iris, 3 classes: the Gram route gives the same fit as the explicit
## feature matrix, without ever holding the D x N matrix.
data(iris)
set.seed(1)
idx <- sample(nrow(iris), 100)
mn <- colMeans(iris[idx, 1:4]); sc <- apply(iris[idx, 1:4], 2, sd)
U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))   # 4 x 100
U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))   # 4 x  50
levs    <- levels(iris$Species)
Y.train <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
rownames(Y.train) <- levs
beta <- nmfkc.kernel.beta.nearest.med(U.train)$beta

## Accumulate S, G0 over blocks of 25 columns (4 blocks here)
g <- nmfkc.signed.rff.gram(Y.train, U.train, beta = beta, D = 50,
                           seed = 1, block.size = 25)
g
#> Gram summary of Random Fourier Features (signed) for nmfkc.signed()
#>   N (samples):        100
#>   D (features):       50
#>   Y rows (P):         3
#>   beta (bandwidth):   5.778
#>   blocks:             4 x 25 columns
#>   S  = A A^T:         50 x 50
#>   G0 = Y A^T:         3 x 50
#> Pass as `A` to nmfkc.signed(); use `$pars` with nmfkc.signed.rff() for new data.
res.g <- nmfkc.signed(Y.train, A = g, rank = 3, verbose = FALSE)

## Same random map, explicit matrix, no warm-start: same fit
Z <- nmfkc.signed.rff(U.train, pars = g$pars)$Z              # 50 x 100
res.m <- nmfkc.signed(Y.train, A = Z, rank = 3, warm.start = FALSE,
                      verbose = FALSE)
all.equal(res.g$C, res.m$C)
#> [1] TRUE

## Test-set prediction: regenerate features with the stored pars
Z.test <- nmfkc.signed.rff(U.test, pars = g$pars)$Z
pred   <- predict(res.g, newA = Z.test, type = "class")
mean(pred == as.character(iris$Species[-idx]))
#> [1] 0.78
# }
```

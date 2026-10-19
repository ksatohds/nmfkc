# Block-wise Gram accumulation of Nystroem kernel covariates for large N

Builds the two summary matrices that
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) needs
from the Nystr\\om kernel covariate matrix \\A = C^\top\\ (\\M \times
N\\, \\C\_{nm} = k(\bm u_n, \bm v_m)\\) without ever forming \\A\\
itself: \$\$S = A A^\top \quad (M \times M), \qquad G_0 = Y A^\top \quad
(P \times M).\$\$ The input `U` is processed in column blocks; for each
block the kernel values \\k(V, U_b)\\ are computed with
[`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md),
added into \\S\\ and \\G_0\\, and discarded. Peak memory is \\O(M^2 +
PM + M \cdot \mathtt{block.size})\\ instead of \\O(MN)\\.

The returned object is passed to
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) as its
`A` argument. The multiplicative updates are unchanged, since they only
ever use \\S\\ and \\G_0\\; the fit equals the one obtained from the
explicit matrix `nmfkc.kernel(V, U, beta = beta)` up to floating-point
summation order. This is the large-\\N\\ form of the NMF-LAB Kernel
design with Nystr\\om approximation (Satoh 2026).

## Usage

``` r
nmfkc.kernel.gram(Y, U, V, beta = NULL, block.size = NULL, ...)
```

## Arguments

- Y:

  Non-negative \\P \times N\\ response matrix (e.g. a one-hot label
  matrix). `NA` is not allowed.

- U:

  A \\p \times N\\ numeric matrix of inputs; columns are samples.

- V:

  Landmarks. Either a \\p \times M\\ numeric matrix of landmark points,
  or a single integer \\M\\, in which case \\M\\ landmarks are chosen
  from a random subsample of the columns of `U` by the method in
  `landmarks` (default k-means++ seeding followed by \\k\\-means
  refinement).

- beta:

  Positive scalar. Gaussian kernel bandwidth \\k(u, v) = \exp(-\beta
  \lVert u - v \rVert^2)\\. If `NULL` (default), the nearest-landmark
  median heuristic
  [`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md)`(U_sub, Uk = V)`
  on the subsample is used, and the value is reported in a message. For
  bandwidth selection build one Gram object per candidate in
  `$beta_candidates` of that helper (see Details).

- block.size:

  Integer. Number of columns of `U` per block. Default: as many columns
  as keep one kernel block near 256 MB (`2^25 / M` columns), capped at
  `N`.

- ...:

  Hidden options:

  - `landmarks`: how to choose \\M\\ landmarks when `V` is an integer:
    `"kmeans++"` (default; Arthur & Vassilvitskii 2007 seeding then
    Lloyd refinement), `"kmeans"`
    ([`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html) with
    `nstart = 5`), or `"random"` (uniform subsample).

  - `sample.size`: size of the random subsample of columns used for
    landmark selection and the bandwidth heuristic (default
    `min(N, 10000)`).

  - `seed`: integer seed for the subsample and the landmark selection
    (default 123; the caller's random stream is restored).

  - `kernel`: kernel name passed to
    [`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md)
    (default `"Gaussian"`); further kernel parameters are passed through
    as well.

## Value

An object of class `"nmfkc.gram"`: a list with

- `S`:

  \\M \times M\\ matrix \\A A^\top\\ (non-negative, positive
  semi-definite).

- `G0`:

  \\P \times M\\ matrix \\Y A^\top\\.

- `N`, `D`:

  Number of samples and of covariates (`D = M`, the number of
  landmarks).

- `type`, `signed`:

  `"nystrom"` and `FALSE` (non-negative features; accepted by
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) and
  [`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md)).

- `landmarks`:

  The \\p \times M\\ landmark matrix \\V\\. Pass it to
  [`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md)
  as its first argument to build covariates for new data.

- `beta`, `kernel`:

  Bandwidth and kernel used.

- `landmarks.method`, `sample.size`:

  How the landmarks were chosen (`"supplied"` when `V` was a matrix).

- `block.size`, `n.blocks`:

  Blocking actually used.

- `A.block`:

  A function `A.block(idx)` returning the \\M \times \|\mathtt{idx}\|\\
  kernel block \\k(V, U\_{idx})\\.
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) uses it
  to rebuild \\B = CA\\ block by block after the fit. It closes over `U`
  and `V`, so the object keeps `U` alive (no copy is made).

- `rownames`:

  Row names of \\A\\ (`colnames(V)` if any, else `NULL`).

## Details

With Gram input the fold-based
[`nmfkc.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.cv.md) is
not available (it must split the columns of \\A\\). Select \\\beta\\ as
in the NMF-LAB paper: fix the landmarks, build one Gram object per
candidate \\\beta\\ on the training columns (pass the landmark matrix
`g$landmarks` as `V` so all candidates share it), fit with
[`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md), and
score validation columns with
`predict(fit, newA = nmfkc.kernel(g$landmarks, U.val, beta = beta))`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Williams, C. K. I., & Seeger, M. (2001). Using the Nystr\\om method to
speed up kernel machines. *Advances in NIPS*, 13.

Zhang, K., Tsang, I. W., & Kwok, J. T. (2008). Improved Nystr\\om
low-rank approximation and error analysis. *ICML*.

Arthur, D., & Vassilvitskii, S. (2007). k-means++: the advantages of
careful seeding. *SODA*.

Satoh, K. (2026). Applying non-negative matrix factorization with
covariates to label matrix for classification. *Japanese Journal of
Statistics and Data Science*.
[doi:10.1007/s42081-026-00349-x](https://doi.org/10.1007/s42081-026-00349-x)

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md),
[`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md),
[`nmfkc.signed.rff.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.gram.md)
(the signed, Random-Fourier counterpart for
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md))

## Examples

``` r
# \donttest{
## Iris, 3 classes: Nystroem covariates with 12 k-means++ landmarks.
## The Gram route gives the same fit as the explicit M x N matrix.
data(iris)
set.seed(1)
idx <- sample(nrow(iris), 100)
mn <- colMeans(iris[idx, 1:4]); sc <- apply(iris[idx, 1:4], 2, sd)
U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))   # 4 x 100
U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))   # 4 x  50
levs    <- levels(iris$Species)
Y.train <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
rownames(Y.train) <- levs

## 12 landmarks by k-means++, bandwidth by the nearest-landmark median
## heuristic, S and G0 accumulated over blocks of 25 columns
g <- nmfkc.kernel.gram(Y.train, U.train, V = 12, block.size = 25, seed = 1)
#> nmfkc.kernel.gram: beta = 2.322 (nearest-landmark median heuristic on 100 samples)
g
#> Gram summary of Nystroem kernel covariates (non-negative) for nmfkc()
#>   N (samples):        100
#>   D (landmarks):      12
#>   Y rows (P):         3
#>   beta (bandwidth):   2.322
#>   landmarks:          kmeans++ on a subsample of 100
#>   blocks:             4 x 25 columns
#>   S  = A A^T:         12 x 12
#>   G0 = Y A^T:         3 x 12
#> Pass as `A` to nmfkc(); use nmfkc.kernel(g$landmarks, U.new, beta = g$beta) for new data.
res.g <- nmfkc(Y.train, A = g, rank = 3, verbose = FALSE)

## Same landmarks and beta, explicit 12 x 100 matrix: same fit
A <- nmfkc.kernel(g$landmarks, U.train, beta = g$beta)
res.m <- nmfkc(Y.train, A = A, rank = 3, verbose = FALSE)
all.equal(res.g$C, res.m$C)
#> [1] TRUE

## Test-set prediction: kernel values against the stored landmarks
A.test <- nmfkc.kernel(g$landmarks, U.test, beta = g$beta)
pred   <- predict(res.g, newA = A.test, type = "class")
mean(pred == as.character(iris$Species[-idx]))
#> [1] 0.96
# }
```

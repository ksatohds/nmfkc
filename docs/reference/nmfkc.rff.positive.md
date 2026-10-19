# Positive (non-negative) random features for the Gaussian kernel

Generates \\\omega_d \sim \mathcal{N}(0, I_p)\\ and applies the
non-negative feature map \$\$\phi_d(u) =
\frac{1}{\sqrt{D}}\exp\\\bigl(\sqrt{2\beta}\\\omega_d^\top u -
2\beta\lVert u\rVert^2 + \kappa\bigr) \\\ge 0,\$\$ for which
\\\mathbb{E}\[\phi(u)^\top\phi(u')\] = e^{2\kappa}\exp(-\beta\lVert
u-u'\rVert^2)\\: the Gaussian kernel with bandwidth `beta`, times a
constant. The \\D \times N\\ matrix \\\Phi\\ therefore satisfies
\\\Phi^\top\Phi \approx e^{2\kappa} K\\ *and is non-negative*, so it can
be used as the covariate matrix of
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) directly
– unlike the cosine Random Fourier Features of
[`nmfkc.signed.rff`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md),
which take both signs and need
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md).
This is the "positive random features" construction of Choromanski et
al. (2021), stated there for the softmax kernel, transported to the
Gaussian kernel by rescaling.

The constant \\\kappa\\ is minus the largest exponent
\\\sqrt{2\beta}\\\omega_d^\top u_n - 2\beta\lVert u_n\rVert^2\\ observed
on the training `U` (stored in `pars`), so every training feature is at
most \\1/\sqrt{D}\\ and [`exp()`](https://rdrr.io/r/base/Log.html)
cannot overflow (points far from the others underflow to zero at large
\\\beta\\, which only means that bandwidth fits poorly). It scales the
kernel estimate by \\e^{2\kappa}\\, which \\\Theta\\ absorbs, so the
fitted model and its predictions do not depend on it. Note that
\\\kappa\\ depends on the realized \\\omega\\: the stabilized features
are an unbiased estimator of the randomly rescaled kernel
\\e^{2\kappa}k\\, not of \\k\\ itself; the unbiasedness statement above
is for the unstabilized map. Inputs should nevertheless be centred and
scaled: the estimator's variance grows with \\\beta\lVert u -
u'\rVert^2\\, so positive features need a larger `D` than cosine
features for the same accuracy.

## Usage

``` r
nmfkc.rff.positive(U, beta = NULL, D = ceiling(ncol(U)/2), seed = NULL, ...)
```

## Arguments

- U:

  A \\p \times N\\ numeric matrix; columns are data points.

- beta:

  Positive scalar. Gaussian kernel bandwidth parameter \\k(u,u') =
  \exp(-\beta\lVert u-u'\rVert^2)\\. Can be obtained via
  [`nmfkc.kernel.beta.nearest.med`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.nearest.med.md).
  Ignored when `pars` is supplied.

- D:

  Integer. Number of random features (default `ceiling(ncol(U) / 2)`;
  choose explicitly for large \\N\\). With `hyperbolic = TRUE` the \\D\\
  features are \\D/2\\ antithetic pairs. Ignored when `pars` is
  supplied.

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing
  \\\omega\\; the caller's random stream is restored.

- ...:

  Hidden options:

  - `pars`: a list from a previous call; when supplied the same random
    map (and the same \\\kappa\\) is reused and `beta`, `D`, `seed` are
    ignored. Use this for test data.

  - `hyperbolic`: logical (default `FALSE`). If `TRUE`, each
    \\\omega_d\\ is paired with \\-\omega_d\\ (the "hyperbolic" /
    antithetic variant of Choromanski et al.).

## Value

A list with `Z`, the non-negative \\D \times N\\ feature matrix, and
`pars = list(omega, D, beta, kappa, hyperbolic, type = "positive")`.
Pass `Z` to
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) (or
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md))
as `A`, and `pars` back to this function to build features for new data.

## Lifecycle

This function is **experimental**. The interface may change in future
versions.

## References

Choromanski, K., Likhosherstov, V., Dohan, D., Song, X., Gane, A.,
Sarlos, T., Hawkins, P., Davis, J., Mohiuddin, A., Kaiser, L., Belanger,
D., Colwell, L., & Weller, A. (2021). Rethinking attention with
Performers. *ICLR*. arXiv:2009.14794.

Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
machines. *Advances in NIPS*, 20.

## See also

[`nmfkc.rff.positive.gram`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rff.positive.gram.md)
(large \\N\\),
[`nmfkc.signed.rff`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.rff.md)
(signed cosine features),
[`nmfkc.kernel`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.md),
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)

## Examples

``` r
set.seed(1)
U <- scale(matrix(stats::rnorm(3 * 40), 3, 40))        # centred inputs
pr <- nmfkc.rff.positive(U, beta = 0.5, D = 20000, seed = 1)
all(pr$Z >= 0)
#> [1] TRUE
## Phi' Phi / exp(2 kappa) approximates the Gaussian kernel
K  <- nmfkc.kernel(U, beta = 0.5)
Kp <- crossprod(pr$Z) / exp(2 * pr$pars$kappa)
max(abs(Kp - K))
#> [1] 0.2354917

## The same random map on new data
U.new <- matrix(stats::rnorm(3 * 5), 3, 5)
dim(nmfkc.rff.positive(U.new, pars = pr$pars)$Z)      # 20000 x 5
#> [1] 20000     5
```

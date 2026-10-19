# Estimate Gaussian/RBF kernel parameter beta from covariates (supports landmarks)

Computes a data-driven reference scale for the Gaussian/RBF kernel from
covariates using a robust "median nearest-neighbor (or nearest-landmark)
distance" heuristic, and returns the corresponding kernel parameter
\\\beta\\.

The Gaussian/RBF kernel is assumed to be written in the form \$\$k(u,v)
= \exp\\-\beta \\u-v\\^2\\ = \exp\\-\\u-v\\^2/(2\sigma^2)\\,\$\$ hence
\\\beta = 1/(2\sigma^2)\\. This function first estimates a typical
distance scale \\\sigma_0\\ by the median of distances, then sets
\\\beta_0 = 1/(2\sigma_0^2)\\.

If `Uk` is `NULL`, \\\sigma_0\\ is estimated as the median of
nearest-neighbor distances within `U` (excluding self-distance). If `Uk`
is provided, \\\sigma_0\\ is estimated as the median of nearest-landmark
distances from each sample in `U` to its closest landmark in `Uk`.

To control memory usage for large `N` (and `M`), distances are computed
in blocks. Optionally, columns of `U` can be randomly subsampled via
`sample.size` to reduce cost.

## Usage

``` r
nmfkc.kernel.beta.nearest.med(
  U,
  Uk = NULL,
  block.size = 1000,
  block.size.Uk = 2000,
  sample.size = NULL,
  ...
)
```

## Arguments

- U:

  A numeric matrix of covariates (`K x N`); columns are samples.

- Uk:

  An optional numeric matrix of landmarks (`K x M`); columns are
  landmark points. If provided, distances are computed from samples in
  `U` to landmarks in `Uk`.

- block.size:

  Integer. Number of columns of `U` processed per block when computing
  distances (controls memory usage). If `N <= 1000`, it is automatically
  set to `N`.

- block.size.Uk:

  Integer. Number of columns of `Uk` processed per block when `Uk` is
  not `NULL` (controls memory usage). If `M <= 2000`, it is
  automatically set to `M`.

- sample.size:

  Integer or `NULL`. If not `NULL`, randomly subsamples this many
  columns of `U` (without replacement) before computing distances, to
  reduce computational cost.

- ...:

  Additional arguments. Hidden option `candidates` controls the
  candidate grid: one of `"7points"` (default), `"4points"`, or a
  numeric vector of \\t\\ values. See Details.

## Value

A list with elements:

- `beta`: Estimated kernel parameter \\\beta_0 = 1/(2\sigma_0^2)\\.

- `beta_candidates`: Numeric vector of candidate \\\beta\\ values
  (logarithmic grid) intended for cross-validation.

- `dist_median`: The estimated distance scale \\\sigma_0\\ (median of
  nearest-neighbor or nearest-landmark distances).

- `block.size.used`: The effective block size(s) used. Either a scalar
  (no `Uk`) or a named vector `c(U=..., Uk=...)` when `Uk` is provided.

- `sample.size.used`: The number of columns of `U` actually used (after
  subsampling).

- `uk_is_u`: Logical flag indicating whether `Uk` was detected as
  identical to `U` (only returned when `Uk` is provided).

## Details

**Candidate grid:** Along with `beta`, the function returns
`beta_candidates`, a logarithmic grid suitable for cross-validation. The
grid is symmetric on the bandwidth scale \\\sigma\\ around \\\sigma_0\\:
\$\$\sigma = \sigma_0 \times 10^{t},\$\$ and since \\\beta =
1/(2\sigma^2)\\, this corresponds to \\\beta = \beta_0 \times
10^{-2t}\\.

The grid of \\t\\ values can be customized through the hidden argument
`candidates` (passed via `...`):

- `"7points"` (default): \\t \in \\-1,-2/3,-1/3,0,1/3,2/3,1\\\\ (7
  candidates spanning one decade, matches the grid used in the RFF-NMF
  research memo).

- `"4points"`: \\t \in \\-1/2, 0, 1/2, 1\\\\ yielding \\\beta_0 \times
  10^{(1,0,-1,-2)}\\ (the legacy short grid).

- A numeric vector: user-specified \\t\\ values. The grid returned is
  \\\beta_0 \times 10^{-2t}\\.

Prior to version 0.6.8, the grid depended on whether `Uk` was supplied
(4 candidates for `Uk = NULL`, 7 for supplied `Uk`). The current
implementation unifies both branches via `candidates`.

**Notes:**

- When `Uk` is identical to `U`, the function detects this case and
  excludes self-distances (distance 0) to avoid \\\sigma_0=0\\.

- `sample.size` performs random subsampling without setting a seed. For
  reproducible results, set
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before calling this
  function.

## See also

[`nmfkc.kernel.gaussian`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.gaussian.md),
[`nmfkc.kernel.beta.cv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.kernel.beta.cv.md)

## Examples

``` r
# Basic (nearest-neighbor within U)
U <- matrix(runif(20), nrow = 2)
beta_info <- nmfkc.kernel.beta.nearest.med(U)
beta0 <- beta_info$beta
betas <- beta_info$beta_candidates

# With landmarks (nearest-landmark distances)
Uk <- matrix(runif(10), nrow = 2)
# \donttest{
beta_info2 <- nmfkc.kernel.beta.nearest.med(U, Uk)
# }
```

# Symmetric NMF for networks (tri / bi / signed)

Single entry point for symmetric NMF of network data with correct
multiplicative updates. Three model types are supported via `type`:

- **tri** (`type="tri"`, default): \\Y \approx X C X^\top\\ with \\X, C
  \ge 0\\ (both non-negative; \\C\\ symmetric by design). Uses
  Frobenius-full bilateral gradient.

- **bi** (`type="bi"`): \\Y \approx X X^\top\\ (\\C\\ fixed to \\I_Q\\),
  cube-root damping (He et al. 2011).

- **signed** (`type="signed"`): \\Y \approx X (C\_{+} - C\_{-}) X^\top\\
  with \\X \ge 0\\ and signed \\C = C\_{+} - C\_{-}\\. Preserves the
  soft-clustering interpretation of \\X\\ while allowing negative
  off-diagonals of \\C\\ (inter-cluster *repulsion*).

**Non-negative adjacency matrix assumption.** All three types assume \\Y
\ge 0\\ (a non-negative adjacency/affinity matrix). The qualifier
“signed” in `type = "signed"` refers to the middle coefficient \\C\\,
not to \\Y\\ itself. The underlying Ding, Li & Jordan (2010)
sign-splitting updates require \\Y \ge 0\\ to guarantee monotone
descent; supplying a signed \\Y\\ triggers an error. For a signed data
matrix, see
[`nmfkc.signed`](https://ksatohds.github.io/nmfkc/reference/nmfkc.signed.md).

## Usage

``` r
nmfkc.net(
  Y,
  rank = 2,
  type = c("tri", "bi", "signed"),
  epsilon = 1e-04,
  maxit = 5000,
  verbose = FALSE,
  ...
)
```

## Arguments

- Y:

  Symmetric (N x N) **non-negative** adjacency matrix. `NA` entries are
  automatically treated as masked edges (equivalent to supplying
  `Y.weights` with 0 at those positions); see the note on `Y.weights`
  below.

- rank:

  Integer Q.

- type:

  `"tri"` (default), `"bi"`, or `"signed"`.

- epsilon, maxit, verbose:

  Standard. **Note**: the convergence test uses the *unpenalized*
  squared error, not the penalized objective the updates minimize, so
  with a large `C.L1` or `X.L2.ortho` the iteration can stop while the
  optimized quantity is still moving. The other fitters test the
  penalized value; this is kept as-is so existing fits reproduce.

- ...:

  Hidden options: `nstart` (default 1; see note below), `seed` (default
  123), `X.restriction`, `X.init`, `C.init` (tri only) or
  `Cp.init`/`Cn.init` (signed only), `Y.weights`, `diag.exclude`, `C.L1`
  (tri only), `X.L2.ortho`, `prefix`.

  **`diag.exclude`** (default `FALSE`): when `TRUE` the diagonal of
  \\Y\\ (self-loops) is neither fitted nor counted in the fit statistics
  – for a network without self-loops, whose zero diagonal is not an
  observation. It is the same as a `Y.weights` with a zero diagonal (and
  is combined with a `Y.weights` or an `NA` mask when one is given); the
  fit then carries `diag.exclude = TRUE`.

  **`Y.weights`** is an optional non-negative N x N weight matrix
  (symmetric, same shape as `Y`). When supplied, the loss becomes \\\sum
  W\_{ij} \\ (Y\_{ij} - \hat Y\_{ij})^2\\
  ([`lm`](https://rdrr.io/r/stats/lm.html)()-style, **linear** in
  \\W\\). Logical matrices (`TRUE` / `FALSE`) are also accepted. Typical
  usage by
  [`nmfkc.net.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.ecv.md)
  is a binary mask (\\W \in \\0,1\\\\) holding out test edges on the
  upper triangle; real-valued weights for edge-level importance
  weighting are also supported. If `Y.weights` is `NULL` (default) and
  `Y` contains `NA`, a binary mask is auto-generated (0 at `NA`
  positions, 1 elsewhere), and the `NA` entries in `Y` are replaced by 0
  so the multiplicative updates can proceed.

  **`X.init`** controls the initialization of the N x Q basis matrix
  \\X\\. Accepted values:

  - `"kmeans"` (default): k-means on the rows of \\Y\\ (equivalently
    columns, since \\Y\\ is symmetric); the Q cluster centers become the
    columns of \\X\\. Each node is treated as an N-dimensional
    connectivity profile, so clusters correspond to nodes with similar
    neighborhood structure – essentially a fast proxy for spectral
    clustering (Kuang, Yun & Park 2015, *SymNMF*). Scales well and is
    the recommended default for network data.

  - `"kmeansar"`: `"kmeans"` followed by filling zero entries of \\X\\
    with \\\mathrm{Uniform}(0, \bar Y / 100)\\ to escape trivial
    stationary points.

  - `"nndsvd"`: Non-negative Double SVD with additive randomness
    (NNDSVDar). Requires a full SVD of \\Y\\, so for very large networks
    (N \> a few thousand) `"kmeans"` is preferable.

  - `"runif"`: Uniform random entries in \\\[0, 1\]\\.

  - `"random"`: Legacy default (pre-v0.6.8), equivalent to
    `abs(rnorm(N * Q)) * 0.1`. Kept for backward compatibility.

  - A numeric N x Q matrix supplied by the user (used as-is).

  The anchor options of
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md)
  (`X.init = "spa"`, `X.anchor`) are not available for the symmetric
  model; passing them is an error rather than being ignored. When
  `nstart > 1`, each restart uses a distinct seed so that k-means /
  runif / NNDSVDar produce different candidate initial values across the
  multi-start loop.

  **Multi-start recommendation.** For `type = "signed"` the \\C =
  C\_{+} - C\_{-}\\ bottleneck can take both positive and negative
  values, so the objective has more local minima than for `"tri"` or
  `"bi"`. A larger `nstart` (e.g., 10-50) is recommended during
  exploration to reduce the chance of being trapped at a suboptimal
  stationary point. The default 1 is intended for fast development;
  raise for publication-grade runs.

## Value

Object of class `c("nmfkc.net.<type>", "nmfkc.net", "nmfkc")`. For
`type = "signed"` the return also carries `$Cp, $Cn`.

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmfkc.net.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.ecv.md),
[`nmfkc.net.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.DOT.md)

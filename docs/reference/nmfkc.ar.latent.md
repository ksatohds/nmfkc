# Latent transition matrices of an NMF-VAR model

This function is **experimental**. The interface may change in future
versions: argument names, defaults and the contents of the returned
object are not yet stable.

Reduces a fitted NMF-VAR to the \\Q\times Q\\ operators that drive its
coefficient vector. Writing \\\bm b_t=\Theta\bm a_t\\ for the fitted
scores, the residual-free recursion is \$\$\bm
b_t=\sum\_{d=1}^{D}G_d\\\bm b\_{t-d}+\bm\theta,\qquad G_d=\Theta_d X \\
\\ (Q\times Q),\$\$ and the \\G_d\\ carry the model's whole
deterministic structure — transition weights, propagation, roots, the
long-run mean. They stay readable where the \\P\times P\\ matrices
\\\Xi_d = X\Theta_d\\ do not: \\4\times4\\ instead of \\47\times47\\.

Use
[`plot.nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.ar.latent.md)
to draw them,
[`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md)
to attach standard errors, and
[`nmfkc.ar.stationarity`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.stationarity.md)
if the stationarity verdict is all that is wanted.

## Usage

``` r
nmfkc.ar.latent(x)
```

## Arguments

- x:

  A fitted `"nmfkc"` object obtained from an
  [`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md)
  design (`nmfkc(ar$Y, ar$A, rank = Q)`).

## Value

An object of class `"nmfkc.ar.latent"`: a list with

- G:

  List of \\D\\ latent transition matrices \\G_d\\ (\\Q\times Q\\),
  named `lag1`, ..., with `dimnames` = basis labels (row = effect,
  column = cause).

- G.sum:

  \\\sum_d G_d\\.

- eigenvalues:

  The \\DQ\\ eigenvalues of the latent companion matrix. The \\PD\times
  PD\\ companion of the \\\Xi_d\\ has the same non-zero eigenvalues, so
  these are all of the model's roots.

- spectral.radius:

  \\\rho\\, the largest modulus among them.

- spectral.radius.sum:

  \\\rho(\sum_d G_d)\\; below 1 iff stationary.

- stationary:

  Logical; `TRUE` when `spectral.radius < 1`.

- cycle.period:

  Period \\2\pi/\|\arg\lambda\|\\ implied by the dominant complex root,
  or `NA` when no root is complex.

- alternating:

  `TRUE` when some root is real and negative, i.e. the decay alternates
  in sign (period 2).

- non.oscillatory:

  `TRUE` only when every root is real and positive. Named after the
  roots, not the trajectory — overshoot is still possible.

- theta0:

  Latent intercept \\\bm\theta\\ (or `NULL`).

- mu.b, mu.y:

  Long-run means \\(I-\sum_d G_d)^{-1}\bm\theta\\ and \\X\bm\mu_b\\;
  `NA` when the fit is not stationary or has no intercept.

- p.star:

  Composition of the long-run mean, \\\bm\mu_b/(\bm 1'\bm\mu_b)\\ —
  equivalently the fixed point of the *residual-free* composition
  dynamics. It is **not** in general the limit of \\\bm b_t/(\bm 1'\bm
  b_t)\\ for the stochastic process, nor its expectation: the ratio of
  expectations is not the expectation of the ratio.

- colsum:

  Column sums \\c_j\\ of \\\sum_d\Xi_d = X\sum_d\Theta_d\\, computed as
  `colSums(X) %*% Lambda` (length \\P\\). Being column sums of a
  non-negative matrix whose radius is \\\rho(\sum_d G_d)\\, they bracket
  it. The unweighted `colSums(Lambda)` would only be correct when \\\bm
  1'X=\bm 1'\\.

- colsum.max:

  \\\max_j c_j\\; below 1 it certifies stationarity.

- df:

  Effective dimension \\Q(P+m-Q)\\ of the fit, where \\m\\ is the number
  of columns of \\\Theta\\ (\\PD+1\\ with an intercept, \\PD\\ without).
  These are the degrees of freedom to use when comparing \\R^2\\ or
  information criteria across \\Q\\: the naive count \\Q(P+m)\\ ignores
  the \\Q^2\\ rotational indeterminacy.

- separability:

  Per-basis maximum of the row-normalised \\X\\; a value near 1 means
  that basis owns a variable loading almost only on it. This is a
  *necessary*-looking but not sufficient condition for uniqueness — see
  the counterexample under Identification.

- X, Theta:

  The basis and the full coefficient matrix the quantities above were
  computed from, carried along so that
  [`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md)
  does not need the fitted object a second time.

- X.restriction:

  The constraint
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) placed
  on the columns of \\X\\, which decides how much of \\T\\ is pinned
  down.

- dims:

  Named vector of `P`, `Q`, `D`.

## Details

**The latent process is VARMA, not VAR.** The display above drops the
observation residual. With \\\bm y_t = X\bm b_t + \bm e_t\\ the exact
identity is \$\$\bm b_t=\sum\_{d=1}^{D}G_d\\\bm b\_{t-d}+\bm\theta
+\sum\_{d=1}^{D}\Theta_d\\\bm e\_{t-d},\$\$ i.e. a VARMA\\(D,D)\\ in
general (verified to machine precision; dropping the moving-average term
leaves an error of order the residuals themselves). \\G_d\\ is therefore
the **autoregressive operator of the fitted VAR reduced to the latent
space**, not the coefficient matrix of a latent VAR that the scores obey
exactly. In particular an entry of \\G_d\\ should not be read as Granger
causality between latent conditions: the moving-average term is common
to all coordinates and is not conditioned on.

**Direction convention**: \\(G_d)\_{q q'}\\ is the effect of
\\b\_{q',\\t-d}\\ on \\b\_{q,\\t}\\ — **row = effect** (at \\t\\),
**column = cause** (at \\t-d\\).

**Identification.** The factorization is unique only up to
\\(X,\Theta)\to(XT,T^{-1}\Theta)\\, under which \\G_d\to T^{-1}G_dT\\.
The eigenvalues, \\\rho\\, \\\Xi_d\\ and \\\bm\mu_y\\ are therefore
invariant, while individual entries of \\G_d\\ are not. Column
normalization removes only the scale part of \\T\\, and a separable
fitted \\X\\ does not certify the rest: with \\X=I_2\\,
\\\Theta=\left(\begin{smallmatrix}2&3\\1&1\end{smallmatrix}\right)\\ and
\\T=\left(\begin{smallmatrix}0.9&0.1\\0.1&0.9\end{smallmatrix}\right)\\
the basis is perfectly separable and column-normalized, yet \\XT\ge0\\,
\\T^{-1}\Theta\ge0\\, the product is unchanged and \\T\\ is no
permutation. Read the entries of \\G_d\\ as descriptive and the
invariants as inferential; see
[`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md).

**Oscillation.** A complex pair means a cycle, but a *negative real*
root also oscillates — it alternates in sign each period, i.e. with
period 2. `cycle.period` covers the first case and `alternating` the
second; `non.oscillatory` is `TRUE` only when neither occurs. Positive
real roots rule out oscillation but **not overshoot**: a non-normal
\\G\\ can send a coordinate up before it comes down, which is exactly
what the off-diagonal impulse responses of the Canada fit do. Note also
that for \\Q=2\\, \\D=1\\ and \\G\ge0\\ the roots are real by
construction (the discriminant is \\(a-d)^2+4bc\ge0\\), so in that
setting the absence of a cycle is a property of the design, not a
finding.

## References

Satoh, K. (2025). Applying non-negative matrix factorization with
covariates to multivariate time series data as a vector autoregression
model. *Japanese Journal of Statistics and Data Science*.
arXiv:2501.17446.
[doi:10.1007/s42081-025-00314-0](https://doi.org/10.1007/s42081-025-00314-0)

## See also

[`plot.nmfkc.ar.latent`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.ar.latent.md),
[`nmfkc.ar.latent.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.latent.inference.md),
[`nmfkc.ar.stationarity`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.stationarity.md),
[`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md),
[`nmfkc.ar.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.DOT.md)

## Examples

``` r
d <- AirPassengers
ar_data <- nmfkc.ar(d, degree = 2)
result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#> Y(1,142)~X(1,1)C(1,3)A(3,142)=XB(1,142)...
#> 0sec
lat <- nmfkc.ar.latent(result)
lat$G$lag1          # latent transition matrix of lag 1
#>           Basis1
#> Basis1 0.9376883
lat$eigenvalues     # all roots of the model
#> [1]  0.96846796 -0.03077969
lat                 # full report
#> Latent transition matrices  G_d = Theta_d %*% X   (row = effect at t, column = cause at t-d)
#> 
#> lag1:
#>        Basis1
#> Basis1  0.938
#> 
#> lag2:
#>        Basis1
#> Basis1   0.03
#> 
#> roots: 2 (companion 2x2)   rho = 0.9685  ->  stationary
#> all roots real but one is negative: alternating decay (period 2)
#> 
#> Long-run mean of the latent scores (mu.b):
#>  Basis1 
#> 337.397 
#> Long-run composition (p*):
#> Basis1 
#>      1 
#> 
#> Effective dimension df = Q(P + PD + 1 - Q) = 3
#> For the stationarity verdict see nmfkc.ar.stationarity(), for standard errors nmfkc.ar.latent.inference().
```

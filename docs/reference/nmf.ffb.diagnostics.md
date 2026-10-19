# Structural Diagnostics of a Fitted NMF-FFB Feedback Matrix

Two quantities that describe *what* the selected feedback is, as opposed
to whether it is there (which is the job of the calibrated
likelihood-ratio test in
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)).

**Cycle structure.** The equilibrium operator is \\A = X\Theta_1\\ on
the indicators and \\A_Q = \Theta_1X\\ on the factors, with the same
spectral radius \\\rho\\. The strongly connected components of the
directed graph of \\A_Q\\ are the cycles. \\\rho = 0\\ with a non-empty
support means \\A\\ is nilpotent: the selected feedback is a cascade of
directed cross-factor paths without a closed loop. Since \\A \ge 0\\,
Perron–Frobenius gives a non-negative leading eigenvector, and its
normalized entries say how much of the loop each factor carries.

**Factor-level alignment.** Feedback of the form \\\Theta_1 =
aX^\top\Psi^{-1}\\ enters each factor through the factor scores of the
factors and is observationally equivalent, under the Gaussian working
model, to a change of \\\Theta_2\\ and of the factor covariance
\\\Phi\\; only the departure from that family is identified, through the
uniquenesses \\\Psi\\. On the free set \\F\\ of the exclusion
restriction the family is a subspace of dimension at most \\Q(Q-1)\\,
and `omega` is the share of the squared norm of \\\Theta_1\\ that lies
in it. `omega` has no natural zero: for an isotropic direction its
expectation is `omega0 = dim / n.free`, so it should be read against
`omega0` and against a simulation in which the true feedback is
indicator-level (there it averages about `1.3 * omega0`).

## Usage

``` r
nmf.ffb.diagnostics(object, which = c("selected", "full"))
```

## Arguments

- object:

  A fitted `nmf.ffb` object from
  [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md)
  with `method = "fiml"`.

- which:

  Which feedback matrix to describe: `"selected"` (default, the
  BIC-selected support) or `"full"` (the unpenalized fit).

## Value

A list with `cycles` (`A.factor`, `rho`, `cycle`, `n.cycles`, `perron`)
and `omega` (`omega`, `omega0`, `dim`, `n.free`).

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y2 <- matrix(runif(2 * 200), 2, 200)
X <- matrix(0.02, 6, 2); X[1:3, 1] <- 1; X[4:6, 2] <- 1
X <- sweep(X, 2, colSums(X), "/")
T1 <- matrix(0, 2, 6); T1[1, 5] <- 0.6; T1[2, 2] <- 0.5
L <- solve(diag(6) - X %*% T1)
Y1 <- pmax(L %*% (X %*% (matrix(c(1, 0.5, 0.3, 1), 2, 2) %*% Y2) +
                  matrix(rnorm(6 * 200, 0, 0.05), 6, 200)), 0)
fit <- nmf.ffb(Y1, Y2, Q = 2)
d <- nmf.ffb.diagnostics(fit)
d$cycles$rho          # spectral radius
#> [1] 0.1782386
d$cycles$cycle        # which factors lie on a cycle
#> Factor1 Factor2 
#>       1       1 
c(d$omega$omega, d$omega$omega0)
#> [1] 0.2960767 0.3333333
# }
```

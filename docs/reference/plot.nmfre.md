# Plot convergence diagnostics for NMF models

Plots the objective function value over iterations for `nmfre` and
`nmf.ffb` objects. (For `nmfkc` and `nmfae`, plot methods are defined in
their respective source files.)

## Usage

``` r
# S3 method for class 'nmfre'
plot(x, ...)

# S3 method for class 'nmf.ffb'
plot(x, ..., which = c("penalized", "reconstruction", "both"))
```

## Arguments

- x:

  A fitted model object.

- ...:

  Additional graphical arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

- which:

  For `plot.nmf.ffb`: which objective to plot. One of `"full"` (default;
  `loss + penalties`, the actual monotonically-decreasing quantity that
  the multiplicative updates minimize), `"reconstruction"` (Frobenius
  distance only, \\\\ Y_1 - X B \\\_F^2\\), or `"both"` (overlay both
  with a legend). `"both"` is useful for diagnosing whether
  regularization is actively shaping the solution: if the two curves
  diverge, the penalties are pulling the optimizer away from the pure
  least-squares minimum. Ignored for `nmf.ffb(method = "fiml")` objects,
  which carry no iteration trace: for them the BIC of the L1 path
  (`x$path`) is drawn against `log(lambda1)`: the line follows the
  smallest BIC at each penalty (over the starts), every individual
  proposal is a small grey point, the number of selected paths is
  printed at each point and the BIC-selected penalty is circled.

## Value

Invisible `NULL`.

## See also

[`nmfre`](https://ksatohds.github.io/nmfkc/reference/nmfre.md),
[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(runif(20), nrow = 4)
A <- diag(5)
res <- nmfre(Y, A, rank = 2)
plot(res)

# }
```

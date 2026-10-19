# Compute model selection criteria for a fitted nmfkc model

`nmfkc.criterion` computes the effective rank (raw and as the
\\\[0,1\]\\ broken-stick index) and the sample-clustering quality
measures (silhouette, CPCC, dist.cor) from a fitted `nmfkc` model.

This function can be called on a model that was fitted with
`detail = "fast"` or `detail = "minimal"` to compute the full set of
criteria afterwards.

## Usage

``` r
nmfkc.criterion(object, Y, detail = c("full", "fast", "minimal"), ...)
```

## Arguments

- object:

  An object of class `"nmfkc"` returned by
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md).

- Y:

  The original observation matrix (P x N) used for fitting.

- detail:

  Character string controlling the level of computation: `"full"`
  (default) computes all criteria including silhouette, CPCC and
  dist.cor; `"fast"` skips the expensive distance-based criteria;
  `"minimal"` skips the fit and clustering statistics.

- ...:

  Additional arguments: `Y.weights` (non-negative weight matrix;
  [`lm()`](https://rdrr.io/r/stats/lm.html)-style loss \\\sum W \\
  r^2\\; default: all ones). See
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) for
  full details.

## Value

A list with components:

- r.squared:

  R-squared between Y and XB.

- sigma:

  Residual standard deviation.

- mae:

  Mean absolute error.

- B.prob:

  Column-normalized coefficient matrix (soft-clustering probabilities).

- B.cluster:

  Hard clustering labels (argmax of B.prob per column).

- X.prob:

  Row-normalized basis matrix.

- X.cluster:

  Hard clustering labels per row of X.

- criterion:

  Named list: effective.rank, effective.rank.index, silhouette, CPCC,
  dist.cor.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmfkc.rank`](https://ksatohds.github.io/nmfkc/reference/nmfkc.rank.md)

## Examples

``` r
Y <- t(iris[, -5])
res <- nmfkc(Y, rank = 3, detail = "fast")
#> Y(4,150)~X(4,3)B(3,150)...
#> 0sec
crit <- nmfkc.criterion(res, Y)
crit$criterion$silhouette
#> [1] 0.3440654
```

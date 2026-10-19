# Summary method for nmfae objects

`summary.nmfae` produces a summary of a fitted NMF-AE model, including
dimensions, convergence status, goodness-of-fit statistics, and
structure diagnostics (sparsity of factor matrices).

## Usage

``` r
# S3 method for class 'nmfae'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"nmfae"` returned by
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md).

- ...:

  Additional arguments (currently unused).

## Value

An object of class `"summary.nmfae"`, a list with components:

- call:

  The matched call.

- dims:

  Named vector `c(P1, P2, N)`.

- Q:

  Decoder rank.

- R:

  Encoder rank.

- n.params:

  Total number of parameters (P1*Q + Q*R + R\*P2).

- autoencoder:

  Logical; TRUE if P1 == P2 and Y1 was used as Y2.

- niter:

  Number of iterations.

- runtime:

  Elapsed time.

- objfunc:

  Final objective value.

- r.squared:

  R-squared.

- sigma:

  Residual standard error (RMSE).

- mae:

  Mean absolute error.

- n.missing:

  Number of missing elements.

- prop.missing:

  Percentage of missing elements.

- X1.sparsity:

  Proportion of near-zero elements in X1.

- C.sparsity:

  Proportion of near-zero elements in C.

- X2.sparsity:

  Proportion of near-zero elements in X2.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md),
[`print.summary.nmfae`](https://ksatohds.github.io/nmfkc/reference/print.summary.nmfae.md)

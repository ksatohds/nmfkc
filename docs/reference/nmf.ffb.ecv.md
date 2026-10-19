# Choose the number of factors for NMF-FFB by element-wise cross-validation

Step 1 of the NMF-FFB workflow: select \\Q\\ by element-wise
cross-validation of the *feed-forward* fit. Held-out entries of \\Y_1\\
are predicted from the remaining ones, so the criterion measures how
well a \\Q\\-factor non-negative basis describes the outcomes, which is
what \\Q\\ is for.

## Usage

``` r
nmf.ffb.ecv(
  Y1,
  Y2,
  rank = NULL,
  X.init = "nndsvd",
  X.L2.ortho = 100,
  epsilon = 1e-06,
  maxit = 5000,
  ...
)
```

## Arguments

- Y1, Y2:

  Endogenous and exogenous matrices (variables in rows, units in
  columns).

- rank:

  Candidate values of \\Q\\. Defaults to `seq_len(min(3, nrow(Y1)))`.

- X.init:

  Stage-1 initialization, passed to
  [`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md).

- X.L2.ortho:

  Orthogonality penalty on the basis in stage 1. Keep the value that
  will be used in
  [`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md):
  \\Q\\ and the penalty are not separable.

- epsilon, maxit:

  Convergence tolerance and iteration cap for stage 1.

- ...:

  Passed to
  [`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md):
  `nfolds`, `seed`, and `cores` (or `ncores`) to evaluate the rank x
  fold grid in parallel, defaulting to `getOption("mc.cores", 1L)` as
  elsewhere in the package.

## Value

The object returned by
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md):
the held-out error by rank, with the selected `rank`.

## Details

Feedback cannot be cross-validated, and that is not a limitation of the
implementation. By the reduced-form equivalence, the equilibrium mapping
\\M = (I - X\Theta_1)^{-1}X\Theta_2\\ is reproduced exactly by a
feed-forward model with the exogenous matrix \\L_Q\Theta_2\\: the two
models predict \\Y_1\\ from \\Y_2\\ identically, so no prediction
criterion – column-wise CV, element-wise CV, or a test-sample error –
can prefer one over the other. Feedback is identified by the
*conditional covariance*, not by prediction, which is why it is
addressed by
[`nmf.ffb.test`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md)
and not here.

## Workflow

    ecv <- nmf.ffb.ecv(Y1, Y2, rank = 1:5)   # 1. choose Q            <- this function
    fit <- nmf.ffb(Y1, Y2, rank = ecv$rank)  # 2. estimate; BIC selects the support
    tst <- nmf.ffb.test(fit, Y1, Y2)         # 3. test the feed-forward null
    dgn <- nmf.ffb.diagnostics(fit)          # 4. cycles, spectral radius, identifiability
    inf <- nmf.ffb.inference(fit, Y1, Y2)    # 5. intervals for the retained entries

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`nmf.ffb.test`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.test.md),
[`nmfkc.ecv`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ecv.md)

## Examples

``` r
set.seed(1)
Y <- t(as.matrix(iris[, 1:4]))
Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
ecv <- nmf.ffb.ecv(Y1, Y2, rank = 1:2, nfolds = 3)
#> Performing Element-wise CV for Q = 1,2 (3-fold)...
ecv$rank
#> [1] 1 2
```

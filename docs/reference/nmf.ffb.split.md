# Heuristic Variable Splitting for NMF-FFB

Infers a heuristic partition of observed variables into exogenous
(\\Y_2\\) and endogenous (\\Y_1\\) blocks for use in NMF-FFB. The method
is based on positive-SEM logic, causal ordering, and optional sign
alignment using the first principal component (PC1).

The procedure:

- internally standardizes variables (mean 0, sd 1),

- optionally flips signs so that most variables align positively with
  PC1,

- infers a causal ordering by repeatedly regressing each variable on the
  remaining ones and selecting the variable with the largest minimum
  standardized coefficient,

- determines an exogenous block by scanning the ordering from upstream
  and stopping at the first variable whose strongest parent coefficient
  exceeds `threshold`.

If `n.exogenous` is supplied, it overrides the automatic threshold rule.

## Usage

``` r
nmf.ffb.split(
  x,
  n.exogenous = NULL,
  threshold = 0.1,
  auto.flipped = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A numeric matrix or data frame with **rows = samples** and **columns =
  observed variables**.

- n.exogenous:

  Optional integer specifying the number of exogenous variables
  (\\Y_2\\). If `NULL`, the number is inferred automatically by the
  coefficient cut-off rule.

- threshold:

  Standardized regression-coefficient threshold used in the automatic
  exogenous–endogenous split. A variable is treated as endogenous once
  its maximum standardized parent coefficient exceeds this value.
  (Default: `0.1`)

- auto.flipped:

  Logical; if `TRUE`, applies PC1-based automatic sign flipping after
  standardization to ensure consistent orientation. (Default: `TRUE`)

- verbose:

  Logical; if `TRUE`, prints progress messages and the resulting
  variable split. (Default: `FALSE`)

- ...:

  Reserved for future use; currently unused (also accepted by the
  `nmf.ffb.split` alias for argument forwarding).

## Value

A list with:

- endogenous.variables:

  Character vector of variables selected as endogenous (\\Y_1\\).

- exogenous.variables:

  Character vector of variables selected as exogenous (\\Y_2\\).

- ordered.variables:

  Variables in inferred causal order (from exogenous to endogenous).

- is.flipped:

  Logical vector indicating which variables were sign-flipped during
  processing.

- n.exogenous:

  Integer giving the number of exogenous variables.

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md)

## Examples

``` r
# Infer exogenous/endogenous split from iris
sp <- nmf.ffb.split(iris[, -5], n.exogenous = 2)
sp$endogenous.variables
#> [1] "Petal.Width"  "Petal.Length"
sp$exogenous.variables
#> [1] "Sepal.Length" "Sepal.Width" 
```

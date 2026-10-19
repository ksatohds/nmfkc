# DOT graph visualization for nmfae objects

`nmfae.DOT` generates a DOT language string for visualizing the
structure of a three-layer NMF model. Two graph types are supported:
`"XCX"` shows encoder factors, \\\Theta\\, and decoder factors;
`"YXCXY"` shows the full structure from \\Y_2\\ through \\X_2\\,
\\\Theta\\, \\X_1\\ to \\Y_1\\.

Edge widths are proportional to matrix element values, and edges below
`threshold` are omitted for clarity.

## Usage

``` r
nmf.rrr.DOT(
  result,
  type = c("XCX", "YXCXY"),
  threshold = 0.01,
  sig.level = 0.1,
  rankdir = "LR",
  fill = TRUE,
  weight_scale = 5,
  weight_scale_x1 = weight_scale,
  weight_scale_theta = weight_scale,
  weight_scale_x2 = weight_scale,
  Y1.label = NULL,
  X1.label = NULL,
  X2.label = NULL,
  Y2.label = NULL,
  Y1.title = "Output (Y1)",
  X1.title = "Response (X1)",
  X2.title = "Covariate (X2)",
  Y2.title = "Input (Y2)",
  hide.isolated = TRUE
)
```

## Arguments

- result:

  An object of class `"nmfae"` returned by
  [`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md).

- type:

  Character. Graph type: `"XCX"` (default) or `"YXCXY"`.

- threshold:

  Numeric. Edges with values below this are omitted. Default is 0.01.

- sig.level:

  Numeric or `NULL`. Significance level for filtering C edges when
  inference results are available. Only edges with p-value below
  `sig.level` are shown, with significance stars. Set to `NULL` to
  disable. Default is 0.1.

- rankdir:

  Character. Graph direction for DOT layout. Default is `"LR"` (left to
  right).

- fill:

  Logical. If `TRUE`, nodes are filled with color. Default is `TRUE`.

- weight_scale:

  Numeric. Base scale factor for edge widths. Default is 5.

- weight_scale_x1:

  Numeric. Scale factor for \\X_1\\ edges.

- weight_scale_theta:

  Numeric. Scale factor for \\\Theta\\ edges.

- weight_scale_x2:

  Numeric. Scale factor for \\X_2\\ edges.

- Y1.label:

  Character vector of output variable labels.

- X1.label:

  Character vector of decoder basis labels.

- X2.label:

  Character vector of encoder basis labels.

- Y2.label:

  Character vector of input variable labels.

- Y1.title:

  Character. Title for output node group. Default is `"Output (Y1)"`.

- X1.title:

  Character. Title for the response-basis node group. Default is
  `"Response (X1)"`.

- X2.title:

  Character. Title for the covariate-basis node group. Default is
  `"Covariate (X2)"`.

- Y2.title:

  Character. Title for input node group. Default is `"Input (Y2)"`.

- hide.isolated:

  Logical. If `TRUE` (default), Y1 and Y2 nodes that have no edges at or
  above `threshold` are excluded from the graph. Only applies when
  `type = "YXCXY"`.

## Value

A character string containing the DOT graph specification.

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmf.rrr`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.md)

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(runif(20), nrow = 4)
res <- nmf.rrr(Y, rank1 = 2)
dot <- nmf.rrr.DOT(res)
# }
```

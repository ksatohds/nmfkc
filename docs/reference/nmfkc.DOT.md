# Generate Graphviz DOT Scripts for NMF or NMF-with-Covariates Models

Produces a Graphviz DOT script visualizing the structure of an NMF model
(\\Y \approx X C A\\) or its simplified forms.

Supported visualization types:

- `"YX"` — Standard NMF view: latent factors \\X\\ map to observations
  \\Y\\.

- `"YA"` — Direct regression view: covariates \\A\\ map directly to
  \\Y\\ using the combined coefficient matrix \\X C\\.

- `"YXA"` — Full tri-factorization: \\A \rightarrow C \rightarrow X
  \rightarrow Y\\.

Edge widths are scaled by coefficient magnitude, and nodes with no edges
above the threshold are omitted from the visualization.

## Usage

``` r
nmfkc.DOT(
  result,
  type = c("YX", "YA", "YXA"),
  threshold = 0.01,
  C.signed = NULL,
  sig.level = 0.1,
  rankdir = "LR",
  fill = TRUE,
  weight_scale = 5,
  weight_scale_ax = weight_scale,
  weight_scale_xy = weight_scale,
  weight_scale_ay = weight_scale,
  Y.label = NULL,
  X.label = NULL,
  A.label = NULL,
  Y.title = "Observation (Y)",
  X.title = "Basis (X)",
  A.title = "Covariates (A)",
  hide.isolated = TRUE
)
```

## Arguments

- result:

  The return value from `nmfkc`, containing matrices `X`, `B`, and
  optionally `C`.

- type:

  Character string specifying the visualization style: one of `"YX"`,
  `"YA"`, `"YXA"`.

- threshold:

  Minimum coefficient magnitude to display an edge. When
  `C.signed = TRUE` (signed \\\Theta\\) this is an absolute-value cut,
  i.e. an edge is drawn when \\\|coef\| \ge\\ `threshold`.

- C.signed:

  Logical or `NULL`. Whether the coefficient matrix \\\Theta\\ (`= C`)
  may be signed (real-valued). When `TRUE`, the `threshold` is applied
  to \\\|coef\|\\, edge widths are scaled by \\\|coef\|\\, and negative
  edges are drawn as black **dashed** lines (positive edges remain
  solid) to distinguish them; the numeric edge labels keep their sign.
  When `FALSE`, the historical non-negative behaviour is used (negative
  coefficients fall below the threshold and are hidden). When `NULL`
  (default), the mode is auto-detected from the fit (`result$C.signed`)
  or from the presence of negative entries in \\C\\ / \\XC\\. The basis
  \\X\\ is always non-negative, so \\X \rightarrow Y\\ edges are
  unaffected.

- sig.level:

  Significance level for filtering C edges when inference results are
  available (i.e., `x` is of class `"nmfkc.inference"`). Only edges with
  p-value below `sig.level` are shown, decorated with significance stars
  (`*`, `**`, `***`). Set to `NULL` to disable filtering and show all
  edges above `threshold`. Default is 0.1.

- rankdir:

  Graphviz rank direction (e.g., `"LR"`, `"TB"`).

- fill:

  Logical; whether nodes should be drawn with filled shapes.

- weight_scale:

  Base scaling factor for edge widths.

- weight_scale_ax:

  Scaling factor for edges \\A \rightarrow X\\ (type `"YXA"`).

- weight_scale_xy:

  Scaling factor for edges \\X \rightarrow Y\\.

- weight_scale_ay:

  Scaling factor for edges \\A \rightarrow Y\\ (type `"YA"`).

- Y.label:

  Optional character vector for labels of Y nodes.

- X.label:

  Optional character vector for labels of X (latent factor) nodes.

- A.label:

  Optional character vector for labels of A (covariate) nodes.

- Y.title:

  Cluster title for Y nodes.

- X.title:

  Cluster title for X nodes.

- A.title:

  Cluster title for A nodes.

- hide.isolated:

  Logical. If `TRUE` (default), Y and A nodes that have no edges at or
  above `threshold` are excluded from the graph.

## Value

A character string representing a Graphviz DOT script.

## See also

[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`nmf.rrr.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.rrr.DOT.md),
[`nmf.ffb.DOT`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.DOT.md),
[`nmfkc.ar.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.DOT.md),
[`plot.nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.DOT.md)

## Examples

``` r
Y <- matrix(cars$dist, nrow = 1)
A <- rbind(1, cars$speed)
result <- nmfkc(Y, A, rank = 1)
#> Y(1,50)~X(1,1)C(1,2)A(2,50)=XB(1,50)...
#> 0sec
dot <- nmfkc.DOT(result)
cat(dot)
#> digraph NMF {
#>   graph [rankdir=LR compound=true];
#>   splines=true; nodesep=0.4; ranksep=0.7; fontname="Arial";
#> 
#>   // Output variables (Y)
#>   subgraph cluster_Y{label="Observation (Y)" style="rounded" color="black" penwidth=1.0;
#>   node [shape=box, style="filled,rounded", fillcolor="lightblue", color=black, penwidth=1.5];
#>     Y_1 [label="Y1"];
#>   }
#> 
#>   // Latent factors (X)
#>   subgraph cluster_X{label="Basis (X)" style="rounded" color="black" penwidth=1.0;
#>   node [shape=ellipse, style="filled,rounded", fillcolor="wheat", color=black, penwidth=1.0];
#>     X_1 [label="Basis1"];
#>   }
#> 
#>   edge [fontname="Arial", fontsize=8, arrowhead=open];
#> 
#>   // X -> Y edges (X)
#>   edge [color="gray0", fontcolor="gray0", style=solid];
#>   X_1 -> Y_1 [label="1.00", penwidth=5.00];
#> }
```

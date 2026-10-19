# Generate a Graphviz DOT Diagram for NMF-AR / NMF-VAR Models

Produces a Graphviz DOT script for visualizing autoregressive
NMF-with-covariates models constructed via `nmfkc.ar` + `nmfkc`.

The diagram displays three types of directed relationships:

- Lagged predictors: \\T\_{t-k} \rightarrow X\\,

- Current latent factors: \\X \rightarrow T_t\\,

- Optional intercept effects: `Const -> X`.

Importantly, *no direct edges from lagged variables to current outputs*
(\\T\_{t-k} \rightarrow T_t\\) are drawn, in accordance with the NMF-AR
formulation.

Each block of lagged variables is displayed in its own DOT subgraph
(e.g., “T-1”, “T-2”, ...), while latent factor nodes and current-time
outputs are arranged in separate clusters.

## Usage

``` r
nmfkc.ar.DOT(
  result,
  degree = 1,
  intercept = any(colnames(result$C) == "(Intercept)"),
  threshold = 0.1,
  rankdir = "RL",
  fill = TRUE,
  weight_scale_xy = 5,
  weight_scale_lag = 5,
  weight_scale_int = 3,
  hide.isolated = TRUE
)
```

## Arguments

- result:

  A fitted `nmfkc` object representing the AR model. Must contain
  matrices `X` and `C`.

- degree:

  Maximum AR lag to visualize.

- intercept:

  Logical; if `TRUE`, draws intercept nodes for columns named
  "(Intercept)" in matrix `C`. Default is `TRUE` when an intercept
  column is detected in `C`, `FALSE` otherwise (auto-detected).

- threshold:

  Minimum coefficient magnitude required to draw an edge.

- rankdir:

  Graphviz rank direction (e.g., `"RL"`, `"LR"`, `"TB"`).

- fill:

  Logical; whether nodes are filled with color.

- weight_scale_xy:

  Scaling factor for edges \\X \rightarrow T\\.

- weight_scale_lag:

  Scaling factor for lagged edges \\T-k \rightarrow X\\.

- weight_scale_int:

  Scaling factor for intercept edges.

- hide.isolated:

  Logical. If `TRUE` (default), Y nodes that have no edges at or above
  `threshold` are excluded from the graph.

## Value

A character string representing a Graphviz DOT file.

## See also

[`nmfkc.ar`](https://ksatohds.github.io/nmfkc/reference/nmfkc.ar.md),
[`nmfkc`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md),
[`plot.nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.DOT.md)

## Examples

``` r
d <- AirPassengers
ar_data <- nmfkc.ar(d, degree = 2)
result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#> Y(1,142)~X(1,1)C(1,3)A(3,142)=XB(1,142)...
#> 0sec
dot <- nmfkc.ar.DOT(result, degree = 2)
cat(dot)
#> digraph NMF_AR {
#>   graph [rankdir=RL compound=true];
#>   splines=true; nodesep=0.4; ranksep=0.7; fontname="Arial";
#> 
#>   // Current-time outputs (T)
#>   subgraph cluster_Y{label="T" style="rounded" color="black" penwidth=1.0;
#>   node [shape=box, style="filled,rounded", fillcolor="lightblue", color=black, penwidth=1.5];
#>     1 [label="1"];
#>   }
#> 
#>   // Latent variables (X)
#>   subgraph cluster_X{label="Latent Variables" style="rounded" color="black" penwidth=1.0;
#>   node [shape=ellipse, style="filled,rounded", fillcolor="wheat", color=black, penwidth=1.0];
#>     Basis1 [label="Basis1"];
#>   }
#> 
#>   edge [fontname="Arial", fontsize=8, arrowhead=open];
#> 
#>   // X -> T edges (factor loadings)
#>   edge [color="gray0", fontcolor="gray0", style=solid];
#>   Basis1 -> 1 [label="1.0", penwidth=5.00];
#> 
#>   // Intercept nodes
#>   Const1 [shape=circle, label="11.0"];
#>   Const1 -> Basis1 [penwidth=3.00];
#> 
#>   // Lag block T-1
#>   subgraph cluster_C1 {label="T-1" style="rounded";
#>     node [shape=box, style="filled,rounded", fillcolor="lightcoral", color=black, penwidth=1.5];
#>     1_1 [label="1"];
#>   }
#>   // T-1 -> X edges
#>   edge [color=black, fontcolor=black, style=solid];
#>   1_1 -> Basis1 [label="0.9", penwidth=5.00];
#> }
```

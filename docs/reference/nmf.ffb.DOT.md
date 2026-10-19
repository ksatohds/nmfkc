# Generate a Graphviz DOT Diagram for an NMF-FFB Model

Creates a Graphviz DOT script that visualizes the structural network
estimated by `nmf.ffb`. The resulting diagram displays:

- endogenous observed variables (\\Y_1\\),

- exogenous observed variables (\\Y_2\\),

- latent factors (\\F_1\\, ..., \\F_Q\\),

together with the non-negative path coefficients whose magnitudes exceed
a user-specified threshold.

Directed edges represent estimated relationships:

- \\Y_2 \rightarrow F_q\\: entries of `C2` (exogenous loadings),

- \\F_q \rightarrow Y_1\\: rows of `X` (factor-to-endogenous mappings),

- \\Y_1 \rightarrow F_q\\: entries of `C1` (feedback paths).

Edge widths are scaled by coefficient magnitude, and nodes are placed in
optional visual clusters. Only variables participating in edges above
the threshold are displayed, while latent factors are always shown.

## Usage

``` r
nmf.ffb.DOT(
  result,
  weight_scale = 5,
  weight_scale_c2 = weight_scale,
  weight_scale_x1 = weight_scale,
  weight_scale_feedback = weight_scale,
  threshold = 0.01,
  sig.level = 0.1,
  rankdir = "LR",
  fill = TRUE,
  cluster.box = c("normal", "faint", "invisible", "none"),
  cluster.labels = NULL,
  hide.isolated = TRUE,
  model = c("selected", "null", "full"),
  ...
)
```

## Arguments

- result:

  A list returned by `nmf.ffb`, containing matrices `X`, `C1`, and `C2`.

- weight_scale:

  Base scaling factor for edge widths.

- weight_scale_c2:

  Scaling factor for edges \\Y_2 \rightarrow F_q\\ (C2 matrix). Defaults
  to `weight_scale`.

- weight_scale_x1:

  Scaling factor for edges \\F_q \rightarrow Y_1\\ (X matrix). Defaults
  to `weight_scale`.

- weight_scale_feedback:

  Scaling factor for feedback edges \\Y_1 \rightarrow F_q\\ (C1 matrix).
  Defaults to `weight_scale`.

- threshold:

  Minimum coefficient value needed for an edge to be drawn.

- sig.level:

  Significance level for filtering structural edges (\\C_1\\ feedback
  and \\C_2\\ exogenous loadings) when inference results are present. If
  `result` contains a `coefficients` data frame from
  [`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md),
  only edges with `prob.unsupported < sig.level` are drawn, with
  significance stars (`*` `**` `***`) appended to the edge label. The
  \\X\\ (factor-to-\\Y_1\\) edges are never starred since the basis is
  not the inference target. Set to `NULL` to disable significance
  filtering and fall back to the `threshold` magnitude filter for both
  \\C_1\\ and \\C_2\\. Default is `0.1`.

- rankdir:

  Graphviz rank direction (e.g., `"LR"`, `"TB"`).

- fill:

  Logical; whether to use filled node shapes.

- cluster.box:

  Character string controlling the visibility and style of cluster
  frames around Y2, factors, and Y1 blocks. One of `"normal"`,
  `"faint"`, `"invisible"`, `"none"`.

- cluster.labels:

  Optional character vector of length 3 giving custom labels for the Y2,
  factor, and Y1 clusters.

- hide.isolated:

  Logical. If `TRUE` (default), Y1 and Y2 nodes that have no edges at or
  above `threshold` are excluded from the graph.

- model:

  Which fit of a `method = "fiml"` object to draw: `"selected"`
  (default; the BIC-selected model, with significance stars when
  inference results are present), `"null"` (the feed-forward null:
  \\\Theta_1 = 0\\, \\\Theta_2\\ from `result$null`) or `"full"` (the
  unpenalized feedback fit from `result$full`). The latter two ignore
  `result$coefficients`, which refer to the selected model. Ignored for
  `method = "mu"` objects, which carry a single fit.

- ...:

  For backward compatibility: accepts deprecated names
  `weight_scale_y2f` (use `weight_scale_c2`) and `weight_scale_fy1` (use
  `weight_scale_x1`).

## Value

A character string representing a valid Graphviz DOT script.

## See also

[`nmf.ffb`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.md),
[`nmf.ffb.inference`](https://ksatohds.github.io/nmfkc/reference/nmf.ffb.inference.md),
[`plot.nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.DOT.md)

## Examples

``` r
Y <- t(iris[, -5])
Y1 <- Y[1:2, ]
Y2 <- Y[3:4, ]
result <- nmf.ffb(Y1, Y2, rank = 2, maxit = 500)
dot <- nmf.ffb.DOT(result)
cat(dot)
#> digraph NMF_SEM_Full_Mechanism {
#>   graph [rankdir=LR compound=true];
#>   splines=true; nodesep=0.4; ranksep=0.7; fontname="Arial";
#> 
#>   // Exogenous variables (Y2)
#>   subgraph cluster_Y2{label="Exogenous (Y2)" style="rounded" color="black" penwidth=1.0;
#>   node [shape=box, style="filled,rounded", fillcolor="lightcoral", color=black, penwidth=1.5];
#>     Y2_1 [label="Petal.Length"];
#>     Y2_2 [label="Petal.Width"];
#>   }
#> 
#>   // Endogenous variables (Y1)
#>   subgraph cluster_Y1{label="Endogenous (Y1)" style="rounded" color="black" penwidth=1.0;
#>   node [shape=box, style="filled,rounded", fillcolor="lightblue", color=black, penwidth=1.5];
#>     Y1_1 [label="Sepal.Length"];
#>     Y1_2 [label="Sepal.Width"];
#>   }
#> 
#>   // Latent Factors (F)
#>   subgraph cluster_F{label="Latent Factors" style="rounded" color="black" penwidth=1.0;
#>   node [shape=ellipse, style="filled,rounded", fillcolor="wheat", color=black, penwidth=1.0];
#>     F_1 [label="Factor 1"];
#>     F_2 [label="Factor 2"];
#>   }
#> 
#>   edge [fontname="Arial", fontsize=8, arrowhead=open];
#> 
#>   // 1. External Driving (Y2 -> Factor) [C2]
#>   edge [color=black, fontcolor=black, style=solid];
#>   Y2_1 -> F_1 [label="1.41", penwidth=5.00];
#>   Y2_1 -> F_2 [label="0.40", penwidth=1.43];
#>   Y2_2 -> F_2 [label="0.56", penwidth=1.99];
#> 
#>   // 2. Generation (Factor -> Y1) [X]
#>   edge [color="gray0", fontcolor="gray0", style=solid];
#>   F_1 -> Y1_1 [label="0.96", penwidth=4.81];
#>   F_1 -> Y1_2 [label="0.04", penwidth=0.50];
#>   F_2 -> Y1_2 [label="0.99", penwidth=5.00];
#> 
#>   // 3. Internal Feedback (Y1 -> Factor) [C1]
#>   edge [style=dashed, color="gray0", fontcolor="gray0"];
#> }
```

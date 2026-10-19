# Generate a Graphviz DOT Diagram for a Symmetric NMF Network

Creates a Graphviz DOT script that visualizes the two-layer structure of
a symmetric NMF model (\\Y \approx X C X^\top\\).

The resulting diagram displays:

- outer nodes: the original network nodes (rows/columns of Y),

- inner nodes: the latent basis/group nodes (columns of X),

- directed edges from basis to node: membership weights from X,

- undirected edges between basis nodes: inter-group interactions from C
  matrix (for tri-symmetric models only).

## Usage

``` r
nmfkc.net.DOT(
  result,
  threshold = 0.01,
  sig.level = 0.1,
  weight_scale = 5,
  weight_scale_xy = 1,
  weight_scale_xx = weight_scale,
  rankdir = "TB",
  fill = TRUE,
  hide.isolated = TRUE,
  Y.label = NULL,
  X.label = NULL,
  Y.title = "Nodes (Y)",
  X.title = "Basis (X)",
  show.theta = NULL,
  signed = inherits(result, "nmfkc.net.signed"),
  cluster.box = c("none", "normal", "faint", "invisible"),
  layout = c("neato", "fdp", "twopi", "circo", "dot"),
  X.color = NULL,
  Y.cluster = c("soft", "hard")
)
```

## Arguments

- result:

  A list returned by
  [`nmfkc()`](https://ksatohds.github.io/nmfkc/reference/nmfkc.md) with
  `Y.symmetric = "bi"` or `"tri"`, or the newer
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  with `type = "bi"` or `"tri"`. If inference results are present (from
  [`nmfkc.net.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.inference.md)),
  C edges are decorated with significance stars.

- threshold:

  Minimum coefficient value to display an edge. It is read on the
  membership scale for X edges (`X.prob`, rows summing to one) and on
  \\C\\ for the basis-to-basis edges, and it does two jobs at once,
  switching over at \\1/Q\\:

  - **At or below \\1/Q\\** it only thins edges. No node can be hidden,
    because a row of `X.prob` sums to one and so has a maximum of at
    least \\1/Q\\; every node genuinely keeps an edge.

  - **Above \\1/Q\\** it may also start hiding nodes (when
    `hide.isolated = TRUE`). Not immediately: nodes go only once the
    threshold passes the smallest row maximum of `X.prob`, which is
    data-dependent and can be well above \\1/Q\\. Inspect
    `apply(fit$X.prob, 1, max)` to see where that is.

  The default `0.01` therefore draws every node — guaranteed for any \\Q
  \le 100\\, since only then is \\0.01 \le 1/Q\\ — and most edges. For a
  readable community graph, 0.2–0.4 usually thins the edges enough while
  keeping all nodes.

- sig.level:

  Significance level for filtering C edges (if inference results are
  present). Set to `NULL` to show all edges above threshold.

- weight_scale:

  Base scaling factor for edge widths.

- weight_scale_xy:

  Scaling factor for X edges (basis -\> node).

- weight_scale_xx:

  Scaling factor for C edges (basis \<-\> basis).

- rankdir:

  Graphviz rank direction (`"TB"` or `"LR"`).

- fill:

  Logical; whether nodes are filled with color.

- hide.isolated:

  Logical; if TRUE, omit outer nodes with no X edge above `threshold`.
  The test is on the membership scale (`X.prob`), the same quantity the
  X edges are drawn from, so it does not shift with `type` or
  `X.restriction`. See `threshold` for when this actually removes
  anything — at the default it never does for any \\Q \le 100\\, which
  is correct rather than a limitation.

- Y.label:

  Character vector of labels for outer nodes.

- X.label:

  Character vector of labels for basis nodes.

- Y.title:

  Cluster title for outer nodes.

- X.title:

  Cluster title for basis nodes.

- show.theta:

  Logical or NULL. Whether to draw C edges between basis nodes. NULL =
  auto-detect (TRUE for tri, FALSE for bi).

- signed:

  Logical. If `TRUE`, \\C\\ is treated as a signed matrix: positive
  entries rendered as solid edges and negative as dashed, with edge
  visibility threshold on \\\|C\|\\. Default is
  `inherits(result, "nmfkc.net.signed")` so that results from
  [`nmfkc.net`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.md)
  are auto-detected.

- cluster.box:

  Style of cluster box: `"none"`, `"normal"`, `"faint"`, `"invisible"`.

- layout:

  Graphviz layout engine, in recommended order: `"neato"` (default;
  spring model, clearest for small/medium community graphs), `"fdp"`
  (force-directed, scales to larger graphs), `"twopi"` (radial),
  `"circo"` (circular), `"dot"` (hierarchical). For community networks,
  `"neato"` or `"fdp"` with a raised `threshold` (e.g.\\ 0.2–0.3)
  separate the groups best.

- X.color:

  Color palette for basis nodes (length Q).

- Y.cluster:

  Coloring mode for outer nodes: `"soft"` (weighted mix) or `"hard"`
  (most probable basis color).

## Value

A character string of class `c("nmfkc.net.DOT", "nmfkc.DOT")`
representing a valid Graphviz DOT script. Use
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to render
(requires DiagrammeR).

## Lifecycle

This function is **experimental**. The interface may change in future
versions; details are to be described in an upcoming paper.

## See also

[`nmfkc.net.inference`](https://ksatohds.github.io/nmfkc/reference/nmfkc.net.inference.md),
[`nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/nmfkc.DOT.md),
[`plot.nmfkc.DOT`](https://ksatohds.github.io/nmfkc/reference/plot.nmfkc.DOT.md)

## Examples

``` r
# \donttest{
library(nmfkc)
Y <- matrix(c(0,1,1,0,0,0,
              1,0,1,0,0,0,
              1,1,0,1,0,0,
              0,0,1,0,1,1,
              0,0,0,1,0,1,
              0,0,0,1,1,0), 6, 6)
res <- nmfkc.net(Y, rank = 2, type = "tri", nstart = 20)
dot <- nmfkc.net.DOT(res)
plot(dot)
# }
```

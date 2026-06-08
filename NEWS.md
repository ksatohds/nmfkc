# nmfkc 0.7.4 (development)

### **`nmfkc.net.DOT()`: default layout is now `"neato"`**
- The `layout` choices are reordered by recommendation
  (`neato`, `fdp`, `twopi`, `circo`, `dot`), so the default changes from
  `"fdp"` to `"neato"`, which separates community graphs more clearly.  Raising
  `threshold` (e.g.\ 0.2--0.3) further declutters weak membership edges.

### **Bug fix: `nmfkc.net.DOT()` mis-detected `type = "bi"` as `"tri"`**
- The bi-vs-tri auto-detection ignored the result's `$type` field and fell back
  to `all.equal(C, diag(Q))`, which fails when `C` carries dimnames (it reports a
  names mismatch).  A `type = "bi"` fit was therefore treated as `"tri"`, drawing
  the inter-class interaction layer that the bi model (with \eqn{C = I}) should
  not have.  Detection now uses `$type` first (falling back to the dimnames-safe
  identity check), so `"bi"` correctly draws no inter-class edges.

### **`nmfkc.ard()`: better default prior scale**
- The default `b` is now the initial per-component energy scale
  `(nrow(Y) + ncol(Y)) / K * mean(Y)` instead of a fixed `0.001 * mean(Y)`.
  The old fixed fraction over-pruned (winner-take-all collapse onto one
  dominant component) when `(F + N)/K` was large; the new scale-aware
  default recovers genuine low-rank structure stably (e.g. a clean rank-3
  signal: relevance `1, 0.99, 0.87, 0, ...`, all restarts agree).

### **New `nmfkc.ard()`: ARD rank determination (Tan & Fevotte 2013, prototype)**
- Automatic Relevance Determination for the NMF rank (Euclidean).  Fits
  NMF *once* at an over-complete rank and prunes automatically: each
  component carries a relevance weight with an inverse-gamma prior and the
  multiplicative updates gain a penalty (`L2` half-normal / `L1`
  exponential) that drives unsupported components to zero.  The number of
  surviving components is the estimated rank -- no rank scan.  Returns an
  `"nmfkc.ard"` object with `print` and a relevance-bar `plot`.  Plain NMF
  only; a sensitive point estimate (depends on prior / start / init), so a
  complement to the CV / consensus engines, not a sole criterion.

### **New `nmfkc.consensus()`: consensus-clustering rank selection (Brunet 2004)**
- The bioinformatics-standard stability approach, as a lightweight engine
  like `nmfkc.ecv` / `nmfkc.bicv`.  For each rank it runs NMF `nrun` times
  from random initializations (`X.init = "runif"`), builds the consensus
  matrix from the per-run hard clusterings, and returns two stability
  scores per rank: `cophenetic` (cophenetic correlation coefficient,
  Brunet et al. 2004) and `dispersion` (Kim & Park 2007, in `[0,1]`).
  Unlike the CV engines, a good rank *maximizes* stability.  Optional
  `keep.consensus = TRUE` returns the consensus matrices.
- Also reports `pac`, the Proportion of Ambiguous Clustering
  (Senbabaoglu et al. 2014; fraction of consensus entries in the
  ambiguous interval `pac.range`, default `(0.1, 0.9)`).  Lower is better
  and it is more sensitive than the often-saturated `cophenetic`.  The
  `print`/criteria-`plot` show all three metrics.
- Returns an `"nmfkc.consensus"` object with `print` and `plot` methods:
  `plot(cs)` (`type = "criteria"`) draws the stability curves;
  `plot(cs, type = "heatmap", rank = ...)` draws the consensus matrix
  heatmap(s) reordered by hierarchical clustering (default = all ranks in
  a `n2mfrow` grid; `mfrow` overridable).

### **New `nmfkc.bicv()`: bi-cross-validation for rank selection**
- Owen & Perry's (2009) bi-cross-validation (BCV), a **lightweight CV
  engine in the spirit of `nmfkc.ecv`**: it returns the held-out error
  per rank (`objfunc`, `sigma`) and nothing more.  Holds out a row-block
  *and* a column-block at once, fits NMF only on the retained block, and
  predicts the held-out block by folding the held-out rows/columns onto
  the fixed factors via non-negative regression (no information leakage,
  unlike element-wise `nmfkc.ecv`).  `nfolds = 2` (leave out half rows /
  half columns) per Owen & Perry's recommendation.

### **`*.rank`: eff.rank.idx shown for context (no best marker)**
- The broken-stick-corrected effective-rank index (`eff.rank.idx`,
  green) is drawn for context only and no longer carries a "Best (Max)"
  marker: it is a factor-utilization diagnostic (most even relative to
  the random null), not a predictive rank optimum.  The recommended rank
  is driven solely by the ECV minimum and the R-squared elbow.

### **`*.rank`: broken-stick-corrected effective-rank index**
- The `*.rank` criteria table gains `effective.rank.expected` (the
  broken-stick / uniform-Dirichlet null `exp(H_Q - 1)`, `H_Q` = the
  `Q`-th harmonic number) and `effective.rank.index`, the \[0, 1] index
  `(effective.rank - expected) / (Q - expected)` (clamped).  The index
  anchors 0 at the random null and 1 at perfect evenness, removing the
  small-rank inflation of the raw `effective.rank / Q`.  Its maximum is
  a meaningful rank, so the diagnostics plot now draws this corrected
  index (green, `eff.rank.idx`) with a restored "Best (Max)" marker in
  place of the raw ratio.

### **`*.rank` results gain `plot()` / `print()` methods**
- The rank-selection functions (`nmfkc.rank()`, `nmfkc.net.rank()`,
  `nmfkc.signed.rank()`, `nmfae.rank()`, `nmfae.signed.rank()`) now
  return a classed object (`"nmf.rank"`).  `plot()` redraws the
  three-criterion diagnostics plot (honouring `main`, `xlab`, `ylab`,
  `lwd`) and `print()` shows the recommended rank, the per-criterion
  best ranks, and the criteria table.  As before the constructor draws
  immediately when `plot = TRUE`; the `$rank.best` and `$criteria`
  fields are unchanged, so existing code keeps working.

### **New `nmf.cluster.flow()`: cluster-flow diagram across ranks**
- `nmf.cluster.flow()` and `nmf.cluster.criteria()` now treat the
  supplied `fits` as a generic \emph{sequence of results} (kept in the
  given order, \strong{not} sorted by rank), so the same rank fitted as
  different models is also supported.  Both gain a `names` argument for
  the x-axis tick labels (default: each result's `$rank`), and in
  `nmf.cluster.flow()` the `reference` argument is now the \strong{index}
  (1-based position) of the result that defines the colours -- not a rank
  value -- defaulting to the central result
  `floor(length(fits) / 2) + 1` (e.g.\ the 2nd of 2 or 3 results).
- The adjusted Rand index (ARI) between each pair of adjacent ranks is
  now computed and printed along the top of the figure (and returned in
  `$ARI`, length \eqn{R - 1}), summarizing how much the hard clustering
  changes from one rank to the next.
- Each cluster box is now tinted by the \strong{majority} reference
  colour among the individuals it contains (the colour shared by the
  most member lines); ties are broken in favour of the earliest palette
  entry (the smallest reference-cluster id).  This shows at a glance
  which reference cluster dominates each box at each rank.
- `nmf.cluster.flow()` now inserts a gap of one average cluster
  (\eqn{N / k}) between clusters in the per-rank layout and sizes each
  grey box exactly to the minimum/maximum position of its members, so
  the cluster boxes are clearly separated with the gaps maximized.  Each
  rank is normalized to the full height independently.
- The cluster number is the dominant-factor index (argmax of the
  coefficient) of each fit, kept as-is so it matches the factor/basis
  numbering of the supplied models.  A factor that never dominates any
  individual leaves an empty, unused cluster number (a gap, e.g.\ labels
  `2, 3` with no `1`) -- this is correct and consistent with the fit,
  and the labels are not renumbered.
- `nmf.cluster.flow()` now returns a classed object with a dedicated
  `plot()` method, so the diagram can be (re)drawn with
  `plot(fl, col = , lwd = , xlab = , ylab = , main = )` -- the colour
  vector (indexed by reference cluster), line width, axis labels and
  title are all honoured.  The constructor still draws immediately by
  default (`plot = TRUE`) and forwards graphical arguments to the plot
  method; use `plot = FALSE` to build the object and plot it later.
  Its `print()` method shows the adjacent-rank ARI and the full
  \eqn{N \times R} cluster table.
- `nmf.cluster.flow(fits, reference = )` takes a list of models fitted
  at different ranks (any non-negative MU family) and draws an
  alluvial / Sankey-style diagram of how the hard sample clustering
  changes with the rank \eqn{Q}: each individual flows left-to-right
  across the ranks (x-axis), its vertical position is set by its cluster
  (clusters reordered per rank by a barycenter heuristic to reduce
  crossings), and lines are coloured by the cluster at the
  \code{reference} rank -- so one can watch the reference clusters split
  or merge.  At every rank a translucent grey box is drawn \emph{in
  front} of each cluster's members with the cluster number centred
  inside, so the grouping and labels are visible at all ranks (not only
  the reference).  The default line palette is now a strong,
  well-separated qualitative set (ColorBrewer \dQuote{Dark 2}, no pale
  colours) and can be overridden with \code{col}.  Returns (invisibly)
  the \eqn{N \times R} table with rows = individuals, columns = rank,
  entries = cluster number.

### **New `nmf.cluster.criteria()`: sample-clustering quality across ranks**
- `nmf.cluster.criteria(fits, Y)` takes a \strong{list of fits} (one per
  rank; a single fit is also accepted) and reports the clustering-quality
  criteria `silhouette`, `CPCC`, and `dist.cor` for each rank, returning
  a per-rank `$criteria` table (mirroring `nmf.cluster.flow()`).  It has
  `plot()` (line plot of the three criteria vs rank) and `print()` (the
  table) methods, and draws immediately when `plot = TRUE`.  Works for
  any family (`nmfkc`, `nmfkc.signed`, `nmfae`, `nmfae.signed`,
  `nmfkc.net`, `nmfre`, `nmf.sem`/`nmf.ffb`; the last needs the exogenous
  block via `Y2`).  These are **clustering-stability** diagnostics,
  deliberately separate from the rank-selection `*.rank` functions
  (r.squared / effective rank / ECV).
- Hard sample clustering needs a non-negative coefficient/score matrix
  (a valid membership simplex).  `nmf.cluster.criteria()` detects this from the
  actual coefficient: when it is non-negative the hard-label
  `silhouette` (and cluster sizes) are returned; when it is signed
  `silhouette` is `NA` while the distance-based `CPCC` and `dist.cor`
  are still computed.  (ARI is not reported here -- it compares two
  clusterings, e.g.\ across ranks or resamples, so it is not a
  single-fit quantity.)
- `nmfkc.rank()` no longer carries `ARI`, `silhouette`, `CPCC`, or
  `dist.cor` in its `criteria` table -- those clustering-stability
  metrics now live in `nmf.cluster.criteria()`.  All five `*.rank` functions
  return the **same five columns** (`rank`, `effective.rank`,
  `effective.rank.ratio`, `r.squared`, `sigma.ecv`).  Per-rank fits use
  `detail = "fast"`, so the expensive O(N^2) distance computations are
  skipped during rank selection.  `rank.best` is unchanged.  The
  `*.rank` functions now emit a one-line message pointing to
  `nmf.cluster.criteria()` for clustering quality.

### **Rank-selection functions for the other NMF families**
- New `nmfkc.net.rank()`, `nmfkc.signed.rank()`, `nmfae.rank()`
  (paired \eqn{Q = R}) and `nmfae.signed.rank()` (paired) bring
  `nmfkc.rank`-style rank selection to the other multiplicative-update
  models.  Each reports the three criteria that are well defined for
  every family -- `r.squared`, the effective rank (utilization), and the
  element-wise CV error `sigma.ecv` -- and returns
  `list(rank.best, criteria)`.  (`nmf.ffb` / `nmfre` are not covered:
  they do not support the element masking that ECV needs.)
- **`nmfkc.rank()` plot simplified and unified.**  All `*.rank`
  functions now share one back-end `.rank.finish()` and draw the same
  concise three-criterion figure: `r.squared` (red), `eff.rank` (green),
  and `sigma.ecv` (blue, right axis), each as a line with points,
  rank-number labels, and a highlighted best marker -- "Best (Elbow)"
  for the R-squared knee, "Best (Peak)" for the effective-rank
  utilization, and "Best (Min)" for the CV minimum.  `nmfkc.rank()`
  still computes `ARI`, `silhouette`, `CPCC`, and `dist.cor` into its
  `criteria` table, but no longer plots them.
- The four new `*.rank` functions gain a `detail` argument matching
  `nmfkc.rank`: `"full"` (default) runs the element-wise CV and reports
  `sigma.ecv`; `"fast"` skips the (expensive) CV, so the plot shows only
  `r.squared` and `eff.rank` and the recommended rank falls back to the
  R-squared elbow.

### **Internal: shared element-wise CV helpers**
- The four element-wise cross-validation functions (`nmfkc.ecv()`,
  `nmfae.ecv()`, `nmfkc.signed.ecv()`, `nmfae.signed.ecv()`) now build
  their folds through a single internal helper `.ecv.make.folds()`,
  removing four near-identical copies of the fold-partitioning loop.
  `nmfkc.net.ecv()` keeps its symmetric upper-triangle folds.
- \strong{All five} element-wise CV functions now share one
  config-indexed loop driver `.ecv.run(labels, nfolds, run_one,
  progress)`: the single-rank ones (`nmfkc.ecv()`, `nmfkc.net.ecv()`,
  `nmfkc.signed.ecv()`) and the \eqn{(Q, R)}-grid ones (`nmfae.ecv()`,
  `nmfae.signed.ecv()`).  Each supplies a model-specific
  `run_one(i, k)` closure (mask fold, refit config `i`, return held-out
  loss) and an optional progress callback; `.ecv.run()` handles the
  config-by-fold loop, the `objfunc`/`sigma`/`objfunc.fold` aggregation,
  and naming.  This removes the last copies of the CV-loop machinery,
  including the per-grid reshaping in `nmfae.ecv()`.
- The refactor is behaviour-preserving: for the same seed the folds and
  all CV values (`objfunc`, `sigma`, `objfunc.fold`, names/labels) are
  byte-for-byte identical to before, verified across EU and KL losses,
  the symmetric (upper-triangle) case, and both paired and full
  \eqn{(Q, R)} grids.

### **Unified summary print blocks**
- New shared internal helpers `.print.fit.statistics()` and
  `.print.structure.diagnostics()` render the "Statistics" /
  "Goodness of fit" and "Structure Diagnostics" blocks for
  `summary.nmfkc()`, `summary.nmfae()`, and `summary.nmfkc.net()`
  (incl. the signed variant).  Labels are padded to a common width so
  values are column-aligned, fields absent from a given model are
  skipped automatically (e.g.\ `nmfkc.net` has no residual SE), and any
  future fit statistic or sparsity row is now added in one place
  instead of per-summary.

### **Effective Rank in all five MU-family summaries**
- `summary()` now reports the **Effective Rank** as `x.xx / Q  (NN.N%)`
  -- the absolute value, the nominal rank, and the utilization ratio
  `effective.rank / Q` as a percentage -- for
  `nmfkc()`, `nmfkc.net()`, `nmfae()`, `nmf.ffb()` / `nmf.sem()`, and
  `nmfre()` — previously only `nmfkc()` showed it.  Each is computed by
  the new shared internal helper `.effective.rank(B)` from the model's
  natural \eqn{Q \times N} coefficient/score matrix: the coefficients
  \eqn{B} (`nmfkc`), the latent encoding \eqn{H} (`nmfae`), the node
  membership \eqn{X^\top} (`nmfkc.net`), the latent scores
  \eqn{C_1 Y_1 + C_2 Y_2} (`nmf.ffb`), and the BLUP scores
  \eqn{\Theta A + U} (`nmfre`).  `NA` at \eqn{Q = 1}.

### **Rank-selection diagnostics: silhouette / CPCC fixed, IC removed**
- **`silhouette` is now computed in the original data space.**  It used
  to be evaluated on the rank-\eqn{Q} `B.prob` simplex, whose dimension
  changes with \eqn{Q}; that made it monotone in \eqn{Q} (always
  favouring the smallest rank) and hid genuine cluster structure.  It is
  now the standard mean silhouette width over `dist(t(Y))` (the fixed
  original-data sample distances) with the per-sample hard labels — the
  k-means convention.  On data with real clusters it now shows an
  interior optimum (e.g. the road-OD network peaks at the same rank as
  the cross-validation minimum).
- **`CPCC` is now the classic cophenetic correlation of `dist(t(B))`.**
  It used to be computed from the soft co-membership `t(B.prob) %*%
  B.prob`, which was nearly flat across \eqn{Q}.  It is now
  `cor(dist(t(B)), cophenetic(hclust(dist(t(B)))))` — how well a
  hierarchical clustering of the rank-\eqn{Q} coefficient distances
  reproduces those distances (Sokal & Rohlf).  It now varies with
  \eqn{Q} and recovers an interior optimum.
- **Removed `ICp`, `AIC`, and `BIC`** from `nmfkc()`'s `criterion` list,
  from `summary.nmfkc()`, and from `nmfkc.rank()`'s table.  Empirically
  (across three real datasets) `ICp` was monotone increasing (always
  selecting \eqn{Q=1}) and `AIC` monotone decreasing (always selecting
  the largest \eqn{Q}); for NMF, where the parameter count grows as
  \eqn{Q(P+N)}, these information criteria do not have a usable interior
  optimum, so they were misleading rather than informative.
- The internal helper `.silhouette.simple()` (centroid-approximate, took
  a `B.prob` matrix) was replaced by `.silhouette.mean(D, labels)`,
  which returns the exact mean silhouette width from a distance matrix
  and labels.

### **Breaking change: symmetric NMF removed from `nmfkc()`**
- The `Y.symmetric = "bi" / "tri"` option (deprecated in v0.7.x) has been
  **removed** from `nmfkc()` and `nmfkc.ecv()`.  Symmetric NMF of network
  data now lives exclusively in the dedicated `nmfkc.net()` /
  `nmfkc.net.ecv()` functions, which use the correct Frobenius
  bilateral-gradient updates.  Passing `Y.symmetric` to `nmfkc()` or
  `nmfkc.ecv()` now stops with a message pointing to the replacement:
  `nmfkc.net(Y, rank, type = "tri")` (types `"tri"`, `"bi"`, `"signed"`).
  This also removes the bi/tri code branches (cube-root damping, fixed
  `C = I`, tri C-update, upper-triangle CV folds) from `nmfkc()`,
  simplifying the core function.

### **New diagnostic: effective rank**
- `nmfkc()` now reports `criterion$effective.rank`, the **effective
  rank** of the fit: `exp` of the Shannon entropy of the
  explained-variance distribution
  `p_k = var(B[k, ]) / sum_j var(B[j, ])`.  By the trace identity
  `sum_k var(B[k, ]) = tr(Cov(B))`, each `p_k` is the exact fraction
  of the total coefficient variance carried by factor `k`, so the
  entropy is a genuine additive decomposition (variances add;
  standard deviations do not, which is why variance — not sd — is the
  natural partner for the entropy here).  It ranges in `[1, Q]` and
  counts how many latent factors actively shape across-sample
  variation (dead, zero-variance factors drop out).  This is the
  PCA-style explained-variance / effective-dimensionality measure and
  reuses the `exp(entropy)` functional form of Roy & Vetterli (2007).
- `summary.nmfkc()` prints `Effective Rank: x.xx / Q`.
- `nmfkc.rank()` adds an `effective.rank` column to its criteria table.
  When effective rank plateaus well below the nominal rank, the extra
  factors are not carrying additional coefficient variance — a signal
  that the rank is over-specified.
- `nmfkc.rank(plot = TRUE)` overlays an `eff.rank` curve (effective
  rank divided by nominal rank, in `[0, 1]`, solid green line) on the
  diagnostics plot.  A peak in this utilization curve marks the rank at
  which the latent factors carry the most evenly distributed variance.

### **Diagnostics cleanup: B.prob crispness metrics**
- Removed `B.prob.sd.min` and `B.prob.entropy.mean` from `nmfkc()`'s
  `criterion` list, from `summary.nmfkc()`, and from `nmfkc.rank()`'s
  criteria table and plot.  All three `B.prob.*` peakedness metrics are
  monotone in the rank `Q`, so they carry no peak/elbow signal for rank
  selection (verified empirically); the principled rank signals are
  ECV, the R-squared elbow, and the new `effective.rank` utilization.
- `B.prob.max.mean` (clustering crispness) is retained, but only in
  `summary.nmfkc()` ("Clustering Crispness") and the `criterion` list.
  At a fixed `Q` it remains a useful confidence check — the mean
  dominant-cluster membership — before treating `B.cluster` as hard
  labels.  It is no longer shown in `nmfkc.rank()` (cross-`Q`), where
  its `1/Q` baseline shift makes it misleading.
- `summary.nmfkc()` no longer prints "Clustering Entropy" (it duplicated
  the crispness information).

### **Improvements**
- **Unified three-variant R² across all NMF functions.**  Every NMF
  variant (`nmfkc()`, `nmfae()`, `nmfae.signed()`, `nmfkc.net()`,
  `nmfkc.signed()`, `nmfre()`) now returns three goodness-of-fit
  summaries on the same scale, computed by the new internal helper
  `.r.squared.all()`:
  - `r.squared`: Pearson \eqn{\mathrm{cor}(Y, \widehat Y)^2}
    (scale-invariant, in \eqn{[0, 1]}).  Unchanged from before.
  - `r.squared.uncentered`: \eqn{1 - \|Y - \widehat Y\|_F^2 / \|Y\|_F^2}.
    Baseline = the zero matrix (natural for non-negative factorizations
    without an intercept); matches the "uncentered R²" of intercept-free
    regression.
  - `r.squared.centered`: \eqn{1 - \|Y - \widehat Y\|_F^2 / \|Y - \bar Y_{p\cdot}\|_F^2}.
    Baseline = per-row mean; the standard ("centered") multivariate-
    regression \eqn{R^2}; equals 0 when the model predicts the row mean.

  The two suffixed variants differ only in their baseline (denominator);
  both use the Frobenius norm in the numerator.  Naming follows the
  centered/uncentered \eqn{R^2} distinction used by statistics software
  (e.g. statsmodels).
  All three respect `Y.weights == 0` masking (the standard NA-hold-out
  convention).  For `nmfre()` the same three variants are also
  reported on the fixed-only prediction as `r.squared.fixed.*`.
  Displayed by all `summary.*` methods.

### **Bug Fixes**
- `nmfkc.net()`: `r.squared` now correctly excludes weight-zero (NA-masked)
  entries when `Y.weights` is supplied or auto-masking is in effect,
  matching the convention used by `nmfkc()`, `nmfae()`, `nmfae.signed()`,
  and `nmfkc.signed()`.  Previously the correlation was computed over the
  full matrix including replaced-NA cells, giving a distorted r.squared.

### **Documentation**
- `nmfkc()`: removed Examples 3 & 4 (deprecated `Y.symmetric = "bi"/"tri"`);
  the documentation now points users to `\link{nmfkc.net}()` for symmetric
  NMF.
- `summary.nmf.sem()`: example code, `@param`, and `@seealso` updated to
  use the canonical `nmf.ffb` name (the S3 method continues to dispatch
  correctly via `c("nmf.ffb", "nmf.sem")` inheritance).

# nmfkc 0.7.3

### **Documentation**
- README and `nmf-sem-with-nmfkc.Rmd` vignette code now reference the
  canonical `nmf.ffb.*` aliases (`nmf.ffb()`, `nmf.ffb.cv()`,
  `nmf.ffb.DOT()`) instead of the legacy `nmf.sem.*` names.  Both
  names continue to work; the change only affects what users see on
  the GitHub Pages homepage and in the vignette source.

# nmfkc 0.7.2

### **Headline: NMF-FFB rebrand and full bootstrap inference**
- `nmf.ffb*` family added as the canonical alias for `nmf.sem*` (Satoh
  2025, arXiv:2512.18250 adopts "NMF-FFB" — Non-negative Matrix
  Factorization with Feed-Forward + Feedback — as the model's
  canonical name).  `nmf.sem*` continues to work and shares the same
  return classes (`c("nmf.ffb", "nmf.sem")` and
  `c("nmf.ffb.inference", "nmf.sem.inference", ...)`), so existing
  scripts are unaffected.
- `nmf.sem.inference()` / `nmf.ffb.inference()`: replaced the legacy
  1-step Newton wild bootstrap with a **full X-fixed pair bootstrap**.
  Resamples columns of (Y1, Y2), refits (C1, C2) with X held at the
  original fit, and reports per-element `support_rate = mean(|c_b| >
  threshold)` together with percentile CIs.  Significance markers
  (`*` / `**` / `***` at sup > 0.95 / 0.99 / 0.999) follow the lavaan
  convention.  Both Theta_1 (feedback) and Theta_2 (exogenous) are
  inference targets (previous version covered only Theta_2).
- `nmf.sem()` / `nmf.ffb()`: now runs `nmfkc(Y1, A = Y2)`
  **internally by default** when `X.init` is a string method,
  forwarding `X.init`, `X.L2.ortho`, `epsilon`, `maxit`, `seed`.  The
  feedforward fit is used both as the X warm-start and as the
  baseline for `SC.map`.  `nmfkc.baseline = FALSE` opts out.

### **Bug Fixes**
- `nmf.sem.inference()`: fixed dimension bug in the Leontief identity
  matrix (`I_mat <- diag(Q)` should have been `diag(P1)`); previously
  every replicate was silently marked invalid when `P1 != Q`.
- `nmfkc.net()`: now auto-masks NA entries of `Y` (parity with the
  other four NMF variants); previously errored at the `min(Y) < 0`
  check when `Y` contained NA.
- `nmfkc()`: Fixed C matrix asymmetry in tri-symmetric NMF (`Y.symmetric = "tri"`). The C update was using stale B and XB computed from the old X; now B and XB are recomputed after X is updated. Also fixed column reordering to permute both rows and columns of C. Previously the relative asymmetry could reach ~46%; now it is at machine precision (~1e-14).

### **Improvements**
- `Y.weights` semantics unified to `lm()`-style weighted least squares
  across `nmfkc()`, `nmfae()`, `nmfkc.net()`, `nmfkc.signed()`,
  `nmfae.signed()`: loss is now `sum(W * (Y - Yhat)^2)` (linear in W,
  matching `lm()`'s `weights` argument).  Binary masks (W ∈ {0, 1};
  the standard ECV / NA-mask case) are unaffected since W = W^2.
- All MU functions now emit a `"maximum iterations (N) reached..."`
  warning when `maxit` is exhausted without meeting the relative-
  tolerance criterion (previously silent in `nmfae`, `nmfae.signed`,
  `nmfkc.net`, `nmfkc.signed`, `nmfre`, and `nmf.sem`).
- All MU functions now share `maxit = 5000` as the default (was 5000
  / 20000 / 50000 inconsistently).  Together with the maxit warning
  above, users see explicit feedback when 5000 is insufficient and
  can opt into a larger cap.
- New shared internal helper `.init_X_method()` for X initialization
  via `"nndsvd"` / `"kmeans"` / `"kmeansar"` / `"runif"` / numeric
  matrix.  All NMF families now use the same dispatch logic; previous
  ad-hoc inline implementations are removed.
- `nmf.sem()` returns `SC.map` (input-output structural fidelity:
  correlation between the equilibrium operator and the feedforward
  baseline mapping; Satoh 2025 §4.SC.map) automatically when
  `nmfkc.baseline` is supplied or computed internally.
- `summary.nmf.sem()`: rewritten to display the full-bootstrap
  inference output — separate Theta_1 / Theta_2 blocks with
  `Estimate | CI_low | CI_high | support | Pr(>0) | sig`, plus a
  bootstrap meta-info header.
- `coef.nmf.sem()`: now returns a long-format data frame with rows
  for every entry of both C1 and C2 (`Type | Basis | Covariate |
  Estimate`); previously returned only the C2 matrix when no
  inference had been run.  Schema matches the inference-augmented
  output for uniformity.
- `plot.nmf.sem()`: default trace is now `objfunc.full` (loss +
  penalties — the actual monotonically-decreasing quantity that the
  multiplicative updates minimize) instead of `objfunc` (reconstruction
  only).  New argument `which = "full" | "reconstruction" | "both"`.
- `nmf.sem.DOT()`: significance stars now appear on Theta_1 (feedback
  Y1 → F) edges in addition to Theta_2 (exogenous Y2 → F); X (F → Y1)
  edges remain unstarred since the basis is not the inference target.
- `plot.nmfae.ecv()`: Heatmap cell text color is now always black for better readability on light-colored cells.
- `nmfkc()`: `X.init = "runif"` now supports `nstart > 1` for multi-start initialization. Multiple random starting points are evaluated with 10 standard NMF iterations, and the best (lowest Frobenius error) is selected.
- `nmfae()`, `nmfre()`: `r.squared` is now computed as `cor(Y, fitted)^2` (squared correlation between observed and fitted values), consistent with `nmfkc()`. Previously `nmfae()` used `1 - SS_res/SS_tot` and `nmfre()` used the same regression-style R-squared, which can behave unexpectedly for intercept-free non-negative models.
- `nmfkc.kernel.beta.nearest.med()`: added a `candidates` argument controlling the bandwidth grid. Options: `"7points"` (new default, `t = {-1,-2/3,-1/3,0,1/3,2/3,1}`), `"4points"` (`t = {-1/2, 0, 1/2, 1}`), or a user-supplied numeric vector of \eqn{t} values. Previously the grid silently differed between the no-landmark (`Uk = NULL`; 4 points) and landmark (7 points) branches.

### **New Functions (Signed NMF family)**
- `nmfkc.signed()`: NMF-KC with signed covariate/coefficient.  Model \eqn{Y \approx X \Theta A} with \eqn{X \ge 0}, \eqn{\Theta = C_{+} - C_{-}} (signed), \eqn{A} real-valued.  Uses Ding et al. (2010) sign-splitting + Direct MU; \eqn{Y} may also contain negative entries (semi-NMF regression).  Supports `Y.weights` for element-wise masking.
- `nmfkc.signed.cv()`, `nmfkc.signed.ecv()`: column-wise and element-wise k-fold CV for rank selection on signed data.
- `nmfae.signed()`: Three-layer autoencoder with \strong{signed bottleneck} \eqn{Y_1 \approx X_1 (C_{+} - C_{-}) X_2 Y_2}.  \eqn{X_1, X_2 \ge 0} preserve soft clustering on both decoder and encoder sides while the bottleneck \eqn{\Theta} can carry negative weights (e.g., anti-correlated properties).  Hybrid warm-start (from `nmfae()`) + Direct MU with multi-restart.
- `nmfae.signed.ecv()`: element-wise CV for (decoder-rank, encoder-rank) selection.
- `nmfae.signed.inference()`: sandwich SE + wild bootstrap for \eqn{\Theta} (no non-negativity projection on \eqn{\Theta} since it is signed).
- S3 methods `predict.*.signed()`, `plot.*.signed()`, `summary.*.signed()`, and `nmfae.signed.rename()` helper.

### **New Functions (Network NMF family)**
- `nmfkc.net()`: Single unified entry point for symmetric NMF of network data, with `type = "tri" | "bi" | "signed"`.  All three variants use the Frobenius-full bilateral gradient (supersedes the one-sided approximation in `nmfkc(Y.symmetric = ...)`).  `type = "signed"` supports signed \eqn{C = C_{+} - C_{-}} via Ding et al. (2010) sign-splitting, preserving \eqn{X \ge 0} for soft clustering while allowing inter-cluster repulsion.  The returned object's fields are uniform across types: \code{$Cp} and \code{$Cn} are \code{NULL} for tri/bi, and populated matrices for signed.  \code{$C} is always populated (identity for bi, non-negative for tri, signed for signed).
- `nmfkc.net.ecv()`: Element-wise cross-validation with upper-triangle folds (mirrored to the lower triangle to prevent symmetry leakage). Unified entry point for `type = "tri" | "bi" | "signed"` (calls `nmfkc.net()` with the matching `type` for each fold).
- `nmfkc.net.DOT()`: Graphviz DOT visualization for symmetric NMF networks. Displays basis-to-node membership edges and inter-basis interaction edges (C matrix) with significance stars. Now has `signed` parameter (auto-detected from class) to render negative `C` entries as dashed edges.
- `nmfkc.net.inference()`: Statistical inference for symmetric NMF. Wrapper around `nmfkc.inference()` with `A = t(X)`. Returns off-diagonal C coefficients with sandwich SE and wild bootstrap.

### **Deprecations**
- `nmfkc(Y, Y.symmetric = "bi"|"tri")`: **Deprecated** in favor of `nmfkc.net(Y, type = "bi"|"tri")`.  The old implementation uses a one-sided gradient approximation that empirically converges for \eqn{C \ge 0} but is theoretically incorrect and does not extend to signed \eqn{C}.  The deprecated branch still works in v0.6.8 (with a deprecation warning) and will be removed in a future release.

### **Parameter Renames** (old names remain usable for backward compatibility)
- `nmf.sem.DOT()`: `weight_scale_y2f` → `weight_scale_c2`, `weight_scale_fy1` → `weight_scale_x1` (matrix-name-based naming, consistent with `nmfae.DOT()` and `nmfkc.DOT()`).
- `nmf.sem.DOT()`: `sig.level` moved to after `threshold` for consistency with other `.DOT` functions.

### **Documentation**
- README, vignettes, and roxygen `@title` / `@description` updated to
  use **NMF-FFB** as the canonical model name (with "(formerly
  NMF-SEM)" attached on first mention for discoverability of the
  legacy term).  File names (`R/nmf.sem.R`, `vignettes/nmf-sem-with-
  nmfkc.Rmd`, `man/nmf.sem.Rd`), function names (`nmf.sem*`), and S3
  classes (`"nmf.sem"`) are unchanged so URLs and existing scripts
  continue to work.

# nmfkc 0.6.7

### **Bug Fixes**
- Added `fitted.nmfae()` and `residuals.nmfae()` S3 methods; previously `fitted()` on an `nmfae` object silently returned `NULL` because the wrong field name (`$XB` instead of `$Y1hat`) was used.

### **Naming Unification** (old names remain usable for backward compatibility)
- Coefficient tables: all inference functions now use `Basis` / `Covariate` columns (was `Factor`/`Exogenous` in `nmf.sem.inference()`, `Decoder`/`Encoder` in `nmfae.inference()`).
- Wild bootstrap defaults unified: `wild.B = 500`, `wild.seed = 123` across all inference functions.
- First argument of all `.DOT` functions renamed to `result` for consistency.
- CV tuning parameters (`nfolds`, `seed`, `shuffle`) moved to `...` in `nmfkc.ecv()`, `nmfae.ecv()`, `nmfae.cv()`, `nmf.sem.cv()`; `div` also accepted for backward compatibility.

# nmfkc 0.6.6

### **New Functions**
- `nmfkc.criterion()`: Extracted criterion computation from `nmfkc()` as a standalone exported function. Supports `detail = "full"` / `"fast"` / `"minimal"` to control computation cost.
- `nmfre.inference()`: Separated statistical inference from `nmfre()` optimization. Returns coefficient table with SE, z-values, and p-values via wild bootstrap.
- `nmf.sem.inference()`: Statistical inference for the C2 parameter matrix in NMF-SEM. Uses sandwich SE and wild bootstrap.
- S3 methods `coef()`, `fitted()`, `residuals()` for all model classes (`nmfkc`, `nmfae`, `nmfre`, `nmf.sem`).
- S3 methods `plot()` for `nmfre` and `nmf.sem` (convergence diagnostics).
- `summary.nmf.sem()`: Stability diagnostics, fit statistics, and C2 coefficient table.

### **Parameter Renames** (old names remain usable for backward compatibility)
- `nmfkc()`, `nmfkc.rank()`: `save.time` / `save.memory` → `detail`
- `nmfae()`: `Q` → `rank`, `R` → `rank.encoder`
- `nmfre()`: `Q` → `rank`, `dfU.cap.rate` → `df.rate`
- `nmfre.dfU.scan()`, `nmfkc.ar.degree.cv()`: `Q` → `rank`
- `nmfkc.residual.plot()`: `Y_XB_palette` → `fitted.palette`, `E_palette` → `residual.palette`
- `nmfkc.kernel.beta.nearest.med()`: `block_size` → `block.size`, `sample_size` → `sample.size`

### **Other Improvements**
- `hide.isolated` option added to all `.DOT` functions (default `TRUE`).
- `nmf.sem.DOT()`: Added `sig.level` parameter; C2 edges decorated with significance stars.
- `nmfkc()`: Added `X.restriction = "none"` option and `X.init = "kmeansar"` initialization.
- Added arXiv/DOI references to roxygen documentation for all main functions.
- `@section Lifecycle: Experimental` added to `nmfae()`.
- Removed `mc.cores` parallel option from `nmfae.ecv()` for CRAN compliance.

# nmfkc 0.6.0
### **Bug Fixes**
- Fixed variable `T` shadowing `TRUE` in information criterion computation.
- Fixed `nmfkc.ecv()` to use KL divergence for evaluation when `method="KL"`.
- Added performance flags (`save.time=TRUE`) to `nmfkc.ecv()` inner calls.
- Fixed zero-division in `nmfkc.rank()` elbow normalization when R-squared values are identical.
- Fixed parameter name mismatch (`rank` → `Q`) in `nmfkc.rank()` call to `nmfkc.ecv()`.
- Fixed descending loop in `nmf.sem.split()` when P=2.
- Added input validation for `n.exogenous` in `nmf.sem.split()`.

### **Documentation**
- Added roxygen documentation for `summary.nmfkc()` and `print.summary.nmfkc()`.
- Added `@return` for `plot.nmfkc()` and `predict.nmfkc()`.
- Added missing `@return` items (`method`, `n.missing`, `n.total`, `rank`, `mae`) to `nmfkc()`.

### **Code Quality**
- Replaced `T`/`F` with `TRUE`/`FALSE`.
- Replaced `1:length()` with `seq_along()`.
- Changed default font from Meiryo to Arial in DOT functions.
- Aligned `nmf.sem.cv()` defaults with `nmf.sem()`.

# nmfkc 0.5.8
### **Graphviz DOT Output Consolidation and Cleanup**
- Harmonized all DOT-generating functions (`nmf.sem.DOT`, `nmfkc.DOT`, `nmfkc.ar.DOT`) for consistent structure, naming conventions, and visualization logic.
- Standardized node and edge formatting rules, including unified cluster behavior, color schemes, and edge-scaling conventions.
- Implemented threshold-aware coefficient labeling so that displayed numerical precision aligns with the visualization threshold, preventing misleadingly detailed labels.
- Removed unused or redundant DOT fragments and improved compatibility across Graphviz engines.
- Enhanced layout readability through consistent indentation, node grouping, and suppression of isolated nodes in specific visualization modes (e.g., `type = "YA"` in `nmfkc.DOT`).
- Refactored and expanded internal DOT helper functions (`.nmfkc_dot_format_coef`, `.nmfkc_dot_digits_from_threshold`, `.nmfkc_dot_cluster_nodes`, etc.) for better maintainability and uniform behavior.

* **New Function:** Implemented `nmfkc.ecv()` for Element-wise Cross-Validation (Wold's CV).
  - This function randomly masks elements of the observation matrix to evaluate structural reconstruction error.
  - It provides a statistically robust criterion for rank selection, avoiding the monotonic error decrease often seen in standard column-wise CV.
  - Supports vector input for `rank` to evaluate multiple ranks simultaneously.
* **Missing Value & Weight Support:**
  - `nmfkc()` and `nmfkc.cv()` now fully support missing values (`NA`) and observation weights via the hidden argument `Y.weights` (passed through `...`).
  - If `Y` contains `NA`s, they are automatically detected and masked (assigned a weight of 0) during optimization.
* **Rank Selection Diagnostics (`nmfkc.rank`):**
  - **Dual-Axis Visualization:** The plot now displays fitting metrics ($R^2$, etc.) on the left axis and ECV Sigma (RMSE) on the right axis (blue line).
  - **Automatic Best Rank labeling:** The plot explicitly marks the "Best" rank based on two criteria:
    - **Elbow:** Geometric elbow point of the $R^2$ curve.
    - **Min:** Minimum error point of the Element-wise CV.
  - `save.time` defaults to `FALSE`, enabling the robust Element-wise CV calculation by default.
* **Argument Standardization:**
  - Unified the rank argument name to `rank` across all functions (`nmfkc`, `nmfkc.cv`, `nmfkc.ecv`, `nmfkc.rank`).
  - The legacy argument `Q` is still supported for backward compatibility but internally mapped to `rank`.
* **Summary Improvements:**
  - Updated `summary()` and `print()` methods to report:
    - Sparsity of Basis ($X$) and Coefficients ($B$).
    - Clustering Entropy (indicating "Crisp" vs "Ambiguous" clustering).
    - Clustering Crispness (Mean Max Probability).
    - Number and percentage of missing values in $Y$.
* **Other Improvements:**
  - Added a validation check in `nmfkc.ar()` to ensure the input `Y` has no missing values (as they cannot be propagated to the covariate matrix `A` in VAR models).
  - Refined `nmfkc.residual.plot()` layout margins for better visibility of titles.
  - Updated documentation to reflect all changes.
* Regularization Update:  
  The regularization scheme has been revised from **L2 (ridge)** to **L1 (lasso-type)** penalties.  
  - `gamma` now controls the **L1 penalty on the coefficient matrix** \( B = C A \), promoting sparsity in sample-wise coefficients.  
  - A new argument `lambda` has been added to control the **L1 penalty on the parameter matrix** \( C \), encouraging sparsity in the shared template structure.  
  Both parameters can be passed through the ellipsis (`...`) to `nmfkc()` and related functions.  
* Function Signature Simplification:** Many less-frequently used arguments in `nmfkc()` (e.g., `gamma`, `X.restriction`, `X.init`) and in `nmfkc.cv()` (e.g., `div`, `seed`) have been moved into the ellipsis (`...`) for a cleaner function signature.
* Performance Improvement: The internal function `.silhouette.simple` was vectorized and optimized to reduce computational cost, particularly for the calculation of `a(i)` and `b(i)`.
* Removed the `fast.calc` option from the `nmfkc()` function.
* Added the `X.init` argument to the `nmfkc()` function, allowing selection between `'kmeans'` and `'nndsvd'` initialization methods.
* The penalty term has been changed from `tr(CC')` to `tr(BB')` = `tr(CAA'C')`.
* Implemented the internal `.z` and `xnorm` functions.
* Added the fast.calc option to the `nmfkc()` function.
* Optimized internal calculations for improved performance.
* Updated `citation("nmfkc")` and added AIC/BIC to the output.
* Implemented the `nmfkc.ar.stationarity()` function.
* Modified the `z()` function.
* Used `crossprod()` for faster matrix multiplication.
* Implemented the `nmfkc.ar.DOT()` function.
* Added logic to sort the columns of `X` to form a unit matrix in special cases.
* Implemented `nmfkc.kernel.beta.cv()` and `nmfkc.ar.degree.cv()` functions.
* Set the default column names of `X` to `Basis1`, `Basis2`, etc.
* Added `X.prob` and `X.cluster` to the return object.
* Skipped CPCC and silhouette calculations when `save.time = TRUE`.
* Added a prototype for the `nmfkc.ar()` function.
* Added the `criterion` argument to the `nmfkc()` function to support multiple criteria.
* Updated the `nmfkc.rank()` function.
* Added the `criterion` argument to the `nmfkc.rank()` function.
* Implemented the `save.time` argument.
* Implemented the `nmfkc.rank()` function.
* Implemented the `nstart` option from the `kmeans()` function.
* Added an experimental implementation of the `nmfkc.rank()` function.
* Removed zero-variance columns and rows with a warning.
* Added source and references to the documentation.
* Renamed several components for clarity:
    * `nmfkcreg` to `nmfkc`
    * `create.kernel` to `nmfkc.kernel`
    * `nmfkcreg.cv` to `nmfkc.cv`
    * `P` to `B.prob`
    * `cluster` to `B.cluster`
    * `unit` to `X.column`
    * `trace` to `print.trace`
    * `dims` to `print.dims`
* Added the `r.squared` argument to the `nmfkcreg.cv()` function.
* In `nmfkcreg()`:
    * Added the `dims` argument to check matrix sizes.
    * Added the `unit` argument to normalize the basis matrix columns.
* Modified the `create.kernel()` function to support prediction.
* Updated examples on GitHub.
* Removed the `YHAT` return value; use `XB` instead.
* Added the `cluster` return value for hard clustering.

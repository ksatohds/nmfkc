# nmfkc 0.6.8

### **Bug Fixes**
- `nmfkc()`: Fixed C matrix asymmetry in tri-symmetric NMF (`Y.symmetric = "tri"`). The C update was using stale B and XB computed from the old X; now B and XB are recomputed after X is updated. Also fixed column reordering to permute both rows and columns of C. Previously the relative asymmetry could reach ~46%; now it is at machine precision (~1e-14).

### **Improvements**
- `plot.nmfae.ecv()`: Heatmap cell text color is now always black for better readability on light-colored cells.
- `nmfkc()`: `X.init = "runif"` now supports `nstart > 1` for multi-start initialization. Multiple random starting points are evaluated with 10 standard NMF iterations, and the best (lowest Frobenius error) is selected.
- `nmfae()`, `nmfre()`: `r.squared` is now computed as `cor(Y, fitted)^2` (squared correlation between observed and fitted values), consistent with `nmfkc()`. Previously `nmfae()` used `1 - SS_res/SS_tot` and `nmfre()` used the same regression-style R-squared, which can behave unexpectedly for intercept-free non-negative models.
- `nmfkc.kernel.beta.nearest.med()`: added a `candidates` argument controlling the bandwidth grid. Options: `"7points"` (new default, `t = {-1,-2/3,-1/3,0,1/3,2/3,1}`), `"4points"` (`t = {-1/2, 0, 1/2, 1}`), or a user-supplied numeric vector of \eqn{t} values. Previously the grid silently differed between the no-landmark (`Uk = NULL`; 4 points) and landmark (7 points) branches.

### **New Functions**
- `nmfkc.net()`: Dedicated entry point for symmetric NMF of network data (\eqn{Y \approx X C X^\top}, \eqn{X, C \ge 0}).  Uses the correct Frobenius-full bilateral gradient (supersedes the one-sided approximation in `nmfkc(Y.symmetric = ...)`).  Supports `type = "tri"` (default) and `type = "bi"` (with cube-root damping, He et al. 2011).  `C` is symmetric by design; initialization is symmetrized so the bilateral gradient is exact throughout iteration.
- `nmfkc.net.signed()`: Symmetric NMF with signed `C = Cp - Cn` (\eqn{X \ge 0} preserves soft clustering; \eqn{C} allows inter-cluster repulsion).  Uses bilateral gradient + Ding et al. (2010) sign-splitting.
- `nmfkc.net.ecv()`, `nmfkc.net.signed.ecv()`: Element-wise cross-validation with upper-triangle folds (mirrored to the lower triangle to prevent symmetry leakage).
- `nmfkc.net.DOT()`: Graphviz DOT visualization for symmetric NMF networks. Displays basis-to-node membership edges and inter-basis interaction edges (C matrix) with significance stars. Now has `signed` parameter (auto-detected from class) to render negative `C` entries as dashed edges.
- `nmfkc.net.inference()`: Statistical inference for symmetric NMF. Wrapper around `nmfkc.inference()` with `A = t(X)`. Returns off-diagonal C coefficients with sandwich SE and wild bootstrap.

### **Deprecations**
- `nmfkc(Y, Y.symmetric = "bi"|"tri")`: **Deprecated** in favor of `nmfkc.net(Y, type = "bi"|"tri")`.  The old implementation uses a one-sided gradient approximation that empirically converges for \eqn{C \ge 0} but is theoretically incorrect and does not extend to signed \eqn{C}.  The deprecated branch still works in v0.6.8 (with a deprecation warning) and will be removed in a future release.

### **Parameter Renames** (old names remain usable for backward compatibility)
- `nmf.sem.DOT()`: `weight_scale_y2f` → `weight_scale_c2`, `weight_scale_fy1` → `weight_scale_x1` (matrix-name-based naming, consistent with `nmfae.DOT()` and `nmfkc.DOT()`).
- `nmf.sem.DOT()`: `sig.level` moved to after `threshold` for consistency with other `.DOT` functions.

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

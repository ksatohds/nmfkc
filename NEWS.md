# nmfkc 0.5.5
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

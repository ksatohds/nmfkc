# nmfkc 0.5.4
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

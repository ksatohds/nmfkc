# nmfkc 0.5.0
0.5.0: Set values smaller than `.Machine$double.eps` to 0 to avoid negative results in `z()`.
0.4.9: Optimized internal calculations for improved performance.
0.4.8: Updated `citation("nmfkc")` and added AIC/BIC to the output.
0.4.7: Implemented the `nmfkc.ar.stationarity()` function.
0.4.6: Modified the `z()` function.
0.4.5: Used `crossprod()` for faster matrix multiplication.
0.4.4: Implemented the `nmfkc.ar.DOT()` function.
0.4.3: Added logic to sort the columns of `X` to form a unit matrix in special cases.
0.4.2: Implemented `nmfkc.kernel.beta.cv()` and `nmfkc.ar.degree.cv()` functions.
0.4.1: Set the default column names of `X` to `Basis1`, `Basis2`, etc.
0.3.9: Added `X.prob` and `X.cluster` to the return object.
0.3.8: Skipped CPCC and silhouette calculations when `save.time = TRUE`.
0.3.7: Added a prototype for the `nmfkc.ar()` function.
0.3.5: Added the `criterion` argument to the `nmfkc()` function to support multiple criteria.
0.3.4: Updated the `nmfkc.rank()` function.
0.3.2: Added the `criterion` argument to the `nmfkc.rank()` function.
0.3.0: Implemented the `save.time` argument.
0.2.9: Implemented the `nmfkc.rank()` function.
0.2.8: Implemented the `nstart` option from the `kmeans()` function.
0.2.7: Added an experimental implementation of the `nmfkc.rank()` function.
0.2.6: Removed zero-variance columns and rows with a warning.
0.2.5: Added source and references to the documentation.
0.2.4: Renamed several components for clarity:
     `nmfkcreg` to `nmfkc`
     `create.kernel` to `nmfkc.kernel`
     `nmfkcreg.cv` to `nmfkc.cv`
     `P` to `B.prob`
     `cluster` to `B.cluster`
     `unit` to `X.column`
     `trace` to `print.trace`
     `dims` to `print.dims`
0.2.3: Added the `r.squared` argument to the `nmfkcreg.cv()` function.
0.2.2: In `nmfkcreg()`:
     Added the `dims` argument to check matrix sizes.
     Added the `unit` argument to normalize the basis matrix columns.
0.2.1: Modified the `create.kernel()` function to support prediction.
0.2.0: Updated examples on GitHub.
0.1.6:
     Removed the `YHAT` return value; use `XB` instead.
     Added the `cluster` return value for hard clustering.

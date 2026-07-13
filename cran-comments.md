## R CMD check results

0 errors | 0 warnings | 1 note

## Test environments

* Windows 11 (local), R 4.4.1
* win-builder (R-devel, R-release)
* R-hub v2 (Linux, macOS, Windows)

## Notes

* "checking for future file timestamps ... NOTE" — appears on environments
  that cannot reach the CRAN time server to verify the current time; not a
  package issue.

## This is an update

This is a feature/maintenance update from v0.8.2 (currently on CRAN,
published 2026-06-14) to v0.8.8.  All changes since v0.8.2 are listed in
NEWS.md; the major items are:

* New `nmfre()` family: non-negative matrix factorization as a linear
  mixed model (random-effect basis coefficients), estimated by an
  outer-inner ECM with automatic ridge selection from the variance
  components, plus `nmfre.inference()` / `nmfre.ecv()` and S3 methods.
* Optional MAP penalties rolled out consistently across the multiplicative-
  update fitters: `X.L2.ortho` (basis orthogonality), `X.L2.smooth`
  (row smoothness), `C.L1` (sparsity of the parameter matrix), and `C.L2`
  (ridge on the signed coefficient matrix).  All default to 0 (off).
* New `by` argument on the coefficient-table `print` methods to choose the
  grouping order (by covariate or by basis).
* Canonical names `nmf.rrr*` (three-layer NMF as reduced-rank regression)
  and `nmf.ffb*` (feed-forward + feedback) are now the primary interface;
  the previous `nmfae*` and `nmf.sem*` names are retained as deprecated
  aliases that forward with a `.Deprecated()` note, so existing code and
  saved objects continue to work unchanged.
* Removed the `B.L1` penalty (an L1 on the fitted field `C A` that acted as a
  degenerate global shrinkage); use `C.L1` for sparsity / variable selection.

## Additional checks

* All tests pass (testthat).
* All examples run without errors (including --run-donttest).
* All vignettes build without errors.
* No reverse dependencies on CRAN.

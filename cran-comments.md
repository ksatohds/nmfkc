## R CMD check results

0 errors | 0 warnings | 0 notes (win-builder R-devel and R-release)

## Test environments

* win-builder (R-devel 2026-06-13 r90149, R-release 4.6.0): Status OK
* Windows 11 (local), R 4.4.1: 1 NOTE (see below)
* R-hub v2 (Linux, macOS, Windows)

## Notes

* The only NOTE seen anywhere is "checking for future file timestamps ...
  NOTE" on the local machine; it appears when the environment cannot reach the
  CRAN time server to verify the current time, and is not a package issue.  It
  does not appear on win-builder, where the check is fully clean.

## This is an update

This is a feature/maintenance update from v0.7.3 (currently on CRAN,
published 2026-05-14) to v0.8.2.  All changes since v0.7.3 are listed in
NEWS.md; the major items are:

* New rank-selection engines, complementing the existing element-wise
  cross-validation: `nmfkc.bicv()` (bi-cross-validation, Owen & Perry 2009),
  `nmfkc.consensus()` (consensus clustering with cophenetic / dispersion / PAC,
  Brunet et al. 2004; Kim & Park 2007; Senbabaoglu et al. 2014), and
  `nmfkc.ard()` (automatic relevance determination, Tan & Fevotte 2013).
* `nmfkc.rank()` now reports the broken-stick-corrected effective rank
  (Roy & Vetterli 2007) together with R-squared and the element-wise CV error
  (Wold 1978); the heavier clustering criteria moved to `nmfkc.consensus()`.
* Leaner, safer interfaces for `nmfkc.ard()`, `nmfkc.bicv()` and
  `nmfkc.consensus()` (fine-tuning arguments moved into `...` with the same
  safe defaults; existing named-argument calls are unaffected).
* Bug fix in `nmfkc.net.DOT()`: a `type = "bi"` fit was mis-detected as
  `"tri"` when the C matrix carried dimnames; detection now uses the result's
  `$type`.  The default Graphviz layout is now `"neato"`.
* Two new vignettes ("Choosing the NMF rank on data with a known true rank"
  and "Soft community detection in networks with nmfkc.net").

## Additional checks

* All tests pass (testthat).
* All examples run without errors (including --run-donttest).
* All vignettes build without errors.
* No reverse dependencies on CRAN.

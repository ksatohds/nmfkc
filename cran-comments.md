## R CMD check results

0 errors | 0 warnings | 1 note

## Test environments

* Windows 11 (local), R 4.4.1
* Ubuntu 24.04 LTS (remote server), R 4.5.3
* win-builder (R-devel, R-release)
* R-hub v2 (Linux, macOS, Windows)

## Notes

* "unable to verify current time" — appears on local/server environments
  without internet access to CRAN servers; not a package issue.

## This is an update

This is a feature/maintenance update from v0.6.7 (currently on CRAN,
published 2026-04-15) to v0.7.3.  Highlights of the changes since v0.6.7
are listed in NEWS.md; major items include:

* New `nmf.ffb()` family (canonical alias for `nmf.sem()`; the rebrand
  matches the published manuscript, Satoh 2025 arXiv:2512.18250).  Legacy
  `nmf.sem*` names continue to work and share the same return classes.
* Full X-fixed pair bootstrap inference for the feedback/exogenous
  coefficient matrices (replacing the legacy 1-step Newton wild bootstrap).
* New `nmfkc.signed()` / `nmfae.signed()` family for signed coefficient
  NMF (Ding et al. 2010 sign-splitting; admits negative entries in
  coefficient matrices and observed data).
* `Y.weights` semantics unified to `lm()`-style weighted least squares
  across all five NMF variants.
* All MU functions now share `maxit = 5000` default and emit a warning
  when `maxit` is exhausted without convergence.
* Several bug fixes (NA handling in `nmfkc.net()`, C matrix asymmetry in
  tri-symmetric NMF, dimension bug in `nmf.sem.inference()`).

## Additional checks

* All tests pass (testthat).
* All examples run without errors (including --run-donttest).
* All vignettes build without errors.
* No reverse dependencies on CRAN.

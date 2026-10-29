## This is an update

An update from v0.9.8 (currently on CRAN, published 2026-09-23) to v1.0.2.
All changes are listed in NEWS.md.  (A v1.0.1 was tagged on GitHub on
2026-10-20 but never submitted; its changes are included here.)

### No exported function was added or removed

The 72 exports, the S3 methods and the data sets are those of v0.9.8, and no
formal argument was added or removed.  The new controls are optional and
passed through `...`; with them absent a call gives the result it gave
before (checked bit-identical on several hundred existing calls), with the
single exception listed below.  The main additions:

* `X.anchor` and `X.init = "spa"` (anchor rows of the basis, chosen by the
  successive projection algorithm or given) in `nmfkc()` and the other
  fitters, and a separability diagnostic in `summary()`;
* an identification diagnostic (`identified.approx`, `identification`) in
  the five inference functions: an algebraic check on the point estimate,
  no re-fitting;
* `wild.unit = "block"`, a block wild bootstrap over time points, in
  `nmfkc.ar.latent.inference()` and `nmfkc.inference()`;
* `class.effects` and `X.init = "nmf"` in `nmf.gmm()`;
* `diag.exclude`, `folds`, `seeds` and `pred` in `nmfkc.net.ecv()`.

### Calls whose result changes

* **One deliberate change of a default.**  The re-fitting wild bootstrap of
  `nmfkc.inference(method = "refit")` now draws one multiplier per column of
  `Y` (`wild.unit = "column"`) instead of one per cell (`"element"`), so that
  it treats the columns as the independent units, as the sandwich standard
  errors and the default `"onestep"` engine of the same function already do.
  This changes the bootstrap SEs, intervals and p-values of such calls when
  `Y` has more than one row and `wild.unit` is not set;
  `wild.unit = "element"` reproduces the earlier results exactly, and with
  one row the two units give identical results.  Point estimates, sandwich
  SEs and the `"onestep"` engine are unaffected.  NEWS.md states this.
* Otherwise only calls that previously failed or silently ignored an
  argument:
  * `nmfkc()` with two or more all-zero columns of `Y` stopped with
    "missing value where TRUE/FALSE needed"; the initialization is now
    repaired.
  * A `Y.weights` given to `nmfkc.ecv()`, `nmfkc.signed.ecv()`,
    `nmf.rrr.ecv()`, `nmf.rrr.signed.ecv()` or `nmfkc.net.ecv()` was dropped
    without a word; it is now used.
  * `X.anchor` given to `nmfre()`, `nmf.ffb()`, `nmf.rrr()` or `nmf.gmm()`
    was dropped without a word; it is now imposed.

The `Description` field no longer names `nmf.sem()`, which v0.9.8 removed.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* win-builder, R-devel: uploaded twice (2026-10-29, 06:26 and 16:42 JST),
  no result was returned; the same tarball passed R-release below
* win-builder, R-release (4.6.1 ucrt): Status OK, install 23 s, check 225 s
* Ubuntu 26.04.1 LTS, R 4.6.1, `--as-cran` with
  `_R_CHECK_CRAN_INCOMING_USE_ASPELL_=true`: Status OK, 278 s
* the same with `_R_CHECK_DEPENDS_ONLY_=true` (no Suggests, as CRAN's
  noSuggests flavour): Status OK, 269 s

The word "FFB" (feed-forward + feedback, the model fitted by `nmf.ffb()`) is
in the package dictionary `.aspell/nmfkc.rds`.

## On the check time

Unchanged since v0.9.6: only `tests/testthat/test-cran-smoke.R` and a few
toy-sized tests run by default; the remaining blocks begin with `skip_unless_full()`
and run when `NMFKC_FULL_TESTS` is set, which is what we run before every
release (4510 assertions, no failures).  Nothing was deleted to achieve
this.

win-builder reports 225 s (R-release) for this submission, against 126 s for
v0.9.8.

## Additional checks

* All tests pass (testthat), with and without `NMFKC_FULL_TESTS`.
* All examples run without errors, including `--run-donttest`.
* All vignettes build without errors.
* No reverse dependencies on CRAN.

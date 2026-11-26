## This is an update

An update from v1.0.2 (currently on CRAN, published 2026-10-30) to v1.0.6.
All changes are listed in NEWS.md.  (Versions 1.0.3 to 1.0.5 were
development versions on GitHub and were never submitted; their changes are
included here.)

### Exports

Two functions were added and none was removed:

* `nmfkc.ident()`, with `print()`, `summary()` and `plot()` methods: whether
  the factorization of an `nmfkc()` fit is unique, its identification region
  and profile-likelihood intervals;
* `nmfkc.net.consensus()`: consensus clustering of the nodes of `nmfkc.net()`
  fits.

No formal argument of an existing function was added or removed; the new
controls are passed through `...`.

### Calls whose result changes

Deliberate changes of defaults, each stated in NEWS.md with the argument
that reproduces the earlier result:

* `nmfre()`: the starts (`X.init = "mix2"`, variances started from the
  data), the exact X-step (`X.update = "exact"`), exact updates for
  `C.signed = FALSE`, and, for `rank = 3`, the geomin member of the
  identification region (`orient = "geomin"`).  These change the fits where
  the earlier defaults stopped short of the optimum or at a worse start.
* `nmfkc.ar.latent.inference()`: the p-values of the entries of the latent
  transition matrices are one-sided by default (`G.p.side = "one.sided"`;
  `"two.sided"` gives the earlier values), since the entries are
  non-negative.
* `summary.nmfkc()` prints anchor purity and an identification verdict in
  place of the separability line; the components `separability`,
  `anchor.purity` and `anchor.row` of `fit$criterion` are replaced by
  `X.anchor.purity` and `X.anchor.row`.

The other changes are additions (options through `...`, new components of
results, `print()` / `plot()` details) and corrections of help pages; where
an existing function was changed internally, its earlier results were checked
to be reproduced bit for bit, as recorded in NEWS.md.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* win-builder, R-devel (2026-11-17 r90926 ucrt): Status OK, install 25 s,
  check 309 s
* win-builder, R-release (4.6.2 ucrt): Status OK, install 23 s, check 326 s
  (examples 21 s, tests 72 s, PDF manual 62 s)
* Ubuntu 26.04.1 LTS, R 4.6.1, `--as-cran` with
  `_R_CHECK_CRAN_INCOMING_USE_ASPELL_=true`: Status OK, 335 s
* the same with `_R_CHECK_DEPENDS_ONLY_=true` (no Suggests, as CRAN's
  noSuggests flavour): Status OK, 300 s

## On the check time

Unchanged since v0.9.6: only `tests/testthat/test-cran-smoke.R` and a few
toy-sized tests run by default; the remaining blocks begin with
`skip_unless_full()` and run when `NMFKC_FULL_TESTS` is set, which is what we
run before every release (5079 assertions, no failures).

## Additional checks

* All tests pass (testthat), with and without `NMFKC_FULL_TESTS`.
* All examples run without errors, including `--run-donttest`.
* All vignettes build without errors.
* No reverse dependencies on CRAN (`tools::package_dependencies(reverse = TRUE)`,
  all dependency types, 2026-11-26).

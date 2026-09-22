## This is an update

An update from v0.9.6 (currently on CRAN, published 2026-08-25) to v0.9.8.
All changes are listed in NEWS.md.  Two of them are worth naming here because
they are visible from outside the package.

### Twenty exported functions were removed

| removed (20) | use instead |
|:--|:--|
| `nmf.sem()`, `nmf.sem.inference()`, `nmf.sem.cv()`, `nmf.sem.split()`, `nmf.sem.DOT()` | `nmf.ffb*` |
| `nmfae()`, `nmfae.inference()`, `nmfae.ecv()`, `nmfae.cv()`, `nmfae.rank()`, `nmfae.DOT()`, `nmfae.heatmap()`, `nmfae.kernel.beta.cv()`, `nmfae.rename()` and the six `nmfae.signed*` counterparts | `nmf.rrr*` / `nmf.rrr.signed*` |

Every one of them was a pure forwarder, e.g.

```r
nmf.sem <- function(...) { .Deprecated("nmf.ffb"); nmf.ffb(...) }
```

They were deprecated in v0.8.8, announced in that release's NEWS, and have
emitted a `.Deprecated()` note through two CRAN releases (0.8.8 on 2026-07-14
and 0.9.6 on 2026-08-25).  The replacement names have been exported since
v0.8.8, so code that moved to them works on the currently published version as
well.  There are no reverse dependencies on CRAN.

### Twelve functions were added

`nmf.gmm()`, `nmf.gmm.inference()`, `nmf.gmm.select()`, `nmf.gmm.twostage()`
(a Gaussian-mixture latent-class extension, new to CRAN with this release and
documented as experimental on every one of its help pages: its interface is
still settling);
`nmf.ffb.ecv()`, `nmf.ffb.test()`, `nmf.ffb.diagnostics()`;
`nmfkc.kernel.gram()`, `nmfkc.rff.beta.cv()`, `nmfkc.rff.positive()`,
`nmfkc.rff.positive.gram()`, `nmfkc.signed.rff.gram()`.

### Other changes a user could notice

* `nmf.ffb()` was renamed internally to the package's house field names.  This
  is a breaking change confined to the likelihood branch, which was added
  after v0.9.6 and has never been on CRAN.
* `X.rowSums.min` and `X.restriction = "rowSums"` are removed from `nmfkc()` /
  `nmfkc.signed()`.  Both acted on the rows of the basis, which is not a gauge
  fix, so the objective was not monotone and a fit could oscillate to `maxit`.
  Both options were added after v0.9.6 and have never been on CRAN.
* Both fitters now report `epsilon.iter` and `objfunc.increases`, so a fit that
  is oscillating rather than descending is visible in `print()` / `summary()`.

## R CMD check results

0 errors | 0 warnings | 0-1 notes (the one note is environmental; see below)

## Test environments

* win-builder, R-devel (2026-09-21 r90579 ucrt): **Status OK**, install 13 s,
  check 126 s
* Windows 11 (local), R 4.4.1, `--as-cran`: 1 note (see below)
* Windows 11 (local), `_R_CHECK_DEPENDS_ONLY_=true` (no Suggests): 1 note (the
  same one)

## Notes

* "checking for future file timestamps ... NOTE" — appears on environments that
  cannot reach the CRAN time server to verify the current time; not a package
  issue.  win-builder reports Status OK with no note at all.

## On the check time

0.9.4 and 0.9.5 were archived by the pretest for exceeding the ten-minute
budget.  The environment-variable arrangement introduced for 0.9.6 is
unchanged: only `tests/testthat/test-cran-smoke.R` runs by default (55
assertions, about 3 seconds), and the remaining blocks begin with
`skip_unless_full()` and run when `NMFKC_FULL_TESTS` is set, which is what we
run before every release (4140 assertions, 566 seconds locally, no failures).
Nothing was deleted to achieve this.

win-builder reports 126 seconds for this submission, against 232 for 0.9.6.

## Additional checks

* All tests pass (testthat): 55 in the default smoke suite, 4140 with
  `NMFKC_FULL_TESTS` set.
* All examples run without errors, including `--run-donttest`.
* All vignettes build without errors.
* No reverse dependencies on CRAN.

## On the submission interval

v0.9.6 was published on 2026-08-25.

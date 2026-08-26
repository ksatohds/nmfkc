## Resubmission (second)

0.9.5 was archived by the pretest for "Overall checktime 11 min > 10 min",
almost all of it `checking tests ... [442s]`.  Thank you for the suggestion to
run the less important tests conditionally on an environment variable; 0.9.6
does exactly that, and nothing else changed.

Only `tests/testthat/test-cran-smoke.R` now runs by default: it exercises
every exported fitter and its S3 methods on 6 x 20 toy matrices, 35
assertions in under a second, with no bootstrap, no cross-validation and no
restarts.  The other 145 blocks begin with `skip_unless_full()` and run when
`NMFKC_FULL_TESTS` is set, which is what we run before every release (3173
assertions, 144 seconds locally).  Nothing was deleted.

Measured here: 1.7 s in CRAN mode, 144 s in full mode, no failures in either.

## Resubmission (first)

This resubmits the 0.9.4 update, which the incoming pretest rejected for its
overall check time (48 minutes on the pretest Windows machine).  0.9.5 is
0.9.4 plus check-time reductions only — no computed value changes: the
expensive bootstrap regression tests now carry skip_on_cran() with CRAN-sized
copies still running there; the timeseries example sweeps fewer lag orders;
and five of the nine vignettes moved to the package website
(<https://ksatohds.github.io/nmfkc/articles/>), keeping four in the tarball.
The pretest's other note, "Possibly misspelled ... Tokuda", is a co-author's
surname (see Notes below).

## R CMD check results

0 errors | 0 warnings | 1-2 notes (see below; both are environmental or a proper name)

## Test environments

* Windows 11 (local), R 4.4.1, `--as-cran`
* Ubuntu Linux (local server), R 4.5.3
* win-builder, R-devel (2026-08-24 r90445): 232 s total, tests 15 s
* R-hub v2 (Linux, macOS, macOS-arm64, Windows, nosuggests)

## Notes

* "checking for future file timestamps ... NOTE" — appears on environments
  that cannot reach the CRAN time server to verify the current time; not a
  package issue.
* "Possibly misspelled words in DESCRIPTION: Tokuda" — a co-author's surname,
  in the citation for <doi:10.48550/arXiv.2607.27474>. Spelled correctly.

## This is an update

This is a maintenance update from v0.8.8 (currently on CRAN, published
2026-07-14) to v0.9.6.  All changes are listed in NEWS.md.  The update is
mostly corrections; the items a user could notice are:

* Correctness fixes in the inference and RNG paths: the refit bootstrap
  returned p-values that were always 0; several inference and
  cross-validation wrappers reset the caller's random stream, so a user loop
  drawing random numbers after a fit silently repeated itself; the weighted
  path of `nmfkc.signed()` could diverge to `Inf`/`NaN`.
* `nmfkc.net.DOT()` judged which nodes to draw on the raw basis rather than
  on the membership scale the edges come from, so with `type = "tri"` the
  whole outer layer could vanish from the graph.
* Uniform `print()` / `summary()` across the fitters, including a single
  convergence line, so a run that merely exhausted `maxit` is no longer
  displayed as if it had converged.
* `nmf.rrr()` now returns the two score matrices `B1` (decoder side) and `B2`
  (encoder side) with their memberships, replacing the single `B.prob` /
  `B.cluster`.  `H` is retained as a deprecated alias of `B1`.
* The `nmf.rrr` family drops its `rank` / `rank.encoder` argument aliases in
  favour of `rank1` / `rank2` (`Q` / `R` still work).  The aliases had to be
  declared as formals *after* `...`, because `rank` is a prefix of both
  `rank1` and `rank2` and R's partial matching would otherwise reject
  `rank = 3` outright; that put deprecated names in every signature.  Passing
  a removed name now raises an error naming its replacement rather than being
  silently absorbed by `...`.
* Speed: multiplicative-update loops hoist loop invariants, and the
  cross-validation / restart wrappers take an opt-in `cores` argument.  Both
  were verified to leave results unchanged.

The last two items are breaking changes, and both are announced in NEWS.md.
They are confined to the `nmf.rrr` family, which is documented as
experimental and was introduced only in the current release cycle.

## Additional checks

* All tests pass (testthat): 35 in the default smoke suite, 3173 with
  NMFKC_FULL_TESTS set.
* All examples run without errors, including `--run-donttest`.
* All vignettes build without errors.
* No reverse dependencies on CRAN.

## On the submission interval

v0.8.8 was published on 2026-07-14, five weeks ago.  The update is offered now
because it is largely bug fixes — in particular the always-zero bootstrap
p-values and the RNG-stream pollution, both of which can silently affect a
user's results — rather than new features.  We are happy to hold the
submission if a longer interval is preferred.

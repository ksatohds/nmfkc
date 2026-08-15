## R CMD check results

0 errors | 0 warnings | 1 note

## Test environments

* Windows 11 (local), R 4.4.1, `--as-cran`
* Ubuntu Linux (local server), R 4.5.3
* win-builder, R-release 4.6.1 and R-devel (2026-08-14 r90407)
* R-hub v2 (Linux, macOS, macOS-arm64, Windows, nosuggests)

## Notes

* "checking for future file timestamps ... NOTE" — appears on environments
  that cannot reach the CRAN time server to verify the current time; not a
  package issue.
* "Possibly misspelled words in DESCRIPTION: Tokuda" — a co-author's surname,
  in the citation for <doi:10.48550/arXiv.2607.27474>. Spelled correctly.

## This is an update

This is a maintenance update from v0.8.8 (currently on CRAN, published
2026-07-14) to v0.9.4.  All changes are listed in NEWS.md.  The update is
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

Both changes are breaking, and both are announced in NEWS.md.  They are
confined to the `nmf.rrr` family, which is documented as experimental and was
introduced only in the current release cycle.
* Speed: multiplicative-update loops hoist loop invariants, and the
  cross-validation / restart wrappers take an opt-in `cores` argument.  Both
  were verified to leave results unchanged.

## Additional checks

* All tests pass (testthat).
* All examples run without errors, including `--run-donttest`.
* All vignettes build without errors.
* No reverse dependencies on CRAN.

## On the submission interval

v0.8.8 was published on 2026-07-14, one month ago.  The update is offered now
because it is largely bug fixes — in particular the always-zero bootstrap
p-values and the RNG-stream pollution, both of which can silently affect a
user's results — rather than new features.  We are happy to hold the
submission if a longer interval is preferred.

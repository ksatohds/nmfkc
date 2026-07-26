# =====================================================================
#  Internal cross-platform parallel lapply
#
#  Used by the multi-fit wrappers (nmf.gmm.select K-sweep,
#  nmf.cox.cv / nmf.cox.cf fold x grid, ...) to evaluate independent,
#  self-seeded fits concurrently.
#
#  Results are IDENTICAL to lapply(): parLapply() and mclapply() both
#  return their results in INPUT order (not completion order), and each
#  task is a deterministic fit that sets its own RNG seed, so the
#  downstream aggregation and model selection are unchanged regardless
#  of `cores` or platform.
#
#    cores <= 1  : sequential lapply()          (no overhead, all platforms)
#    Windows     : PSOCK cluster + parLapply()  (fork is unavailable)
#    Unix / macOS: mclapply() (fork)            (cheaper; no data copy)
#
#  The default number of cores throughout the package is taken from
#  getOption("mc.cores", 1L): the CRAN-safe default is 1 (sequential),
#  and a user can enable parallelism everywhere at once with, e.g.,
#  options(mc.cores = 4).
# =====================================================================
#' Protect the caller's random stream from a self-seeding fit (Internal)
#'
#' The fitters seed themselves (each has an internal default \code{seed}) so
#' that a given call is reproducible.  Without protection that leaks: every
#' call leaves \code{.Random.seed} at the same state, so a user loop drawing
#' random numbers after a fit silently repeats itself and an outer
#' \code{set.seed()} has no effect on anything downstream.
#'
#' Usage inside a fitter, right after \code{seed} is resolved:
#' \preformatted{
#'   .rng <- .nmfkc.rng.save(seed)
#'   on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
#' }
#' When \code{seed} is \code{NULL} nothing is saved and nothing is restored, so
#' that path keeps whatever RNG behaviour it had.  Note that most fitters take
#' \code{seed} through \code{\dots} and resolve it with
#' \code{if (!is.null(extra$seed)) extra$seed else <default>}, which cannot
#' distinguish "\code{seed = NULL} was passed" from "\code{seed} was not
#' passed" -- so for those an explicit \code{NULL} yields the default, not a
#' free-running stream.  Only where \code{seed} is a real formal (e.g.
#' \code{nmf.ffb}, \code{nmfkc.signed.rff}) does \code{NULL} mean "do not
#' seed".
#'
#' @param seed The seed the fit is about to set, or \code{NULL}.
#' @return Opaque state for \code{.nmfkc.rng.restore()}: \code{NULL} when there
#'   is nothing to protect.
#' @keywords internal
#' @noRd
.nmfkc.rng.save <- function(seed) {
  if (base::is.null(seed)) return(NULL)
  if (base::exists(".Random.seed", envir = base::globalenv()))
    base::get(".Random.seed", envir = base::globalenv())
  else
    base::structure(base::list(), class = "nmfkc.rng.absent")
}

#' Restore a random stream saved by .nmfkc.rng.save() (Internal)
#' @param state Value returned by \code{\link{.nmfkc.rng.save}}.
#' @return \code{NULL}, invisibly.
#' @keywords internal
#' @noRd
.nmfkc.rng.restore <- function(state) {
  if (base::is.null(state)) return(base::invisible(NULL))
  if (base::inherits(state, "nmfkc.rng.absent")) {
    ## There was no stream before the call; leave the session as we found it.
    if (base::exists(".Random.seed", envir = base::globalenv()))
      base::rm(".Random.seed", envir = base::globalenv())
  } else {
    base::assign(".Random.seed", state, envir = base::globalenv())
  }
  base::invisible(NULL)
}


.nmfkc.parlapply <- function(X, FUN, cores = 1L, envir = parent.frame(),
                             packages = "nmfkc") {
  cores <- suppressWarnings(as.integer(cores))
  if (length(cores) != 1L || is.na(cores) || cores <= 1L) return(lapply(X, FUN))
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    ## Attach the required packages in each PSOCK worker. A forked worker
    ## (non-Windows) inherits the parent's attached search path, but a fresh
    ## PSOCK worker does not, so e.g. formula terms like survival::Surv() would
    ## otherwise be unresolvable. Callers pass every package whose exported
    ## symbols the worker closures (or their formulas) reference.
    parallel::clusterCall(cl, function(pkgs)
      for (p in pkgs) suppressMessages(library(p, character.only = TRUE)), packages)
    ## Export the caller's locals so the worker closures resolve their free
    ## variables. clusterExport also forces any lazy promises (e.g. the Y/A
    ## function arguments) in `envir`, so the subsequently serialized closure
    ## carries values rather than un-evaluated promises bound to a missing env.
    ## Skip any binding that errors when forced -- an un-evaluated missing-arg
    ## promise (e.g. a `data` argument unused in matrix mode) is not needed by
    ## the task closure (the sequential path never forced it either), and
    ## clusterExport would otherwise abort on get()-ing it.
    vars <- ls(envir, all.names = TRUE)
    vars <- vars[vapply(vars, function(v)
      tryCatch({ force(get(v, envir = envir)); TRUE }, error = function(e) FALSE),
      logical(1))]
    if (length(vars)) parallel::clusterExport(cl, varlist = vars, envir = envir)
    parallel::parLapply(cl, X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = cores)
  }
}

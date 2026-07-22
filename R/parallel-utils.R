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

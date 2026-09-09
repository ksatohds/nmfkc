# ============================================================
#  R/nmf.sem-deprecated.R -- deprecated NMF-SEM aliases
#
#  The model formerly exposed as "NMF-SEM" is now the canonical
#  NMF-FFB (Non-negative Matrix Factorization with Feed-Forward +
#  Feedback; Satoh 2025, arXiv:2512.18250).  The implementation now
#  lives under the nmf.ffb* names; the nmf.sem* functions remain as
#  thin deprecated wrappers so existing code keeps working.  Each
#  emits a .Deprecated() note pointing at its nmf.ffb* replacement.
#
#  Fitted objects carry class c("nmf.ffb", "nmf.sem", "nmf").  The S3
#  methods are defined on nmf.ffb; the nmf.sem methods at the end of
#  this file are one-line aliases, kept so that an object saved by a
#  version that wrote only c("nmf.sem", "nmf") still dispatches.  When
#  nmf.sem goes, this whole file goes with it.
# ============================================================

#' Deprecated NMF-SEM aliases
#'
#' @description
#' These functions are deprecated aliases retained for backward
#' compatibility.  Use the \code{nmf.ffb*} names instead:
#' \code{\link{nmf.ffb}}, \code{\link{nmf.ffb.inference}},
#' \code{\link{nmf.ffb.cv}}, \code{\link{nmf.ffb.split}} and
#' \code{\link{nmf.ffb.DOT}}.
#'
#' @param ... Arguments passed on to the corresponding \code{nmf.ffb*}
#'   function.
#' @return As the corresponding \code{nmf.ffb*} function.
#' @name nmf.sem-deprecated
#' @keywords internal
NULL

#' @rdname nmf.sem-deprecated
#' @export
nmf.sem <- function(...) { .Deprecated("nmf.ffb"); nmf.ffb(...) }

#' @rdname nmf.sem-deprecated
#' @export
nmf.sem.inference <- function(...) { .Deprecated("nmf.ffb.inference"); nmf.ffb.inference(...) }

#' @rdname nmf.sem-deprecated
#' @export
nmf.sem.cv <- function(...) { .Deprecated("nmf.ffb.cv"); nmf.ffb.cv(...) }

#' @rdname nmf.sem-deprecated
#' @export
nmf.sem.split <- function(...) { .Deprecated("nmf.ffb.split"); nmf.ffb.split(...) }

#' @rdname nmf.sem-deprecated
#' @export
nmf.sem.DOT <- function(...) { .Deprecated("nmf.ffb.DOT"); nmf.ffb.DOT(...) }


# ------------------------------------------------------------------
#  Deprecated S3 aliases.  These do not call .Deprecated(): a method is
#  reached by dispatch, not by name, so the note would fire on ordinary
#  use of an old object and there would be nothing the user could do
#  about it.  The functions above are the ones worth warning on.
# ------------------------------------------------------------------

#' @rdname plot.nmfre
#' @export
plot.nmf.sem <- function(x, ...) plot.nmf.ffb(x, ...)

#' @rdname summary.nmf.ffb
#' @export
summary.nmf.sem <- function(object, ...) summary.nmf.ffb(object, ...)

#' @rdname print.summary.nmf.ffb
#' @export
print.summary.nmf.sem <- function(x, ...) print.summary.nmf.ffb(x, ...)

#' @rdname coef.nmf
#' @export
coef.nmf.sem <- function(object, ...) coef.nmf.ffb(object, ...)

#' @rdname fitted.nmf
#' @export
fitted.nmf.sem <- function(object, ...) fitted.nmf.ffb(object, ...)

#' @rdname residuals.nmf
#' @export
residuals.nmf.sem <- function(object, Y, ...) residuals.nmf.ffb(object, Y, ...)

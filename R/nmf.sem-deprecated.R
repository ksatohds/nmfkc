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
#  Fitted objects still carry class c("nmf.ffb", "nmf.sem", "nmf"),
#  so all S3 methods (summary.nmf.sem, plot.nmf.sem, ...) are shared
#  by both names via inheritance.
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

# ============================================================
#  R/nmfae-deprecated.R -- deprecated NMF-AE / NMF-RRR aliases
#
#  The three-layer model formerly exposed as "NMF-AE" is now the
#  canonical NMF-RRR (Non-negative reduced-rank regression).  The
#  implementation lives under the nmf.rrr* / nmf.rrr.signed* names
#  (in nmfae.R / nmfae.signed.R); the nmfae* / nmfae.signed* names
#  remain as thin deprecated wrappers so existing code keeps working.
#  Each emits a .Deprecated() note pointing at its nmf.rrr* name.
#
#  Fitted objects still carry the "nmfae" / "nmfae.signed" S3 classes
#  (e.g. class(fit) = c("nmf.rrr", "nmfae", "nmf")), so every existing
#  S3 method (summary.nmfae, plot.nmfae, predict.nmfae, ...) is shared
#  by both names via inheritance.
# ============================================================

#' Deprecated NMF-AE aliases
#'
#' @description
#' These functions are deprecated aliases retained for backward
#' compatibility.  Use the \code{nmf.rrr*} names instead, e.g.
#' \code{\link{nmf.rrr}}, \code{\link{nmf.rrr.inference}},
#' \code{\link{nmf.rrr.ecv}}, \code{\link{nmf.rrr.cv}},
#' \code{\link{nmf.rrr.rank}}, \code{\link{nmf.rrr.DOT}},
#' \code{\link{nmf.rrr.heatmap}}, \code{\link{nmf.rrr.kernel.beta.cv}},
#' \code{\link{nmf.rrr.rename}}, and the corresponding
#' \code{nmf.rrr.signed*} functions.
#'
#' @param ... Arguments passed on to the corresponding \code{nmf.rrr*}
#'   function.
#' @return As the corresponding \code{nmf.rrr*} function.
#' @name nmfae-deprecated
#' @keywords internal
NULL

## ---- main NMF-AE family ----

#' @rdname nmfae-deprecated
#' @export
nmfae <- function(...) { .Deprecated("nmf.rrr"); nmf.rrr(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.inference <- function(...) { .Deprecated("nmf.rrr.inference"); nmf.rrr.inference(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.ecv <- function(...) { .Deprecated("nmf.rrr.ecv"); nmf.rrr.ecv(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.cv <- function(...) { .Deprecated("nmf.rrr.cv"); nmf.rrr.cv(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.rank <- function(...) { .Deprecated("nmf.rrr.rank"); nmf.rrr.rank(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.DOT <- function(...) { .Deprecated("nmf.rrr.DOT"); nmf.rrr.DOT(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.heatmap <- function(...) { .Deprecated("nmf.rrr.heatmap"); nmf.rrr.heatmap(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.kernel.beta.cv <- function(...) { .Deprecated("nmf.rrr.kernel.beta.cv"); nmf.rrr.kernel.beta.cv(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.rename <- function(...) { .Deprecated("nmf.rrr.rename"); nmf.rrr.rename(...) }

## ---- signed-bottleneck NMF-AE family ----

#' @rdname nmfae-deprecated
#' @export
nmfae.signed <- function(...) { .Deprecated("nmf.rrr.signed"); nmf.rrr.signed(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.signed.inference <- function(...) { .Deprecated("nmf.rrr.signed.inference"); nmf.rrr.signed.inference(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.signed.ecv <- function(...) { .Deprecated("nmf.rrr.signed.ecv"); nmf.rrr.signed.ecv(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.signed.rank <- function(...) { .Deprecated("nmf.rrr.signed.rank"); nmf.rrr.signed.rank(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.signed.heatmap <- function(...) { .Deprecated("nmf.rrr.signed.heatmap"); nmf.rrr.signed.heatmap(...) }

#' @rdname nmfae-deprecated
#' @export
nmfae.signed.rename <- function(...) { .Deprecated("nmf.rrr.signed.rename"); nmf.rrr.signed.rename(...) }

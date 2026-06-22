# ============================================================
#  R/nmf.rrr.R -- NMF-RRR aliases for the nmfae* / nmfae.signed* families
#
#  "NMF-RRR" (a tri-factorized Non-negative Matrix Factorization read as a
#  Reduced-Rank Regression) is the name adopted in Satoh & Tokuda.  It is the
#  same implementation as nmfae* -- only the name differs -- fitting
#  Y1 ~ X1 Theta X2 Y2, i.e. tri-factorizing the non-negative regression
#  coefficient M = X1 Theta X2 to co-cluster the response variables (X1,
#  "Resp") and the covariate variables (X2, "Cov"), with Theta the block
#  correspondence.
#
#  Thin wrappers (the legacy nmfae* / nmfae.signed* remain fully functional;
#  no deprecation):
#      nmf.rrr*         = nmfae*
#      nmf.rrr.signed*  = nmfae.signed*
#  The two fitters prepend "nmf.rrr" / "nmf.rrr.signed" to the result class;
#  existing S3 methods (summary/plot/predict/...) are reused by inheritance
#  (results keep the "nmfae" / "nmfae.signed" classes).
# ============================================================


#' NMF-RRR: tri-factorized non-negative reduced-rank regression
#'
#' @description
#' \code{nmf.rrr} and the \code{nmf.rrr.*} functions are the NMF-RRR names for
#' the \code{\link{nmfae}} and \code{\link{nmfae.signed}} families.  NMF-RRR
#' fits \eqn{Y_1 \approx X_1 \Theta X_2 Y_2}: tri-factorizing the non-negative
#' regression coefficient \eqn{M = X_1 \Theta X_2} co-clusters the response
#' variables (\eqn{X_1}, labelled \dQuote{Resp}) and the covariate variables
#' (\eqn{X_2}, \dQuote{Cov}), with \eqn{\Theta} the block correspondence whose
#' entries are estimated and tested.  These are thin aliases: see the wrapped
#' \code{nmfae*} / \code{nmfae.signed*} function for the arguments and full
#' details.  \code{nmf.rrr()} and \code{nmf.rrr.signed()} additionally prepend
#' \code{"nmf.rrr"} / \code{"nmf.rrr.signed"} to the fitted object's class.
#'
#' @param ... Arguments passed to the wrapped \code{nmfae*} /
#'   \code{nmfae.signed*} function.
#' @return As the wrapped function (the two fitters additionally prepend the
#'   NMF-RRR class to the result).
#' @seealso \code{\link{nmfae}}, \code{\link{nmfae.inference}},
#'   \code{\link{nmfae.ecv}}, \code{\link{nmfae.signed}}
#' @name nmf.rrr
#' @export
nmf.rrr <- function(...) {
  res <- nmfae(...)
  class(res) <- base::unique(c("nmf.rrr", class(res)))
  res
}

#' @rdname nmf.rrr
#' @export
nmf.rrr.inference <- function(...) nmfae.inference(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.ecv <- function(...) nmfae.ecv(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.cv <- function(...) nmfae.cv(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.rank <- function(...) nmfae.rank(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.DOT <- function(...) nmfae.DOT(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.heatmap <- function(...) nmfae.heatmap(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.kernel.beta.cv <- function(...) nmfae.kernel.beta.cv(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.rename <- function(...) nmfae.rename(...)


# ---- signed family ----

#' @rdname nmf.rrr
#' @export
nmf.rrr.signed <- function(...) {
  res <- nmfae.signed(...)
  class(res) <- base::unique(c("nmf.rrr.signed", class(res)))
  res
}

#' @rdname nmf.rrr
#' @export
nmf.rrr.signed.inference <- function(...) nmfae.signed.inference(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.signed.ecv <- function(...) nmfae.signed.ecv(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.signed.rank <- function(...) nmfae.signed.rank(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.signed.heatmap <- function(...) nmfae.signed.heatmap(...)

#' @rdname nmf.rrr
#' @export
nmf.rrr.signed.rename <- function(...) nmfae.signed.rename(...)

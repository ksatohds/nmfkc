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
#' NMF-RRR fits \eqn{Y_1 \approx X_1 \Theta X_2 Y_2}: tri-factorizing the
#' non-negative regression coefficient \eqn{M = X_1 \Theta X_2} co-clusters the
#' response variables (\eqn{X_1}, labelled \dQuote{Resp}) and the covariate
#' variables (\eqn{X_2}, \dQuote{Cov}), with \eqn{\Theta} the block
#' correspondence whose entries are estimated and tested.
#'
#' \code{nmf.rrr()} fits the model; \code{nmf.rrr.signed()} allows a signed
#' bottleneck \eqn{\Theta}. The \code{.inference}, \code{.ecv}, \code{.cv},
#' \code{.rank}, \code{.DOT}, \code{.heatmap}, \code{.kernel.beta.cv} and
#' \code{.rename} helpers provide inference, cross-validation / rank selection,
#' visualization and relabelling. The two fitters prepend
#' \code{"nmf.rrr"} / \code{"nmf.rrr.signed"} to the fitted object's class.
#'
#' @param Y1 Response matrix \eqn{Y_1} (P1 x N), non-negative.
#' @param Y2 Covariate matrix \eqn{Y_2} (P2 x N). Default \code{Y1} (autoencoder).
#' @param rank1 Integer. Rank of the response basis \eqn{X_1} (default 2).
#' @param rank2 Integer. Rank of the covariate basis \eqn{X_2}. Default
#'   (\code{NULL}) sets \code{rank2 = rank1}.
#' @param epsilon Positive convergence tolerance (default \code{1e-4}).
#' @param maxit Maximum number of multiplicative-update iterations (default 5000).
#' @param verbose Logical; print progress (default \code{FALSE}).
#' @param ... Further arguments (e.g. \code{X.init}, \code{C.L1},
#'   \code{Y1.weights}, \code{method}, \code{nstart}, and the helper-specific
#'   arguments). For the helpers (\code{.ecv}, \code{.cv}, \code{.rank}, ...),
#'   pass their arguments here (see the examples). The legacy rank aliases
#'   \code{rank.encoder} / \code{Q} / \code{R} are also accepted via \code{...}
#'   (use \code{rank1} / \code{rank2} in new code).
#' @return A fitted NMF-RRR object (the two fitters additionally prepend the
#'   NMF-RRR class); the helpers return as their respective function.
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.DOT}}
#' @examples
#' \donttest{
#' Y1 <- t(iris[, 1:2])
#' Y2 <- t(iris[, 3:4])
#' fit <- nmf.rrr(Y1, Y2, rank1 = 2, rank2 = 2, maxit = 500)
#' summary(fit)
#' # rank selection and inference
#' # ecv <- nmf.rrr.ecv(Y1, Y2, rank1 = 1:3, rank2 = 1:3)
#' # inf <- nmf.rrr.inference(fit, Y1, Y2)
#' }
#' @name nmf.rrr
#' @export
nmf.rrr <- function(Y1, Y2 = Y1, rank1 = 2, rank2 = NULL,
                    epsilon = 1e-4, maxit = 5000, verbose = FALSE, ...) {
  res <- nmfae(Y1, Y2, rank1 = rank1, rank2 = rank2,
               epsilon = epsilon, maxit = maxit, verbose = verbose, ...)
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
nmf.rrr.signed <- function(Y1, Y2 = Y1, rank1 = 2, rank2 = NULL,
                           epsilon = 1e-4, maxit = 5000, verbose = FALSE, ...) {
  res <- nmfae.signed(Y1, Y2, rank1 = rank1, rank2 = rank2,
                      epsilon = epsilon, maxit = maxit, verbose = verbose, ...)
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

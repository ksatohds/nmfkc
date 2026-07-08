# ============================================================
#  R/nmf.ffb.R — NMF-FFB aliases for the nmf.sem* family
#
#  "NMF-FFB" (Non-negative Matrix Factorization with Feed-Forward +
#  Feedback) is the canonical name adopted in Satoh (2025), arXiv:
#  2512.18250.  It is identical to "NMF-SEM" -- only the name differs.
#
#  This file provides thin wrappers
#      nmf.ffb()           = nmf.sem()
#      nmf.ffb.inference() = nmf.sem.inference()
#      nmf.ffb.cv()        = nmf.sem.cv()
#      nmf.ffb.DOT()       = nmf.sem.DOT()
#      nmf.ffb.split()     = nmf.sem.split()
#  so that the package's user-facing API matches the paper.  The
#  legacy nmf.sem* functions remain fully functional (no deprecation
#  is planned at this time).
#
#  Class hierarchy: result objects carry both new and legacy classes,
#      class(res) = c("nmf.ffb", "nmf.sem")
#  so existing S3 methods (summary.nmf.sem, plot.nmf.sem, ...) are
#  automatically reused via S3 inheritance.
# ============================================================


#' NMF-FFB: NMF with Feed-Forward and Feedback structure
#'
#' @description
#' \code{nmf.ffb} and the \code{nmf.ffb.*} functions fit NMF-FFB
#' (Non-negative Matrix Factorization with Feed-Forward + Feedback), the
#' canonical name adopted in Satoh (2025), arXiv:2512.18250. They are thin
#' wrappers over the internal \code{nmf.sem*} implementation (an earlier name
#' for the same model); result objects carry \code{class = c("nmf.ffb",
#' "nmf.sem")} so all S3 methods are reused by inheritance.
#' \itemize{
#'   \item \code{nmf.ffb()} fits the model.
#'   \item \code{nmf.ffb.inference()} adds bootstrap inference on the feedback /
#'     feed-forward coefficients.
#'   \item \code{nmf.ffb.cv()} performs cross-validation for rank selection.
#'   \item \code{nmf.ffb.DOT()} renders the fitted structure as a Graphviz graph.
#'   \item \code{nmf.ffb.split()} heuristically splits variables into
#'     exogenous / endogenous groups.
#' }
#'
#' @param Y1,Y2 Non-negative data matrices (endogenous and exogenous blocks).
#' @param rank Integer basis rank (\code{NULL} lets the model choose).
#' @param X.init Basis initialization method (default \code{"nndsvd"}); accepts
#'   the same menu as \code{\link{nmfkc}} (e.g. \code{"kmeans"}, \code{"kmeans++"}).
#' @param X.L2.ortho,C1.L1,C2.L1 Non-negative penalty parameters (basis
#'   orthogonality, and L1 on the two coefficient blocks).
#' @param epsilon,maxit Convergence tolerance and maximum iterations.
#' @param seed Integer RNG seed.
#' @param object A fitted \code{nmf.ffb} object (for \code{nmf.ffb.inference}).
#' @param B Number of bootstrap replicates.
#' @param threshold Minimum coefficient magnitude to display / test.
#' @param ci.level Confidence level for bootstrap intervals.
#' @param result A fitted object (for \code{nmf.ffb.DOT}).
#' @param x A fitted object (for \code{nmf.ffb.split}).
#' @param ... Further arguments passed to the underlying implementation.
#' @return As the corresponding fitter / helper (fits carry the \code{"nmf.ffb"}
#'   class in addition to \code{"nmf.sem"}).
#' @references Satoh, K. (2025). arXiv:2512.18250.
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.DOT}}
#' @name nmf.ffb
#' @export
nmf.ffb <- function(Y1, Y2, rank = NULL, X.init = "nndsvd",
                    X.L2.ortho = 100.0, C1.L1 = 1.0, C2.L1 = 0.1,
                    epsilon = 1e-6, maxit = 5000, seed = 123, ...) {
  nmf.sem(Y1, Y2, rank = rank, X.init = X.init,
          X.L2.ortho = X.L2.ortho, C1.L1 = C1.L1, C2.L1 = C2.L1,
          epsilon = epsilon, maxit = maxit, seed = seed, ...)
}


#' @rdname nmf.ffb
#' @export
nmf.ffb.inference <- function(object, Y1, Y2,
                              B = 1000L, threshold = 0.01,
                              ci.level = 0.95,
                              C1.L1 = 1.0, C2.L1 = 0.1,
                              seed = 123L, ...) {
  nmf.sem.inference(object, Y1, Y2,
                    B = B, threshold = threshold,
                    ci.level = ci.level,
                    C1.L1 = C1.L1, C2.L1 = C2.L1,
                    seed = seed, ...)
}


#' @rdname nmf.ffb
#' @export
nmf.ffb.cv <- function(...) nmf.sem.cv(...)


#' @rdname nmf.ffb
#' @export
nmf.ffb.DOT <- function(result, ...) nmf.sem.DOT(result, ...)


#' @rdname nmf.ffb
#' @export
nmf.ffb.split <- function(x, ...) nmf.sem.split(x, ...)

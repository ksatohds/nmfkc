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


#' @rdname nmf.sem
#' @export
nmf.ffb <- function(Y1, Y2, rank = NULL, X.init = "nndsvd",
                    X.L2.ortho = 100.0, C1.L1 = 1.0, C2.L1 = 0.1,
                    epsilon = 1e-6, maxit = 5000, seed = 123, ...) {
  nmf.sem(Y1, Y2, rank = rank, X.init = X.init,
          X.L2.ortho = X.L2.ortho, C1.L1 = C1.L1, C2.L1 = C2.L1,
          epsilon = epsilon, maxit = maxit, seed = seed, ...)
}


#' @rdname nmf.sem.inference
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


#' @rdname nmf.sem.cv
#' @export
nmf.ffb.cv <- function(...) nmf.sem.cv(...)


#' @rdname nmf.sem.DOT
#' @export
nmf.ffb.DOT <- function(result, ...) nmf.sem.DOT(result, ...)


#' @rdname nmf.sem.split
#' @export
nmf.ffb.split <- function(x, ...) nmf.sem.split(x, ...)

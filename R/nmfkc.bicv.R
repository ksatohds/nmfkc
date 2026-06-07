## ============================================================
## nmfkc.bicv(): Bi-cross-validation for NMF rank selection
## (Owen & Perry 2009, AoAS 3(2):564-594).  PROTOTYPE / experimental.
##
## Idea: split rows into folds and columns into folds.  For a held-out
## row-block I and column-block J, the matrix is partitioned into
##   [ A B ]   A = Y[I, J]   (held-out block, predicted)
##   [ C D ]   B = Y[I, Jc]  C = Y[Ic, J]  D = Y[Ic, Jc]  (retained)
## Fit NMF only on D (D ~ L_D R_D), then "fold in" the held-out rows and
## columns by NON-NEGATIVE REGRESSION onto the fixed D-factors:
##   L_I = argmin_{>=0} ||B - L_I R_D||   (held rows' loadings)
##   R_J = argmin_{>=0} ||C - L_D R_J||   (held cols' scores)
## and predict A_hat = L_I R_J.  The held-out rows AND columns never
## enter the NMF fit, so there is no information leakage.
## ============================================================

#' @title Non-negative least squares by multiplicative updates (Internal)
#' @description Solve \eqn{\min_{L\ge 0}\|M - L F\|^2_F} with \eqn{F}
#'   fixed (\code{side = "left"}), or \eqn{\min_{R\ge 0}\|M - F R\|^2_F}
#'   with \eqn{F} fixed (\code{side = "right"}), via a few multiplicative
#'   updates.  Used by \code{\link{nmfkc.bicv}} to fold held-out rows /
#'   columns onto fixed factors.
#' @param M Target matrix.
#' @param F Fixed factor.
#' @param side \code{"left"} solves for \eqn{L} (\code{nrow(M) x k});
#'   \code{"right"} solves for \eqn{R} (\code{k x ncol(M)}).
#' @param maxit Number of multiplicative updates.
#' @return The non-negative factor.
#' @keywords internal
#' @noRd
.nnls.mu <- function(M, F, side = c("left", "right"), maxit = 100) {
  side <- base::match.arg(side)
  eps <- 1e-10
  if (side == "left") {
    k <- base::nrow(F)
    L <- base::matrix(base::mean(M) / (k * base::mean(F) + eps) + eps,
                      base::nrow(M), k)
    FFt <- F %*% base::t(F); MFt <- M %*% base::t(F)
    for (it in base::seq_len(maxit)) L <- L * MFt / (L %*% FFt + eps)
    L
  } else {
    k <- base::ncol(F)
    R <- base::matrix(base::mean(M) / (k * base::mean(F) + eps) + eps,
                      k, base::ncol(M))
    FtF <- base::t(F) %*% F; FtM <- base::t(F) %*% M
    for (it in base::seq_len(maxit)) R <- R * FtM / (FtF %*% R + eps)
    R
  }
}


#' @title Bi-cross-validation for NMF rank selection
#' @description
#' Owen & Perry's (2009) bi-cross-validation (BCV) for choosing the NMF
#' rank.  A lightweight CV engine in the spirit of \code{\link{nmfkc.ecv}}:
#' it returns the held-out error per rank and nothing more (no plot, no
#' table) -- pass the result to \code{which.min(sigma)}, or build your own
#' diagnostics.
#'
#' Unlike the element-wise CV of \code{\link{nmfkc.ecv}} (which holds out
#' scattered \emph{entries} and refits with weights), BCV holds out a whole
#' \strong{row-block and column-block} simultaneously: the model is fit only
#' on the retained block \eqn{D}, and the held-out block \eqn{A} is predicted
#' by folding the held-out rows/columns onto the fixed \eqn{D}-factors via
#' non-negative regression (\eqn{\hat A = L_I R_J}).  Because the held-out
#' rows and columns never enter the fit, there is no information leakage.
#' Covariates are ignored (plain NMF).  The recommended setting is to leave
#' out roughly half the rows and half the columns (\code{nfolds = 2}).
#'
#' @param Y Observation matrix (\eqn{P \times N}), non-negative.
#' @param rank Integer vector of ranks to evaluate.
#' @param nfolds Number of row folds and column folds (the grid is
#'   \code{nfolds x nfolds}).  \code{2} = leave out half rows / half
#'   columns (Owen & Perry's recommendation).
#' @param seed Optional integer seed for the fold assignment.
#' @param nnls.maxit Multiplicative-update iterations for the fold-in
#'   non-negative regressions.
#' @details Each fold keeps about \eqn{(1 - 1/\text{nfolds})} of the rows
#'   and columns, so the retained block \eqn{D} must have more than
#'   \code{rank} rows \strong{and} columns; ranks that violate this for a
#'   given fold are skipped (the held-out error is \code{NA}).  With
#'   \code{nfolds = 2} this needs roughly \eqn{P/2 > \text{rank}} and
#'   \eqn{N/2 > \text{rank}}.
#' @param ... Passed to \code{\link{nmfkc}} for the per-block fits
#'   (e.g.\ \code{maxit}, \code{method}).
#' @return A list (cf.\ \code{\link{nmfkc.ecv}}) with:
#'   \item{objfunc}{Held-out mean squared error for each rank.}
#'   \item{sigma}{Its square root (RMSE) for each rank.}
#'   \item{rank}{The evaluated rank vector.}
#'   \item{nfolds}{The number of folds used.}
#' @references A. B. Owen and P. O. Perry (2009).  Bi-cross-validation of
#'   the SVD and the nonnegative matrix factorization.
#'   \emph{Ann. Appl. Stat.} 3(2):564--594. \doi{10.1214/08-AOAS227}.
#' @seealso \code{\link{nmfkc.ecv}}, \code{\link{nmfkc.rank}}
#' @export
#' @examples
#' \donttest{
#' ## rank-3 non-negative data; bi-CV needs enough kept rows/cols per
#' ## fold (> rank), so use a matrix with ample dimensions.
#' set.seed(1)
#' X <- matrix(abs(rnorm(30 * 3)), 30, 3)
#' B <- matrix(abs(rnorm(3 * 40)), 3, 40)
#' bv <- nmfkc.bicv(X %*% B, rank = 1:6, nfolds = 2)
#' bv$sigma                  # held-out RMSE per rank
#' bv$rank[which.min(bv$sigma)]
#' }
nmfkc.bicv <- function(Y, rank = 1:3, nfolds = 2, seed = 123,
                       nnls.maxit = 100, ...) {
  Y <- base::as.matrix(Y)
  P <- base::nrow(Y); N <- base::ncol(Y)
  if (!base::is.null(seed)) base::set.seed(seed)
  row.fold <- base::sample(base::rep_len(1:nfolds, P))
  col.fold <- base::sample(base::rep_len(1:nfolds, N))

  ## nmfkc args for the (plain-NMF) block fits: drop covariates, skip the
  ## O(N^2) clustering criteria, stay quiet about dimensions.
  ea <- base::list(...); ea$A <- NULL; ea$Q <- NULL
  ea$detail <- "fast"; ea$print.dims <- FALSE

  objfunc <- stats::setNames(base::numeric(base::length(rank)),
                             base::paste0("rank=", rank))
  base::message(base::sprintf(
    "bi-CV: ranks %s, %dx%d fold grid (Owen-Perry 2009)...",
    base::paste(rank, collapse = ","), nfolds, nfolds))

  for (ki in base::seq_along(rank)) {
    k <- rank[ki]
    sse <- 0; cnt <- 0L
    for (fi in 1:nfolds) for (fj in 1:nfolds) {
      I  <- base::which(row.fold == fi); J  <- base::which(col.fold == fj)
      Ic <- base::which(row.fold != fi); Jc <- base::which(col.fold != fj)
      if (base::length(Ic) <= k || base::length(Jc) <= k) next
      D  <- Y[Ic, Jc, drop = FALSE]
      Bm <- Y[I,  Jc, drop = FALSE]   # held rows x kept cols
      Cm <- Y[Ic, J,  drop = FALSE]   # kept rows x held cols
      Am <- Y[I,  J,  drop = FALSE]   # held-out block
      fit <- base::suppressMessages(
        base::do.call("nmfkc", c(base::list(Y = D, rank = k), ea)))
      L_D <- fit$X; R_D <- fit$B
      L_I <- .nnls.mu(Bm, R_D, side = "left",  maxit = nnls.maxit)
      R_J <- .nnls.mu(Cm, L_D, side = "right", maxit = nnls.maxit)
      A_hat <- L_I %*% R_J
      sse <- sse + base::sum((Am - A_hat)^2)
      cnt <- cnt + base::length(Am)
    }
    objfunc[ki] <- if (cnt > 0) sse / cnt else NA_real_
  }

  base::list(objfunc = objfunc, sigma = base::sqrt(objfunc),
             rank = rank, nfolds = nfolds)
}

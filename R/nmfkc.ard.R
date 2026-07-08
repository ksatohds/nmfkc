## ============================================================
## nmfkc.ard(): Automatic Relevance Determination for NMF rank
## (Tan & Fevotte 2013, beta = 2 / Euclidean).  PROTOTYPE.
##
## Fit NMF at an over-complete rank K and let a per-component relevance
## weight lambda_k (inverse-gamma prior) automatically prune unneeded
## components: the multiplicative updates gain a penalty term, and
## components whose energy collapses to ~0 are dropped.  The number of
## surviving components is the estimated rank -- a single fit, no scan.
## ============================================================

## Internal: one ARD-NMF fit from a single random initialization.
#' @keywords internal
#' @noRd
.ard.fit <- function(Y, K, prior, a, b, maxit, epsilon, tol, seed) {
  Fr <- base::nrow(Y); Nc <- base::ncol(Y); eps <- 1e-10; m <- base::mean(Y)
  if (!base::is.null(seed)) base::set.seed(seed)
  W <- base::matrix(stats::runif(Fr * K), Fr, K) * base::sqrt(m / K)
  H <- base::matrix(stats::runif(K * Nc), K, Nc) * base::sqrt(m / K)
  clam <- if (prior == "L2") (Fr + Nc) / 2 + a + 1 else (Fr + Nc) + a + 1

  objfunc <- base::numeric(maxit)
  for (it in base::seq_len(maxit)) {
    lam <- if (prior == "L2")
      (b + 0.5 * (base::colSums(W^2) + base::rowSums(H^2))) / clam
    else
      (b + base::colSums(W) + base::rowSums(H)) / clam
    num <- Y %*% base::t(H); den <- W %*% (H %*% base::t(H))
    pen <- if (prior == "L2") base::sweep(W, 2, 1 / lam, "*")
           else base::matrix(1 / lam, Fr, K, byrow = TRUE)
    W <- W * num / (den + pen + eps)
    num <- base::t(W) %*% Y; den <- (base::t(W) %*% W) %*% H
    pen <- if (prior == "L2") base::sweep(H, 1, 1 / lam, "*")
           else base::matrix(1 / lam, K, Nc)
    H <- H * num / (den + pen + eps)
    objfunc[it] <- 0.5 * base::sum((Y - W %*% H)^2)
    if (it > 1 && base::abs(objfunc[it - 1] - objfunc[it]) <
        epsilon * objfunc[it - 1]) { objfunc <- objfunc[1:it]; break }
  }

  energy <- base::sqrt(base::colSums(W^2) * base::rowSums(H^2))
  ord <- base::order(energy, decreasing = TRUE)
  energy <- energy[ord]; W <- W[, ord, drop = FALSE]; H <- H[ord, , drop = FALSE]
  lam <- lam[ord]
  relevance <- if (base::max(energy) > 0) energy / base::max(energy) else energy
  base::list(rank = base::sum(relevance > tol), relevance = relevance,
             lambda = lam, W = W, H = H, objfunc = objfunc)
}

#' @title Automatic relevance determination for NMF rank (experimental)
#' @description
#' \strong{Prototype} of Tan & Fevotte's (2013) ARD-NMF (Euclidean /
#' \eqn{\beta = 2}).  Unlike the cross-validation (\code{\link{nmfkc.ecv}},
#' \code{\link{nmfkc.bicv}}) and stability (\code{\link{nmfkc.consensus}})
#' engines that \emph{scan} a range of ranks, ARD fits NMF \strong{once} at
#' an over-complete rank and prunes automatically: each component \eqn{k}
#' carries a relevance weight \eqn{\lambda_k} with an inverse-gamma prior;
#' the multiplicative updates gain a penalty \eqn{+ w_{fk}/\lambda_k} (L2 /
#' half-normal prior) or \eqn{+ 1/\lambda_k} (L1 / exponential), which
#' drives unsupported components to zero.  The number of surviving
#' components is the estimated rank.  Covariates are ignored (plain NMF).
#'
#' This is a model-based point estimate: the result depends on the prior, the
#' starting rank and the random initialization, and can vary run to run.  Use
#' it as a complement to the CV / consensus engines, not a sole criterion.
#'
#' @details
#' \strong{Relation to Tan & Fevotte (2013).}  The update equations reproduce
#' the paper's \eqn{\ell_2} / \eqn{\ell_1} ARD-NMF for the \strong{Euclidean
#' (\eqn{\beta = 2})} case exactly: the multiplicative penalties, the
#' closed-form \eqn{\lambda_k} update \eqn{(f(w_k) + f(h_k) + b)/c} and the
#' constant \eqn{c} (their Eq. 33).  Two deliberate simplifications keep it
#' practical: only \eqn{\beta = 2} is implemented (the paper covers the general
#' \eqn{\beta}-divergence), and the default \code{b} is an \emph{empirical}
#' per-component energy scale \eqn{(F + N)/K \cdot \bar{Y}} that reliably avoids
#' the winner-take-all collapse, rather than the paper's method-of-moments value
#' (their Eq. 38).
#'
#' @param Y Observation matrix (\eqn{F \times N}), non-negative.
#' @param rank Over-complete starting rank \eqn{K} (must exceed the true
#'   rank).  \code{NULL} (default) uses \code{min(F, N, 20)}.
#' @param nrun Number of random-initialization restarts (default \code{10}).
#'   ARD is a sensitive point estimate, so several restarts are advisable; the
#'   reported rank is the \strong{mode} of the per-run estimates
#'   (\code{rank.runs}), with a representative modal fit kept for
#'   \code{plot}/\code{W}/\code{H}.
#' @param plot Logical; draw the relevance bar plot.
#' @param ... Advanced options, rarely needed (defaults in parentheses):
#'   \code{prior} (\code{"L2"}: half-normal / squared-energy group
#'   shrinkage, \code{"L1"}: exponential / sparser); \code{seed} (\code{123},
#'   random-initialization seed); \code{a} (\code{1}) and
#'   \code{b} (\code{(F + N)/K * mean(Y)}), the inverse-gamma prior (a smaller
#'   \code{b} over-prunes, a larger one prunes nothing); \code{maxit}
#'   (\code{3000}) and \code{epsilon} (\code{1e-6}) for optimisation control;
#'   \code{tol} (\code{1e-3}), the relevance threshold below which a component is
#'   counted as pruned.
#' @return An object of class \code{"nmfkc.ard"}: a list with
#'   \code{rank} (estimated = mode over restarts), \code{rank.runs} (the
#'   per-run estimates), \code{relevance} (representative run, descending),
#'   \code{lambda}, \code{W}, \code{H} (ordered by relevance),
#'   \code{rank.init}, \code{prior}, \code{nrun}, \code{objfunc} (final
#'   objective value of the representative run) and \code{objfunc.iter} (its
#'   per-iteration trajectory).
#' @references V. Y. F. Tan and C. Fevotte (2013).  Automatic relevance
#'   determination in nonnegative matrix factorization with the
#'   beta-divergence.  \emph{IEEE TPAMI} 35(7):1592--1605.
#'   \doi{10.1109/TPAMI.2012.240}.
#' @seealso \code{\link{nmfkc.ecv}}, \code{\link{nmfkc.bicv}},
#'   \code{\link{nmfkc.consensus}}, \code{\link{nmfkc.rank}}
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(abs(rnorm(40 * 3)), 40, 3)
#' B <- matrix(abs(rnorm(3 * 60)), 3, 60)
#' ar <- nmfkc.ard(X %*% B, rank = 10)   # over-complete start
#' ar$rank                                # ~ 3 surviving components
#' plot(ar)
#' }
nmfkc.ard <- function(Y, rank = NULL, nrun = 10, plot = FALSE, ...) {
  Y <- base::as.matrix(Y); Fr <- base::nrow(Y); Nc <- base::ncol(Y)
  K <- if (base::is.null(rank)) base::min(Fr, Nc, 20L) else rank

  ## Advanced options live in `...` with safe defaults (rarely changed).
  ## Default prior scale ties b to the initial per-component energy scale
  ## (F + N)/K * mean(Y); a fixed small fraction of mean(Y) over-prunes
  ## (winner-take-all collapse) when (F + N)/K is large.
  dots    <- base::list(...)
  prior   <- base::match.arg(if (base::is.null(dots$prior)) "L2" else dots$prior,
                             c("L2", "L1"))
  seed    <- if (base::is.null(dots$seed))    123     else dots$seed
  a       <- if (base::is.null(dots$a))       1       else dots$a
  b       <- if (base::is.null(dots$b))       (Fr + Nc) / K * base::mean(Y) else dots$b
  maxit   <- if (base::is.null(dots$maxit))   3000    else dots$maxit
  epsilon <- if (base::is.null(dots$epsilon)) 1e-6    else dots$epsilon
  tol     <- if (base::is.null(dots$tol))     1e-3    else dots$tol

  ## nrun random-initialization restarts (ARD is a sensitive point
  ## estimate); aggregate the integer rank by its MODE.
  fits <- base::lapply(base::seq_len(nrun), function(r)
    .ard.fit(Y, K, prior, a, b, maxit, epsilon, tol,
             seed = if (base::is.null(seed)) NULL else seed + r))
  rank.runs <- base::vapply(fits, function(f) f$rank, base::integer(1))
  tab <- base::sort(base::table(rank.runs), decreasing = TRUE)
  rank.mode <- base::as.integer(base::names(tab)[1])
  ## representative fit = a modal-rank run with the smallest residual
  modal <- base::which(rank.runs == rank.mode)
  rep.i <- modal[base::which.min(base::vapply(modal,
            function(i) fits[[i]]$objfunc[base::length(fits[[i]]$objfunc)],
            base::numeric(1)))]
  rep <- fits[[rep.i]]

  out <- base::structure(
    base::list(rank = rank.mode, rank.runs = rank.runs,
               relevance = rep$relevance, lambda = rep$lambda,
               W = rep$W, H = rep$H, rank.init = K, prior = prior,
               nrun = nrun,
               ## house style: objfunc = final scalar, objfunc.iter = trajectory
               objfunc = rep$objfunc[base::length(rep$objfunc)],
               objfunc.iter = rep$objfunc, tol = tol),
    class = "nmfkc.ard")
  if (plot) graphics::plot(out)
  out
}


#' @title Print method for nmfkc.ard objects
#' @param x An object of class \code{"nmfkc.ard"}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.nmfkc.ard <- function(x, ...) {
  base::cat(base::sprintf(
    "ARD-NMF rank determination (Tan-Fevotte 2013, %s prior, %d run%s)\n",
    x$prior, x$nrun, if (x$nrun > 1) "s" else ""))
  base::cat(base::sprintf("  starting rank: %d  ->  estimated rank (mode): %d\n",
                          x$rank.init, x$rank))
  if (!base::is.null(x$rank.runs) && x$nrun > 1) {
    tab <- base::table(x$rank.runs)
    base::cat("  rank over runs:",
              base::paste(base::sprintf("%s:%d", base::names(tab), tab),
                          collapse = "  "),
              base::sprintf("(median %g)\n", stats::median(x$rank.runs)))
  }
  base::cat("  relevance (representative run):",
            base::paste(base::round(x$relevance, 3), collapse = " "), "\n")
  base::invisible(x)
}


#' @title Plot method for nmfkc.ard objects
#' @description Bar plot of the per-component relevance (descending), with a
#'   line at the pruning threshold; bars above it (the estimated rank) are
#'   highlighted.
#' @param x An object of class \code{"nmfkc.ard"}.
#' @param main Plot title.
#' @param ... Passed to \code{\link[graphics]{barplot}}.
#' @return \code{x}, invisibly.
#' @export
plot.nmfkc.ard <- function(x, main = NULL, ...) {
  if (base::is.null(main))
    main <- base::sprintf("ARD-NMF: estimated rank = %d (of %d)",
                          x$rank, x$rank.init)
  keep <- x$relevance > x$tol
  cols <- base::ifelse(keep, "steelblue", "grey80")
  graphics::barplot(x$relevance, col = cols, border = NA,
                    names.arg = base::seq_along(x$relevance),
                    xlab = "component (by relevance)", ylab = "relevance",
                    main = main, ...)
  graphics::abline(h = x$tol, col = "red", lty = 2)
  base::invisible(x)
}

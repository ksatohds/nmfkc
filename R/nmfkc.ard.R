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
#' This is a model-based point estimate: the result depends on the prior
#' (\code{a}, \code{b}), the starting rank and the random initialization,
#' and can vary run to run.  Use it as a complement to the CV / consensus
#' engines, not a sole criterion.
#'
#' @param Y Observation matrix (\eqn{F \times N}), non-negative.
#' @param rank Over-complete starting rank \eqn{K} (must exceed the true
#'   rank).  \code{NULL} (default) uses \code{min(F, N, 20)}.
#' @param prior \code{"L2"} (half-normal; group-shrinks by squared energy)
#'   or \code{"L1"} (exponential).
#' @param a,b Inverse-gamma hyperparameters for \eqn{\lambda_k}.  \code{b}
#'   defaults to a small data-scaled value \code{0.001 * mean(Y)}.
#' @param maxit,epsilon Maximum iterations and relative-objective tolerance.
#' @param tol Relevance threshold (relative to the largest component) below
#'   which a component is counted as pruned.
#' @param seed Random-initialization seed.
#' @param plot Logical; draw the relevance bar plot.
#' @return An object of class \code{"nmfkc.ard"}: a list with
#'   \code{rank} (estimated), \code{relevance} (per-component, descending),
#'   \code{lambda}, \code{W}, \code{H} (ordered by relevance),
#'   \code{rank.init}, \code{prior} and \code{objfunc}.
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
nmfkc.ard <- function(Y, rank = NULL, prior = c("L2", "L1"),
                      a = 1, b = NULL, maxit = 3000, epsilon = 1e-6,
                      tol = 1e-3, seed = 123, plot = FALSE) {
  prior <- base::match.arg(prior)
  Y <- base::as.matrix(Y); Fr <- base::nrow(Y); Nc <- base::ncol(Y)
  eps <- 1e-10
  K <- if (base::is.null(rank)) base::min(Fr, Nc, 20L) else rank
  m <- base::mean(Y)
  if (base::is.null(b)) b <- 0.001 * m

  if (!base::is.null(seed)) base::set.seed(seed)
  W <- base::matrix(stats::runif(Fr * K), Fr, K) * base::sqrt(m / K)
  H <- base::matrix(stats::runif(K * Nc), K, Nc) * base::sqrt(m / K)
  clam <- if (prior == "L2") (Fr + Nc) / 2 + a + 1 else (Fr + Nc) + a + 1

  objfunc <- base::numeric(maxit)
  for (it in base::seq_len(maxit)) {
    ## relevance weights (closed form given W, H)
    lam <- if (prior == "L2")
      (b + 0.5 * (base::colSums(W^2) + base::rowSums(H^2))) / clam
    else
      (b + base::colSums(W) + base::rowSums(H)) / clam
    ## W update (penalised MU)
    num <- Y %*% base::t(H); den <- W %*% (H %*% base::t(H))
    pen <- if (prior == "L2") base::sweep(W, 2, 1 / lam, "*")
           else base::matrix(1 / lam, Fr, K, byrow = TRUE)
    W <- W * num / (den + pen + eps)
    ## H update (penalised MU)
    num <- base::t(W) %*% Y; den <- (base::t(W) %*% W) %*% H
    pen <- if (prior == "L2") base::sweep(H, 1, 1 / lam, "*")
           else base::matrix(1 / lam, K, Nc)
    H <- H * num / (den + pen + eps)

    objfunc[it] <- 0.5 * base::sum((Y - W %*% H)^2)
    if (it > 1 && base::abs(objfunc[it - 1] - objfunc[it]) <
        epsilon * objfunc[it - 1]) { objfunc <- objfunc[1:it]; break }
  }

  ## component energy -> relevance; surviving components = estimated rank
  energy <- base::sqrt(base::colSums(W^2) * base::rowSums(H^2))
  ord <- base::order(energy, decreasing = TRUE)
  energy <- energy[ord]; W <- W[, ord, drop = FALSE]; H <- H[ord, , drop = FALSE]
  lam <- lam[ord]
  relevance <- if (base::max(energy) > 0) energy / base::max(energy)
               else energy
  rank.ard <- base::sum(relevance > tol)

  out <- base::structure(
    base::list(rank = rank.ard, relevance = relevance, lambda = lam,
               W = W, H = H, rank.init = K, prior = prior,
               objfunc = objfunc, tol = tol),
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
    "ARD-NMF rank determination (Tan-Fevotte 2013, %s prior)\n", x$prior))
  base::cat(base::sprintf("  starting rank: %d  ->  estimated rank: %d\n",
                          x$rank.init, x$rank))
  base::cat("  relevance (descending):",
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

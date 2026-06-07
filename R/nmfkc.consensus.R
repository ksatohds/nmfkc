## ============================================================
## nmfkc.consensus(): consensus-clustering rank selection for NMF
## (Brunet et al. 2004).  A lightweight engine in the spirit of
## nmfkc.ecv / nmfkc.bicv: returns a stability score per rank.
## ============================================================

#' @title Consensus-clustering rank selection for NMF (Brunet 2004)
#' @description
#' The bioinformatics-standard stability approach to choosing the NMF rank.
#' A lightweight engine like \code{\link{nmfkc.ecv}} / \code{\link{nmfkc.bicv}}:
#' it returns one stability score per rank and nothing more.
#'
#' For each rank, NMF is run \code{nrun} times from different random
#' initializations (\code{X.init = "runif"}).  Each run gives a hard
#' clustering of the samples (the column \eqn{\arg\max} of the coefficient
#' matrix).  Averaging the \eqn{N \times N} connectivity matrices (1 if two
#' samples share a cluster) over the runs yields the \strong{consensus
#' matrix}; its crispness measures how reproducible the clustering is.  Two
#' summaries are reported per rank:
#' \itemize{
#'   \item \code{cophenetic}: the cophenetic correlation coefficient (CPCC)
#'     of the consensus matrix (Brunet et al. 2004).  Close to 1 = stable.
#'   \item \code{dispersion}: the Kim & Park (2007) dispersion
#'     \eqn{\frac{1}{N^2}\sum_{ij} 4 (C_{ij} - 1/2)^2 \in [0, 1]}; 1 when
#'     every consensus entry is exactly 0 or 1 (perfectly crisp).
#' }
#' Unlike the cross-validation engines (where the rank \emph{minimizes} the
#' error), here a good rank \strong{maximizes} stability, or is the largest
#' rank before it drops.
#'
#' @param Y Observation matrix (\eqn{P \times N}), non-negative.
#' @param A Optional covariate matrix passed to \code{\link{nmfkc}}.
#' @param rank Integer vector of ranks to evaluate (\eqn{\ge 2}).
#' @param nrun Number of random-initialization runs per rank (default 30).
#' @param seed Base integer seed; run \eqn{r} of rank index \eqn{i} uses
#'   \code{seed + 1000 * i + r}.
#' @param keep.consensus Logical; if \code{TRUE} also return the list of
#'   consensus matrices (one \eqn{N \times N} matrix per rank).
#' @param ... Passed to \code{\link{nmfkc}} (e.g.\ \code{maxit},
#'   \code{method}).  \code{X.init} is forced to \code{"runif"}.
#' @return A list with:
#'   \item{cophenetic}{Cophenetic correlation coefficient for each rank.}
#'   \item{dispersion}{Dispersion coefficient (\eqn{[0, 1]}) for each rank.}
#'   \item{rank}{The evaluated rank vector.}
#'   \item{nrun}{Number of runs per rank.}
#'   \item{consensus}{List of consensus matrices, or \code{NULL}.}
#' @references
#' Brunet, J.-P., Tamayo, P., Golub, T. R., Mesirov, J. P. (2004).
#' Metagenes and molecular pattern discovery using matrix factorization.
#' \emph{PNAS} 101(12):4164--4169. \doi{10.1073/pnas.0308531101}.
#' Kim, H., Park, H. (2007).  Sparse non-negative matrix factorizations
#' \dots \emph{Bioinformatics} 23(12):1495--1502.
#' @seealso \code{\link{nmfkc.ecv}}, \code{\link{nmfkc.bicv}},
#'   \code{\link{nmfkc.rank}}, \code{\link{nmf.cluster.criteria}}
#' @export
#' @examples
#' \donttest{
#' Y <- t(as.matrix(iris[, 1:4]))
#' cs <- nmfkc.consensus(Y, rank = 2:5, nrun = 20)
#' cs$cophenetic     # stability per rank (closer to 1 = more reproducible)
#' }
nmfkc.consensus <- function(Y, A = NULL, rank = 2:4, nrun = 30, seed = 123,
                            keep.consensus = FALSE, ...) {
  Y <- base::as.matrix(Y)
  N <- base::ncol(Y)

  ## nmfkc args: random init (required for stability), plain fast fit.
  ea <- base::list(...); ea$Q <- NULL; ea$rank <- NULL; ea$seed <- NULL
  ea$X.init <- "runif"; ea$detail <- "fast"; ea$print.dims <- FALSE

  cophenetic <- stats::setNames(base::numeric(base::length(rank)),
                                base::paste0("rank=", rank))
  dispersion <- cophenetic
  consensus.list <- if (keep.consensus)
    base::vector("list", base::length(rank)) else NULL

  base::message(base::sprintf(
    "consensus: ranks %s, %d runs/rank (Brunet 2004); %d fits total...",
    base::paste(rank, collapse = ","), nrun, nrun * base::length(rank)))

  for (ki in base::seq_along(rank)) {
    k <- rank[ki]
    Cmat <- base::matrix(0, N, N)
    for (r in base::seq_len(nrun)) {
      fit <- base::suppressMessages(base::do.call("nmfkc",
        c(base::list(Y = Y, A = A, rank = k,
                     seed = seed + 1000L * ki + r), ea)))
      cl <- base::apply(fit$B, 2, base::which.max)   # hard cluster per sample
      Cmat <- Cmat + (base::outer(cl, cl, "==") + 0) # connectivity
    }
    Cmat <- Cmat / nrun                              # consensus in [0, 1]
    if (keep.consensus) consensus.list[[ki]] <- Cmat

    ## Kim-Park dispersion: 1 = perfectly crisp (all 0/1).
    dispersion[ki] <- base::sum(4 * (Cmat - 0.5)^2) / (N * N)

    ## Cophenetic correlation of the consensus (Brunet).
    d <- stats::as.dist(1 - Cmat)
    if (k >= 2 && base::any(d > 0)) {
      hc <- stats::hclust(d, method = "average")
      cophenetic[ki] <- stats::cor(base::as.vector(d),
                                   base::as.vector(stats::cophenetic(hc)))
    } else {
      cophenetic[ki] <- NA_real_
    }
  }

  base::list(cophenetic = cophenetic, dispersion = dispersion,
             rank = rank, nrun = nrun, consensus = consensus.list)
}

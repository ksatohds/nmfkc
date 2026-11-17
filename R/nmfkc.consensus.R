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
#'   \item \code{pac}: the Proportion of Ambiguous Clustering
#'     (Senbabaoglu et al. 2014) -- the fraction of off-diagonal consensus
#'     entries falling in the ambiguous interval \code{pac.range}
#'     (default \eqn{(0.1, 0.9)}).  \strong{Lower is better} (less
#'     ambiguity); a more sensitive readout than \code{cophenetic}, which
#'     tends to saturate.
#' }
#' Unlike the cross-validation engines (where the rank \emph{minimizes} the
#' error), here a good rank \strong{maximizes} stability, or is the largest
#' rank before it drops.
#'
#' @param Y Observation matrix (\eqn{P \times N}), non-negative.
#' @param A Optional covariate matrix passed to \code{\link{nmfkc}}.
#' @param rank Integer vector of ranks to evaluate (\eqn{\ge 2}).
#' @param nrun Number of random-initialization runs per rank (default 30).
#' @param keep.consensus Logical; if \code{TRUE} also return the list of
#'   consensus matrices (one \eqn{N \times N} matrix per rank).
#' @param ... Advanced options, rarely needed (defaults in parentheses):
#'   \code{seed} (\code{123}, base seed; run \eqn{r} of rank index \eqn{i} uses
#'   \code{seed + 1000 * i + r}) and \code{pac.range} (\code{c(0.1, 0.9)}, the
#'   ambiguous interval \eqn{(u_1, u_2)} for the PAC measure).  Any other
#'   arguments are passed to \code{\link{nmfkc}} (e.g.\ \code{maxit});
#'   \code{X.init} is forced to \code{"runif"}.  \code{cores}
#'   (\code{getOption("mc.cores", 1L)}) runs the \code{nrun} restarts of each
#'   rank in parallel; results are identical to the sequential loop for any
#'   \code{cores} (PSOCK cluster on Windows, forking elsewhere).
#' @return An object of class \code{"nmfkc.consensus"} (a list) with:
#'   \item{cophenetic}{Cophenetic correlation coefficient for each rank.}
#'   \item{dispersion}{Dispersion coefficient (\eqn{[0, 1]}) for each rank.}
#'   \item{pac}{Proportion of Ambiguous Clustering (\eqn{[0, 1]}, lower is
#'     better) for each rank.}
#'   \item{rank}{The evaluated rank vector.}
#'   \item{nrun}{Number of runs per rank.}
#'   \item{consensus}{List of consensus matrices, or \code{NULL}.}
#'   It has \code{\link{print.nmfkc.consensus}} and
#'   \code{\link{plot.nmfkc.consensus}} (\code{type = "criteria"} /
#'   \code{"heatmap"}) methods.
#' @references
#' Brunet, J.-P., Tamayo, P., Golub, T. R., Mesirov, J. P. (2004).
#' Metagenes and molecular pattern discovery using matrix factorization.
#' \emph{PNAS} 101(12):4164--4169. \doi{10.1073/pnas.0308531101}.
#' Kim, H., Park, H. (2007).  Sparse non-negative matrix factorizations
#' \dots \emph{Bioinformatics} 23(12):1495--1502.
#' Senbabaoglu, Y., Michailidis, G., Li, J. Z. (2014).  Critical
#' limitations of consensus clustering in class discovery.
#' \emph{Sci. Rep.} 4:6207. \doi{10.1038/srep06207}.
#' @seealso \code{\link{nmfkc.rank}}, \code{\link{nmfkc.ecv}},
#'   \code{\link{nmfkc.bicv}}, \code{\link{nmfkc.ard}},
#'   \code{\link{nmf.cluster.criteria}}
#' @export
#' @examples
#' \donttest{
#' Y <- t(as.matrix(iris[, 1:4]))
#' cs <- nmfkc.consensus(Y, rank = 2:5, nrun = 20, keep.consensus = TRUE)
#' cs                       # stability table per rank
#' plot(cs)                 # type = "criteria": stability curves
#' plot(cs, type = "heatmap")            # all ranks, n2mfrow grid
#' plot(cs, type = "heatmap", rank = 3)  # one rank, with labels
#' }
nmfkc.consensus <- function(Y, A = NULL, rank = 2:4, nrun = 30,
                            keep.consensus = FALSE, ...) {
  Y <- base::as.matrix(Y)
  N <- base::ncol(Y)

  ## Advanced options live in `...` (the rest passes to nmfkc per run).
  dots      <- base::list(...)
  seed      <- if (base::is.null(dots$seed))      123         else dots$seed
  pac.range <- if (base::is.null(dots$pac.range)) c(0.1, 0.9) else dots$pac.range
  cores     <- if (base::is.null(dots$cores)) base::getOption("mc.cores", 1L) else dots$cores
  dots$seed <- NULL; dots$pac.range <- NULL; dots$cores <- NULL

  ## nmfkc args: random init (required for stability), plain fast fit.
  ea <- dots; ea$Q <- NULL; ea$rank <- NULL; ea$seed <- NULL
  ea$X.init <- "runif"; ea$detail <- "fast"; ea$print.dims <- FALSE

  cophenetic <- stats::setNames(base::numeric(base::length(rank)),
                                base::paste0("rank=", rank))
  dispersion <- cophenetic
  pac        <- cophenetic
  consensus.list <- if (keep.consensus)
    base::vector("list", base::length(rank)) else NULL

  base::message(base::sprintf(
    "consensus: ranks %s, %d runs/rank (Brunet 2004); %d fits total...",
    base::paste(rank, collapse = ","), nrun, nrun * base::length(rank)))

  for (ki in base::seq_along(rank)) {
    k <- rank[ki]
    ## Each of the nrun restarts is an independent, self-seeded nmfkc fit
    ## (seed + 1000*ki + r). .nmfkc.parlapply preserves input order, so the
    ## per-run connectivity matrices are summed in the same r order as the
    ## sequential loop; the integer-valued sum is bit-identical for any cores.
    run_one <- function(r) {
      fit <- base::suppressMessages(base::do.call("nmfkc",
        c(base::list(Y = Y, A = A, rank = k,
                     seed = seed + 1000L * ki + r), ea)))
      cl <- base::apply(fit$B, 2, base::which.max)   # hard cluster per sample
      base::outer(cl, cl, "==") + 0                  # connectivity
    }
    conn <- .nmfkc.parlapply(base::seq_len(nrun), run_one,
                             cores = if (nrun > 1L) cores else 1L)
    Cmat <- base::Reduce("+", conn) / nrun           # consensus in [0, 1]
    if (keep.consensus) {
      if (!base::is.null(base::colnames(Y)))
        base::dimnames(Cmat) <- base::list(base::colnames(Y), base::colnames(Y))
      consensus.list[[ki]] <- Cmat
    }
    cr <- .consensus.criteria(Cmat, k, pac.range)
    dispersion[ki] <- cr$dispersion
    pac[ki] <- cr$pac
    cophenetic[ki] <- cr$cophenetic
  }

  out <- base::list(cophenetic = cophenetic, dispersion = dispersion,
                    pac = pac, rank = rank, nrun = nrun,
                    consensus = consensus.list)
  base::class(out) <- "nmfkc.consensus"
  out
}


## The three stability summaries of one consensus matrix (shared by
## nmfkc.consensus() and nmfkc.net.consensus()).
.consensus.criteria <- function(Cmat, k, pac.range) {
  N <- base::nrow(Cmat)
  ## Kim-Park dispersion: 1 = perfectly crisp (all 0/1).
  dispersion <- base::sum(4 * (Cmat - 0.5)^2) / (N * N)

  ## PAC (Senbabaoglu 2014): proportion of off-diagonal consensus entries
  ## in the ambiguous interval (pac.range); LOWER = crisper / better.
  lt <- Cmat[base::lower.tri(Cmat)]
  pac <- base::mean(lt > pac.range[1] & lt < pac.range[2])

  ## Cophenetic correlation of the consensus (Brunet).
  d <- stats::as.dist(1 - Cmat)
  cophenetic <- if (k >= 2 && base::any(d > 0)) {
    hc <- stats::hclust(d, method = "average")
    stats::cor(base::as.vector(d), base::as.vector(stats::cophenetic(hc)))
  } else NA_real_
  base::list(dispersion = dispersion, pac = pac, cophenetic = cophenetic)
}


#' @title Consensus-clustering rank selection for nmfkc.net (Brunet 2004)
#' @description
#' The consensus-clustering stability of \code{\link{nmfkc.consensus}}, for
#' the network models of \code{\link{nmfkc.net}}
#' (\eqn{Y \approx X C X^\top}, \eqn{X X^\top} or \eqn{X (C_+ - C_-) X^\top}).
#' The clustered objects are the \eqn{N} nodes.
#'
#' For each rank, \code{\link{nmfkc.net}} is run \code{nrun} times with
#' different seeds, from a k-means++ start by default (see below).  Each run
#' gives a hard clustering of the nodes, the \code{X.cluster} of the fit (the
#' largest entry of each row of \eqn{X}).  The \eqn{N \times N} connectivity
#' matrices (1 if two nodes share a cluster) are averaged over the runs into
#' the consensus matrix, which is summarized by \code{cophenetic},
#' \code{dispersion} and \code{pac} exactly as in
#' \code{\link{nmfkc.consensus}}.
#'
#' The clustering is taken from \eqn{X}, not from the column maxima of
#' \code{fit$B} \eqn{= C X^\top} that \code{\link{nmfkc.consensus}} uses:
#' with off-diagonal entries in \eqn{C} (or negative ones, for
#' \code{type = "signed"}), the columns of \eqn{C X^\top} mix the bases.
#'
#' \strong{Starting values.}  The consensus needs starts that differ from run
#' to run: at too small a rank, which blocks are merged then changes from run
#' to run, and the rank shows as unstable.  The default is
#' \code{X.init = "kmeans++"} (\eqn{D^2} seeding, then Lloyd's iterations),
#' which differs with the seed.  Two other starts fail in opposite ways.
#' \itemize{
#'   \item The default of \code{\link{nmfkc.net}}, \code{"kmeans"}, keeps the
#'     best of 10 k-means restarts.  The runs then nearly always merge the
#'     same blocks, and a too small rank looks as stable as the true one.
#'   \item The uniform start of \code{\link{nmfkc.consensus}},
#'     \code{"runif"}, stops in poor local minima for \code{"tri"} and
#'     \code{"signed"}, so the consensus measures the optimizer rather than
#'     the data.
#' }
#' On 20 networks of three planted blocks (10 nodes each, links with
#' probability 0.8 inside a block and 0.05 between blocks), rank 3 had the
#' largest dispersion among ranks 2 to 6, with no tie, in 19 or 20 networks
#' for each type with \code{"kmeans++"}.  With \code{"kmeans"} it did in 3;
#' in most of the others ranks 2 and 3 were both fully stable.  With
#' \code{"runif"} it did in 16 for \code{"bi"} and in none for \code{"tri"}
#' and \code{"signed"}.  With \code{"nndsvd"} (NNDSVDar, whose only
#' randomness is the small fill of its zeros) it did in 11 for \code{"tri"},
#' 12 for \code{"signed"} and none for \code{"bi"}.  With weaker blocks (0.4
#' and 0.15) or blocks of unequal size, no start picked rank 3 reliably.  Stability favours small
#' ranks, so read it together with \code{\link{nmfkc.net.rank}}.  Other
#' starts can be given through \code{X.init}.
#'
#' @param Y Symmetric (\eqn{N \times N}) non-negative adjacency matrix, as
#'   for \code{\link{nmfkc.net}}.
#' @param rank Integer vector of ranks to evaluate.
#' @param type \code{"tri"} (default), \code{"bi"} or \code{"signed"}.
#' @param nrun Number of runs (seeds) per rank (default 30).
#' @param keep.consensus Logical; if \code{TRUE} also return the list of
#'   consensus matrices (one \eqn{N \times N} matrix per rank, named by the
#'   nodes).
#' @param ... Advanced options, as for \code{\link{nmfkc.consensus}}:
#'   \code{seed} (\code{123}; run \eqn{r} of rank index \eqn{i} uses
#'   \code{seed + 1000 * i + r}), \code{pac.range} (\code{c(0.1, 0.9)}) and
#'   \code{cores} (\code{getOption("mc.cores", 1L)}; results are identical for
#'   any \code{cores}).  Any other arguments are passed to
#'   \code{\link{nmfkc.net}} (e.g.\ \code{maxit}, \code{epsilon},
#'   \code{diag.exclude}, \code{Y.weights}, \code{X.init}, \code{nstart}).
#' @return An object of class \code{"nmfkc.consensus"} with the components of
#'   \code{\link{nmfkc.consensus}} and \code{type}; its
#'   \code{\link{print.nmfkc.consensus}} and \code{\link{plot.nmfkc.consensus}}
#'   methods apply.
#' @references
#' Brunet, J.-P., Tamayo, P., Golub, T. R., Mesirov, J. P. (2004).
#' Metagenes and molecular pattern discovery using matrix factorization.
#' \emph{PNAS} 101(12):4164--4169. \doi{10.1073/pnas.0308531101}.
#' @seealso \code{\link{nmfkc.net}}, \code{\link{nmfkc.net.rank}},
#'   \code{\link{nmfkc.net.ecv}}, \code{\link{nmfkc.consensus}}
#' @export
#' @examples
#' \donttest{
#' ## three planted blocks of 10 nodes, linked with probability 0.8 inside a
#' ## block and 0.05 between blocks
#' set.seed(1)
#' g <- rep(1:3, each = 10)
#' P <- ifelse(outer(g, g, "=="), 0.8, 0.05)
#' Y <- matrix(rbinom(length(P), 1, P), nrow(P))
#' Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]; diag(Y) <- 0
#' cs <- nmfkc.net.consensus(Y, rank = 2:6, nrun = 20, keep.consensus = TRUE)
#' cs                                    # the most stable rank is 3
#' plot(cs)
#' plot(cs, type = "heatmap", rank = 3)
#' }
nmfkc.net.consensus <- function(Y, rank = 2:4, type = c("tri", "bi", "signed"),
                                nrun = 30, keep.consensus = FALSE, ...) {
  type <- base::match.arg(type)
  Y <- base::as.matrix(Y)
  if (base::nrow(Y) != base::ncol(Y))
    base::stop("Y must be a square (N x N) adjacency matrix.")

  dots      <- base::list(...)
  seed      <- if (base::is.null(dots$seed))      123         else dots$seed
  pac.range <- if (base::is.null(dots$pac.range)) c(0.1, 0.9) else dots$pac.range
  cores     <- if (base::is.null(dots$cores)) base::getOption("mc.cores", 1L) else dots$cores
  dots$seed <- NULL; dots$pac.range <- NULL; dots$cores <- NULL

  ## nmfkc.net args; k-means++ starts unless X.init is given (see the help).
  ea <- dots; ea$rank <- NULL; ea$Q <- NULL; ea$type <- NULL; ea$verbose <- NULL
  if (base::is.null(ea$X.init)) ea$X.init <- "kmeans++"

  nodes <- if (!base::is.null(base::rownames(Y))) base::rownames(Y) else base::colnames(Y)
  cophenetic <- stats::setNames(base::numeric(base::length(rank)),
                                base::paste0("rank=", rank))
  dispersion <- cophenetic
  pac        <- cophenetic
  consensus.list <- if (keep.consensus)
    base::vector("list", base::length(rank)) else NULL

  base::message(base::sprintf(
    "consensus (nmfkc.net, type = %s): ranks %s, %d runs/rank (Brunet 2004); %d fits total...",
    type, base::paste(rank, collapse = ","), nrun, nrun * base::length(rank)))

  for (ki in base::seq_along(rank)) {
    k <- rank[ki]
    ## Independent, self-seeded fits; .nmfkc.parlapply keeps the order, so the
    ## integer-valued sum is identical for any cores.
    run_one <- function(r) {
      fit <- base::suppressMessages(base::do.call("nmfkc.net",
        c(base::list(Y = Y, rank = k, type = type, verbose = FALSE,
                     seed = seed + 1000L * ki + r), ea)))
      cl <- fit$X.cluster                              # hard cluster per node
      base::outer(cl, cl, "==") + 0                    # connectivity
    }
    conn <- .nmfkc.parlapply(base::seq_len(nrun), run_one,
                             cores = if (nrun > 1L) cores else 1L)
    Cmat <- base::Reduce("+", conn) / nrun
    base::dimnames(Cmat) <- if (!base::is.null(nodes)) base::list(nodes, nodes) else NULL
    if (keep.consensus) consensus.list[[ki]] <- Cmat
    cr <- .consensus.criteria(Cmat, k, pac.range)
    dispersion[ki] <- cr$dispersion
    pac[ki] <- cr$pac
    cophenetic[ki] <- cr$cophenetic
  }

  out <- base::list(cophenetic = cophenetic, dispersion = dispersion,
                    pac = pac, rank = rank, nrun = nrun,
                    consensus = consensus.list, type = type)
  base::class(out) <- "nmfkc.consensus"
  out
}


#' @title Print method for nmfkc.consensus objects
#' @param x An object of class \code{"nmfkc.consensus"}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.nmfkc.consensus <- function(x, ...) {
  base::cat(base::sprintf("Consensus rank selection (Brunet 2004), %d runs/rank%s\n",
                          x$nrun,
                          if (base::is.null(x$type)) ""
                          else base::sprintf(", nmfkc.net(type = \"%s\")", x$type)))
  base::cat(base::sprintf("  best: dispersion max = rank %d | PAC min = rank %d\n",
                          x$rank[base::which.max(x$dispersion)],
                          x$rank[base::which.min(x$pac)]))
  df <- base::data.frame(rank = x$rank,
                         cophenetic = base::round(x$cophenetic, 3),
                         dispersion = base::round(x$dispersion, 3))
  if (!base::is.null(x$pac)) df$pac <- base::round(x$pac, 3)
  base::print(df, row.names = FALSE)
  if (base::is.null(x$consensus))
    base::cat("  (consensus matrices not stored; use keep.consensus = TRUE for heatmaps)\n")
  base::invisible(x)
}


#' @title Plot a consensus rank-selection (nmfkc.consensus) object
#' @description
#' Two views of a \code{\link{nmfkc.consensus}} result:
#' \itemize{
#'   \item \code{type = "criteria"} (default): the stability curves
#'     \code{cophenetic} (blue) and \code{dispersion} (red) against rank.
#'     No "best" marker is drawn -- a stable rank tends to be a coarse one,
#'     so the curves are shown for inspection rather than as an optimum.
#'   \item \code{type = "heatmap"}: the consensus matrix of each requested
#'     \code{rank}, reordered by average-linkage hierarchical clustering of
#'     \eqn{1 - \bar C} (blue = 0 / different cluster, red = 1 / same
#'     cluster).  Requires \code{keep.consensus = TRUE} at compute time.
#' }
#' @param x An object of class \code{"nmfkc.consensus"}.
#' @param type \code{"criteria"} or \code{"heatmap"}.
#' @param rank For \code{type = "heatmap"}, the rank(s) to display.
#'   \code{NULL} (default) shows every rank stored in \code{x}.
#' @param mfrow Panel layout for multiple heatmaps, as
#'   \code{c(nrow, ncol)}.  \code{NULL} (default) uses
#'   \code{\link[grDevices]{n2mfrow}} for a near-square grid.
#' @param mar Per-panel margins \code{c(b, l, t, r)} for the heatmap
#'   view.  \code{NULL} (default) uses tight margins for a grid and wider
#'   margins (for sample labels) when a single rank is shown.
#' @param col Heatmap colour palette (length-50 blue-to-red by default).
#' @param main Plot title for the \code{"criteria"} view (heatmap panels
#'   are titled per rank).
#' @param ... Further arguments passed to the underlying
#'   \code{\link[graphics]{plot}} / \code{\link[graphics]{image}}.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{nmfkc.consensus}}
#' @export
plot.nmfkc.consensus <- function(x, type = c("criteria", "heatmap"),
                                 rank = NULL, mfrow = NULL, mar = NULL,
                                 col = grDevices::hcl.colors(50, "Blue-Red"),
                                 main = NULL, ...) {
  type <- base::match.arg(type)

  if (type == "criteria") {
    graphics::plot(x$rank, x$cophenetic, type = "o", col = "blue", pch = 16,
                   lwd = 2, ylim = c(0, 1), xlab = "rank (Q)",
                   ylab = "stability (0-1)",
                   main = if (base::is.null(main)) "Consensus stability" else main,
                   ...)
    graphics::lines(x$rank, x$dispersion, type = "o", col = "red", pch = 17,
                    lwd = 2)
    leg <- c("cophenetic (Brunet, high=good)", "dispersion (Kim-Park, high=good)")
    lcol <- c("blue", "red"); lpch <- c(16, 17)
    if (!base::is.null(x$pac)) {
      graphics::lines(x$rank, x$pac, type = "o", col = "darkgreen", pch = 15,
                      lwd = 2)
      leg <- c(leg, "PAC (Senbabaoglu, low=good)")
      lcol <- c(lcol, "darkgreen"); lpch <- c(lpch, 15)
    }
    graphics::legend("left", leg, col = lcol, pch = lpch, lwd = 2,
                     bg = "white", cex = 0.78)
    return(base::invisible(x))
  }

  ## type == "heatmap"
  if (base::is.null(x$consensus))
    base::stop("No consensus matrices stored. Re-run with keep.consensus = TRUE.")
  ranks <- if (base::is.null(rank)) x$rank else base::intersect(rank, x$rank)
  if (!base::length(ranks))
    base::stop("Requested rank(s) are not among the evaluated ranks.")
  n <- base::length(ranks)
  if (base::is.null(mfrow)) mfrow <- grDevices::n2mfrow(n)
  show.labels <- (n == 1)
  if (base::is.null(mar))
    mar <- if (show.labels) c(6, 6, 3, 1) else c(1.5, 1.5, 3, 1)
  old <- graphics::par(mfrow = mfrow, mar = mar, oma = c(2, 0, 0, 0))
  base::on.exit(graphics::par(old))

  for (k in ranks) {
    i <- base::which(x$rank == k)
    C <- x$consensus[[i]]
    hc <- stats::hclust(stats::as.dist(1 - C), method = "average")
    Cr <- C[hc$order, hc$order]; N <- base::nrow(Cr)
    graphics::image(1:N, 1:N, base::t(Cr[N:1, ]), col = col, zlim = c(0, 1),
                    axes = FALSE, xlab = "", ylab = "",
                    main = base::sprintf("rank=%d  coph=%.3f  disp=%.3f",
                                         k, x$cophenetic[i], x$dispersion[i]),
                    ...)
    if (show.labels && !base::is.null(base::rownames(Cr))) {
      cax <- base::min(0.7, 28 / N)
      graphics::axis(1, 1:N, base::colnames(Cr), las = 2, cex.axis = cax)
      graphics::axis(2, 1:N, base::rev(base::rownames(Cr)), las = 2, cex.axis = cax)
    } else {
      graphics::box()
    }
  }
  graphics::mtext("blue = 0 (different cluster)  ->  red = 1 (same cluster)",
                  side = 1, outer = TRUE, cex = 0.7)
  base::invisible(x)
}

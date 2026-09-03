# =====================================================
# nmfkc.rff.beta.cv.R -- bandwidth (and feature-count) selection for random
#                        feature covariates by cross-validation
#
# The counterpart of nmfkc.kernel.beta.cv() for Random Fourier Features:
# for each candidate beta (and each D) the random feature matrix is
# regenerated on the training columns and the fit is cross-validated with
# nmfkc.signed.cv() (signed cosine features) or nmfkc.cv() (positive
# features).  For large N the selection can run on a random subsample of
# the columns (the bandwidth choice is stable on a subsample) while the
# final fit uses the block-wise Gram route on all N.
# =====================================================


#' Select the bandwidth of random-feature covariates by cross-validation
#'
#' @description
#' Chooses the Gaussian-kernel bandwidth \eqn{\beta} (and optionally the
#' number of random features \eqn{D}) for Random Fourier Feature
#' covariates by cross-validation, in the same way that
#' \code{\link{nmfkc.kernel.beta.cv}} does for exact kernel covariates.
#' For every candidate the random features are regenerated on the columns
#' of \code{U} and the model is cross-validated column-wise:
#' \itemize{
#'   \item \code{type = "signed"}: cosine features from
#'     \code{\link{nmfkc.signed.rff}} (signed), fitted with
#'     \code{\link{nmfkc.signed}} through \code{\link{nmfkc.signed.cv}};
#'   \item \code{type = "positive"}: positive random features from
#'     \code{\link{nmfkc.rff.positive}} (non-negative), fitted with
#'     \code{\link{nmfkc}} through \code{\link{nmfkc.cv}}.
#' }
#' The candidate with the smallest cross-validated prediction error is
#' returned.  Since the same random map cannot be shared across folds
#' without leaking (the map is data-independent, so it can: only the fit
#' is cross-validated), one draw per candidate is used and its
#' \code{seed} is fixed so the comparison across \eqn{\beta} is paired.
#'
#' For large \eqn{N} set \code{sample.size}: the selection then runs on a
#' random subsample of the columns (default: all columns), after which the
#' chosen \eqn{\beta} is used with \code{\link{nmfkc.signed.rff.gram}} /
#' \code{\link{nmfkc.rff.positive.gram}} on the full data.
#'
#' @param Y Response matrix \eqn{P \times N} (non-negative for
#'   \code{type = "positive"}).
#' @param rank Number of bases \eqn{Q}.
#' @param U Input matrix \eqn{p \times N}; columns are samples.
#' @param beta Numeric vector of candidate bandwidths.  Default
#'   \code{NULL}: the seven-point grid of
#'   \code{\link{nmfkc.kernel.beta.nearest.med}} around the median
#'   nearest-neighbour heuristic (computed on the subsample when
#'   \code{sample.size} is set).
#' @param D Integer vector of candidate numbers of random features.  A
#'   single value (default \code{500}) selects \eqn{\beta} only; several
#'   values select over the \eqn{(\beta, D)} grid.
#' @param type \code{"signed"} (cosine RFF, default) or \code{"positive"}.
#' @param plot Logical.  If \code{TRUE} (default), plot the CV objective
#'   against \code{beta} (one line per \code{D}).
#' @param ... Passed to the CV function (\code{\link{nmfkc.signed.cv}} or
#'   \code{\link{nmfkc.cv}}: e.g. \code{nfolds}, \code{epsilon},
#'   \code{maxit}).  Also accepts:
#'   \itemize{
#'     \item \code{sample.size}: number of columns used for the selection
#'       (default all).
#'     \item \code{seed}: integer seed for the subsample and the random
#'       maps (default 123; the caller's stream is restored).
#'     \item \code{cores}: evaluate the candidates in parallel with
#'       \code{.nmfkc.parlapply} (default \code{getOption("mc.cores", 1L)}).
#'       Each candidate is deterministic given its index, so the result is
#'       identical for any \code{cores}.
#'     \item \code{hyperbolic}: passed to \code{\link{nmfkc.rff.positive}}.
#'   }
#'
#' @return A list with
#'   \describe{
#'     \item{\code{beta}, \code{D}}{The selected bandwidth and feature count.}
#'     \item{\code{objfunc}}{Matrix of CV objective values,
#'       \code{length(beta)} rows by \code{length(D)} columns.}
#'     \item{\code{beta.candidates}, \code{D.candidates}}{The grids.}
#'     \item{\code{type}, \code{sample.size}, \code{n.used}}{Settings used.}
#'   }
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.  The interface may change in
#' future versions.
#'
#' @seealso \code{\link{nmfkc.kernel.beta.cv}} (exact kernel; also works for
#'   Nystr\"om covariates with \code{U = landmarks, V = data}),
#'   \code{\link{nmfkc.signed.rff.gram}}, \code{\link{nmfkc.rff.positive.gram}}
#'
#' @examples
#' \donttest{
#' data(iris)
#' set.seed(1)
#' idx <- sample(nrow(iris), 100)
#' U <- t(scale(iris[idx, 1:4]))
#' levs <- levels(iris$Species)
#' Y <- sapply(iris$Species[idx], function(s) as.integer(levs == s)); rownames(Y) <- levs
#'
#' ## signed RFF: choose beta over the median-heuristic grid at D = 50
#' sel <- nmfkc.rff.beta.cv(Y, rank = 3, U, D = 50, plot = FALSE)
#' sel$beta
#' g   <- nmfkc.signed.rff.gram(Y, U, beta = sel$beta, D = 50, seed = 1)
#' fit <- nmfkc.signed(Y, A = g, rank = 3, verbose = FALSE)
#'
#' ## positive RFF: select over a (beta, D) grid
#' sel2 <- nmfkc.rff.beta.cv(Y, rank = 3, U, D = c(100, 400), type = "positive",
#'                           plot = FALSE, verbose = FALSE)
#' sel2$objfunc
#' }
#'
#' @export
nmfkc.rff.beta.cv <- function(Y, rank = 2, U, beta = NULL, D = 500L,
                              type = c("signed", "positive"), plot = TRUE, ...) {
  type <- match.arg(type)
  extra <- list(...)
  if (!is.null(extra$Q)) rank <- extra$Q
  sample.size <- extra$sample.size
  seed        <- if (!is.null(extra$seed))  extra$seed  else 123L
  cores       <- if (!is.null(extra$cores)) extra$cores else getOption("mc.cores", 1L)
  hyperbolic  <- if (!is.null(extra$hyperbolic)) isTRUE(extra$hyperbolic) else FALSE
  cv_args <- extra[!names(extra) %in% c("Q", "sample.size", "seed", "cores", "hyperbolic")]

  if (is.vector(Y) && !is.matrix(Y)) Y <- matrix(Y, nrow = 1)
  Y <- as.matrix(Y); U <- as.matrix(U); storage.mode(U) <- "double"
  N <- ncol(U)
  if (ncol(Y) != N) stop("ncol(Y) must equal ncol(U).")
  D <- as.integer(D)
  if (any(is.na(D)) || any(D < 1L)) stop("'D' must be positive integers.")

  .rng <- .nmfkc.rng.save(seed)
  on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
  set.seed(seed)

  ## subsample of columns for the selection
  if (!is.null(sample.size) && sample.size < N) {
    sub <- sort(sample.int(N, as.integer(sample.size)))
    Ys <- Y[, sub, drop = FALSE]; Us <- U[, sub, drop = FALSE]
  } else {
    Ys <- Y; Us <- U
  }
  n_used <- ncol(Us)

  if (is.null(beta)) {
    beta <- nmfkc.kernel.beta.nearest.med(Us)$beta_candidates
    if (is.null(beta) || length(beta) == 0)
      stop("Failed to determine beta candidates from the nearest-neighbour median.")
  }
  beta <- as.numeric(beta)

  grid <- expand.grid(ib = seq_along(beta), id = seq_along(D))
  run_one <- function(k) {
    b <- beta[grid$ib[k]]; d <- D[grid$id[k]]
    map_seed <- seed + 1000L * grid$id[k] + grid$ib[k]     # paired across candidates
    v <- tryCatch({
      if (type == "signed") {
        Z  <- nmfkc.signed.rff(Us, beta = b, D = d, seed = map_seed)$Z
        cv <- do.call(nmfkc.signed.cv, c(list(Y = Ys, A = Z, rank = rank), cv_args))
      } else {
        Z  <- nmfkc.rff.positive(Us, beta = b, D = d, seed = map_seed, hyperbolic = hyperbolic)$Z
        cv <- do.call(nmfkc.cv, c(list(Y = Ys, A = Z, rank = rank), cv_args))
      }
      as.numeric(cv$objfunc)[1]
    }, error = function(e) {
      warning("candidate beta = ", format(b, digits = 4), ", D = ", d, " failed (",
              conditionMessage(e), "); treated as Inf.", call. = FALSE)
      Inf
    })
    ## a candidate whose fit blew up (non-finite objective) simply loses
    if (!is.finite(v)) v <- Inf
    v
  }
  vals <- unlist(.nmfkc.parlapply(seq_len(nrow(grid)), run_one, cores = cores))
  objfunc <- matrix(vals, nrow = length(beta), ncol = length(D),
                    dimnames = list(format(beta, digits = 4), paste0("D=", D)))
  if (all(!is.finite(objfunc))) stop("every candidate produced a non-finite CV objective.")
  best <- which(objfunc == min(objfunc, na.rm = TRUE), arr.ind = TRUE)[1, ]

  if (isTRUE(plot)) {
    graphics::matplot(log10(beta), objfunc, type = "b", pch = 16, lty = 1,
                      xlab = expression(log[10](beta)), ylab = "CV objective",
                      main = sprintf("nmfkc.rff.beta.cv (%s)", type))
    if (length(D) > 1)
      graphics::legend("topright", legend = colnames(objfunc), col = seq_along(D),
                       lty = 1, pch = 16, bty = "n")
    graphics::abline(v = log10(beta[best[1]]), lty = 2)
  }

  list(beta = beta[best[1]], D = D[best[2]], objfunc = objfunc,
       beta.candidates = beta, D.candidates = D, type = type,
       sample.size = sample.size, n.used = n_used)
}

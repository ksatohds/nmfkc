# =====================================================
# nmfkc.kernel.gram.R -- block-wise Gram accumulation of Nystroem kernel
#                        covariates for nmfkc(), plus the S3 print method
#                        and the guard shared by every "nmfkc.gram" object
#
# NMF-LAB Kernel with M landmarks v_1..v_M (Satoh 2026, JJSD) uses the
# covariate matrix
#     A = C^T,   C[n, m] = k(u_n, v_m),          A is M x N, non-negative,
# in place of the N x N kernel matrix.  The Euclidean multiplicative
# updates of nmfkc() touch A only through
#     S  = A A^T   (M x M)      and      G0 = Y A^T   (P x M),
# so for large N the M x N matrix never has to exist.  This file walks the
# columns of U in blocks, computes the kernel block k(V, U_b) with
# nmfkc.kernel(), accumulates S and G0, and returns them as a "nmfkc.gram"
# object that nmfkc() accepts in place of A.  Landmarks can be supplied or
# chosen here by k-means++ / k-means on a random subsample (the paper's
# recipe: k-means centroids of a 10,000-image subsample).
#
# The same class is produced for signed Random Fourier Features by
# nmfkc.signed.rff.gram(); the `signed` field tells the fitters apart
# (nmfkc() refuses signed objects, nmfkc.signed() accepts both).
#
# References:
#   Williams, C. K. I. & Seeger, M. (2001). Using the Nystroem method to
#     speed up kernel machines. NIPS 13.
#   Zhang, K., Tsang, I. W. & Kwok, J. T. (2008). Improved Nystroem low-rank
#     approximation and error analysis. ICML.
#   Arthur, D. & Vassilvitskii, S. (2007). k-means++: the advantages of
#     careful seeding. SODA.
# =====================================================


#' Block-wise Gram accumulation of Nystroem kernel covariates for large N
#'
#' @description
#' Builds the two summary matrices that \code{\link{nmfkc}} needs from
#' the Nystr\"om kernel covariate matrix
#' \eqn{A = C^\top} (\eqn{M \times N}, \eqn{C_{nm} = k(\bm u_n, \bm v_m)})
#' without ever forming \eqn{A} itself:
#' \deqn{S = A A^\top \quad (M \times M), \qquad
#'       G_0 = Y A^\top \quad (P \times M).}
#' The input \code{U} is processed in column blocks; for each block the
#' kernel values \eqn{k(V, U_b)} are computed with
#' \code{\link{nmfkc.kernel}}, added into \eqn{S} and \eqn{G_0}, and
#' discarded.  Peak memory is
#' \eqn{O(M^2 + PM + M \cdot \mathtt{block.size})} instead of
#' \eqn{O(MN)}.
#'
#' The returned object is passed to \code{\link{nmfkc}} as its \code{A}
#' argument.  The multiplicative updates are unchanged, since they only
#' ever use \eqn{S} and \eqn{G_0}; the fit equals the one obtained from the
#' explicit matrix \code{nmfkc.kernel(V, U, beta = beta)} up to
#' floating-point summation order.  This is the large-\eqn{N} form of the
#' NMF-LAB Kernel design with Nystr\"om approximation (Satoh 2026).
#'
#' @param Y Non-negative \eqn{P \times N} response matrix (e.g. a one-hot
#'   label matrix).  \code{NA} is not allowed.
#' @param U A \eqn{p \times N} numeric matrix of inputs; columns are
#'   samples.
#' @param V Landmarks.  Either a \eqn{p \times M} numeric matrix of
#'   landmark points, or a single integer \eqn{M}, in which case \eqn{M}
#'   landmarks are chosen from a random subsample of the columns of
#'   \code{U} by the method in \code{landmarks} (default k-means++
#'   seeding followed by \eqn{k}-means refinement).
#' @param beta Positive scalar.  Gaussian kernel bandwidth
#'   \eqn{k(u, v) = \exp(-\beta \lVert u - v \rVert^2)}.  If \code{NULL}
#'   (default), the nearest-landmark median heuristic
#'   \code{\link{nmfkc.kernel.beta.nearest.med}(U_sub, Uk = V)} on the
#'   subsample is used, and the value is reported in a message.  For
#'   bandwidth selection build one Gram object per candidate in
#'   \code{$beta_candidates} of that helper (see Details).
#' @param block.size Integer.  Number of columns of \code{U} per block.
#'   Default: as many columns as keep one kernel block near 256 MB
#'   (\code{2^25 / M} columns), capped at \code{N}.
#' @param ... Hidden options:
#'   \itemize{
#'     \item \code{landmarks}: how to choose \eqn{M} landmarks when
#'       \code{V} is an integer: \code{"kmeans++"} (default; Arthur &
#'       Vassilvitskii 2007 seeding then Lloyd refinement),
#'       \code{"kmeans"} (\code{stats::kmeans} with \code{nstart = 5}), or
#'       \code{"random"} (uniform subsample).
#'     \item \code{sample.size}: size of the random subsample of columns
#'       used for landmark selection and the bandwidth heuristic
#'       (default \code{min(N, 10000)}).
#'     \item \code{seed}: integer seed for the subsample and the landmark
#'       selection (default 123; the caller's random stream is restored).
#'     \item \code{kernel}: kernel name passed to \code{\link{nmfkc.kernel}}
#'       (default \code{"Gaussian"}); further kernel parameters are passed
#'       through as well.
#'   }
#'
#' @details
#' With Gram input the fold-based \code{\link{nmfkc.cv}} is not available
#' (it must split the columns of \eqn{A}).  Select \eqn{\beta} as in the
#' NMF-LAB paper: fix the landmarks, build one Gram object per candidate
#' \eqn{\beta} on the training columns (pass the landmark matrix
#' \code{g$landmarks} as \code{V} so all candidates share it), fit with
#' \code{nmfkc()}, and score validation columns with
#' \code{predict(fit, newA = nmfkc.kernel(g$landmarks, U.val, beta = beta))}.
#'
#' @return An object of class \code{"nmfkc.gram"}: a list with
#'   \describe{
#'     \item{\code{S}}{\eqn{M \times M} matrix \eqn{A A^\top}
#'       (non-negative, positive semi-definite).}
#'     \item{\code{G0}}{\eqn{P \times M} matrix \eqn{Y A^\top}.}
#'     \item{\code{N}, \code{D}}{Number of samples and of covariates
#'       (\code{D = M}, the number of landmarks).}
#'     \item{\code{type}, \code{signed}}{\code{"nystrom"} and \code{FALSE}
#'       (non-negative features; accepted by \code{\link{nmfkc}} and
#'       \code{\link{nmfkc.signed}}).}
#'     \item{\code{landmarks}}{The \eqn{p \times M} landmark matrix
#'       \eqn{V}.  Pass it to \code{\link{nmfkc.kernel}} as its first
#'       argument to build covariates for new data.}
#'     \item{\code{beta}, \code{kernel}}{Bandwidth and kernel used.}
#'     \item{\code{landmarks.method}, \code{sample.size}}{How the landmarks
#'       were chosen (\code{"supplied"} when \code{V} was a matrix).}
#'     \item{\code{block.size}, \code{n.blocks}}{Blocking actually used.}
#'     \item{\code{A.block}}{A function \code{A.block(idx)} returning the
#'       \eqn{M \times |\mathtt{idx}|} kernel block \eqn{k(V, U_{idx})}.
#'       \code{\link{nmfkc}} uses it to rebuild \eqn{B = CA} block by
#'       block after the fit.  It closes over \code{U} and \code{V}, so the
#'       object keeps \code{U} alive (no copy is made).}
#'     \item{\code{rownames}}{Row names of \eqn{A} (\code{colnames(V)} if
#'       any, else \code{NULL}).}
#'   }
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.  The interface may change in
#' future versions.
#'
#' @references
#' Williams, C. K. I., & Seeger, M. (2001).  Using the Nystr\"om method to
#' speed up kernel machines.  \emph{Advances in NIPS}, 13.
#'
#' Zhang, K., Tsang, I. W., & Kwok, J. T. (2008).  Improved Nystr\"om
#' low-rank approximation and error analysis.  \emph{ICML}.
#'
#' Arthur, D., & Vassilvitskii, S. (2007).  k-means++: the advantages of
#' careful seeding.  \emph{SODA}.
#'
#' Satoh, K. (2026).  Applying non-negative matrix factorization with
#' covariates to label matrix for classification.  \emph{Japanese Journal
#' of Statistics and Data Science}.  \doi{10.1007/s42081-026-00349-x}
#'
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.kernel}},
#'   \code{\link{nmfkc.kernel.beta.nearest.med}},
#'   \code{\link{nmfkc.signed.rff.gram}} (the signed, Random-Fourier
#'   counterpart for \code{\link{nmfkc.signed}})
#'
#' @examples
#' \donttest{
#' ## Iris, 3 classes: Nystroem covariates with 12 k-means++ landmarks.
#' ## The Gram route gives the same fit as the explicit M x N matrix.
#' data(iris)
#' set.seed(1)
#' idx <- sample(nrow(iris), 100)
#' mn <- colMeans(iris[idx, 1:4]); sc <- apply(iris[idx, 1:4], 2, sd)
#' U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))   # 4 x 100
#' U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))   # 4 x  50
#' levs    <- levels(iris$Species)
#' Y.train <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
#' rownames(Y.train) <- levs
#'
#' ## 12 landmarks by k-means++, bandwidth by the nearest-landmark median
#' ## heuristic, S and G0 accumulated over blocks of 25 columns
#' g <- nmfkc.kernel.gram(Y.train, U.train, V = 12, block.size = 25, seed = 1)
#' g
#' res.g <- nmfkc(Y.train, A = g, rank = 3, verbose = FALSE)
#'
#' ## Same landmarks and beta, explicit 12 x 100 matrix: same fit
#' A <- nmfkc.kernel(g$landmarks, U.train, beta = g$beta)
#' res.m <- nmfkc(Y.train, A = A, rank = 3, verbose = FALSE)
#' all.equal(res.g$C, res.m$C)
#'
#' ## Test-set prediction: kernel values against the stored landmarks
#' A.test <- nmfkc.kernel(g$landmarks, U.test, beta = g$beta)
#' pred   <- predict(res.g, newA = A.test, type = "class")
#' mean(pred == as.character(iris$Species[-idx]))
#' }
#'
#' @export
nmfkc.kernel.gram <- function(Y, U, V, beta = NULL, block.size = NULL, ...) {
  extra <- list(...)
  landmarks   <- if (!is.null(extra$landmarks))   extra$landmarks   else "kmeans++"
  sample.size <- if (!is.null(extra$sample.size)) extra$sample.size else NULL
  seed        <- if (!is.null(extra$seed))        extra$seed        else 123L
  kernel      <- if (!is.null(extra$kernel))      extra$kernel      else "Gaussian"
  kern_extra  <- extra[setdiff(names(extra),
                               c("landmarks", "sample.size", "seed", "kernel"))]
  landmarks <- match.arg(landmarks, c("kmeans++", "kmeans", "random"))

  if (is.vector(Y) && !is.matrix(Y)) Y <- matrix(Y, nrow = 1)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(U)) U <- as.matrix(U)
  storage.mode(U) <- "double"
  N <- ncol(U); p <- nrow(U)
  if (ncol(Y) != N) stop("ncol(Y) must equal ncol(U).")
  if (any(is.na(Y)))
    stop("Y contains NA; a Gram object cannot mask missing entries.")
  if (min(Y) < 0) stop("The matrix Y should be non-negative.")
  if (any(is.na(U))) stop("U contains NA; please impute or remove.")

  ## Keep our own seeding (subsample, landmark selection) out of the
  ## caller's random stream (CONVENTIONS.md 2).
  .rng <- .nmfkc.rng.save(seed)
  on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
  set.seed(seed)

  ## ---- subsample: used for landmark selection and the beta heuristic ----
  if (is.null(sample.size)) sample.size <- min(N, 10000L)
  sample.size <- as.integer(min(sample.size, N))
  sub <- if (sample.size < N) sort(sample.int(N, sample.size)) else seq_len(N)
  U_sub <- U[, sub, drop = FALSE]

  ## ---- landmarks ----
  if (is.matrix(V) || (is.numeric(V) && length(V) > 1L)) {
    V <- as.matrix(V); storage.mode(V) <- "double"
    if (nrow(V) != p) stop("'V' must have nrow(U) = ", p, " rows (one landmark per column).")
    landmarks.method <- "supplied"
  } else {
    M <- as.integer(V)
    if (length(M) != 1L || is.na(M) || M < 1L)
      stop("'V' must be a p x M landmark matrix or a single positive integer M.")
    if (M > ncol(U_sub))
      stop("M = ", M, " landmarks requested but the subsample has only ",
           ncol(U_sub), " columns; raise sample.size.")
    pts <- t(U_sub)                                   # sample.size x p
    V <- switch(landmarks,
      "kmeans++" = {
        ## D^2-weighted seeding then Lloyd refinement (shared helper, the
        ## same one behind X.init = "kmeans++").  Fall back to the seeds if
        ## stats::kmeans() fails (e.g. duplicated points).
        seeds <- .kmeanspp.seed(pts, M)
        km <- tryCatch(stats::kmeans(pts, centers = seeds, iter.max = 100L),
                       error = function(e) NULL)
        if (is.null(km)) t(seeds) else t(km$centers)
      },
      "kmeans" = {
        km <- tryCatch(stats::kmeans(pts, centers = M, iter.max = 100L, nstart = 5L),
                       error = function(e) NULL)
        if (is.null(km)) pts_sel(pts, M) else t(km$centers)
      },
      "random" = pts_sel(pts, M))
    landmarks.method <- landmarks
  }
  M <- ncol(V)
  colnames(V) <- if (!is.null(colnames(V))) colnames(V) else NULL

  ## ---- bandwidth ----
  if (is.null(beta)) {
    if (kernel != "Gaussian")
      stop("'beta' must be specified for kernel = \"", kernel, "\".")
    beta <- nmfkc.kernel.beta.nearest.med(U_sub, Uk = V)$beta
    message("nmfkc.kernel.gram: beta = ", format(beta, digits = 4),
            " (nearest-landmark median heuristic on ", ncol(U_sub), " samples)")
  }

  ## ---- block-wise accumulation of S and G0 ----
  if (is.null(block.size)) block.size <- max(1L, min(N, as.integer(2^25 %/% M)))
  block.size <- as.integer(block.size)
  if (length(block.size) != 1L || is.na(block.size) || block.size < 1L)
    stop("'block.size' must be a positive integer.")
  block.size <- min(block.size, N)

  kernel_block <- function(Ub)
    unclass_matrix(do.call(nmfkc.kernel,
                           c(list(U = V, V = Ub, kernel = kernel, beta = beta), kern_extra)))

  P  <- nrow(Y)
  S  <- matrix(0, M, M)
  G0 <- matrix(0, P, M)
  starts <- seq.int(1L, N, by = block.size)
  for (s in starts) {
    idx <- s:min(s + block.size - 1L, N)
    Ab <- kernel_block(U[, idx, drop = FALSE])       # M x |idx|, >= 0
    S  <- S  + tcrossprod(Ab)                        # += A_b A_b^T
    G0 <- G0 + tcrossprod(Y[, idx, drop = FALSE], Ab) # += Y_b A_b^T
  }

  ## Block generator for nmfkc()'s post-fit B = C A.  Built in its own
  ## environment so the closure references U, V, beta, kernel only.
  A.block <- local({
    U <- U; kernel_block <- kernel_block
    function(idx) kernel_block(U[, idx, drop = FALSE])
  })

  structure(list(
    S                = S,
    G0               = G0,
    N                = N,
    D                = M,
    type             = "nystrom",
    signed           = FALSE,
    landmarks        = V,
    landmarks.method = landmarks.method,
    sample.size      = sample.size,
    beta             = beta,
    kernel           = kernel,
    block.size       = block.size,
    n.blocks         = length(starts),
    A.block          = A.block,
    rownames         = colnames(V)
  ), class = "nmfkc.gram")
}

## Uniform landmark subsample (rows of pts -> p x M matrix).
pts_sel <- function(pts, M) t(pts[sample.int(nrow(pts), M), , drop = FALSE])

## nmfkc.kernel() attaches bookkeeping attributes to its result; strip them
## so the accumulation and the block generator hand back plain matrices.
unclass_matrix <- function(K) {
  K <- unclass(K)
  attributes(K) <- list(dim = dim(K))
  K
}


#' Print method for nmfkc.gram
#'
#' @param x Object of class \code{"nmfkc.gram"} (from
#'   \code{\link{nmfkc.kernel.gram}} or \code{\link{nmfkc.signed.rff.gram}}).
#' @param ... Unused.
#' @return Invisible \code{x}.
#' @section Lifecycle:
#' This function is \strong{experimental}.
#' @seealso \code{\link{nmfkc.kernel.gram}}, \code{\link{nmfkc.signed.rff.gram}}
#' @export
print.nmfkc.gram <- function(x, ...) {
  type <- if (is.null(x$type)) "nystrom" else x$type
  is_nys <- identical(type, "nystrom")
  cat(switch(type,
    rff     = "Gram summary of Random Fourier Features (signed) for nmfkc.signed()\n",
    prf     = "Gram summary of positive random features (non-negative) for nmfkc()\n",
              "Gram summary of Nystroem kernel covariates (non-negative) for nmfkc()\n"))
  cat(sprintf("  N (samples):        %d\n", x$N))
  cat(sprintf("  D (%s):%s%d\n", if (is_nys) "landmarks" else "features",
              if (is_nys) "      " else "       ", x$D))
  cat(sprintf("  Y rows (P):         %d\n", nrow(x$G0)))
  if (!is.null(x$beta))
    cat(sprintf("  beta (bandwidth):   %s\n", format(x$beta, digits = 4)))
  if (identical(type, "prf") && isTRUE(x$pars$hyperbolic))
    cat("  variant:            hyperbolic (antithetic pairs)\n")
  if (is_nys && !is.null(x$landmarks.method))
    cat(sprintf("  landmarks:          %s%s\n", x$landmarks.method,
                if (x$landmarks.method == "supplied") ""
                else sprintf(" on a subsample of %d", x$sample.size)))
  cat(sprintf("  blocks:             %d x %d columns\n", x$n.blocks, x$block.size))
  cat(sprintf("  S  = A A^T:         %d x %d\n", nrow(x$S), ncol(x$S)))
  cat(sprintf("  G0 = Y A^T:         %d x %d\n", nrow(x$G0), ncol(x$G0)))
  cat(switch(type,
    rff = "Pass as `A` to nmfkc.signed(); use `$pars` with nmfkc.signed.rff() for new data.\n",
    prf = "Pass as `A` to nmfkc(); use `$pars` with nmfkc.rff.positive() for new data.\n",
          "Pass as `A` to nmfkc(); use nmfkc.kernel(g$landmarks, U.new, beta = g$beta) for new data.\n"))
  invisible(x)
}


## Fold-based helpers split the columns of A, which a Gram object cannot
## provide.  Fail with a pointer to the validation-set procedure rather than
## with as.matrix() on a list.
.nmfkc.no.gram <- function(A, fname) {
  if (inherits(A, "nmfkc.gram"))
    stop(fname, "() needs the covariate matrix itself: a Gram object ",
         "cannot be split into folds. For large N select tuning parameters ",
         "on a held-out validation set, scoring with predict() on covariates ",
         "built for the validation columns (nmfkc.kernel(g$landmarks, U.val, ",
         "beta = g$beta) or nmfkc.signed.rff(U.val, pars = g$pars)).",
         call. = FALSE)
  invisible(TRUE)
}

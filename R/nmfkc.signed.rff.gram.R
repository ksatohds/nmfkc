# =====================================================
# nmfkc.signed.rff.gram.R -- block-wise Gram accumulation of Random
#                            Fourier Features for nmfkc.signed()
#
# The unweighted Direct MU loop in nmfkc.signed() touches the covariate
# matrix A (D x N) only through
#     S  = A A^T   (D x D)      and      G0 = Y A^T   (Q_obs x D),
# so for large N the D x N matrix never has to exist.  With RFF the
# features are a deterministic function of (omega, b) and the input U, so
# they can be regenerated for any column block at will.  This file
# accumulates S and G0 over column blocks of U and returns them as a
# "nmfkc.gram" object that nmfkc.signed() accepts in place of A.
#
# The object also carries a block generator A.block(idx) so that
# nmfkc.signed() can rebuild B = C A (Q x N, small) after the fit, and so
# that predict() on new data is a matter of calling nmfkc.signed.rff()
# with the stored pars.
#
# References:
#   Rahimi, A. & Recht, B. (2007). Random features for large-scale
#     kernel machines. NIPS 20.
# =====================================================


#' Block-wise Gram accumulation of Random Fourier Features for large N
#'
#' @description
#' Builds the two summary matrices that \code{\link{nmfkc.signed}} needs
#' from a Random Fourier Feature (RFF) covariate matrix
#' \eqn{Z} (\eqn{D \times N}) without ever forming \eqn{Z} itself:
#' \deqn{S = Z Z^\top \quad (D \times D), \qquad
#'       G_0 = Y Z^\top \quad (Q_{\mathrm{obs}} \times D).}
#' The input \code{U} is processed in column blocks; for each block the
#' features \eqn{Z_b = \sqrt{2/D}\cos(\omega U_b + b)} are generated with
#' \code{\link{nmfkc.signed.rff}}, added into \eqn{S} and \eqn{G_0}, and
#' discarded.  Peak memory is therefore
#' \eqn{O(D^2 + Q_{\mathrm{obs}} D + D \cdot \mathtt{block.size})}
#' instead of \eqn{O(DN)}, which is what makes the RFF route usable when
#' \eqn{N} is in the hundreds of thousands and \eqn{D} exceeds the input
#' dimension \eqn{p} by a large factor (the \eqn{D \times N} matrix is
#' \eqn{D/p} times the size of the data itself).
#'
#' The returned object is passed to \code{\link{nmfkc.signed}} as its
#' \code{A} argument.  The multiplicative updates are unchanged, since
#' they only ever use \eqn{S} and \eqn{G_0}; the fit equals the one
#' obtained from the explicit matrix \code{nmfkc.signed.rff(U, ...)$Z}
#' with \code{warm.start = FALSE} up to floating-point summation order.
#'
#' @param Y Real-valued \eqn{Q_{\mathrm{obs}} \times N} response matrix
#'   (e.g. a one-hot label matrix).  \code{NA} is not allowed.
#' @param U A \eqn{p \times N} numeric matrix of inputs; columns are
#'   samples.
#' @param beta Positive scalar.  Gaussian kernel bandwidth
#'   \eqn{k(u, u') = \exp(-\beta \lVert u - u' \rVert^2)}.  Can be obtained
#'   with \code{\link{nmfkc.kernel.beta.nearest.med}}, which itself works
#'   on a random subsample.  Ignored when \code{pars} is supplied.
#' @param D Integer.  Number of random features.  Required (there is no
#'   \code{N/2} default here: on the large-\eqn{N} problems this function
#'   is for, that default would be enormous, and the per-iteration cost of
#'   \code{nmfkc.signed()} is \eqn{O(QD^2)}).  Ignored when \code{pars} is
#'   supplied.
#' @param seed Optional integer passed to \code{set.seed()} before
#'   generating \eqn{\omega, b}; the caller's random stream is restored
#'   afterwards.  Ignored when \code{pars} is supplied.
#' @param block.size Integer.  Number of columns of \code{U} per block.
#'   Default: as many columns as keep one block of features near 256 MB
#'   (\code{2^25 / D} columns), capped at \code{N}.
#' @param ... Hidden option \code{pars}: a list
#'   \code{list(omega, b, D, beta)} from a previous
#'   \code{\link{nmfkc.signed.rff}} or \code{nmfkc.signed.rff.gram} call.
#'   When supplied, the same random map is reused and \code{beta},
#'   \code{D}, \code{seed} are ignored.
#'
#' @return An object of class \code{"nmfkc.gram"}: a list with
#'   \describe{
#'     \item{\code{S}}{\eqn{D \times D} matrix \eqn{Z Z^\top} (signed,
#'       positive semi-definite).}
#'     \item{\code{G0}}{\eqn{Q_{\mathrm{obs}} \times D} matrix
#'       \eqn{Y Z^\top} (signed).}
#'     \item{\code{N}, \code{D}}{Number of samples and of features.}
#'     \item{\code{pars}}{The RFF parameters \code{list(omega, b, D, beta)}.
#'       Pass them to \code{\link{nmfkc.signed.rff}} to generate features
#'       for new data, e.g. for \code{\link{predict.nmfkc.signed}}.}
#'     \item{\code{block.size}, \code{n.blocks}}{Blocking actually used.}
#'     \item{\code{A.block}}{A function \code{A.block(idx)} returning the
#'       \eqn{D \times |\mathtt{idx}|} feature block for the columns
#'       \code{idx} of \code{U}.  \code{\link{nmfkc.signed}} uses it to
#'       rebuild \eqn{B = CZ} block by block after the fit.  It closes over
#'       \code{U} and \code{pars}, so the object keeps \code{U} alive
#'       (no copy is made).}
#'     \item{\code{rownames}}{Always \code{NULL} (RFF features are
#'       unnamed); present so that \code{nmfkc.signed()} can treat the
#'       object like a matrix.}
#'   }
#'
#' @section Choosing the bandwidth:
#' With Gram input the fold-based \code{\link{nmfkc.signed.cv}} is not
#' available (it must split the columns of \eqn{A}).  Select \eqn{\beta}
#' as in the NMF-LAB paper: build one Gram object per candidate on the
#' training columns, fit, and score the validation columns with
#' \code{predict()} on features from \code{nmfkc.signed.rff(U.val,
#' pars = g$pars)}.
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.  The interface may change in
#' future versions.
#'
#' @references
#' Rahimi, A., & Recht, B. (2007).  Random features for large-scale
#' kernel machines.  \emph{Advances in NIPS}, 20.
#'
#' Satoh, K. (2026).  Applying non-negative matrix factorization with
#' covariates to label matrix for classification.  \emph{Japanese Journal
#' of Statistics and Data Science}.  \doi{10.1007/s42081-026-00349-x}
#'
#' @seealso \code{\link{nmfkc.signed}}, \code{\link{nmfkc.signed.rff}},
#'   \code{\link{predict.nmfkc.signed}},
#'   \code{\link{nmfkc.kernel.beta.nearest.med}}
#'
#' @examples
#' \donttest{
#' ## Iris, 3 classes: the Gram route gives the same fit as the explicit
#' ## feature matrix, without ever holding the D x N matrix.
#' data(iris)
#' set.seed(1)
#' idx <- sample(nrow(iris), 100)
#' mn <- colMeans(iris[idx, 1:4]); sc <- apply(iris[idx, 1:4], 2, sd)
#' U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))   # 4 x 100
#' U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))   # 4 x  50
#' levs    <- levels(iris$Species)
#' Y.train <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
#' rownames(Y.train) <- levs
#' beta <- nmfkc.kernel.beta.nearest.med(U.train)$beta
#'
#' ## Accumulate S, G0 over blocks of 25 columns (4 blocks here)
#' g <- nmfkc.signed.rff.gram(Y.train, U.train, beta = beta, D = 50,
#'                            seed = 1, block.size = 25)
#' g
#' res.g <- nmfkc.signed(Y.train, A = g, rank = 3, verbose = FALSE)
#'
#' ## Same random map, explicit matrix, no warm-start: same fit
#' Z <- nmfkc.signed.rff(U.train, pars = g$pars)$Z              # 50 x 100
#' res.m <- nmfkc.signed(Y.train, A = Z, rank = 3, warm.start = FALSE,
#'                       verbose = FALSE)
#' all.equal(res.g$C, res.m$C)
#'
#' ## Test-set prediction: regenerate features with the stored pars
#' Z.test <- nmfkc.signed.rff(U.test, pars = g$pars)$Z
#' pred   <- predict(res.g, newA = Z.test, type = "class")
#' mean(pred == as.character(iris$Species[-idx]))
#' }
#'
#' @export
nmfkc.signed.rff.gram <- function(Y, U, beta = NULL, D = NULL,
                                  seed = NULL, block.size = NULL, ...) {
  extra <- list(...)
  pars  <- extra$pars

  if (is.vector(Y) && !is.matrix(Y)) Y <- matrix(Y, nrow = 1)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(U)) U <- as.matrix(U)
  N <- ncol(U)
  if (ncol(Y) != N) stop("ncol(Y) must equal ncol(U).")
  if (any(is.na(Y)))
    stop("Y contains NA; a Gram object cannot mask missing entries.")
  if (any(is.na(U))) stop("U contains NA; please impute or remove.")

  if (is.null(pars)) {
    if (is.null(beta)) stop("'beta' must be specified (or supply 'pars' via ...).")
    if (is.null(D))
      stop("'D' (number of random features) must be specified explicitly ",
           "(or supply 'pars' via ...).")
    D <- as.integer(D)
    if (length(D) != 1L || is.na(D) || D < 1L) stop("'D' must be a positive integer.")
  } else {
    D <- as.integer(pars$D)
  }

  if (is.null(block.size)) {
    ## One block of features Z_b (D x block.size doubles) near 256 MB.
    block.size <- max(1L, min(N, as.integer(2^25 %/% D)))
  }
  block.size <- as.integer(block.size)
  if (length(block.size) != 1L || is.na(block.size) || block.size < 1L)
    stop("'block.size' must be a positive integer.")
  block.size <- min(block.size, N)

  Q_obs <- nrow(Y)
  S  <- matrix(0, D, D)
  G0 <- matrix(0, Q_obs, D)
  starts <- seq.int(1L, N, by = block.size)
  for (s in starts) {
    idx <- s:min(s + block.size - 1L, N)
    if (is.null(pars)) {
      ## First block also draws (omega, b); nmfkc.signed.rff() seeds and
      ## restores the caller's stream itself.
      r    <- nmfkc.signed.rff(U[, idx, drop = FALSE], beta = beta, D = D,
                               seed = seed)
      pars <- r$pars
    } else {
      r <- nmfkc.signed.rff(U[, idx, drop = FALSE], pars = pars)
    }
    Zb <- r$Z                                   # D x |idx|, signed
    S  <- S  + tcrossprod(Zb)                   # += Z_b Z_b^T
    G0 <- G0 + tcrossprod(Y[, idx, drop = FALSE], Zb)   # += Y_b Z_b^T
  }

  ## Block generator for nmfkc.signed()'s post-fit B = C Z.  Built in its
  ## own environment so the closure references U and pars only (not Y, S,
  ## G0); R's copy-on-modify means U is shared, not duplicated.
  A.block <- local({
    U <- U; pars <- pars
    function(idx) nmfkc.signed.rff(U[, idx, drop = FALSE], pars = pars)$Z
  })

  structure(list(
    S          = S,
    G0         = G0,
    N          = N,
    D          = D,
    type       = "rff",
    signed     = TRUE,      # RFF features take both signs: nmfkc.signed() only
    pars       = pars,
    beta       = pars$beta,
    block.size = block.size,
    n.blocks   = length(starts),
    A.block    = A.block,
    rownames   = NULL
  ), class = "nmfkc.gram")
}


## print.nmfkc.gram() and the .nmfkc.no.gram() guard are shared with the
## Nystroem constructor and live in R/nmfkc.kernel.gram.R.

# =====================================================
# nmfkc.rff.positive.R -- positive (non-negative) random features for the
#                         Gaussian kernel, for use with nmfkc()
#
# Ordinary Random Fourier Features z(u) = sqrt(2/D) cos(omega'u + b) take
# both signs, so they force the signed solver nmfkc.signed() and break the
# NMF-LAB reading of B = Theta A as memberships.  Positive random features
# (Choromanski et al. 2021, "FAVOR+", stated there for the softmax kernel)
# give an unbiased non-negative feature map for the Gaussian kernel: with
# omega ~ N(0, I_p),
#     phi(u) = exp( sqrt(2 beta) omega'u - 2 beta ||u||^2 ),
#     E[ phi(u) phi(u') ] = exp( -beta ||u - u'||^2 ).
# Proof: E exp(sqrt(2 beta) omega'(u + u')) = exp(beta ||u + u'||^2), and
# -2b||u||^2 - 2b||u'||^2 + b||u||^2 + b||u'||^2 + 2b u'u' = -b||u - u'||^2.
#
# Because phi(u) >= 0, A = Phi works with nmfkc() (standard NMF-LAB, with
# the membership interpretation intact) and with the block-wise Gram route.
# The price is variance: the estimator is heavy-tailed for far-apart points,
# so more features are needed than for the cos map, and inputs should be
# centred / scaled.  A constant kappa = -max_{d,n} L_{dn}, minus the largest
# exponent L_{dn} = sqrt(2 beta) omega_d'u_n - 2 beta ||u_n||^2 observed on
# the training data, is added to every exponent so that every training
# feature is <= 1/sqrt(D) and exp() cannot overflow (far points underflow
# to 0 instead, harmlessly).  It multiplies the kernel estimate by
# exp(2 kappa), a global scale that Theta absorbs.  Because kappa depends
# on the realized omega, the stabilized features are not an unbiased
# estimator of the fixed kernel k but of the randomly rescaled exp(2 kappa) k;
# the unbiasedness / variance statements refer to the unstabilized map.
#
# References:
#   Choromanski, K. et al. (2021). Rethinking attention with Performers.
#     ICLR.  arXiv:2009.14794.
#   Rahimi, A. & Recht, B. (2007). Random features for large-scale kernel
#     machines. NIPS 20.
# =====================================================


#' Positive (non-negative) random features for the Gaussian kernel
#'
#' @description
#' Generates \eqn{\omega_d \sim \mathcal{N}(0, I_p)} and applies the
#' non-negative feature map
#' \deqn{\phi_d(u) = \frac{1}{\sqrt{D}}\exp\!\bigl(\sqrt{2\beta}\,\omega_d^\top u
#'   - 2\beta\lVert u\rVert^2 + \kappa\bigr) \;\ge 0,}
#' for which
#' \eqn{\mathbb{E}[\phi(u)^\top\phi(u')] = e^{2\kappa}\exp(-\beta\lVert u-u'\rVert^2)}:
#' the Gaussian kernel with bandwidth \code{beta}, times a constant.  The
#' \eqn{D \times N} matrix \eqn{\Phi} therefore satisfies
#' \eqn{\Phi^\top\Phi \approx e^{2\kappa} K} \emph{and is non-negative}, so it
#' can be used as the covariate matrix of \code{\link{nmfkc}} directly --
#' unlike the cosine Random Fourier Features of
#' \code{\link{nmfkc.signed.rff}}, which take both signs and need
#' \code{\link{nmfkc.signed}}.  This is the "positive random features"
#' construction of Choromanski et al. (2021), stated there for the softmax
#' kernel, transported to the Gaussian kernel by rescaling.
#'
#' The constant \eqn{\kappa} is minus the largest exponent
#' \eqn{\sqrt{2\beta}\,\omega_d^\top u_n - 2\beta\lVert u_n\rVert^2} observed
#' on the training \code{U} (stored in \code{pars}), so every training
#' feature is at most \eqn{1/\sqrt{D}} and \code{exp()} cannot overflow
#' (points far from the others underflow to zero at large \eqn{\beta}, which
#' only means that bandwidth fits poorly).  It scales the kernel estimate by
#' \eqn{e^{2\kappa}}, which \eqn{\Theta} absorbs, so the fitted model and its
#' predictions do not depend on it.  Note that \eqn{\kappa} depends on the
#' realized \eqn{\omega}: the stabilized features are an unbiased estimator of
#' the randomly rescaled kernel \eqn{e^{2\kappa}k}, not of \eqn{k} itself; the
#' unbiasedness statement above is for the unstabilized map.
#' Inputs should nevertheless be centred and scaled: the
#' estimator's variance grows with \eqn{\beta\lVert u - u'\rVert^2}, so
#' positive features need a larger \code{D} than cosine features for the
#' same accuracy.
#'
#' @param U A \eqn{p \times N} numeric matrix; columns are data points.
#' @param beta Positive scalar.  Gaussian kernel bandwidth parameter
#'   \eqn{k(u,u') = \exp(-\beta\lVert u-u'\rVert^2)}.  Can be obtained via
#'   \code{\link{nmfkc.kernel.beta.nearest.med}}.  Ignored when
#'   \code{pars} is supplied.
#' @param D Integer.  Number of random features (default
#'   \code{ceiling(ncol(U) / 2)}; choose explicitly for large \eqn{N}).
#'   With \code{hyperbolic = TRUE} the \eqn{D} features are \eqn{D/2}
#'   antithetic pairs.  Ignored when \code{pars} is supplied.
#' @param seed Optional integer passed to \code{set.seed()} before drawing
#'   \eqn{\omega}; the caller's random stream is restored.
#' @param ... Hidden options:
#'   \itemize{
#'     \item \code{pars}: a list from a previous call; when supplied the same
#'       random map (and the same \eqn{\kappa}) is reused and \code{beta},
#'       \code{D}, \code{seed} are ignored.  Use this for test data.
#'     \item \code{hyperbolic}: logical (default \code{FALSE}).  If
#'       \code{TRUE}, each \eqn{\omega_d} is paired with \eqn{-\omega_d}
#'       (the "hyperbolic" / antithetic variant of Choromanski et al.).
#'   }
#'
#' @return A list with \code{Z}, the non-negative \eqn{D \times N} feature
#'   matrix, and \code{pars = list(omega, D, beta, kappa, hyperbolic, type =
#'   "positive")}.  Pass \code{Z} to \code{\link{nmfkc}} (or
#'   \code{\link{nmfkc.signed}}) as \code{A}, and \code{pars} back to this
#'   function to build features for new data.
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.  The interface may change in
#' future versions.
#'
#' @references
#' Choromanski, K., Likhosherstov, V., Dohan, D., Song, X., Gane, A.,
#' Sarlos, T., Hawkins, P., Davis, J., Mohiuddin, A., Kaiser, L.,
#' Belanger, D., Colwell, L., & Weller, A. (2021).  Rethinking attention
#' with Performers.  \emph{ICLR}.  arXiv:2009.14794.
#'
#' Rahimi, A., & Recht, B. (2007).  Random features for large-scale kernel
#' machines.  \emph{Advances in NIPS}, 20.
#'
#' @seealso \code{\link{nmfkc.rff.positive.gram}} (large \eqn{N}),
#'   \code{\link{nmfkc.signed.rff}} (signed cosine features),
#'   \code{\link{nmfkc.kernel}}, \code{\link{nmfkc}}
#'
#' @examples
#' set.seed(1)
#' U <- scale(matrix(stats::rnorm(3 * 40), 3, 40))        # centred inputs
#' pr <- nmfkc.rff.positive(U, beta = 0.5, D = 20000, seed = 1)
#' all(pr$Z >= 0)
#' ## Phi' Phi / exp(2 kappa) approximates the Gaussian kernel
#' K  <- nmfkc.kernel(U, beta = 0.5)
#' Kp <- crossprod(pr$Z) / exp(2 * pr$pars$kappa)
#' max(abs(Kp - K))
#'
#' ## The same random map on new data
#' U.new <- matrix(stats::rnorm(3 * 5), 3, 5)
#' dim(nmfkc.rff.positive(U.new, pars = pr$pars)$Z)      # 20000 x 5
#'
#' @export
nmfkc.rff.positive <- function(U, beta = NULL,
                               D = ceiling(ncol(U) / 2),
                               seed = NULL, ...) {
  extra <- list(...)
  pars  <- extra$pars
  hyperbolic <- if (!is.null(extra$hyperbolic)) isTRUE(extra$hyperbolic) else FALSE

  if (!is.matrix(U)) U <- as.matrix(U)
  storage.mode(U) <- "double"
  if (any(is.na(U))) stop("U contains NA; please impute or remove.")
  sq_norm <- colSums(U * U)

  if (is.null(pars)) {
    if (is.null(beta))
      stop("'beta' must be specified (or supply 'pars' via ...).")
    if (!is.numeric(beta) || length(beta) != 1L || beta <= 0)
      stop("'beta' must be a positive scalar.")
    D <- as.integer(D)
    if (length(D) != 1L || is.na(D) || D < 1L) stop("'D' must be a positive integer.")
    .rng <- .nmfkc.rng.save(seed)
    on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    D_draw <- if (hyperbolic) as.integer(ceiling(D / 2)) else D
    omega <- matrix(stats::rnorm(D_draw * nrow(U)), nrow = D_draw, ncol = nrow(U))
    ## Stabilizer: subtract the largest exponent seen on the training data
    ## (Performer-style), so every training feature is <= 1/sqrt(D) and
    ## exp() cannot overflow; far points underflow to 0 instead.  A global
    ## constant, so it only rescales the kernel by exp(2 kappa).
    L0 <- sqrt(2 * beta) * (omega %*% U)
    if (hyperbolic) L0 <- rbind(L0, -L0)
    L0 <- sweep(L0, 2, -2 * beta * sq_norm, "+")
    pars <- list(
      omega      = omega,
      D          = if (hyperbolic) 2L * D_draw else D,
      beta       = beta,
      kappa      = -max(L0),
      hyperbolic = hyperbolic,
      type       = "positive"
    )
  }
  if (ncol(pars$omega) != nrow(U))
    stop("'pars' were generated for p = ", ncol(pars$omega), " input rows, but nrow(U) = ", nrow(U), ".")

  proj  <- sqrt(2 * pars$beta) * (pars$omega %*% U)            # D_draw x N
  shift <- -2 * pars$beta * sq_norm + pars$kappa                # length N
  if (isTRUE(pars$hyperbolic)) {
    Z <- rbind(exp(sweep( proj, 2, shift, "+")),
               exp(sweep(-proj, 2, shift, "+"))) / sqrt(pars$D)
  } else {
    Z <- exp(sweep(proj, 2, shift, "+")) / sqrt(pars$D)
  }
  list(Z = Z, pars = pars)
}


#' Block-wise Gram accumulation of positive random features for large N
#'
#' @description
#' The \code{\link{nmfkc.rff.positive}} counterpart of
#' \code{\link{nmfkc.signed.rff.gram}}: walks the columns of \code{U} in
#' blocks, generates the non-negative feature block
#' \eqn{\Phi_b} with the same random map, accumulates
#' \eqn{S = \Phi\Phi^\top} (\eqn{D \times D}) and \eqn{G_0 = Y\Phi^\top}
#' (\eqn{P \times D}), and discards the block, so the \eqn{D \times N}
#' feature matrix never exists in memory.  Because the features are
#' non-negative, the returned \code{"nmfkc.gram"} object is accepted by
#' \code{\link{nmfkc}} (as well as by \code{\link{nmfkc.signed}}), and the
#' NMF-LAB membership interpretation of \eqn{B = \Theta\Phi} is retained.
#'
#' @param Y Non-negative \eqn{P \times N} response matrix (e.g. one-hot
#'   labels).  \code{NA} is not allowed.
#' @param U A \eqn{p \times N} numeric matrix of inputs; columns are samples.
#' @param beta Positive scalar.  Gaussian kernel bandwidth.  Ignored when
#'   \code{pars} is supplied.
#' @param D Integer.  Number of random features.  Required (no \eqn{N/2}
#'   default here).  Ignored when \code{pars} is supplied.
#' @param seed Optional integer seed for \eqn{\omega}; the caller's stream
#'   is restored.  Ignored when \code{pars} is supplied.
#' @param block.size Integer.  Columns of \code{U} per block.  Default:
#'   about 256 MB of features per block (\code{2^25 / D} columns), capped
#'   at \code{N}.
#' @param ... Hidden options \code{pars} (reuse a random map from
#'   \code{\link{nmfkc.rff.positive}} or a previous call) and
#'   \code{hyperbolic} (see \code{\link{nmfkc.rff.positive}}).
#'
#' @return An object of class \code{"nmfkc.gram"} with \code{S}, \code{G0},
#'   \code{N}, \code{D}, \code{type = "prf"}, \code{signed = FALSE},
#'   \code{pars}, \code{beta}, \code{block.size}, \code{n.blocks},
#'   \code{A.block} (block generator) and \code{rownames = NULL}.  Features
#'   for new data: \code{nmfkc.rff.positive(U.new, pars = g$pars)$Z}.
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.  The interface may change in
#' future versions.
#'
#' @references
#' Choromanski, K. et al. (2021).  Rethinking attention with Performers.
#' \emph{ICLR}.  arXiv:2009.14794.
#'
#' @seealso \code{\link{nmfkc.rff.positive}}, \code{\link{nmfkc}},
#'   \code{\link{nmfkc.kernel.gram}} (Nystr\"om), \code{\link{nmfkc.signed.rff.gram}}
#'
#' @examples
#' \donttest{
#' data(iris)
#' set.seed(1)
#' idx <- sample(nrow(iris), 100)
#' mn <- colMeans(iris[idx, 1:4]); sc <- apply(iris[idx, 1:4], 2, sd)
#' U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))
#' U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))
#' levs    <- levels(iris$Species)
#' Y.train <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
#' rownames(Y.train) <- levs
#' beta <- nmfkc.kernel.beta.nearest.med(U.train)$beta
#'
#' g <- nmfkc.rff.positive.gram(Y.train, U.train, beta = beta, D = 200,
#'                              seed = 1, block.size = 25)
#' g
#' res <- nmfkc(Y.train, A = g, rank = 3, verbose = FALSE)   # standard NMF-LAB
#' Z.test <- nmfkc.rff.positive(U.test, pars = g$pars)$Z
#' pred   <- predict(res, newA = Z.test, type = "class")
#' mean(pred == as.character(iris$Species[-idx]))
#' }
#'
#' @export
nmfkc.rff.positive.gram <- function(Y, U, beta = NULL, D = NULL,
                                    seed = NULL, block.size = NULL, ...) {
  extra <- list(...)
  pars  <- extra$pars
  hyperbolic <- if (!is.null(extra$hyperbolic)) isTRUE(extra$hyperbolic) else FALSE

  if (is.vector(Y) && !is.matrix(Y)) Y <- matrix(Y, nrow = 1)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(U)) U <- as.matrix(U)
  storage.mode(U) <- "double"
  N <- ncol(U)
  if (ncol(Y) != N) stop("ncol(Y) must equal ncol(U).")
  if (any(is.na(Y))) stop("Y contains NA; a Gram object cannot mask missing entries.")
  if (min(Y) < 0) stop("The matrix Y should be non-negative.")
  if (any(is.na(U))) stop("U contains NA; please impute or remove.")

  if (is.null(pars)) {
    if (is.null(beta)) stop("'beta' must be specified (or supply 'pars' via ...).")
    if (is.null(D))
      stop("'D' (number of random features) must be specified explicitly ",
           "(or supply 'pars' via ...).")
    D <- as.integer(D)
    if (length(D) != 1L || is.na(D) || D < 1L) stop("'D' must be a positive integer.")
    ## Draw omega once; kappa (minus the largest exponent over all training
    ## columns) is filled in by a first pass over the blocks below.
    .rng <- .nmfkc.rng.save(seed)
    on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    D_draw <- if (hyperbolic) as.integer(ceiling(D / 2)) else D
    pars <- list(
      omega      = matrix(stats::rnorm(D_draw * nrow(U)), nrow = D_draw, ncol = nrow(U)),
      D          = if (hyperbolic) 2L * D_draw else D,
      beta       = beta,
      kappa      = NA_real_,      # filled below from a first pass over the blocks
      hyperbolic = hyperbolic,
      type       = "positive"
    )
  }
  D <- as.integer(pars$D)

  if (is.null(block.size)) block.size <- max(1L, min(N, as.integer(2^25 %/% D)))
  block.size <- as.integer(block.size)
  if (length(block.size) != 1L || is.na(block.size) || block.size < 1L)
    stop("'block.size' must be a positive integer.")
  block.size <- min(block.size, N)
  starts <- seq.int(1L, N, by = block.size)

  if (is.na(pars$kappa)) {
    ## First pass: the largest exponent over all training columns -- the same
    ## stabilizer nmfkc.rff.positive() uses -- without holding the features.
    mx <- -Inf
    for (s in starts) {
      idx <- s:min(s + block.size - 1L, N)
      Ub <- U[, idx, drop = FALSE]
      L0 <- sqrt(2 * pars$beta) * (pars$omega %*% Ub)
      if (isTRUE(pars$hyperbolic)) L0 <- rbind(L0, -L0)
      L0 <- sweep(L0, 2, -2 * pars$beta * colSums(Ub * Ub), "+")
      mx <- max(mx, max(L0))
    }
    pars$kappa <- -mx
  }

  P  <- nrow(Y)
  S  <- matrix(0, D, D)
  G0 <- matrix(0, P, D)
  for (s in starts) {
    idx <- s:min(s + block.size - 1L, N)
    Zb <- nmfkc.rff.positive(U[, idx, drop = FALSE], pars = pars)$Z   # D x |idx|, >= 0
    S  <- S  + tcrossprod(Zb)
    G0 <- G0 + tcrossprod(Y[, idx, drop = FALSE], Zb)
  }

  A.block <- local({
    U <- U; pars <- pars
    function(idx) nmfkc.rff.positive(U[, idx, drop = FALSE], pars = pars)$Z
  })

  structure(list(
    S          = S,
    G0         = G0,
    N          = N,
    D          = D,
    type       = "prf",
    signed     = FALSE,
    pars       = pars,
    beta       = pars$beta,
    block.size = block.size,
    n.blocks   = length(starts),
    A.block    = A.block,
    rownames   = NULL
  ), class = "nmfkc.gram")
}

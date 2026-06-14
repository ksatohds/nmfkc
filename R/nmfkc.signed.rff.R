# =====================================================
# nmfkc.signed.rff.R — Random Fourier Features generator
#                      for use with nmfkc.signed()
#
# Given training covariates U (p x N), generates the signed RFF
# feature matrix
#     Z = sqrt(2 / D) * cos(omega U + b),       D x N,   Z^T Z ~ K,
# where omega_d ~ N(0, 2 beta I_p) and b_d ~ Uniform(0, 2 pi)
# (Rahimi & Recht, 2007).  The result is a list containing the
# sign-unrestricted feature matrix Z and the generating parameters
# pars = list(omega, b, D, beta), so that the same random map can be
# re-applied to new data and nmfkc.signed() can store pars for
# downstream summary.
#
# Design: kept intentionally minimal.  The kernel-faithful direct MU
# algorithm and all S3 methods (predict, plot, summary) live in
# nmfkc.signed.R; this file only provides the feature map.
#
# References:
#   Rahimi, A. & Recht, B. (2007). Random features for large-scale
#     kernel machines. NIPS 20.
# =====================================================


#' Random Fourier Features for nmfkc.signed()
#'
#' @description
#' Generates RFF random parameters
#' \eqn{\omega_d \sim \mathcal{N}(0, 2\beta I_p)},
#' \eqn{b_d \sim \mathrm{Uniform}(0, 2\pi)} (Rahimi & Recht, 2007) and
#' applies the RFF transform
#' \deqn{z_d(u) = \sqrt{2/D}\, \cos(\omega_d^\top u + b_d)}
#' to each column of \code{U}, yielding a sign-unrestricted
#' \eqn{D \times N} feature matrix \eqn{Z} such that
#' \eqn{Z^\top Z \approx K}, the Gaussian kernel matrix with bandwidth
#' \code{beta}.
#'
#' The return value is a list with the feature matrix \code{Z} and
#' the generating parameters \code{pars = list(omega, b, D, beta)} so
#' that the same random map can be re-applied to new data (by passing
#' \code{pars} back) and \code{\link{nmfkc.signed}} can record the
#' parameters for downstream summary.
#'
#' @param U A \eqn{p \times N} numeric matrix; columns are data points.
#' @param beta Positive scalar. Gaussian kernel bandwidth parameter.
#'   Can be obtained via \code{\link{nmfkc.kernel.beta.nearest.med}}.
#'   May be \code{NULL} only when \code{pars} is supplied through
#'   \code{...}.
#' @param D Integer. Number of random features. Defaults to
#'   \code{ceiling(ncol(U) / 2)}.  This default is intended for the
#'   \strong{training-time fresh generation} only; for test data, always
#'   supply \code{pars} via \code{...} to inherit the training \eqn{D}
#'   together with \eqn{\omega, b}.  For very large \eqn{N} the default
#'   may be excessive (direct MU cost is \eqn{O(QD^2)}); choose a
#'   smaller \eqn{D} manually.  For very small \eqn{N}, RFF is not
#'   recommended; use a full kernel matrix with \code{\link{nmfkc}}
#'   instead.
#' @param seed Optional integer passed to \code{set.seed()} before
#'   generating \eqn{\omega, b}, for reproducibility.  Ignored when
#'   \code{pars} is supplied.
#' @param ... Hidden option \code{pars}: a list
#'   \code{list(omega, b, D, beta)} obtained from a previous call
#'   (\code{Ztrain$pars}).  When supplied, \eqn{\omega, b} are reused
#'   and the \code{beta}, \code{D}, \code{seed} arguments are ignored.
#'   Use this to apply the same random map to test data.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{Z}}{A \eqn{D \times N} sign-unrestricted numeric
#'       matrix.  Pass this to \code{\link{nmfkc.signed}} as its
#'       \code{A} argument.}
#'     \item{\code{pars}}{A list \code{list(omega, b, D, beta)}.  Pass
#'       this to \code{\link{nmfkc.signed}} via its \code{pars}
#'       argument (for summary display) and to subsequent
#'       \code{nmfkc.signed.rff()} calls (to reuse the same random map
#'       on new data).}
#'   }
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @seealso \code{\link{nmfkc.signed}},
#'   \code{\link{nmfkc.kernel.beta.nearest.med}}
#'
#' @examples
#' \donttest{
#' ## Iris 3-class classification with RFF + direct MU (Ding-Li-Jordan)
#' data(iris)
#' set.seed(1)
#' idx <- sample(nrow(iris), 100)    # 100 training, 50 test
#'
#' ## Scale features using TRAINING mean/sd; transpose to p x N layout
#' mn      <- colMeans(iris[idx, 1:4])
#' sc      <- apply(iris[idx, 1:4], 2, sd)
#' U.train <- t(scale(iris[idx,  1:4], center = mn, scale = sc))   # 4 x 100
#' U.test  <- t(scale(iris[-idx, 1:4], center = mn, scale = sc))   # 4 x  50
#'
#' ## One-hot encode training labels as a Q_obs x N target matrix
#' levs      <- levels(iris$Species)
#' Y.train   <- sapply(iris$Species[idx], function(s) as.integer(levs == s))
#' rownames(Y.train) <- levs                                       # 3 x 100
#' lab.train <- iris$Species[idx]
#' lab.test  <- iris$Species[-idx]
#'
#' ## Beta candidates from nearest-neighbour median heuristic
#' beta_info <- nmfkc.kernel.beta.nearest.med(U.train)
#' betas     <- beta_info$beta_candidates
#'
#' ## CV over beta candidates: for each beta, generate RFF, fit, and
#' ## evaluate column-wise CV-MSE on training data
#' cv_mse <- numeric(length(betas))
#' for (i in seq_along(betas)) {
#'   rff_i   <- nmfkc.signed.rff(U.train, beta = betas[i], D = 50, seed = 1)
#'   cv_i    <- nmfkc.signed.cv(Y.train, A = rff_i$Z, rank = 3, seed = 123)
#'   cv_mse[i] <- cv_i$objfunc
#' }
#' beta_best <- betas[which.min(cv_mse)]
#'
#' ## Generate signed RFF features with the best beta
#' rff.train <- nmfkc.signed.rff(U.train, beta = beta_best, D = 50, seed = 1)
#' rff.test  <- nmfkc.signed.rff(U.test, pars = rff.train$pars)
#'
#' ## Fit on training data only
#' res <- nmfkc.signed(Y.train, A = rff.train$Z, rank = 3,
#'                     pars = rff.train$pars, verbose = FALSE)
#'
#' ## Predict on training and test separately
#' pred.train <- predict(res, newA = rff.train$Z, type = "class")
#' pred.test  <- predict(res, newA = rff.test$Z,  type = "class")
#' mean(pred.train == as.character(lab.train))
#' mean(pred.test  == as.character(lab.test))
#' }
#'
#' @export
nmfkc.signed.rff <- function(U, beta = NULL,
                             D = ceiling(ncol(U) / 2),
                             seed = NULL, ...) {
  extra <- list(...)
  pars  <- extra$pars

  if (is.null(pars)) {
    if (is.null(beta))
      stop("'beta' must be specified (or supply 'pars' via ...).")
    if (!is.null(seed)) set.seed(seed)
    p <- nrow(U)
    pars <- list(
      omega = matrix(stats::rnorm(D * p, mean = 0, sd = sqrt(2 * beta)),
                     nrow = D, ncol = p),
      b     = stats::runif(D, 0, 2 * pi),
      D     = D,
      beta  = beta
    )
  }

  Z <- sqrt(2 / pars$D) * cos(pars$omega %*% U + pars$b)
  list(Z = Z, pars = pars)
}

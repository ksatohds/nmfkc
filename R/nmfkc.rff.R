# =====================================================
# nmfkc.rff.R — Kernel-faithful NMF with Random Fourier Features
#
# Implements the Direct Multiplicative Update (Direct MU) algorithm for
#   Y = X H Z + E,   Y, X, H = Hp - Hn, Z = Zp - Zn (RFF features),
# where Z^T Z approximates a Gaussian kernel K (Rahimi & Recht 2007).
#
# References:
#   Rahimi, A. & Recht, B. (2007). Random features for large-scale kernel
#     machines. NIPS.
# =====================================================


############################################################
## nmfkc.rff.random : generate & apply RFF in one call
############################################################

#' Compute Random Fourier Features of input data
#'
#' @description
#' Generates RFF random parameters
#' \eqn{\omega_d \sim \mathcal{N}(0, 2\beta I_p)},
#' \eqn{b_d \sim \mathrm{Uniform}(0, 2\pi)} (Rahimi & Recht, 2007) and
#' applies the RFF transform
#' \deqn{z_d(u) = \sqrt{2/D}\, \cos(\omega_d^\top u + b_d)}
#' to each column of \code{U}, yielding a \eqn{D \times N} feature matrix
#' \eqn{Z} such that \eqn{Z^\top Z \approx K} (the Gaussian kernel matrix
#' with bandwidth \code{beta}).
#'
#' By default, the non-negative posneg split
#' \eqn{(Z_{+}, Z_{-})} is returned, suitable for direct use with
#' \code{\link{nmfkc.rff}}.  The generating parameters (\code{omega},
#' \code{b}) are returned alongside so that the same random map can be
#' reapplied to new data (e.g., a test set).
#'
#' @param U A \eqn{p \times N} numeric matrix; columns are data points.
#' @param beta Positive scalar.  Gaussian kernel bandwidth parameter.
#'   Can be obtained via \code{\link{nmfkc.kernel.beta.nearest.med}}.
#'   May be \code{NULL} only when \code{pars} is supplied through
#'   \code{...}.
#' @param D Integer. Number of random features. Defaults to
#'   \code{ceiling(ncol(U) / 2)}.  This default is intended for the
#'   \strong{training-time fresh generation} only; for test data, always
#'   supply \code{pars} via \code{...} to inherit the training \eqn{D}
#'   together with \eqn{\omega, b}.  For very large \eqn{N} the default
#'   may be excessive (iteration cost is \eqn{O(QD^2)}); choose a
#'   smaller \eqn{D} manually.  For very small \eqn{N}, RFF is not
#'   recommended; use a full kernel matrix with \code{\link{nmfkc}}
#'   instead.
#' @param seed Optional integer passed to \code{set.seed()} before
#'   generating \eqn{\omega, b}, for reproducibility.  Ignored when
#'   \code{pars} is supplied.
#' @param nonneg Logical. If \code{TRUE} (default), returns the posneg
#'   split as \code{list(Zp = max(Z, 0), Zn = max(-Z, 0), pars = ...)}.
#'   If \code{FALSE}, returns the raw sign-unrestricted \eqn{D \times N}
#'   matrix \eqn{Z} with the generating \code{pars} attached as
#'   attribute.
#' @param ... Hidden option \code{pars}: a list
#'   \code{list(omega, b, D, beta)} obtained from a previous call.
#'   When supplied, \eqn{\omega, b} are reused and \code{beta}/\code{D}/
#'   \code{seed} arguments are ignored.  Use this to apply the same
#'   random map to test data.
#'
#' @return When \code{nonneg = TRUE}: a list with elements \code{Zp},
#'   \code{Zn} (non-negative matrices, each \eqn{D \times N}) and
#'   \code{pars}.  When \code{nonneg = FALSE}: a \eqn{D \times N}
#'   sign-unrestricted matrix with \code{attr(., "pars")} attached.
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @seealso \code{\link{nmfkc.rff}},
#'   \code{\link{nmfkc.kernel.beta.nearest.med}}
#'
#' @examples
#' set.seed(1)
#' U.train <- matrix(stats::rnorm(5 * 40), 5, 40)
#' U.test  <- matrix(stats::rnorm(5 * 10), 5, 10)
#'
#' ## Training: only beta required; D defaults to ncol(U)/2
#' Ztr <- nmfkc.rff.random(U.train, beta = 0.5, seed = 1)
#' dim(Ztr$Zp)              # 20 x 40 (D = ceiling(40/2) = 20)
#'
#' ## Test: reuse the same omega/b via pars
#' Zte <- nmfkc.rff.random(U.test, pars = Ztr$pars)
#' dim(Zte$Zp)              # 20 x 10
#'
#' @export
nmfkc.rff.random <- function(U, beta = NULL,
                             D = ceiling(ncol(U) / 2),
                             seed = NULL,
                             nonneg = TRUE, ...) {
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

  proj <- pars$omega %*% U + pars$b
  Z <- sqrt(2 / pars$D) * cos(proj)

  if (isTRUE(nonneg)) {
    list(Zp = pmax(Z, 0), Zn = pmax(-Z, 0), pars = pars)
  } else {
    attr(Z, "pars") <- pars
    Z
  }
}


############################################################
## nmfkc.rff : Direct MU algorithm (Algorithm 2)
############################################################

#' Kernel-faithful NMF with Random Fourier Features
#'
#' @description
#' Solves the kernel-faithful NMF-KC problem
#' \deqn{Y \approx X\,(H_{+} - H_{-})\,(Z_{+} - Z_{-}),
#'   \quad X, H_{+}, H_{-} \ge 0,}
#' by a direct multiplicative-update (MU) algorithm that avoids the
#' \eqn{N \times N} kernel matrix.  The implicit kernel is the true RFF
#' approximation \eqn{Z^\top Z} (not the posneg-split surrogate
#' \eqn{Z_{+}^\top Z_{+} + Z_{-}^\top Z_{-}}), and the objective
#' \eqn{\|Y - X H Z\|_F^2} decreases monotonically per update.
#'
#' All iteration-time matrix operations use only the precomputed
#' \eqn{D \times D} Gram matrix \eqn{S = Z Z^\top} and
#' \eqn{G_0 = Y Z^\top} (\eqn{Q_{\mathrm{obs}} \times D}), so the cost per
#' iteration is \eqn{O(QD^2)} and independent of \eqn{N}.
#'
#' @param Y Non-negative \eqn{Q_{\mathrm{obs}} \times N} response matrix.
#' @param Zp Non-negative \eqn{D \times N} matrix \eqn{Z_{+} = \max(Z, 0)},
#'   where \eqn{Z} is produced by \code{\link{nmfkc.rff.random}}.
#' @param Zn Non-negative \eqn{D \times N} matrix \eqn{Z_{-} = \max(-Z, 0)}.
#' @param rank Integer. Number of latent components \eqn{Q} in \eqn{X}.
#' @param warm.start Logical. If \code{TRUE} (default), internally runs the
#'   posneg-split NMF via \code{\link{nmfkc}} on \code{Y, A = rbind(Zp, Zn)}
#'   and uses it to seed \eqn{X, H_{+}, H_{-}}.  This usually accelerates
#'   convergence.  If \code{FALSE}, \eqn{X, H_{+}, H_{-}} are initialized
#'   from Gaussian random draws (used in the research memo for comparison).
#'   Ignored when both \code{X.init} and \code{C.init} are supplied via
#'   \code{...}.
#' @param maxit Maximum number of iterations (default 5000).
#' @param epsilon Relative convergence tolerance on the objective
#'   (default 1e-4).
#' @param print.trace Logical. If \code{TRUE}, prints summary messages at
#'   start and end.
#' @param ... Additional arguments:
#'   \itemize{
#'     \item \code{Q}: alias for \code{rank} (backward compatibility).
#'     \item \code{pars}: the \code{list(omega, b, D, beta)} obtained from
#'           \code{\link{nmfkc.rff.random}}.  When supplied, the RFF
#'           generating parameters are stored in the returned object
#'           (as \code{$pars}) so that \code{summary()} can report
#'           \code{beta} and the random map can be re-applied to new data.
#'     \item \code{X.init} (\eqn{Q_{\mathrm{obs}} \times Q}): explicit
#'           initial basis matrix.
#'     \item \code{C.init} (\eqn{Q \times 2D}): explicit initial
#'           coefficient matrix (first \eqn{D} columns are \eqn{H_{+}},
#'           last \eqn{D} are \eqn{H_{-}}).  Same layout as
#'           \code{$C} returned by \code{\link{nmfkc}} on
#'           \code{Y, A = rbind(Zp, Zn)}.
#'   }
#'   When both \code{X.init} and \code{C.init} are supplied, they are
#'   used directly; \code{warm.start} has no effect.  When only one is
#'   supplied, the missing one is filled by the posneg warm-start
#'   (mirrors \code{\link{nmfre}}'s behavior with \code{X.init} / \code{C.init}).
#'
#' @return An object of class \code{c("nmfkc.rff", "nmfkc")} with elements
#' \itemize{
#'   \item \code{X}: \eqn{Q_{\mathrm{obs}} \times Q} basis matrix (non-negative).
#'   \item \code{Hp}, \code{Hn}: \eqn{Q \times D} non-negative parts of \eqn{H}.
#'   \item \code{C}: \code{cbind(Hp, Hn)} for compatibility with
#'     \code{nmfkc}-style methods.
#'   \item \code{objfunc.iter}: objective-function values per iteration.
#'   \item \code{objfunc}: final objective value.
#'   \item \code{r.squared}: \eqn{\mathrm{cor}(Y, \widehat Y)^2} of the full
#'     reconstruction \eqn{\widehat Y = X(H_{+} - H_{-})(Z_{+} - Z_{-})}.
#'   \item \code{mae}: mean absolute error of the reconstruction.
#'   \item \code{iter}: number of iterations performed.
#'   \item \code{runtime}: elapsed seconds.
#'   \item \code{call}: the matched call.
#' }
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @references
#' Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
#' machines. \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Ding, C., Li, T., & Jordan, M. I. (2010). Convex and semi-nonnegative
#' matrix factorizations. \emph{IEEE TPAMI}, 32(1), 45-55.
#'
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.rff.random}},
#'   \code{\link{nmfkc.kernel.beta.nearest.med}},
#'   \code{\link{predict.nmfkc.rff}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' U <- matrix(stats::rnorm(5 * 40), 5, 40)
#' Y <- matrix(abs(stats::rnorm(8 * 40)), 8, 40)
#' Z <- nmfkc.rff.random(U, beta = 0.5, seed = 1)
#'
#' ## Default: warm-start via posneg nmfkc() internally
#' res <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = 3, maxit = 200)
#' plot(res); summary(res)
#'
#' ## Random initialization (for research reproduction)
#' res_rand <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = 3, maxit = 200,
#'                        warm.start = FALSE)
#' }
#'
#' @export
nmfkc.rff <- function(Y, Zp, Zn,
                       rank = NULL,
                       warm.start = TRUE,
                       maxit = 5000,
                       epsilon = 1e-4,
                       print.trace = FALSE,
                       ...) {
  cl <- match.call()
  extra_args <- list(...)
  if (is.null(rank) && !is.null(extra_args$Q)) rank <- extra_args$Q
  if (is.null(rank)) stop("'rank' must be specified.")
  Q <- as.integer(rank)

  ## Hidden init options (consistent with nmfre): X.init (Q_obs x Q),
  ## C.init (Q x 2D).  When supplied, they override warm.start.
  X.init <- extra_args$X.init
  C.init <- extra_args$C.init

  ## Hidden: RFF generating parameters (for summary display and
  ## downstream reuse).  If supplied, stored in the returned object.
  pars_rff <- extra_args$pars

  Q_obs <- nrow(Y)
  N     <- ncol(Y)
  D_rff <- nrow(Zp)
  if (!identical(dim(Zp), dim(Zn)))
    stop("'Zp' and 'Zn' must have the same dimensions.")
  if (ncol(Zp) != N)
    stop("ncol(Zp) must equal ncol(Y).")

  t0 <- proc.time()

  ## ---- precomputation (avoids N x N matrices) ----
  Z <- Zp - Zn
  S <- tcrossprod(Z)               # D x D
  S_p <- pmax(S, 0); S_n <- pmax(-S, 0)
  G0 <- Y %*% t(Z)                 # Q_obs x D
  Y_sqnorm <- sum(Y * Y)

  ## ---- initialization ----
  ## If X.init / C.init are supplied, use them directly.
  ## Else if warm.start = TRUE, run posneg NMF to seed the MU iterations.
  ## Else use Gaussian random initialization.
  need_warm <- isTRUE(warm.start) &&
               (is.null(X.init) || is.null(C.init))
  if (need_warm) {
    if (isTRUE(print.trace))
      cat("nmfkc.rff: warm-starting via posneg nmfkc(..) ...\n")
    res0 <- nmfkc(Y, A = rbind(Zp, Zn), rank = Q,
                   maxit = maxit, epsilon = epsilon,
                   print.dims = FALSE)
    if (is.null(X.init)) X.init <- res0$X
    if (is.null(C.init)) C.init <- res0$C
  }

  if (!is.null(X.init) && !is.null(C.init)) {
    if (!identical(dim(X.init), c(Q_obs, Q)))
      stop("X.init must have dimensions (nrow(Y), rank).")
    if (!identical(dim(C.init), c(as.integer(Q), 2L * as.integer(D_rff))))
      stop("C.init must have dimensions (rank, 2 * nrow(Zp)).")
    X  <- X.init
    Hp <- C.init[, 1:D_rff, drop = FALSE]
    Hn <- C.init[, (D_rff + 1):(2 * D_rff), drop = FALSE]
  } else {
    ## warm.start = FALSE and no explicit init supplied -> random
    X  <- matrix(abs(stats::rnorm(Q_obs * Q)) * 0.1, Q_obs, Q)
    Hp <- matrix(abs(stats::rnorm(Q * D_rff)) * 0.1, Q, D_rff)
    Hn <- matrix(abs(stats::rnorm(Q * D_rff)) * 0.01, Q, D_rff)
  }

  small <- 1e-16
  P <- crossprod(X); G <- crossprod(X, G0)
  H <- Hp - Hn; PH <- P %*% H
  obj_prev <- Y_sqnorm - 2 * sum(G * H) + sum(H * (PH %*% S))
  objfunc.iter <- numeric(maxit)

  if (isTRUE(print.trace))
    cat(sprintf("nmfkc.rff: Y(%d,%d) Z(%d,%d) rank=%d\n",
                Q_obs, N, D_rff, N, Q))

  ## ---- iterations ----
  iter <- 0L
  for (iter in seq_len(maxit)) {
    G_p <- pmax(G, 0); G_n <- pmax(-G, 0)
    PHp <- P %*% Hp;   PHn <- P %*% Hn

    ## Hp update (kernel-faithful, true L gradient)
    Hp <- Hp * (G_p + PHp %*% S_n + PHn %*% S_p) /
               (G_n + PHp %*% S_p + PHn %*% S_n + small)

    ## Hn update (Gauss-Seidel: uses updated Hp)
    PHp <- P %*% Hp
    Hn <- Hn * (G_n + PHp %*% S_p + PHn %*% S_n) /
               (G_p + PHp %*% S_n + PHn %*% S_p + small)

    ## X update (avoids N x N via YMt = G0 %*% (Hp-Hn)^T, MMt = H S H^T)
    H   <- Hp - Hn; Ht <- t(H)
    YMt <- G0 %*% Ht
    HS  <- H %*% S; MMt <- HS %*% Ht

    X <- X * (pmax(YMt, 0) + X %*% pmax(-MMt, 0)) /
             (pmax(-YMt, 0) + X %*% pmax(MMt, 0) + small)

    ## refresh dependents & objective
    P  <- crossprod(X); G <- crossprod(X, G0)
    PH <- P %*% H
    obj_cur <- Y_sqnorm - 2 * sum(G * H) + sum(H * (PH %*% S))
    objfunc.iter[iter] <- obj_cur

    if (abs(obj_prev - obj_cur) / max(abs(obj_prev), 1e-12) < epsilon) break
    obj_prev <- obj_cur
  }
  objfunc.iter <- objfunc.iter[seq_len(iter)]

  ## ---- final reconstruction stats (kernel-faithful) ----
  Yhat <- X %*% (Hp %*% Zp + Hn %*% Zn - Hp %*% Zn - Hn %*% Zp)
  r.squared <- tryCatch(
    stats::cor(as.vector(Yhat), as.vector(Y))^2,
    error = function(e) NA_real_
  )
  mae <- mean(abs(Y - Yhat))

  runtime <- as.numeric((proc.time() - t0)[3])

  if (isTRUE(print.trace))
    cat(sprintf("nmfkc.rff: done in %d iterations, %.2fs, R2=%.4f\n",
                iter, runtime, r.squared))

  result <- list(
    X             = X,
    Hp            = Hp,
    Hn            = Hn,
    C             = cbind(Hp, Hn),
    objfunc.iter  = objfunc.iter,
    objfunc       = obj_cur,
    r.squared     = r.squared,
    mae           = mae,
    iter          = iter,
    runtime       = runtime,
    rank          = Q,
    D             = D_rff,
    pars          = pars_rff,   # RFF generating parameters (may be NULL)
    call          = cl
  )
  class(result) <- c("nmfkc.rff", "nmfkc")
  result
}


############################################################
## predict.nmfkc.rff
############################################################

#' Predict method for nmfkc.rff models
#'
#' @description
#' Computes kernel-faithful predictions
#' \deqn{\widehat Y = X\,(H_{+} Z_{+} + H_{-} Z_{-}
#'                        - H_{+} Z_{-} - H_{-} Z_{+})}
#' for new RFF features \code{newZp}, \code{newZn}.  Negative predictions are
#' clipped to zero before being returned or converted to class probabilities.
#'
#' @param object A fitted model of class \code{"nmfkc.rff"}.
#' @param newZp Non-negative \eqn{D \times N_{\mathrm{new}}} matrix
#'   \eqn{Z_{+}^{\mathrm{new}}}.  If \code{NULL}, the training prediction
#'   cannot be returned (training Zp/Zn is not stored by design); supply
#'   both \code{newZp} and \code{newZn} here.
#' @param newZn Non-negative \eqn{D \times N_{\mathrm{new}}} matrix
#'   \eqn{Z_{-}^{\mathrm{new}}}.
#' @param type Output type: \code{"response"} (default), \code{"prob"}
#'   (column-normalized to sum to 1), or \code{"class"} (row name with
#'   maximum probability).
#' @param ... Unused.
#'
#' @return A numeric matrix (\code{"response"} or \code{"prob"}) or a
#'   character vector (\code{"class"}).
#'
#' @seealso \code{\link{nmfkc.rff}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' U <- matrix(stats::rnorm(5 * 40), 5, 40)
#' Y <- matrix(abs(stats::rnorm(8 * 40)), 8, 40)
#' Z <- nmfkc.rff.random(U, beta = 0.5, seed = 1)
#' res <- nmfkc.rff(Y, Z$Zp, Z$Zn, rank = 3, maxit = 100)
#' Yhat <- predict(res, newZp = Z$Zp, newZn = Z$Zn)
#' dim(Yhat)
#' }
#'
#' @export
predict.nmfkc.rff <- function(object, newZp = NULL, newZn = NULL,
                              type = c("response", "prob", "class"),
                              ...) {
  type <- match.arg(type)
  if (is.null(newZp) || is.null(newZn))
    stop("Both 'newZp' and 'newZn' must be supplied.")

  Hp <- object$Hp; Hn <- object$Hn; X <- object$X
  Yhat <- X %*% (Hp %*% newZp + Hn %*% newZn
                 - Hp %*% newZn - Hn %*% newZp)
  Yhat <- pmax(Yhat, 0)
  if (type == "response") return(Yhat)

  col_sums <- colSums(Yhat) + 1e-12
  probs <- sweep(Yhat, 2, col_sums, "/")
  if (type == "prob") return(probs)

  classes <- rownames(Yhat)
  if (is.null(classes)) classes <- as.character(seq_len(nrow(Yhat)))
  classes[apply(probs, 2, which.max)]
}


############################################################
## plot.nmfkc.rff
############################################################

#' Plot method for nmfkc.rff objects (convergence)
#'
#' @description
#' Draws the objective-function trajectory for a fitted
#' \code{\link{nmfkc.rff}} model.
#'
#' @param x An \code{nmfkc.rff} object.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return Called for its side effect (plot). Returns \code{x} invisibly.
#' @seealso \code{\link{nmfkc.rff}}
#' @export
plot.nmfkc.rff <- function(x, ...) {
  extra_args <- list(...)
  args <- list(x = x$objfunc.iter, type = "l")
  if (is.null(extra_args$main))
    args$main <- sprintf("r.squared = %.3f", x$r.squared)
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  all_args <- c(args, extra_args)
  do.call(graphics::plot, all_args)
  invisible(x)
}


############################################################
## summary.nmfkc.rff  +  print.summary.nmfkc.rff
############################################################

#' Summary method for nmfkc.rff objects
#'
#' @description
#' Produces a concise summary of a \code{\link{nmfkc.rff}} fit.
#'
#' @param object A fitted \code{nmfkc.rff} object.
#' @param ... Unused.
#' @return An object of class \code{"summary.nmfkc.rff"}.
#' @seealso \code{\link{nmfkc.rff}}
#' @export
summary.nmfkc.rff <- function(object, ...) {
  ans <- list(
    call       = object$call,
    Q_obs      = nrow(object$X),
    Q          = ncol(object$X),
    D          = object$D,
    iter       = object$iter,
    runtime    = object$runtime,
    objfunc    = object$objfunc,
    r.squared  = object$r.squared,
    mae        = object$mae,
    X.sparsity = if (!is.null(object$X))
                   mean(object$X < 1e-4) else NA_real_,
    Hp.sparsity = if (!is.null(object$Hp))
                    mean(object$Hp < 1e-4) else NA_real_,
    Hn.sparsity = if (!is.null(object$Hn))
                    mean(object$Hn < 1e-4) else NA_real_,
    pars       = object$pars      # RFF generating params (may be NULL)
  )
  class(ans) <- "summary.nmfkc.rff"
  ans
}

#' Print method for summary.nmfkc.rff objects
#'
#' @param x An object of class \code{"summary.nmfkc.rff"}.
#' @param digits Number of significant digits.
#' @param ... Unused.
#' @return Called for its side effect. Returns \code{x} invisibly.
#' @seealso \code{\link{summary.nmfkc.rff}}
#' @export
print.summary.nmfkc.rff <- function(x, digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Kernel-faithful NMF-RFF (Direct MU)\n")
  cat(sprintf("  Y rows (Q_obs):   %d\n", x$Q_obs))
  cat(sprintf("  Rank (Q):         %d\n", x$Q))
  cat(sprintf("  RFF dim (D):      %d\n", x$D))

  ## RFF generating parameters (if captured via pars argument)
  if (!is.null(x$pars)) {
    p_in  <- ncol(x$pars$omega)
    o_dim <- dim(x$pars$omega)
    cat("\nRFF parameters:\n")
    cat(sprintf("  Input dim (p):    %d\n", p_in))
    cat(sprintf("  beta (bandwidth): %s\n",
                format(x$pars$beta, digits = digits)))
    cat(sprintf("  omega:            %d x %d  (range [%s, %s])\n",
                o_dim[1], o_dim[2],
                format(min(x$pars$omega), digits = digits),
                format(max(x$pars$omega), digits = digits)))
    cat(sprintf("  b:                length %d  (range [%s, %s])\n",
                length(x$pars$b),
                format(min(x$pars$b), digits = digits),
                format(max(x$pars$b), digits = digits)))
  }

  cat("\nConvergence:\n")
  cat(sprintf("  Iterations:       %d\n", x$iter))
  cat(sprintf("  Runtime (secs):   %.2f\n", x$runtime))
  cat(sprintf("  Final objfunc:    %s\n", format(x$objfunc, digits = digits)))

  cat("\nGoodness of fit:\n")
  cat(sprintf("  R-squared (cor^2): %s\n",
              format(x$r.squared, digits = digits)))
  cat(sprintf("  MAE:               %s\n",
              format(x$mae, digits = digits)))

  cat("\nSparsity (< 1e-4):\n")
  cat(sprintf("  X:  %.1f%%   Hp: %.1f%%   Hn: %.1f%%\n",
              100 * x$X.sparsity,
              100 * x$Hp.sparsity,
              100 * x$Hn.sparsity))
  cat("\n")
  invisible(x)
}

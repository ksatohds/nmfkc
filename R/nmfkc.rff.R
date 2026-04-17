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
## nmfkc.rff.params : generate RFF random parameters
############################################################

#' Generate Random Fourier Features parameters
#'
#' @description
#' Generates the random parameters \eqn{\omega_d \sim \mathcal{N}(0, 2\beta I_p)}
#' and \eqn{b_d \sim \mathrm{Uniform}(0, 2\pi)} required to compute
#' Random Fourier Features (Rahimi & Recht, 2007) that approximate the
#' Gaussian kernel \eqn{k(u, v) = \exp(-\beta \|u - v\|^2)}.
#'
#' @param p Integer. Dimensionality of the input data.
#' @param D Integer. Number of random features.
#' @param beta Positive scalar. Gaussian kernel bandwidth parameter.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with elements \code{omega} (\eqn{D \times p} matrix),
#'   \code{b} (length-\eqn{D} vector), \code{D}, and \code{beta}.
#'
#' @seealso \code{\link{nmfkc.rff.apply}}, \code{\link{nmfkc.rff.direct}},
#'   \code{\link{nmfkc.kernel.beta.nearest.med}}
#'
#' @examples
#' pars <- nmfkc.rff.params(p = 5, D = 20, beta = 0.5, seed = 1)
#' str(pars)
#'
#' @export
nmfkc.rff.params <- function(p, D, beta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  list(
    omega = matrix(stats::rnorm(D * p, mean = 0, sd = sqrt(2 * beta)),
                   nrow = D, ncol = p),
    b     = stats::runif(D, 0, 2 * pi),
    D     = D,
    beta  = beta
  )
}


############################################################
## nmfkc.rff.apply : compute RFF features
############################################################

#' Compute Random Fourier Features of input data
#'
#' @description
#' Applies the RFF transform
#' \deqn{z_d(u) = \sqrt{2/D}\, \cos(\omega_d^\top u + b_d)}
#' to each column of \code{U}, yielding a \eqn{D \times N} feature matrix
#' \eqn{Z} such that \eqn{Z^\top Z \approx K} (the Gaussian kernel matrix).
#' The matrix \eqn{Z} is sign-unrestricted (contains negative values).
#'
#' @param U A \eqn{p \times N} numeric matrix; columns are data points.
#' @param pars A list returned by \code{\link{nmfkc.rff.params}}.
#'
#' @return A \eqn{D \times N} numeric matrix of RFF features.
#'
#' @seealso \code{\link{nmfkc.rff.params}}, \code{\link{nmfkc.rff.nonneg}}
#'
#' @examples
#' U <- matrix(stats::rnorm(5 * 20), 5, 20)
#' pars <- nmfkc.rff.params(p = 5, D = 50, beta = 0.5, seed = 1)
#' Z <- nmfkc.rff.apply(U, pars)
#' dim(Z)     # 50 x 20
#'
#' @export
nmfkc.rff.apply <- function(U, pars) {
  proj <- pars$omega %*% U + pars$b
  sqrt(2 / pars$D) * cos(proj)
}


############################################################
## nmfkc.rff.nonneg : posneg split
############################################################

#' Positive/negative split of a sign-unrestricted matrix
#'
#' @description
#' Splits a sign-unrestricted matrix \eqn{Z} into its positive and negative
#' parts stacked row-wise,
#' \deqn{[Z_{+}; Z_{-}],\ Z_{+} = \max(Z, 0),\ Z_{-} = \max(-Z, 0),}
#' producing a non-negative matrix suitable as a covariate for standard NMF
#' (the \dQuote{posneg split}; Chen & Plemmons 2009).
#'
#' @param Zraw A \eqn{D \times N} numeric matrix (may contain negative values).
#' @param posneg.split Logical. If \code{TRUE} (default), returns the stacked
#'   \eqn{2D \times N} matrix \eqn{[Z_{+}; Z_{-}]}. If \code{FALSE}, returns
#'   only \eqn{Z_{+} = \max(Z, 0)}.
#'
#' @return A non-negative numeric matrix.
#'
#' @seealso \code{\link{nmfkc.rff.apply}}, \code{\link{nmfkc.rff.direct}}
#'
#' @examples
#' Zraw <- matrix(stats::rnorm(30), 6, 5)
#' Z_aug <- nmfkc.rff.nonneg(Zraw)
#' dim(Z_aug)   # 12 x 5
#'
#' @export
nmfkc.rff.nonneg <- function(Zraw, posneg.split = TRUE) {
  if (isTRUE(posneg.split)) rbind(pmax(Zraw, 0), pmax(-Zraw, 0))
  else pmax(Zraw, 0)
}


############################################################
## nmfkc.rff.direct : Direct MU algorithm (Algorithm 2)
############################################################

#' Kernel-faithful NMF with Random Fourier Features (Direct MU)
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
#'   where \eqn{Z} is produced by \code{\link{nmfkc.rff.apply}}.
#' @param Zn Non-negative \eqn{D \times N} matrix \eqn{Z_{-} = \max(-Z, 0)}.
#' @param rank Integer. Number of latent components \eqn{Q} in \eqn{X}.
#' @param maxit Maximum number of iterations (default 5000).
#' @param epsilon Relative convergence tolerance on the objective
#'   (default 1e-4).
#' @param res.init Optional \code{nmfkc} object for warm-start
#'   initialization.  Its \code{X} seeds the basis and
#'   \code{C[, 1:D]} / \code{C[, (D+1):(2D)]} seed \code{Hp} / \code{Hn}.
#'   Supply the result of running \code{\link{nmfkc}} on
#'   \code{Y, A = rbind(Zp, Zn)} (the posneg model).
#' @param print.trace Logical. If \code{TRUE}, prints summary messages at
#'   start and end.
#' @param ... Additional arguments: \code{Q} is accepted as an alias for
#'   \code{rank} for backward compatibility.
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
#' @references
#' Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel
#' machines. \emph{Advances in Neural Information Processing Systems}, 20.
#'
#' Ding, C., Li, T., & Jordan, M. I. (2010). Convex and semi-nonnegative
#' matrix factorizations. \emph{IEEE TPAMI}, 32(1), 45-55.
#'
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.rff.params}},
#'   \code{\link{nmfkc.rff.apply}}, \code{\link{nmfkc.rff.nonneg}},
#'   \code{\link{nmfkc.kernel.beta.nearest.med}},
#'   \code{\link{predict.nmfkc.rff}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' U <- matrix(stats::rnorm(5 * 40), 5, 40)
#' Y <- matrix(abs(stats::rnorm(8 * 40)), 8, 40)
#' pars <- nmfkc.rff.params(p = 5, D = 50, beta = 0.5, seed = 1)
#' Zraw <- nmfkc.rff.apply(U, pars)
#' Zp <- pmax(Zraw, 0);  Zn <- pmax(-Zraw, 0)
#' res <- nmfkc.rff.direct(Y, Zp, Zn, rank = 3, maxit = 200)
#' plot(res)
#' summary(res)
#' }
#'
#' @export
nmfkc.rff.direct <- function(Y, Zp, Zn,
                              rank = NULL,
                              maxit = 5000,
                              epsilon = 1e-4,
                              res.init = NULL,
                              print.trace = FALSE,
                              ...) {
  cl <- match.call()
  extra_args <- list(...)
  if (is.null(rank) && !is.null(extra_args$Q)) rank <- extra_args$Q
  if (is.null(rank)) stop("'rank' must be specified.")
  Q <- as.integer(rank)

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
  if (!is.null(res.init)) {
    if (is.null(res.init$X) || is.null(res.init$C))
      stop("'res.init' must contain X and C matrices (e.g., from nmfkc()).")
    X  <- res.init$X
    if (ncol(res.init$C) != 2 * D_rff)
      stop("res.init$C must have 2*D columns (posneg-split warm start).")
    Hp <- res.init$C[, 1:D_rff, drop = FALSE]
    Hn <- res.init$C[, (D_rff + 1):(2 * D_rff), drop = FALSE]
  } else {
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
    cat(sprintf("nmfkc.rff.direct: Y(%d,%d) Z(%d,%d) rank=%d\n",
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
    cat(sprintf("nmfkc.rff.direct: done in %d iterations, %.2fs, R2=%.4f\n",
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
#' @seealso \code{\link{nmfkc.rff.direct}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' U <- matrix(stats::rnorm(5 * 40), 5, 40)
#' Y <- matrix(abs(stats::rnorm(8 * 40)), 8, 40)
#' pars <- nmfkc.rff.params(p = 5, D = 50, beta = 0.5, seed = 1)
#' Zraw <- nmfkc.rff.apply(U, pars); Zp <- pmax(Zraw, 0); Zn <- pmax(-Zraw, 0)
#' res <- nmfkc.rff.direct(Y, Zp, Zn, rank = 3, maxit = 100)
#' Yhat <- predict(res, newZp = Zp, newZn = Zn)
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
#' \code{\link{nmfkc.rff.direct}} model.
#'
#' @param x An \code{nmfkc.rff} object.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @return Called for its side effect (plot). Returns \code{x} invisibly.
#' @seealso \code{\link{nmfkc.rff.direct}}
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
#' Produces a concise summary of a \code{\link{nmfkc.rff.direct}} fit.
#'
#' @param object A fitted \code{nmfkc.rff} object.
#' @param ... Unused.
#' @return An object of class \code{"summary.nmfkc.rff"}.
#' @seealso \code{\link{nmfkc.rff.direct}}
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
                    mean(object$Hn < 1e-4) else NA_real_
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

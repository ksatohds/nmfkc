# =====================================================
# nmfkc.signed.R — NMF-KC with signed covariate matrix
#
# Extends nmfkc() to allow the covariate matrix A to contain negative
# entries (and also the response Y).  The model is
#
#   Y ≈ X Θ A,    X ≥ 0,  Θ ∈ R^{Q×D},  A ∈ R^{D×N},
#
# where Θ and A may be signed.  Internally the function uses the
# sign-splitting trick (Ding, Li, & Jordan 2010) to solve
#
#   min_{X, C_+, C_-}  || Y − X (C_+ − C_-)(A_+ − A_-) ||_F^2
#   s.t.  X ≥ 0,  C_+ ≥ 0,  C_- ≥ 0, (optionally colSums(X) = 1)
#
# by a Direct Multiplicative Update (Direct MU) algorithm that avoids
# the N × N Gram matrix A A^T is D × D.  On each iteration, after the
# X update, columns of X are rescaled to satisfy the requested
# X.restriction (default colSums), with the scale absorbed into C_+,
# C_- so that X Θ A is unchanged.
#
# The primary application is kernel approximation via Random Fourier
# Features (RFF; Rahimi & Recht 2007), but A can be any real-valued
# matrix (PCA scores, wavelet coefficients, neural-net activations,
# RFF features, ...).
#
# References:
#   Ding, C. H. Q., Li, T., & Jordan, M. I. (2010).  Convex and
#     semi-nonnegative matrix factorizations.  IEEE TPAMI 32(1), 45–55.
#   Rahimi, A., & Recht, B. (2007).  Random features for large-scale
#     kernel machines.  NIPS.
# =====================================================


#' NMF-KC with signed covariate matrix
#'
#' @description
#' Solves
#' \deqn{Y \approx X\,\Theta\,A,\qquad X \ge 0,\;\Theta\in\R^{Q\times D},
#'   \;A\in\R^{D\times N},}
#' where the covariate matrix \eqn{A} and the coefficient matrix
#' \eqn{\Theta} may be \strong{signed}.  Internally \eqn{A = A_{+} - A_{-}}
#' and \eqn{\Theta = C_{+} - C_{-}} with \eqn{A_{\pm}, C_{\pm} \ge 0}
#' (sign-splitting trick, Ding et al. 2010), and the problem is solved
#' by a Direct Multiplicative Update algorithm whose iteration cost is
#' \eqn{O(Q D^2)}, independent of \eqn{N}.
#'
#' Only \eqn{X} is structurally constrained to be non-negative (Semi-NMF
#' sense of Ding, Li, & Jordan 2010).  In particular, \strong{\eqn{Y} may
#' contain negative entries}, in which case the response is fit in the
#' least-squares sense without any non-negativity requirement on \eqn{Y}.
#'
#' When \eqn{A \ge 0} (so \eqn{A_{-} = 0}), the result reduces to
#' \code{\link{nmfkc}}(Y, A, rank) with Euclidean loss, up to reordering.
#'
#' @param Y Real-valued \eqn{Q_{\mathrm{obs}} \times N} response matrix.
#'   Unlike \code{\link{nmfkc}}, negative entries are allowed.
#' @param A Real-valued \eqn{D \times N} covariate matrix (signed).
#'   A single matrix is passed; its positive and negative parts
#'   \eqn{A_{+} = \max(A, 0)} and \eqn{A_{-} = \max(-A, 0)} are computed
#'   internally.  When using Random Fourier Features (Rahimi & Recht
#'   2007) as \eqn{A}, supply the RFF parameters via the hidden
#'   \code{pars} argument so that \code{predict()} can regenerate
#'   features for new data (see \code{pars} entry in \code{\dots} below).
#' @param rank Integer.  Number of latent components \eqn{Q} in \eqn{X}.
#' @param epsilon Relative convergence tolerance on the objective
#'   (default \code{1e-4}).
#' @param maxit Maximum number of iterations (default \code{5000}).
#' @param verbose Logical.  Print dimensions at start (default \code{TRUE}).
#' @param ... Additional arguments:
#'   \itemize{
#'     \item \code{Q}: alias for \code{rank}.
#'     \item \code{X.restriction}: constraint applied to columns of
#'       \eqn{X} after every update, with the scale absorbed into
#'       \eqn{C_{+}, C_{-}}.  One of \code{"colSums"} (default,
#'       \eqn{\mathrm{colSums}(X) = 1}), \code{"colSqSums"},
#'       \code{"totalSum"}, \code{"none"}, \code{"fixed"}.
#'     \item \code{X.init}: initialization for \eqn{X}.  Either
#'       \code{"random"} (default), \code{"nndsvd"} (requires \eqn{Y \ge 0}),
#'       or an explicit \eqn{Q_{\mathrm{obs}} \times Q} numeric matrix.
#'     \item \code{C.init}: explicit initial \eqn{Q \times D} coefficient
#'       matrix \eqn{\Theta} (signed).  Split internally.
#'     \item \code{warm.start}: logical (default \code{TRUE}).  If
#'       \code{TRUE} \strong{and \eqn{Y \ge 0}}, runs
#'       \code{nmfkc(Y, A = rbind(A_+, A_-), rank = Q)} internally to
#'       seed \eqn{X, C_{+}, C_{-}}.  When the warm-start path is
#'       active, the user's \code{X.init}, \code{seed}, \code{nstart},
#'       and \code{X.restriction} arguments are forwarded to the
#'       internal \code{\link{nmfkc}} call so that initialization
#'       choices propagate (one exception: \code{X.init = "random"} is
#'       a \code{nmfkc.signed()}-specific legacy string meaning
#'       \code{abs(rnorm) * 0.1} and is not recognized by
#'       \code{\link{nmfkc}}; in that case the warm-start falls back to
#'       \code{\link{nmfkc}}'s own default, \code{"kmeans"}).  Ignored
#'       when \eqn{Y} has negative entries (warm-start is disabled).
#'     \item \code{seed}: RNG seed for random initialization (default 123).
#'     \item \code{prefix}: name prefix for rows of \eqn{C} and columns
#'       of \eqn{X} (default \code{"Basis"}).
#'     \item \code{pars}: optional list \code{list(omega, b, D, beta)}
#'       of Random Fourier Feature parameters (Rahimi & Recht 2007;
#'       \code{omega}: frequency matrix, \code{b}: phase offset,
#'       \code{D}: feature dimension, \code{beta}: bandwidth).  When
#'       supplied, it is stored in the returned object so that
#'       \code{summary()} can report \eqn{\beta} and downstream
#'       \code{predict()} calls can regenerate RFF features for new
#'       data.  If \code{A} is not RFF features, leave this \code{NULL}.
#'     \item \code{Y.weights}: Optional non-negative weight matrix
#'       (\eqn{Q_{\mathrm{obs}} \times N}) or vector (length \eqn{N}),
#'       analogous to the \code{weights} argument of
#'       \code{\link[stats]{lm}}.  Loss becomes
#'       \eqn{\sum W_{ij} \, (Y_{ij} - (XCA)_{ij})^2}
#'       (\code{lm()}-style, \strong{linear} in \eqn{W}).  Logical
#'       matrices (\code{TRUE} / \code{FALSE}) are also accepted.
#'       Typical usage by \code{\link{nmfkc.signed.cv}} /
#'       \code{\link{nmfkc.signed.ecv}} passes a binary mask
#'       \eqn{W \in \{0,1\}} to hold out test elements; real-valued
#'       weights for observation-level importance weighting are also
#'       supported.  Default \code{NULL}: if \code{Y} has \code{NA},
#'       a binary mask is auto-constructed (0 for \code{NA}, 1
#'       elsewhere); otherwise no weighting.
#'     \item \code{nstart}: number of random restarts.  \strong{Signed
#'       models have more local minima than non-negative ones} because
#'       \eqn{\Theta = C_{+} - C_{-}} can take both positive and negative
#'       values.  Since \code{nmfkc.signed()} itself does not loop over
#'       restarts (callers control it), set the outer-loop size via e.g.
#'       running the function several times with different \code{seed}
#'       and keeping the fit with the smallest \code{$objfunc}.  A
#'       restart budget of 10-50 is recommended for publication-grade
#'       runs on signed data.
#'   }
#'
#' @return An object of class \code{c("nmfkc.signed", "nmfkc")} with
#' \itemize{
#'   \item \code{X}: \eqn{Q_{\mathrm{obs}} \times Q} basis matrix (non-negative,
#'     column-normalized according to \code{X.restriction}).
#'   \item \code{Cp}, \code{Cn}: \eqn{Q \times D} non-negative parts of
#'     \eqn{\Theta}, so that \eqn{\Theta = C_{+} - C_{-}}.
#'   \item \code{C}: \eqn{C_{+} - C_{-}} (= \eqn{\Theta}), signed.
#'   \item \code{B}: \eqn{C \, A}, \eqn{Q \times N} (signed).
#'   \item \code{objfunc.iter}: objective values per iteration.
#'   \item \code{objfunc}: final objective.
#'   \item \code{r.squared}: \eqn{\mathrm{cor}(Y, \widehat Y)^2}.
#'   \item \code{mae}: mean absolute error.
#'   \item \code{iter}: number of iterations performed.
#'   \item \code{runtime}: elapsed seconds.
#'   \item \code{Y.signed}: logical; whether \eqn{Y} contained negative
#'     entries during fitting.
#'   \item \code{pars}: RFF generating parameters, if supplied.
#'   \item \code{call}: the matched call.
#' }
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.  The interface may change in
#' future versions.
#'
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010).  Convex and
#' semi-nonnegative matrix factorizations.  \emph{IEEE TPAMI}, 32(1), 45–55.
#'
#' Rahimi, A., & Recht, B. (2007).  Random features for large-scale
#' kernel machines.  \emph{Advances in NIPS}, 20.
#'
#' @seealso \code{\link{nmfkc}}, \code{\link{predict.nmfkc.signed}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' ## Example 1: signed A (e.g., hand-built RFF features), non-negative Y
#' ## Build simple signed features Z = sqrt(2/D) * cos(omega^T U + b):
#' U     <- matrix(stats::rnorm(5 * 40), 5, 40)           # raw input
#' D     <- 20                                            # feature dim
#' omega <- matrix(stats::rnorm(5 * D), 5, D)             # random freqs
#' b     <- stats::runif(D, 0, 2 * pi)                    # phase
#' Z     <- sqrt(2 / D) *
#'            cos(t(omega) %*% U + matrix(b, D, 40))      # D x 40, signed
#' Y     <- matrix(abs(stats::rnorm(8 * 40)), 8, 40)
#' res1  <- nmfkc.signed(Y, A = Z, rank = 3, maxit = 200)
#'
#' ## Example 2: signed Y (regression)
#' Y2    <- matrix(stats::rnorm(8 * 40), 8, 40)           # signed response
#' res2  <- nmfkc.signed(Y2, A = Z, rank = 3, maxit = 200,
#'                        warm.start = FALSE)
#' }
#'
#' @export
nmfkc.signed <- function(Y, A, rank = NULL,
                         epsilon = 1e-4, maxit = 5000,
                         verbose = TRUE, ...) {
  cl <- match.call()
  extra_args <- list(...)

  ## --- 1. Parameter extraction ---
  if (is.null(rank) && !is.null(extra_args$Q)) rank <- extra_args$Q
  if (is.null(rank)) stop("'rank' must be specified.")
  Q <- as.integer(rank)

  X.restriction <- if (!is.null(extra_args$X.restriction))
    extra_args$X.restriction else "colSums"
  X.restriction <- match.arg(X.restriction,
    c("colSums", "colSqSums", "totalSum", "none", "fixed"))

  X.init     <- if (!is.null(extra_args$X.init))     extra_args$X.init     else "random"
  C.init     <- if (!is.null(extra_args$C.init))     extra_args$C.init     else NULL
  warm.start <- if (!is.null(extra_args$warm.start)) extra_args$warm.start else TRUE
  seed       <- if (!is.null(extra_args$seed))       extra_args$seed       else 123L
  prefix     <- if (!is.null(extra_args$prefix))     extra_args$prefix     else "Basis"
  pars_rff   <- extra_args$pars
  Y.weights  <- extra_args$Y.weights

  ## --- 2. Input preparation & validation ---
  if (is.vector(Y)) Y <- matrix(Y, nrow = 1)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(A)) A <- as.matrix(A)
  Q_obs <- nrow(Y); N <- ncol(Y); D <- nrow(A)
  if (ncol(A) != N) stop("ncol(A) must equal ncol(Y).")
  if (any(is.na(A))) stop("A contains NA; please impute or remove.")

  ## --- Y.weights handling (supports vector, matrix, NA-autodetect) ---
  ## When supplied, the loss becomes sum(W * (Y - XCA)^2) and held-out
  ## elements (W = 0) are masked.  Y-NA auto-masking matches nmfkc().
  if (!is.null(Y.weights) && is.vector(Y.weights) && !is.matrix(Y.weights)) {
    if (length(Y.weights) == N) {
      Y.weights <- matrix(Y.weights, nrow = Q_obs, ncol = N, byrow = TRUE)
    } else if (length(Y.weights) == 1) {
      Y.weights <- matrix(Y.weights, nrow = Q_obs, ncol = N)
    } else {
      stop("Length of Y.weights vector must match ncol(Y) or be 1.")
    }
  }
  if (is.null(Y.weights) && any(is.na(Y))) {
    Y.weights <- matrix(1, nrow = Q_obs, ncol = N)
    Y.weights[is.na(Y)] <- 0
    Y[is.na(Y)] <- 0
  } else if (!is.null(Y.weights)) {
    Y.weights <- as.matrix(Y.weights)
    if (!all(dim(Y.weights) == dim(Y)))
      stop("Y.weights dimensions must match Y.")
    Y.weights[is.na(Y.weights)] <- 0
    Y[is.na(Y) | Y.weights == 0] <- 0
  } else if (any(is.na(Y))) {
    stop("Y contains NA; please impute, remove, or supply Y.weights.")
  }
  has.weights <- !is.null(Y.weights) && any(Y.weights != 1)
  if (!has.weights) Y.weights <- NULL
  Y_is_nonneg <- all(Y >= 0)

  if (isTRUE(verbose)) {
    msg <- sprintf("Y(%d,%d) ~ X(%d,%d) %%*%% C(%d,%d) %%*%% A(%d,%d)  [signed covariate]%s%s",
                   Q_obs, N, Q_obs, Q, Q, D, D, N,
                   if (Y_is_nonneg) "" else ", Y signed",
                   if (has.weights) ", weighted" else "")
    message(msg)
  }

  t0 <- proc.time()

  ## --- 3. Precomputation (avoids N x N Gram) ---
  Ap <- pmax(A,  0)
  An <- pmax(-A, 0)
  A_diff <- Ap - An                 # = A
  ## For weighted loss, S / G0 are NOT iteration-invariant, so we skip the
  ## precomputation and use a slower weighted path inside the main loop.
  ## For the unweighted loss, the precomputation gives the O(QD^2) speedup.
  if (!has.weights) {
    S   <- tcrossprod(A_diff)         # D x D, signed
    G0  <- Y %*% t(A_diff)            # Q_obs x D, signed (Y signed OK)
    S_p <- pmax(S,  0)
    S_n <- pmax(-S, 0)
  } else {
    S <- NULL; G0 <- NULL; S_p <- NULL; S_n <- NULL
  }
  Y_sqnorm <- sum(Y * Y)            # always >= 0

  ## --- 4. X.restriction helpers ---
  xscale <- switch(X.restriction,
    colSums   = function(X) colSums(X) + 1e-16,
    colSqSums = function(X) sqrt(colSums(X * X)) + 1e-16,
    totalSum  = function(X) rep(sum(X) + 1e-16, ncol(X)),
    none      = function(X) rep(1, ncol(X)),
    fixed     = function(X) rep(1, ncol(X)))

  ## --- 5. Initialization ---
  X <- NULL; Cp <- NULL; Cn <- NULL

  ## 5a. Warm-start via posneg nmfkc() (requires Y >= 0; disabled when weighted
  ## because the warm-start init can produce numerical explosion in the
  ## weighted MU loop — ECV CV callers get random init instead).
  explicit_X_mat <- is.matrix(X.init) || (is.numeric(X.init) && length(X.init) > 1)
  need_warm <- isTRUE(warm.start) && Y_is_nonneg && !has.weights &&
               is.null(C.init) && !explicit_X_mat
  if (need_warm) {
    warm_args <- list(Y, A = rbind(Ap, An), rank = Q,
                      epsilon = epsilon, maxit = maxit, verbose = FALSE,
                      seed = seed, X.restriction = X.restriction)
    if (has.weights) warm_args$Y.weights <- Y.weights
    ## Forward X.init to nmfkc() so that user-selected initialization
    ## (e.g., "nndsvd", "kmeans", or a user-supplied matrix) propagates
    ## into the posneg warm-start.  "random" is a nmfkc.signed()-specific
    ## string (= abs(rnorm)*0.1) not recognized by nmfkc(); in that case
    ## fall back to nmfkc()'s own default ("kmeans") for the warm-start.
    if (!identical(X.init, "random")) warm_args$X.init <- X.init
    ## Forward nstart if the user supplied it (nmfkc() default is 1).
    if (!is.null(extra_args$nstart)) warm_args$nstart <- extra_args$nstart
    res0 <- do.call(nmfkc, warm_args)
    X  <- res0$X
    Cp <- res0$C[, 1:D, drop = FALSE]
    Cn <- res0$C[, (D + 1):(2 * D), drop = FALSE]
  } else {
    ## 5b. Explicit X.init / random
    if (explicit_X_mat) {
      X <- as.matrix(X.init)
      if (!identical(dim(X), c(Q_obs, Q)))
        stop("X.init must have dimensions (nrow(Y), rank).")
    } else if (identical(X.init, "nndsvd") && Y_is_nonneg) {
      ## Simple NNDSVD-like deterministic init via SVD of Y
      sv <- svd(Y, nu = Q, nv = 0)
      X  <- pmax(sv$u, 0)
      bad <- colSums(X) < 1e-10
      if (any(bad)) X[, bad] <- matrix(abs(stats::rnorm(Q_obs * sum(bad))) * 0.1,
                                        Q_obs, sum(bad))
    } else {
      set.seed(seed)
      X <- matrix(abs(stats::rnorm(Q_obs * Q)) * 0.1, Q_obs, Q)
    }
    ## C.init or random
    if (!is.null(C.init)) {
      if (!identical(dim(C.init), c(Q, D)))
        stop("C.init must have dimensions (rank, nrow(A)).")
      Cp <- pmax(C.init,  0)
      Cn <- pmax(-C.init, 0)
    } else {
      set.seed(seed + 1L)
      Cp <- matrix(abs(stats::rnorm(Q * D)) * 0.1,  Q, D)
      Cn <- matrix(abs(stats::rnorm(Q * D)) * 0.01, Q, D)
    }
  }

  ## Apply X.restriction initially (absorb into Cp, Cn so X %*% (Cp-Cn) %*% A is unchanged)
  if (X.restriction != "fixed" && X.restriction != "none") {
    d0 <- xscale(X)
    X  <- sweep(X,  2, d0, "/")
    Cp <- sweep(Cp, 1, d0, "*")
    Cn <- sweep(Cn, 1, d0, "*")
  }

  small <- 1e-16
  Wmat <- Y.weights  # short alias; NULL if no weights

  compute_obj <- function(X, Cp, Cn) {
    Yhat <- X %*% (Cp - Cn) %*% A_diff
    ## lm()-style weighted least squares: L = sum(W * (Y - Yhat)^2).
    ## W = W^2 for binary {0,1} masks so this is unchanged for ECV.
    if (has.weights) sum(Wmat * (Y - Yhat)^2) else sum((Y - Yhat)^2)
  }

  if (!has.weights) {
    P <- crossprod(X); G <- crossprod(X, G0)
    H <- Cp - Cn; PH <- P %*% H
    obj_prev <- Y_sqnorm - 2 * sum(G * H) + sum(H * (PH %*% S))
  } else {
    obj_prev <- compute_obj(X, Cp, Cn)
  }
  objfunc.iter <- numeric(maxit)

  ## --- 6. Main loop ---
  iter <- 0L
  for (iter in seq_len(maxit)) {

    if (!has.weights) {
      ## ---- Fast unweighted path (precomputed S, G0) ----
      ## 6a. Cp update
      G_p <- pmax(G, 0); G_n <- pmax(-G, 0)
      PCp <- P %*% Cp;   PCn <- P %*% Cn
      Cp  <- Cp * (G_p + PCp %*% S_n + PCn %*% S_p) /
                  (G_n + PCp %*% S_p + PCn %*% S_n + small)

      ## 6b. Cn update (Gauss-Seidel)
      PCp <- P %*% Cp
      Cn  <- Cn * (G_n + PCp %*% S_p + PCn %*% S_n) /
                  (G_p + PCp %*% S_n + PCn %*% S_p + small)

      ## 6c. X update
      if (X.restriction != "fixed") {
        H   <- Cp - Cn; Ht <- t(H)
        YMt <- G0 %*% Ht
        HS  <- H %*% S; MMt <- HS %*% Ht
        X <- X * (pmax(YMt, 0) + X %*% pmax(-MMt, 0)) /
                 (pmax(-YMt, 0) + X %*% pmax(MMt, 0) + small)
        if (X.restriction != "none") {
          d <- xscale(X)
          X  <- sweep(X,  2, d, "/")
          Cp <- sweep(Cp, 1, d, "*")
          Cn <- sweep(Cn, 1, d, "*")
        }
      }

      ## 6d. Refresh precomputed quantities & evaluate objective in closed form
      P  <- crossprod(X); G <- crossprod(X, G0)
      H  <- Cp - Cn; PH <- P %*% H
      obj_cur <- Y_sqnorm - 2 * sum(G * H) + sum(H * (PH %*% S))

    } else {
      ## ---- Weighted path (no S/G0 precompute) ----
      ## Let Yhat_+ = X Cp A, Yhat_- = X Cn A; residual under W.
      tX <- t(X)
      Yhat_p <- X %*% Cp %*% A_diff          # signed if A signed
      Yhat_n <- X %*% Cn %*% A_diff
      WY     <- Wmat * Y                      # weighted observations
      WYhp   <- Wmat * Yhat_p
      WYhn   <- Wmat * Yhat_n

      ## 6a'. Cp update via weighted split
      ## grad_{Cp} L = -2 X^T (W*(Y - X(Cp-Cn)A)) A^T
      ##             = -2 ( X^T (W*Y) A^T - X^T (W*Yhat_p) A^T + X^T (W*Yhat_n) A^T )
      G_w   <- tX %*% WY   %*% t(A_diff)      # Q x D, signed
      Hp_w  <- tX %*% WYhp %*% t(A_diff)      # Q x D, signed
      Hn_w  <- tX %*% WYhn %*% t(A_diff)      # Q x D, signed
      Gp <- pmax(G_w, 0); Gn <- pmax(-G_w, 0)
      Hpp <- pmax(Hp_w, 0); Hpn <- pmax(-Hp_w, 0)
      Hnp <- pmax(Hn_w, 0); Hnn <- pmax(-Hn_w, 0)
      Cp <- Cp * (Gp + Hpn + Hnp) / (Gn + Hpp + Hnn + small)
      ## Recompute Hp_w with updated Cp
      Yhat_p <- X %*% Cp %*% A_diff
      WYhp   <- Wmat * Yhat_p
      Hp_w   <- tX %*% WYhp %*% t(A_diff)
      Hpp <- pmax(Hp_w, 0); Hpn <- pmax(-Hp_w, 0)
      ## 6b'. Cn update
      Cn <- Cn * (Gn + Hpp + Hnn) / (Gp + Hpn + Hnp + small)

      ## 6c'. X update (weighted)
      ## -dL/(2) = (W*Y) M^T - (W*(XM)) M^T   where M = H A (signed).
      ## Ding split on the two signed terms A1, A2:
      ##   Pull = [A1]_+ + [A2]_-,  Push = [A1]_- + [A2]_+
      if (X.restriction != "fixed") {
        H   <- Cp - Cn
        M   <- H %*% A_diff        # Q x N, signed
        XM  <- X %*% M             # Q_obs x N, signed
        WY  <- Wmat * Y
        WXM <- Wmat * XM
        A1  <- WY  %*% t(M)        # Q_obs x Q, signed
        A2  <- WXM %*% t(M)        # Q_obs x Q, signed
        num <- pmax(A1, 0) + pmax(-A2, 0)
        den <- pmax(-A1, 0) + pmax(A2, 0)
        X <- X * num / (den + small)
        if (X.restriction != "none") {
          d <- xscale(X)
          X  <- sweep(X,  2, d, "/")
          Cp <- sweep(Cp, 1, d, "*")
          Cn <- sweep(Cn, 1, d, "*")
        }
      }

      obj_cur <- compute_obj(X, Cp, Cn)
    }

    objfunc.iter[iter] <- obj_cur
    ## Safety: bail out on NaN/Inf (numerical explosion in weighted MU)
    if (!is.finite(obj_cur)) {
      warning("nmfkc.signed: objective became non-finite at iter ", iter,
              "; stopping early.")
      break
    }
    if (is.finite(obj_prev) &&
        abs(obj_prev - obj_cur) / max(abs(obj_prev), 1e-12) < epsilon) break
    obj_prev <- obj_cur
  }
  objfunc.iter <- objfunc.iter[seq_len(iter)]

  ## --- 7. Post-processing: sort columns of X (nmfkc-style centroid order) ---
  if (ncol(X) > 1 && X.restriction != "fixed") {
    index <- order(as.vector(matrix(1:nrow(X) / nrow(X), nrow = 1) %*% X))
    X  <- X[, index, drop = FALSE]
    Cp <- Cp[index, , drop = FALSE]
    Cn <- Cn[index, , drop = FALSE]
  }

  ## --- 8. Reconstruction statistics ---
  C  <- Cp - Cn                     # Q x D, signed  (= Theta)
  B  <- C %*% A_diff                # Q x N, signed
  XB <- X %*% B                     # Q_obs x N, Yhat
  resid <- Y - XB
  if (has.weights) {
    ## lm()-style weighted least squares (matches compute_obj in the loop).
    objfunc <- sum(Wmat * resid^2)
    valid <- (Wmat > 0)
    r.squared <- tryCatch(
      stats::cor(as.vector(XB)[valid], as.vector(Y)[valid])^2,
      error = function(e) NA_real_
    )
    mae <- sum(Wmat * abs(resid)) / max(sum(Wmat), small)
  } else {
    objfunc <- sum(resid * resid)
    r.squared <- tryCatch(
      stats::cor(as.vector(XB), as.vector(Y))^2,
      error = function(e) NA_real_
    )
    mae <- mean(abs(resid))
  }

  ## --- 9. Names ---
  rownames(X)  <- rownames(Y)
  colnames(X)  <- paste0(prefix, seq_len(Q))
  rownames(Cp) <- paste0(prefix, seq_len(Q))
  rownames(Cn) <- paste0(prefix, seq_len(Q))
  rownames(C)  <- paste0(prefix, seq_len(Q))
  if (!is.null(rownames(A))) {
    colnames(Cp) <- rownames(A)
    colnames(Cn) <- rownames(A)
    colnames(C)  <- rownames(A)
  }

  runtime <- as.numeric((proc.time() - t0)[3])

  ## Fields matching the main `nmfkc()` object for S3 method inheritance
  dims_str <- sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)",
                      Q_obs, N, Q_obs, Q, Q, D, D, N)
  n.valid   <- if (has.weights) sum(Wmat > 0) else Q_obs * N
  n.total   <- Q_obs * N
  n.missing <- n.total - n.valid
  sigma     <- if (n.valid > 0)
                 sqrt(sum(if (has.weights) Wmat * resid^2 else resid^2) / n.valid)
               else NA_real_

  result <- list(
    X             = X,
    Cp            = Cp,
    Cn            = Cn,
    C             = C,
    B             = B,
    XB            = XB,              # reconstruction (alias for fitted.nmf)
    objfunc.iter  = objfunc.iter,
    objfunc       = objfunc,
    r.squared     = r.squared,
    sigma         = sigma,
    mae           = mae,
    iter          = iter,
    runtime       = runtime,
    rank          = Q,
    D             = D,
    dims          = dims_str,
    method        = "EU",
    n.missing     = n.missing,
    n.total       = n.total,
    X.restriction = X.restriction,
    Y.signed      = !Y_is_nonneg,
    pars          = pars_rff,
    call          = cl
  )
  class(result) <- c("nmfkc.signed", "nmfkc", "nmf")
  result
}


#' Predict method for nmfkc.signed
#'
#' @description
#' Computes \eqn{\widehat Y = X \, C \, A_{\mathrm{new}}}
#' (\eqn{= X (C_{+} - C_{-})(A_{+}^{\mathrm{new}} - A_{-}^{\mathrm{new}})}).
#' For \code{type = "response"} the raw prediction is returned
#' (possibly signed).  For \code{type = "prob"} and \code{"class"},
#' negative entries of \eqn{\widehat Y} are clipped to zero before
#' column normalization, since probabilities must be non-negative.
#'
#' @param object A fitted \code{"nmfkc.signed"} object.
#' @param newA Real-valued \eqn{D \times N_{\mathrm{new}}} covariate matrix.
#' @param type Output: \code{"response"} (raw signed), \code{"prob"},
#'   or \code{"class"}.
#' @param ... Unused.
#'
#' @return A numeric matrix (\code{"response"} or \code{"prob"}) or a
#'   character vector (\code{"class"}).
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions.
#'
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
#' semi-nonnegative matrix factorizations. \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, 32(1), 45--55.
#'
#' @export
predict.nmfkc.signed <- function(object, newA = NULL,
                                  type = c("response", "prob", "class"),
                                  ...) {
  type <- match.arg(type)
  if (is.null(newA)) stop("'newA' must be supplied.")
  if (!is.matrix(newA)) newA <- as.matrix(newA)

  Yhat <- object$X %*% object$C %*% newA    # may be signed
  if (type == "response") return(Yhat)

  ## For prob / class, clip negatives and column-normalize
  Yhat <- pmax(Yhat, 0)
  col_sums <- colSums(Yhat) + 1e-12
  probs <- sweep(Yhat, 2, col_sums, "/")
  if (type == "prob") return(probs)

  classes <- rownames(Yhat)
  if (is.null(classes)) classes <- as.character(seq_len(nrow(Yhat)))
  classes[apply(probs, 2, which.max)]
}


#' Plot method for nmfkc.signed (convergence)
#'
#' @param x An \code{nmfkc.signed} object.
#' @param ... Passed to \code{plot()}.
#' @return Invisible \code{x}.
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions.
#'
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
#' semi-nonnegative matrix factorizations. \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, 32(1), 45--55.
#'
#' @export
plot.nmfkc.signed <- function(x, ...) {
  extra_args <- list(...)
  args <- list(x = x$objfunc.iter, type = "l")
  if (is.null(extra_args$main))
    args$main <- sprintf("r.squared = %.3f", x$r.squared)
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  do.call(graphics::plot, c(args, extra_args))
  invisible(x)
}


#' Summary method for nmfkc.signed
#'
#' @param object An \code{nmfkc.signed} object.
#' @param ... Unused.
#' @return An object of class \code{"summary.nmfkc.signed"}.
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions.
#'
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
#' semi-nonnegative matrix factorizations. \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, 32(1), 45--55.
#'
#' @export
summary.nmfkc.signed <- function(object, ...) {
  .negfrac <- function(M) {
    if (is.null(M)) return(NA_real_)
    sum(pmax(-M, 0)) / (sum(abs(M)) + 1e-16)
  }
  .range   <- function(M) if (is.null(M)) c(NA_real_, NA_real_) else range(M)
  ans <- list(
    call          = object$call,
    Q_obs         = nrow(object$X),
    Q             = ncol(object$X),
    D             = object$D,
    X.restriction = object$X.restriction,
    Y.signed      = isTRUE(object$Y.signed),
    iter          = object$iter,
    runtime       = object$runtime,
    objfunc       = object$objfunc,
    r.squared     = object$r.squared,
    mae           = object$mae,
    ## Sparsity
    X.sparsity    = if (!is.null(object$X))  mean(object$X  < 1e-4) else NA_real_,
    Cp.sparsity   = if (!is.null(object$Cp)) mean(object$Cp < 1e-4) else NA_real_,
    Cn.sparsity   = if (!is.null(object$Cn)) mean(object$Cn < 1e-4) else NA_real_,
    ## Range (min, max)
    X.range       = .range(object$X),
    Cp.range      = .range(object$Cp),
    Cn.range      = .range(object$Cn),
    C.range       = .range(object$C),
    B.range       = .range(object$B),
    ## Negative-mass ratio  sum(max(-M,0)) / sum(|M|)
    X.negfrac     = .negfrac(object$X),       # X >= 0, so always 0
    C.negfrac     = .negfrac(object$C),
    B.negfrac     = .negfrac(object$B),
    pars          = object$pars
  )
  class(ans) <- "summary.nmfkc.signed"
  ans
}


#' Print method for summary.nmfkc.signed
#'
#' @param x Object of class \code{"summary.nmfkc.signed"}.
#' @param digits Number of significant digits.
#' @param ... Unused.
#' @return Invisible \code{x}.
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions.
#'
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
#' semi-nonnegative matrix factorizations. \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, 32(1), 45--55.
#'
#' @export
print.summary.nmfkc.signed <- function(x,
                                        digits = max(3L, getOption("digits") - 3L),
                                        ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Semi-NMF-KC with signed covariate (Direct MU)\n")
  cat(sprintf("  Y rows (Q_obs):    %d%s\n", x$Q_obs,
              if (x$Y.signed) "  (Y has negative entries)" else ""))
  cat(sprintf("  Rank (Q):          %d\n", x$Q))
  cat(sprintf("  Covariate dim (D): %d\n", x$D))
  cat(sprintf("  X.restriction:     %s\n", x$X.restriction))

  if (!is.null(x$pars)) {
    p_in  <- ncol(x$pars$omega)
    o_dim <- dim(x$pars$omega)
    cat("\nRFF parameters (for reference):\n")
    cat(sprintf("  Input dim (p):    %d\n", p_in))
    cat(sprintf("  beta (bandwidth): %s\n",
                format(x$pars$beta, digits = digits)))
    cat(sprintf("  omega:            %d x %d  (range [%s, %s])\n",
                o_dim[1], o_dim[2],
                format(min(x$pars$omega), digits = digits),
                format(max(x$pars$omega), digits = digits)))
  }

  cat("\nConvergence:\n")
  cat(sprintf("  Iterations:        %d\n", x$iter))
  cat(sprintf("  Runtime (secs):    %.2f\n", x$runtime))
  cat(sprintf("  Final objfunc:     %s\n", format(x$objfunc, digits = digits)))

  cat("\nGoodness of fit:\n")
  cat(sprintf("  R-squared (cor^2): %s\n",
              format(x$r.squared, digits = digits)))
  cat(sprintf("  MAE:               %s\n",
              format(x$mae, digits = digits)))

  cat("\nStructure (range / sparsity / negative mass):\n")
  fmt_range <- function(r)
    sprintf("[%s, %s]",
            format(r[1], digits = digits, width = 1),
            format(r[2], digits = digits, width = 1))
  cat(sprintf("  X  (Q_obs x Q):  range %s   sparsity %5.1f%%   neg-mass %5.1f%%\n",
              fmt_range(x$X.range),  100 * x$X.sparsity,  100 * x$X.negfrac))
  cat(sprintf("  Cp (Q x D):      range %s   sparsity %5.1f%%   (non-negative)\n",
              fmt_range(x$Cp.range), 100 * x$Cp.sparsity))
  cat(sprintf("  Cn (Q x D):      range %s   sparsity %5.1f%%   (non-negative)\n",
              fmt_range(x$Cn.range), 100 * x$Cn.sparsity))
  cat(sprintf("  C=Cp-Cn (Q x D): range %s                   neg-mass %5.1f%%\n",
              fmt_range(x$C.range),                          100 * x$C.negfrac))
  cat(sprintf("  B=CA (Q x N):    range %s                   neg-mass %5.1f%%\n",
              fmt_range(x$B.range),                          100 * x$B.negfrac))
  cat("\n")
  invisible(x)
}


## ==============================================================
## Cross-validation helpers for nmfkc.signed
## ==============================================================

#' Column-wise k-fold cross-validation for nmfkc.signed
#'
#' @description Column-wise k-fold CV by held-out samples: for each
#' fold, the model is fit on the training columns and evaluated on the
#' held-out columns by solving for new-sample coefficients via a
#' weighted refit with \eqn{X} fixed.
#'
#' @param Y Real-valued \eqn{Q_{\mathrm{obs}} \times N} response matrix
#'   (signed entries allowed).
#' @param A Real-valued \eqn{D \times N} covariate matrix (signed).
#' @param rank Integer \eqn{Q}.
#' @param ... Passed to \code{\link{nmfkc.signed}}; also accepts
#'   \code{nfolds} (default 5; \code{div} alias), \code{seed}
#'   (default 123), \code{shuffle} (default \code{TRUE}).
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.
#'
#' @return A list with \code{objfunc} (mean squared prediction error),
#'   \code{sigma} (RMSE), \code{objfunc.block} (per-fold MSE vector),
#'   \code{block} (integer fold assignment of length \eqn{N}).  Field
#'   names match \code{\link{nmfkc.cv}}.
#' @seealso \code{\link{nmfkc.signed}}, \code{\link{nmfkc.signed.ecv}}
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
#' semi-nonnegative matrix factorizations. \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, 32(1), 45--55.
#'
#' @export
nmfkc.signed.cv <- function(Y, A, rank = 2, ...) {
  extra <- list(...)
  if (!is.null(extra$Q)) rank <- extra$Q
  nfolds  <- if (!is.null(extra$nfolds)) extra$nfolds
             else if (!is.null(extra$div)) extra$div else 5
  seed    <- if (!is.null(extra$seed))    extra$seed    else 123
  shuffle <- if (!is.null(extra$shuffle)) extra$shuffle else TRUE

  Y <- as.matrix(Y); A <- as.matrix(A)
  N <- ncol(Y)
  if (ncol(A) != N) stop("ncol(A) must equal ncol(Y).")

  ## Strip CV-specific / user-overridable args
  fit_args <- extra
  fit_args$nfolds <- NULL; fit_args$div <- NULL
  fit_args$seed <- NULL; fit_args$shuffle <- NULL
  fit_args$Q <- NULL
  fit_args$verbose <- NULL; fit_args$Y.weights <- NULL

  ## Create folds over columns
  set.seed(seed)
  idx <- if (isTRUE(shuffle)) sample.int(N) else seq_len(N)
  chunk <- N %/% nfolds; rem <- N %% nfolds
  folds <- vector("list", nfolds); s <- 1L
  for (k in 1:nfolds) {
    sz <- chunk + ifelse(k <= rem, 1L, 0L)
    folds[[k]] <- idx[s:(s + sz - 1L)]; s <- s + sz
  }

  ## Per-fold: fit on train columns (via Y.weights = 0 on test columns),
  ## then score on held-out columns using Yhat = X C A
  objfunc.block <- numeric(nfolds)
  ## block: integer vector of length N assigning each column to a fold
  ## (matches main nmfkc.cv's output shape).
  block <- integer(N)
  for (k in 1:nfolds) block[folds[[k]]] <- k
  for (k in 1:nfolds) {
    test_cols <- folds[[k]]
    W <- matrix(1, nrow = nrow(Y), ncol = N)
    W[, test_cols] <- 0
    fit <- suppressMessages(do.call(
      nmfkc.signed,
      c(list(Y = Y, A = A, rank = rank, verbose = FALSE,
             Y.weights = W), fit_args)))
    Yhat_test <- fit$X %*% fit$C %*% A[, test_cols, drop = FALSE]
    objfunc.block[k] <- mean((Y[, test_cols, drop = FALSE] - Yhat_test)^2)
  }
  structure(list(
    objfunc = mean(objfunc.block),
    sigma = sqrt(mean(objfunc.block)),
    objfunc.block = objfunc.block,
    block = block
  ), class = c("nmfkc.signed.cv", "nmfkc.cv"))
}


#' Element-wise cross-validation for nmfkc.signed
#'
#' @description Element-wise k-fold CV (Wold's CV): held-out elements
#' are masked via \code{Y.weights = 0} during fitting, and the RMSE on
#' those elements is reported.  Loops over candidate \code{rank} values.
#'
#' @param Y Real-valued \eqn{Q_{\mathrm{obs}} \times N} response matrix
#'   (signed entries allowed).
#' @param A Real-valued \eqn{D \times N} covariate matrix (signed).
#' @param rank Integer vector of candidate ranks (default \code{1:3}).
#' @param ... Passed to \code{\link{nmfkc.signed}}; also accepts
#'   \code{nfolds} (default 5; \code{div} alias), \code{seed}
#'   (default 123).
#'
#' @section Lifecycle:
#' This function is \strong{experimental}.
#'
#' @return A list with \code{objfunc} (MSE per rank), \code{sigma}
#'   (RMSE), \code{objfunc.fold} (per-fold per-rank), \code{folds},
#'   \code{Q.grid}.
#' @seealso \code{\link{nmfkc.signed}}, \code{\link{nmfkc.signed.cv}}
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
#' semi-nonnegative matrix factorizations. \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, 32(1), 45--55.
#'
#' @export
nmfkc.signed.ecv <- function(Y, A, rank = 1:3, ...) {
  extra <- list(...)
  if (!is.null(extra$Q)) rank <- extra$Q
  nfolds <- if (!is.null(extra$nfolds)) extra$nfolds
            else if (!is.null(extra$div)) extra$div else 5
  seed   <- if (!is.null(extra$seed))   extra$seed   else 123

  Y <- as.matrix(Y); A <- as.matrix(A)
  P <- nrow(Y); N <- ncol(Y)
  if (ncol(A) != N) stop("ncol(A) must equal ncol(Y).")

  fit_args <- extra
  fit_args$nfolds <- NULL; fit_args$div <- NULL
  fit_args$seed <- NULL; fit_args$Q <- NULL
  fit_args$verbose <- NULL; fit_args$Y.weights <- NULL

  ## Create folds over valid elements (non-NA)
  set.seed(seed)
  valid <- which(!is.na(Y))
  perm  <- sample(valid)
  chunk <- length(perm) %/% nfolds; rem <- length(perm) %% nfolds
  folds <- vector("list", nfolds); s <- 1L
  for (k in 1:nfolds) {
    sz <- chunk + ifelse(k <= rem, 1L, 0L)
    folds[[k]] <- perm[s:(s + sz - 1L)]; s <- s + sz
  }

  run_one <- function(q, k) {
    test_idx <- folds[[k]]
    W <- matrix(1, nrow = P, ncol = N)
    if (any(is.na(Y))) W[is.na(Y)] <- 0
    W[test_idx] <- 0
    fit <- suppressMessages(do.call(
      nmfkc.signed,
      c(list(Y = Y, A = A, rank = q, verbose = FALSE,
             Y.weights = W), fit_args)))
    ## Yhat = X C A on held-out entries
    Yhat <- fit$X %*% fit$C %*% A
    mean((Y[test_idx] - Yhat[test_idx])^2)
  }

  result_objfunc <- numeric(length(rank))
  result_sigma   <- numeric(length(rank))
  result_fold    <- vector("list", length(rank))
  names(result_objfunc) <- sprintf("Q=%d", rank)
  names(result_sigma)   <- sprintf("Q=%d", rank)
  names(result_fold)    <- sprintf("Q=%d", rank)
  message(sprintf("nmfkc.signed ECV: %d ranks, %d-fold.",
                  length(rank), nfolds))
  for (i in seq_along(rank)) {
    q <- rank[i]
    objs <- numeric(nfolds)
    for (k in 1:nfolds) objs[k] <- run_one(q, k)
    result_fold[[i]]  <- objs
    result_objfunc[i] <- mean(objs)
    result_sigma[i]   <- sqrt(result_objfunc[i])
    message(sprintf("  Q=%d: MSE=%.6f, sigma=%.4f",
                    q, result_objfunc[i], result_sigma[i]))
  }
  structure(list(
    objfunc = result_objfunc, sigma = result_sigma,
    objfunc.fold = result_fold, folds = folds, Q.grid = rank
  ), class = c("nmfkc.signed.ecv", "nmfkc.ecv"))
}

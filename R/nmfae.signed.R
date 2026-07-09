# nmfae.signed.R
# Signed-Bottleneck NMF-AE: Three-layer NMF-AE with signed bottleneck
# Model:  Y1 ~= X1 (Cp - Cn) X2 Y2,  X1, Cp, Cn, X2, Y2 >= 0,  Y1 signed allowed
# Author: Kenichi Satoh
# Date: 2026-04-18
#
# Depends: nmfkc (for warm.start via nmfae())
# References:
#   Ding, C. H. Q., Li, T., & Jordan, M. I. (2010).
#     Convex and Semi-Nonnegative Matrix Factorizations. IEEE TPAMI 32(1), 45-55.
#   Satoh, K. (2026). Signed-Bottleneck NMF-AE: Signed-Bottleneck 3-Layer NMF
#     (research memo 2026-04-18).
#
# This file provides:
#   nmfae.signed()              -- main Direct-MU fit
#   predict.nmfae.signed()      -- prediction on new Y2
#   plot.nmfae.signed()         -- convergence plot
#   summary.nmfae.signed()      -- summary with Cp/Cn diagnostics
#   print.summary.nmfae.signed()
#   nmfae.signed.rename()       -- rename Resp/Cov labels (propagates to Cp, Cn, C)
# ==============================================================

#' @title Signed-Bottleneck NMF-AE: Three-Layer NMF-AE with Signed Bottleneck
#' @keywords internal
#' @description
#' \code{nmfae.signed} fits a three-layer non-negative matrix factorization
#' autoencoder with a \strong{signed} bottleneck, solving
#' \deqn{Y_1 \approx X_1 (C_{+} - C_{-}) X_2 Y_2,
#'   \quad X_1 \ge 0,\; C_{+} \ge 0,\; C_{-} \ge 0,\; X_2 \ge 0,}
#' where \eqn{\Theta = C_{+} - C_{-}} is the signed bottleneck.  The basis
#' matrices \eqn{X_1} (columns sum to 1) and \eqn{X_2} (rows sum to 1) retain
#' their non-negative "parts-based" interpretability, while \eqn{\Theta} can
#' express anti-correlations (e.g., refractive index up vs. Abbe number down).
#'
#' The algorithm uses \strong{Direct Multiplicative Updates} derived from
#' Ding et al. (2010) sign-splitting technique, applied block-wise to the
#' four non-negative blocks \eqn{(C_{+}, C_{-}, X_1, X_2)}.  Each block update
#' monotonically decreases the true objective
#' \eqn{\|Y_1 - X_1(C_{+} - C_{-})X_2 Y_2\|_F^2} (Lee-Seung auxiliary
#' function method).
#'
#' \strong{Relation to \code{\link{nmfae}}:} When \eqn{\Theta \ge 0} suffices
#' (the \code{nmfae} case), \code{nmfae.signed} reduces to \code{nmfae} up to
#' the \eqn{C_{+} - C_{-}} parameterization.  Use \code{nmfae.signed} when the
#' data exhibit negative cross-property correlations that tri-NMF-AE cannot
#' express (e.g., high refractive index <-> low Abbe number trade-off).
#'
#' @param Y1 Output matrix \eqn{Y_1} (P1 x N).  Negative entries allowed
#'   (\code{Y.signed = TRUE} is auto-detected).
#' @param Y2 Input matrix \eqn{Y_2} (P2 x N).  Must be non-negative.
#'   Default \code{Y1} (autoencoder).
#' @param rank1 Integer. Response-basis rank Q. Default 2.
#' @param rank2 Integer. Covariate-basis rank R. Default (\code{NULL}) = \code{rank1}.
#' @param rank,rank.encoder Deprecated aliases of \code{rank1} / \code{rank2}
#'   (\code{Q} / \code{R} also accepted via \code{...}).
#' @param epsilon Relative convergence tolerance on the objective.
#'   Default \code{1e-4}.
#' @param maxit Maximum iterations.  Default 5000.
#' @param verbose Logical.  Print progress.  Default \code{FALSE}.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{warm.start}}{One of \code{TRUE} (default, \strong{hybrid}:
#'       warm-start \eqn{X_1, X_2} from \code{\link{nmfae}} but initialize
#'       \eqn{C_{+}, C_{-}} randomly), \code{"full"} (warm-start everything
#'       including \eqn{C_{+} = C_{\mathrm{tri}}}, \eqn{C_{-} = \delta}),
#'       or \code{FALSE} (random for all blocks).  The hybrid default avoids
#'       the \eqn{C_{-} = 0} local-minimum trap inherited from tri-NMF-AE
#'       while still benefiting from good \eqn{X_1, X_2} initialization.
#'       Ignored when \eqn{Y_1} has negative entries.}
#'     \item{\code{nstart}}{Integer, default 1 (cf. \code{\link{kmeans}}).
#'       Number of random restarts for \eqn{C_{+}, C_{-}}.  Each restart
#'       uses seed \code{seed + 7919 * (s-1)}.  Returns the best run by
#'       final objective.  \strong{Signed models have more local minima
#'       than non-negative ones} because the bottleneck
#'       \eqn{\Theta = C_{+} - C_{-}} can take both positive and negative
#'       values; during exploration a larger \code{nstart} (e.g., 10-50)
#'       reduces the chance of being trapped at an inferior stationary
#'       point (particularly the \eqn{C_{-} = 0} trap from warm-start
#'       from non-negative tri-NMF-AE).  Use the default 1 for fast
#'       development and raise for publication-grade runs.}
#'     \item{\code{X1.L2.ortho}, \code{X2.L2.ortho}}{Non-negative L2
#'       orthogonality penalties (default 0) on the \strong{columns} of
#'       \eqn{X_1} and the \strong{rows} of \eqn{X_2}, penalizing
#'       \eqn{(\lambda/2)\lVert\mathrm{offdiag}(X_1^\top X_1)\rVert^2} and
#'       \eqn{(\lambda/2)\lVert\mathrm{offdiag}(X_2 X_2^\top)\rVert^2}
#'       respectively.  Same convention as \code{\link{nmfae}}; encourage
#'       more distinct (less overlapping) response / covariate bases.}
#'     \item{\code{Y1.weights}}{Optional non-negative weight matrix
#'       (P1 x N) or vector (length N) for \eqn{Y_1}, analogous to the
#'       \code{weights} argument of \code{\link[stats]{lm}}.  Loss
#'       becomes \eqn{\sum W_{ij} \, (Y_{1,ij} - \hat Y_{1,ij})^2}
#'       (\code{lm()}-style, \strong{linear} in \eqn{W}).  Logical
#'       matrices (\code{TRUE} / \code{FALSE}) are also accepted.
#'       Used by \code{\link{nmfae.signed.ecv}} to hold out test
#'       elements via a binary mask \eqn{W \in \{0,1\}}; real-valued
#'       weights for importance weighting are also supported.  Default:
#'       if \code{Y1} has \code{NA}, a binary mask is auto-generated
#'       (0 for \code{NA}, 1 elsewhere).}
#'     \item{\code{Cp.init}, \code{Cn.init}}{Explicit \eqn{Q \times R}
#'       non-negative matrices for initialization.  Overrides warm.start.}
#'     \item{\code{C.init}}{Explicit signed \eqn{Q \times R} matrix,
#'       internally split into \eqn{(C_{+}, C_{-})}.}
#'     \item{\code{X1.init}, \code{X2.init}}{Explicit basis matrices.}
#'     \item{\code{seed}}{RNG seed.  Default 123.}
#'     \item{\code{print.trace}}{Logical.  Print iteration trace.
#'       Default \code{FALSE}.}
#'     \item{\code{prefix.dec}, \code{prefix.enc}}{Label prefixes for the
#'       response/covariate bases.  Default \code{"Resp"}, \code{"Cov"}.}
#'   }
#'
#' @return An object of class \code{c("nmfae.signed", "nmfae", "nmf")} with:
#' \item{X1}{Decoder basis (P1 x Q), column sum 1.}
#' \item{Cp, Cn}{Non-negative parts of \eqn{\Theta} (each Q x R).}
#' \item{C}{Signed bottleneck \eqn{\Theta = C_{+} - C_{-}} (Q x R).}
#' \item{X2}{Encoder basis (R x P2), row sum 1.}
#' \item{Y1hat}{Fitted values \eqn{X_1 (C_{+} - C_{-}) X_2 Y_2}.}
#' \item{H}{Encoding \eqn{(C_{+} - C_{-}) X_2 Y_2} (Q x N, signed).}
#' \item{rank}{\code{c(Q = Q, R = R)}.}
#' \item{dims}{\code{c(P1, P2, N)}.}
#' \item{objfunc, objfunc.iter}{Final and per-iteration objective values.}
#' \item{r.squared}{\eqn{\mathrm{cor}(Y_1, \widehat Y_1)^2} (Pearson; in \eqn{[0,1]}).}
#' \item{r.squared.uncentered}{Uncentered \eqn{R^2 = 1 - \|Y_1 - \widehat Y_1\|_F^2 / \|Y_1\|_F^2} (baseline = zero matrix).}
#' \item{r.squared.centered}{Row-mean centered \eqn{1 - \|Y_1 - \widehat Y_1\|_F^2 / \|Y_1 - \bar Y_{p\cdot}\|_F^2}.}
#' \item{sigma, mae}{Residual SE and mean absolute error.}
#' \item{niter, runtime}{Iterations and elapsed seconds.}
#' \item{Y.signed}{Logical; whether \eqn{Y_1} contained negative entries.}
#' \item{call}{Matched call.}
#'
#' @section Lifecycle:
#' This function is \strong{experimental}; interface may change.
#'
#' @seealso \code{\link{nmfae}}, \code{\link{predict.nmfae.signed}},
#'   \code{\link{summary.nmfae.signed}}, \code{\link{nmfae.signed.rename}}
#'
#' @references
#' Ding, C.H.Q., Li, T., and Jordan, M.I. (2010).  Convex and
#' Semi-Nonnegative Matrix Factorizations.
#' \emph{IEEE TPAMI}, 32(1), 45-55.
#'
#' Satoh, K. (2026).  Signed-Bottleneck NMF-AE: Signed-Bottleneck 3-Layer NMF
#' (research memo, 2026-04-18).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' Y1 <- matrix(abs(rnorm(12)), 3, 4)
#' Y2 <- matrix(abs(rnorm(20)), 5, 4)
#' res <- nmfae.signed(Y1, Y2, rank1 = 2, rank2 = 2, maxit = 500)
#' summary(res)
#' }
#' @references
#' Ding, C. H. Q., Li, T., & Jordan, M. I. (2010). Convex and
#' semi-nonnegative matrix factorizations. \emph{IEEE Transactions on
#' Pattern Analysis and Machine Intelligence}, 32(1), 45--55.
#'
#' @export
nmfae.signed <- function(Y1, Y2 = Y1, rank1 = 2, rank2 = NULL,
                       epsilon = 1e-4, maxit = 5000,
                       verbose = FALSE, ...,
                       rank = NULL, rank.encoder = NULL) {

  cl <- match.call()
  extra_args <- list(...)

  ## ---- 1. Parameter extraction (compatibility with nmfae) ----
  ## rank1 = response basis X1, rank2 = covariate basis X2 (default rank1);
  ## legacy rank / rank.encoder (formals) and Q / R (via ...)
  if (is.null(rank))          rank <- rank1
  if (!is.null(extra_args$Q)) rank <- extra_args$Q
  if (is.null(rank.encoder))  rank.encoder <- rank2
  if (!is.null(extra_args$R)) rank.encoder <- extra_args$R
  if (is.null(rank.encoder))  rank.encoder <- rank
  Q <- as.integer(rank)
  R <- as.integer(rank.encoder)

  warm.start <- if (!is.null(extra_args$warm.start)) extra_args$warm.start else TRUE
  nstart  <- if (!is.null(extra_args$nstart))  extra_args$nstart  else 1L
  ## Basis-init method forwarded to the nmfae() warm-start step (default
  ## "kmeans"; "kmeans++" etc. accepted). String methods only.
  X.init  <- if (!is.null(extra_args$X.init))  extra_args$X.init  else "kmeans"
  X.init.method <- if (is.character(X.init)) X.init else "kmeans"
  Y1.weights <- if (!is.null(extra_args$Y1.weights)) extra_args$Y1.weights else NULL
  Cp.init    <- extra_args$Cp.init
  Cn.init    <- extra_args$Cn.init
  C.init     <- extra_args$C.init
  X1.init    <- extra_args$X1.init
  X2.init    <- extra_args$X2.init
  seed       <- if (!is.null(extra_args$seed))       extra_args$seed       else 123L
  print.trace <- verbose
  if (!is.null(extra_args$print.trace)) print.trace <- extra_args$print.trace
  prefix.dec <- if (!is.null(extra_args$prefix.dec)) extra_args$prefix.dec else "Resp"
  prefix.enc <- if (!is.null(extra_args$prefix.enc)) extra_args$prefix.enc else "Cov"
  ## Basis orthogonality penalties (same convention as nmfae(): off-diagonal
  ## L2 on X1 columns and X2 rows; both default off).
  X1.L2.ortho <- if (!is.null(extra_args$X1.L2.ortho)) extra_args$X1.L2.ortho else 0
  X2.L2.ortho <- if (!is.null(extra_args$X2.L2.ortho)) extra_args$X2.L2.ortho else 0

  ## Validate warm.start value
  if (is.logical(warm.start)) {
    warm_mode <- if (isTRUE(warm.start)) "hybrid" else "random"
  } else {
    warm_mode <- match.arg(as.character(warm.start),
                            c("hybrid", "full", "random"))
  }

  ## ---- 2. Input preparation & validation ----
  ## NA in Y1 is auto-masked via Y1.weights (matches nmfae() behavior);
  ## NA in Y2 is still an error since Y2 is not weighted.
  Y1 <- as.matrix(Y1); storage.mode(Y1) <- "double"
  Y2 <- as.matrix(Y2); storage.mode(Y2) <- "double"
  if (any(is.na(Y2))) stop("Y2 contains NA; please impute or remove.")
  if (min(Y2) < 0) stop("Y2 must be non-negative.")
  if (ncol(Y1) != ncol(Y2)) stop("Y1 and Y2 must have the same number of columns (N).")

  P1 <- nrow(Y1); P2 <- nrow(Y2); N <- ncol(Y1)
  Y1_signed <- any(Y1 < 0)
  small <- 1e-16
  eps_init <- 1e-6
  start.time <- Sys.time()

  if (print.trace) {
    message(sprintf(
      "Y1(%d,%d) ~ X1(%d,%d) [Cp(%d,%d)-Cn(%d,%d)] X2(%d,%d) Y2(%d,%d)%s",
      P1, N, P1, Q, Q, R, Q, R, R, P2, P2, N,
      if (Y1_signed) "  [Y1 signed]" else ""))
  }

  ## ---- 3a. Y1.weights handling (for ECV / missing-element masking) ----
  ## Vector of length N -> per-column; scalar -> all; NULL -> auto-detect NAs
  if (!is.null(Y1.weights) && is.vector(Y1.weights) && !is.matrix(Y1.weights)) {
    if (length(Y1.weights) == N) {
      Y1.weights <- matrix(Y1.weights, nrow = P1, ncol = N, byrow = TRUE)
    } else if (length(Y1.weights) == 1) {
      Y1.weights <- matrix(Y1.weights, nrow = P1, ncol = N)
    } else {
      stop("Length of Y1.weights vector must match ncol(Y1) or be 1.")
    }
  }
  if (is.null(Y1.weights) && any(is.na(Y1))) {
    Y1.weights <- matrix(1, nrow = P1, ncol = N)
    Y1.weights[is.na(Y1)] <- 0
    Y1[is.na(Y1)] <- 0
  } else if (!is.null(Y1.weights)) {
    Y1.weights <- as.matrix(Y1.weights)
    if (!all(dim(Y1.weights) == dim(Y1)))
      stop("Y1.weights dimensions must match Y1.")
    Y1.weights[is.na(Y1.weights)] <- 0
    Y1[is.na(Y1) | Y1.weights == 0] <- 0
  }
  has.weights <- !is.null(Y1.weights) && any(Y1.weights != 1)
  if (!has.weights) Y1.weights <- NULL
  Wmat <- Y1.weights  # short alias for weighted path

  ## ---- 3b. Pre-computation (unweighted path only; N-dependent, outer-loop) ----
  if (!has.weights) {
    S  <- tcrossprod(Y2)             # P2 x P2, >= 0 (since Y2 >= 0)
    G0 <- tcrossprod(Y1, Y2)         # P1 x P2 (= Y1 Y2^T), signed iff Y1 signed
  } else {
    S <- NULL; G0 <- NULL            # weighted path skips these precomputes
  }
  Y1_sqnorm <- sum(Y1 * Y1)

  ## ---- 4. Initialization ----
  ## Strategy:
  ##   hybrid: X1/X2 from tri-NMF-AE, Cp/Cn random (DEFAULT, best for Abbe-like data)
  ##   full:   X1/X2 from tri-NMF-AE, Cp = tri$C, Cn = eps (old tri-NMF solution)
  ##   random: all blocks random (slowest warm-up, widest exploration)
  X1 <- NULL; X2 <- NULL; Cp <- NULL; Cn <- NULL

  use_explicit_X <- !is.null(X1.init) || !is.null(X2.init)
  use_explicit_C <- !is.null(Cp.init) || !is.null(Cn.init) || !is.null(C.init)

  ## 4a. X1, X2: warm-start from tri-NMF-AE unless mode=random or explicit override
  tri_C_for_full <- NULL
  if (warm_mode %in% c("hybrid", "full") && !Y1_signed && !use_explicit_X) {
    if (print.trace) message("  Init: warm-start X1, X2 from nmfae() ...")
    res0 <- nmfae(Y1, Y2, rank1 = Q, rank2 = R,
                          epsilon = epsilon, maxit = maxit,
                          verbose = FALSE, seed = seed,
                          X.init = X.init.method)
    X1 <- res0$X1
    X2 <- res0$X2
    tri_C_for_full <- res0$C          # for warm_mode == "full"
  } else {
    if (!is.null(X1.init)) {
      X1 <- as.matrix(X1.init)
      if (!identical(dim(X1), c(P1, Q)))
        stop("X1.init must have dimensions (nrow(Y1), rank).")
    } else {
      set.seed(seed)
      X1 <- matrix(abs(stats::rnorm(P1 * Q)) * 0.1, P1, Q)
    }
    if (!is.null(X2.init)) {
      X2 <- as.matrix(X2.init)
      if (!identical(dim(X2), c(R, P2)))
        stop("X2.init must have dimensions (rank.encoder, nrow(Y2)).")
    } else {
      set.seed(seed + 1L)
      X2 <- matrix(abs(stats::rnorm(R * P2)) * 0.1, R, P2)
    }
  }

  ## 4b. Cp, Cn: random (hybrid / random modes), tri (full mode), or explicit
  init_Cp_Cn <- function(s) {
    if (!is.null(Cp.init) || !is.null(Cn.init)) {
      if (is.null(Cp.init) || is.null(Cn.init))
        stop("Provide both Cp.init and Cn.init, or neither.")
      cp <- as.matrix(Cp.init); cn <- as.matrix(Cn.init)
      if (!identical(dim(cp), c(Q, R)) || !identical(dim(cn), c(Q, R)))
        stop("Cp.init and Cn.init must have dimensions (Q, R).")
      if (any(cp < 0) || any(cn < 0))
        stop("Cp.init, Cn.init must be non-negative.")
      return(list(Cp = cp, Cn = cn))
    }
    if (!is.null(C.init)) {
      C0 <- as.matrix(C.init)
      if (!identical(dim(C0), c(Q, R)))
        stop("C.init must have dimensions (Q, R).")
      return(list(Cp = pmax(C0, 0), Cn = pmax(-C0, 0)))
    }
    if (warm_mode == "full" && !is.null(tri_C_for_full)) {
      return(list(Cp = tri_C_for_full,
                  Cn = matrix(eps_init, Q, R)))
    }
    ## Sample signed C from U(-1, 1), split into Cp, Cn.  Scale is arbitrary
    ## because the multiplicative updates auto-rescale; non-zero + sign-
    ## symmetric diversity is what matters.
    set.seed(s)
    C0 <- matrix(stats::runif(Q * R, min = -1, max = 1), Q, R)
    list(Cp = pmax(C0, 0), Cn = pmax(-C0, 0))
  }

  ## ---- 5. Objective (weighted if Wmat present) ----
  compute_obj <- function(X1, Cp, Cn, X2) {
    Y1hat <- X1 %*% (Cp - Cn) %*% X2 %*% Y2
    ## lm()-style weighted least squares: L = sum(W * (Y1 - Y1hat)^2).
    ## W = W^2 for binary {0,1} masks so this is unchanged for ECV.
    if (has.weights) sum(Wmat * (Y1 - Y1hat)^2) else sum((Y1 - Y1hat)^2)
  }

  ## Basis-orthogonality penalty value and MU denominator contributions
  ## (same convention as nmfae()): (lambda/2)||offdiag(X1'X1)||^2 for X1
  ## columns and (lambda/2)||offdiag(X2 X2')||^2 for X2 rows; the positive
  ## gradient X1 offdiag(X1'X1) / offdiag(X2 X2') X2 goes to the denominator.
  pen_X12 <- function(X1, X2) {
    p <- 0
    if (X1.L2.ortho > 0) { G <- crossprod(X1);  diag(G) <- 0; p <- p + (X1.L2.ortho / 2) * sum(G^2) }
    if (X2.L2.ortho > 0) { G <- tcrossprod(X2); diag(G) <- 0; p <- p + (X2.L2.ortho / 2) * sum(G^2) }
    p
  }
  den_X1_pen <- function(X1) {
    if (X1.L2.ortho > 0) { G <- crossprod(X1); diag(G) <- 0; X1.L2.ortho * (X1 %*% G) } else 0
  }
  den_X2_pen <- function(X2) {
    if (X2.L2.ortho > 0) { G <- tcrossprod(X2); diag(G) <- 0; X2.L2.ortho * (G %*% X2) } else 0
  }

  ## ---- 6. Main loop wrapped for multi-restart ----
  run_once <- function(X1, X2, Cp, Cn) {
    ## Normalize
    cs <- colSums(X1) + small
    X1 <- sweep(X1, 2, cs, "/"); Cp <- sweep(Cp, 1, cs, "*"); Cn <- sweep(Cn, 1, cs, "*")
    rs <- rowSums(X2) + small
    X2 <- sweep(X2, 1, rs, "/"); Cp <- sweep(Cp, 2, rs, "*"); Cn <- sweep(Cn, 2, rs, "*")
    obj_prev <- compute_obj(X1, Cp, Cn, X2) + pen_X12(X1, X2)
    objfunc.iter <- numeric(maxit)
    iter <- 0L
    for (iter in seq_len(maxit)) {

    if (!has.weights) {
      ## ---- Fast unweighted path (uses precomputed S, G0) ----
      P1m <- crossprod(X1)
      SX  <- X2 %*% S %*% t(X2)
      GX  <- crossprod(X1, G0) %*% t(X2)
      GX_p <- pmax(GX,  0); GX_n <- pmax(-GX, 0)
      ## Cp update
      Cp <- Cp * (GX_p + P1m %*% Cn %*% SX) /
                 (GX_n + P1m %*% Cp %*% SX + small)
      ## Cn update (Gauss-Seidel)
      Cn <- Cn * (GX_n + P1m %*% Cp %*% SX) /
                 (GX_p + P1m %*% Cn %*% SX + small)
      ## X1 update
      G0X <- G0 %*% t(X2)
      G0X_p <- pmax(G0X, 0); G0X_n <- pmax(-G0X, 0)
      MMt_p <- Cp %*% SX %*% t(Cp) + Cn %*% SX %*% t(Cn)
      MMt_n <- Cp %*% SX %*% t(Cn) + Cn %*% SX %*% t(Cp)
      X1 <- X1 * (G0X_p %*% t(Cp) + G0X_n %*% t(Cn) + X1 %*% MMt_n) /
                 (G0X_n %*% t(Cp) + G0X_p %*% t(Cn) + X1 %*% MMt_p + den_X1_pen(X1) + small)
      ## Normalize X1 cols -> absorb into Cp, Cn rows
      cs <- colSums(X1) + small
      X1 <- sweep(X1, 2, cs, "/")
      Cp <- sweep(Cp, 1, cs, "*"); Cn <- sweep(Cn, 1, cs, "*")
      ## X2 update
      Amat <- crossprod(X1, G0)
      A_p <- pmax(Amat, 0); A_n <- pmax(-Amat, 0)
      HG_p <- crossprod(Cp, A_p) + crossprod(Cn, A_n)
      HG_n <- crossprod(Cp, A_n) + crossprod(Cn, A_p)
      P1m <- crossprod(X1)
      HtH_p <- crossprod(Cp, P1m) %*% Cp + crossprod(Cn, P1m) %*% Cn
      HtH_n <- crossprod(Cp, P1m) %*% Cn + crossprod(Cn, P1m) %*% Cp
      X2 <- X2 * (HG_p + HtH_n %*% X2 %*% S) /
                 (HG_n + HtH_p %*% X2 %*% S + den_X2_pen(X2) + small)
      ## Normalize X2 rows -> absorb into Cp, Cn cols
      rs <- rowSums(X2) + small
      X2 <- sweep(X2, 1, rs, "/")
      Cp <- sweep(Cp, 2, rs, "*"); Cn <- sweep(Cn, 2, rs, "*")

    } else {
      ## ---- Weighted path (no S/G0 precompute; ECV / missing data) ----
      F_mat <- X2 %*% Y2                                 # R x N
      ## Cp update
      XCpF <- X1 %*% Cp %*% F_mat                        # P1 x N, >= 0
      XCnF <- X1 %*% Cn %*% F_mat                        # P1 x N, >= 0
      WY1  <- Wmat * Y1
      WXCpF <- Wmat * XCpF
      WXCnF <- Wmat * XCnF
      G_w  <- crossprod(X1, WY1) %*% t(F_mat)            # Q x R, signed iff Y1 signed
      Hp_w <- crossprod(X1, WXCpF) %*% t(F_mat)          # Q x R, >= 0
      Hn_w <- crossprod(X1, WXCnF) %*% t(F_mat)          # Q x R, >= 0
      Gp <- pmax(G_w, 0); Gn <- pmax(-G_w, 0)
      Cp <- Cp * (Gp + Hn_w) / (Gn + Hp_w + small)
      ## Cn update (recompute Hp_w with new Cp)
      XCpF  <- X1 %*% Cp %*% F_mat
      WXCpF <- Wmat * XCpF
      Hp_w  <- crossprod(X1, WXCpF) %*% t(F_mat)
      Cn <- Cn * (Gn + Hp_w) / (Gp + Hn_w + small)
      ## X1 update
      Mp <- Cp %*% F_mat; Mn <- Cn %*% F_mat
      XMp <- X1 %*% Mp; XMn <- X1 %*% Mn
      WY1  <- Wmat * Y1
      WXMp <- Wmat * XMp; WXMn <- Wmat * XMn
      W1Mp <- WY1 %*% t(Mp); W1Mn <- WY1 %*% t(Mn)
      W1Mp_p <- pmax(W1Mp, 0); W1Mp_n <- pmax(-W1Mp, 0)
      W1Mn_p <- pmax(W1Mn, 0); W1Mn_n <- pmax(-W1Mn, 0)
      P_pp <- WXMp %*% t(Mp); P_pn <- WXMp %*% t(Mn)
      P_np <- WXMn %*% t(Mp); P_nn <- WXMn %*% t(Mn)
      X1 <- X1 * (W1Mp_p + W1Mn_n + P_pn + P_np) /
                 (W1Mp_n + W1Mn_p + P_pp + P_nn + den_X1_pen(X1) + small)
      cs <- colSums(X1) + small
      X1 <- sweep(X1, 2, cs, "/")
      Cp <- sweep(Cp, 1, cs, "*"); Cn <- sweep(Cn, 1, cs, "*")
      ## X2 update: Y1 ~= H X2 Y2, H = X1 Cp - X1 Cn
      Hp <- X1 %*% Cp; Hn <- X1 %*% Cn  # both >= 0
      HpF <- Hp %*% X2 %*% Y2; HnF <- Hn %*% X2 %*% Y2   # >= 0
      WHpF <- Wmat * HpF; WHnF <- Wmat * HnF
      WY1 <- Wmat * Y1
      Hp_tW1 <- crossprod(Hp, WY1) %*% t(Y2)              # R x P2, signed iff Y1 signed
      Hn_tW1 <- crossprod(Hn, WY1) %*% t(Y2)
      Hp_tW1_p <- pmax(Hp_tW1, 0); Hp_tW1_n <- pmax(-Hp_tW1, 0)
      Hn_tW1_p <- pmax(Hn_tW1, 0); Hn_tW1_n <- pmax(-Hn_tW1, 0)
      Q_pp <- crossprod(Hp, WHpF) %*% t(Y2)
      Q_pn <- crossprod(Hp, WHnF) %*% t(Y2)
      Q_np <- crossprod(Hn, WHpF) %*% t(Y2)
      Q_nn <- crossprod(Hn, WHnF) %*% t(Y2)
      X2 <- X2 * (Hp_tW1_p + Hn_tW1_n + Q_pn + Q_np) /
                 (Hp_tW1_n + Hn_tW1_p + Q_pp + Q_nn + den_X2_pen(X2) + small)
      rs <- rowSums(X2) + small
      X2 <- sweep(X2, 1, rs, "/")
      Cp <- sweep(Cp, 2, rs, "*"); Cn <- sweep(Cn, 2, rs, "*")
    }

    ## Objective and convergence
    obj_cur <- compute_obj(X1, Cp, Cn, X2) + pen_X12(X1, X2)
    objfunc.iter[iter] <- obj_cur

    if (print.trace && (iter %% 100 == 0 || iter == 1)) {
      message(sprintf("  iter %5d: objfunc = %.6f", iter, obj_cur))
    }
    if (iter > 1) {
      rel <- abs(obj_prev - obj_cur) / (abs(obj_prev) + small)
      if (rel < epsilon) {
        if (print.trace) message(sprintf("  Converged at iter %d", iter))
        break
      }
    }
    obj_prev <- obj_cur
    }  # end for iter
    ## Warn when the MU loop exhausts maxit without convergence
    ## (matches nmfkc() / nmf.sem() convention).  Each restart inside
    ## the multi-start loop emits its own warning if it fails to
    ## converge, so users see how many restarts hit the cap.
    if (iter == maxit && exists("rel") && rel >= epsilon)
      warning(paste0("maximum iterations (", maxit, ") reached..."))
    list(X1 = X1, X2 = X2, Cp = Cp, Cn = Cn,
         objfunc.iter = objfunc.iter[seq_len(iter)], niter = iter,
         objfunc = obj_prev)
  }  # end run_once

  ## Execute run(s)
  best <- NULL
  for (s in seq_len(nstart)) {
    ## Wide seed spacing via a large prime (7919) to maximize init diversity
    s_seed <- seed + 7919L * (s - 1L)
    cc <- init_Cp_Cn(s_seed + 2L)
    if (print.trace && nstart > 1)
      message(sprintf("  Restart %d/%d (seed=%d) ...", s, nstart, s_seed))
    ## For random mode, also regenerate X1, X2 per restart
    Xs1 <- X1; Xs2 <- X2
    if (warm_mode == "random" && !use_explicit_X && nstart > 1) {
      set.seed(s_seed)
      Xs1 <- matrix(abs(stats::rnorm(P1 * Q)) * 0.1, P1, Q)
      set.seed(s_seed + 1L)
      Xs2 <- matrix(abs(stats::rnorm(R * P2)) * 0.1, R, P2)
    }
    out <- run_once(Xs1, Xs2, cc$Cp, cc$Cn)
    if (is.null(best) || out$objfunc < best$objfunc) best <- out
  }
  X1 <- best$X1; X2 <- best$X2; Cp <- best$Cp; Cn <- best$Cn
  objfunc.iter <- best$objfunc.iter
  iter <- best$niter
  niter <- iter
  diff.time <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))

  ## ---- 7. Post-processing ----
  ## Reorder X1 columns by centroid of row indices
  if (Q > 1) {
    idx1 <- order(matrix(seq_len(P1) / P1, nrow = 1) %*% X1)
    X1 <- X1[, idx1, drop = FALSE]
    Cp <- Cp[idx1, , drop = FALSE]
    Cn <- Cn[idx1, , drop = FALSE]
  }
  ## Reorder X2 rows by centroid of column indices
  if (R > 1) {
    idx2 <- order(X2 %*% matrix(seq_len(P2) / P2, ncol = 1))
    X2 <- X2[idx2, , drop = FALSE]
    Cp <- Cp[, idx2, drop = FALSE]
    Cn <- Cn[, idx2, drop = FALSE]
  }

  ## Names
  rownames(X1) <- rownames(Y1)
  colnames(X2) <- rownames(Y2)
  colnames(X1) <- paste0(prefix.dec, seq_len(Q))
  rownames(Cp) <- paste0(prefix.dec, seq_len(Q))
  rownames(Cn) <- paste0(prefix.dec, seq_len(Q))
  colnames(Cp) <- paste0(prefix.enc, seq_len(R))
  colnames(Cn) <- paste0(prefix.enc, seq_len(R))
  rownames(X2) <- paste0(prefix.enc, seq_len(R))

  C <- Cp - Cn                                  # Q x R, signed
  H <- C %*% X2 %*% Y2                          # Q x N, signed encoding
  Y1hat <- X1 %*% H
  resid <- Y1 - Y1hat

  ## Goodness-of-fit statistics.  lm()-style weighted least squares:
  ## the reported objfunc / sigma / mae use sum(W * resid^2) to match the
  ## in-loop objective compute_obj().  For binary {0,1} masks (the
  ## standard ECV / CV / NA-mask case) this restricts statistics to valid
  ## elements; for real-valued weights it returns weighted averages.
  if (has.weights) {
    objfunc  <- sum(Wmat * resid^2)
    valid    <- (Wmat > 0)
    n.valid  <- sum(valid)
    r2_all   <- .r.squared.all(Y1, Y1hat, Y.weights = Wmat)
    sigma    <- if (n.valid > 0) sqrt(objfunc / n.valid) else NA_real_
    mae      <- if (sum(Wmat) > 0) sum(Wmat * abs(resid)) / sum(Wmat) else NA_real_
  } else {
    objfunc  <- sum(resid * resid)
    r2_all   <- .r.squared.all(Y1, Y1hat)
    sigma    <- sqrt(objfunc / (P1 * N))
    mae      <- mean(abs(resid))
  }
  r.squared          <- r2_all$r.squared
  r.squared.uncentered     <- r2_all$r.squared.uncentered
  r.squared.centered <- r2_all$r.squared.centered

  ## Soft/hard clustering of encoding (only meaningful when H has interpretable sign)
  eps_bp <- 1e-16
  ## Clip to non-negative for clustering probabilities
  Hp <- pmax(H, 0)
  col_mass <- colSums(Hp) + eps_bp
  B.prob <- sweep(Hp, 2, col_mass, "/")
  B.cluster <- apply(B.prob, 2, which.max)
  B.cluster[colSums(Hp) == 0] <- NA

  if (print.trace) {
    message(sprintf("  Done: %d iterations, %.1f sec, R2 = %.4f",
                    niter, diff.time, r.squared))
  }

  result <- list(
    call = cl,
    X1 = X1,
    Cp = Cp,
    Cn = Cn,
    C  = C,
    X2 = X2,
    Y1hat = Y1hat,
    H = H,
    B.prob = B.prob,
    B.cluster = B.cluster,
    rank = c(Q = Q, R = R),
    dims = c(P1 = P1, P2 = P2, N = N),
    objfunc = objfunc,
    objfunc.iter = objfunc.iter,
    r.squared          = r.squared,
    r.squared.uncentered     = r.squared.uncentered,
    r.squared.centered = r.squared.centered,
    sigma = sigma,
    mae = mae,
    niter = niter,
    iter = niter,          # house-style alias (matches nmfre/nmf.sem/nmfkc.net)
    runtime = diff.time,
    n.missing = if (has.weights) sum(Wmat == 0) else 0L,
    n.total = P1 * N,
    Y.signed = Y1_signed
  )
  class(result) <- c("nmfae.signed", "nmfae", "nmf")
  result
}


## ==============================================================
#' @title Statistical Inference for Signed-Bottleneck NMF-AE Signed Bottleneck
#' @keywords internal
#' @description
#' Post-estimation inference for the \strong{signed} bottleneck
#' \eqn{\Theta = C_{+} - C_{-}} in the Signed-Bottleneck NMF-AE model
#' \eqn{Y_1 \approx X_1 \Theta X_2 Y_2}, conditional on
#' \eqn{(\hat X_1, \hat X_2)}.  Uses sandwich covariance and wild bootstrap
#' \strong{without} the non-negativity projection that \code{\link{nmfae.inference}}
#' applies (because \eqn{\Theta} is unconstrained in sign here).
#'
#' @param object A fitted \code{"nmfae.signed"} object.
#' @param Y1 Output matrix used during fitting.
#' @param Y2 Input matrix used during fitting.  Default \code{Y1}.
#' @param wild.bootstrap Logical.  Default \code{TRUE}.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{wild.B}}{Bootstrap replicates.  Default 500.}
#'     \item{\code{wild.seed}}{RNG seed.  Default 123.}
#'     \item{\code{wild.level}}{CI confidence level.  Default 0.95.}
#'     \item{\code{sandwich}}{Use sandwich covariance.  Default \code{TRUE}.}
#'     \item{\code{C.p.side}}{P-value type: \code{"two.sided"} (default for
#'       Signed-Bottleneck NMF-AE) or \code{"one.sided"}.}
#'     \item{\code{cov.ridge}}{Ridge stabilization.  Default \code{1e-8}.}
#'     \item{\code{print.trace}}{Logical.  Default \code{FALSE}.}
#'   }
#'
#' @return An object of class \code{c("nmfae.signed.inference",
#'   "nmfae.inference", "nmfae.signed", "nmfae", "nmf")} with added fields:
#' \item{sigma2.used}{Estimated \eqn{\sigma^2}.}
#' \item{C.se, C.se.boot}{Sandwich / bootstrap SEs for \eqn{\Theta} (Q x R).}
#' \item{C.ci.lower, C.ci.upper}{Bootstrap CIs.}
#' \item{coefficients}{Data frame with Estimate, SE, BSE, z, p-value, CI.}
#' \item{C.p.side}{P-value side used.}
#'
#' @seealso \code{\link{nmfae.signed}}, \code{\link{nmfae.inference}}
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
nmfae.signed.inference <- function(object, Y1, Y2 = Y1,
                                 wild.bootstrap = TRUE, ...) {
  if (!inherits(object, "nmfae.signed"))
    stop("object must be of class 'nmfae.signed'")

  extra_args <- list(...)
  wild.B     <- if (!is.null(extra_args$wild.B))     extra_args$wild.B     else 500
  wild.seed  <- if (!is.null(extra_args$wild.seed))  extra_args$wild.seed  else 123
  wild.level <- if (!is.null(extra_args$wild.level)) extra_args$wild.level else 0.95
  sandwich   <- if (!is.null(extra_args$sandwich))   extra_args$sandwich   else TRUE
  C.p.side   <- if (!is.null(extra_args$C.p.side))   extra_args$C.p.side   else "two.sided"
  cov.ridge  <- if (!is.null(extra_args$cov.ridge))  extra_args$cov.ridge  else 1e-8
  print.trace <- if (!is.null(extra_args$print.trace)) extra_args$print.trace else FALSE

  X1 <- object$X1
  C  <- object$C                      # signed Theta = Cp - Cn
  X2 <- object$X2
  Y1hat <- object$Y1hat
  Q  <- unname(object$rank["Q"])
  R  <- unname(object$rank["R"])
  P1 <- nrow(Y1)
  N  <- ncol(Y1)

  Z   <- X2 %*% Y2
  R_C <- Y1 - Y1hat

  denom <- max(P1 * N - Q * R, 1)
  sigma2.used <- sum(R_C^2) / denom

  X1tX1 <- crossprod(X1)
  ZZt   <- tcrossprod(Z)
  Info_core <- kronecker(ZZt, X1tX1)
  Info <- Info_core / max(sigma2.used, 1e-12)
  Info <- Info + diag(cov.ridge, nrow(Info))

  Hinv <- tryCatch(solve(Info), error = function(e) {
    if (requireNamespace("MASS", quietly = TRUE)) MASS::ginv(Info)
    else stop("Information matrix singular; install MASS package.")
  })

  V_sand <- NULL
  if (isTRUE(sandwich)) {
    X1t <- t(X1)
    J <- matrix(0, Q * R, Q * R)
    for (n in 1:N) {
      z_n <- Z[, n, drop = FALSE]
      r_n <- R_C[, n, drop = FALSE]
      g_n <- X1t %*% r_n
      S_n <- -(g_n %*% t(z_n)) / max(sigma2.used, 1e-12)
      s_n <- as.vector(S_n)
      J <- J + tcrossprod(s_n)
    }
    if (N > 1) J <- (N / (N - 1)) * J
    V_sand <- Hinv %*% J %*% Hinv
  }
  C.vec.cov <- if (!is.null(V_sand)) V_sand else Hinv
  se_vec <- sqrt(pmax(diag(C.vec.cov), 0))
  C.se <- matrix(se_vec, nrow = Q, ncol = R, byrow = FALSE)

  ## Wild bootstrap (NO non-negative projection; Theta is signed)
  C.se.boot <- NULL; C.ci.lower <- NULL; C.ci.upper <- NULL
  if (isTRUE(wild.bootstrap)) {
    set.seed(wild.seed)
    X1t <- t(X1)
    score_mat <- matrix(0, Q * R, N)
    for (n in 1:N) {
      z_n <- Z[, n, drop = FALSE]
      r_n <- R_C[, n, drop = FALSE]
      g_n <- X1t %*% r_n
      G_n <- -(g_n %*% t(z_n)) / max(sigma2.used, 1e-12)
      score_mat[, n] <- as.vector(G_n)
    }
    ## NOTE: project = FALSE -- Theta is signed (no non-negative projection)
    C_boot <- .boot.onestep(as.vector(C), score_mat, Hinv, wild.B,
                            dist = "exp", seed = wild.seed, project = FALSE)
    sd_vec <- apply(C_boot, 1, stats::sd, na.rm = TRUE)
    C.se.boot <- matrix(sd_vec, nrow = Q, ncol = R, byrow = FALSE)
    alpha <- 1 - wild.level
    lo <- apply(C_boot, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
    hi <- apply(C_boot, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
    C.ci.lower <- matrix(lo, nrow = Q, ncol = R, byrow = FALSE)
    C.ci.upper <- matrix(hi, nrow = Q, ncol = R, byrow = FALSE)
  }

  Estimate <- as.vector(C)
  SE <- as.vector(C.se)
  BSE <- if (!is.null(C.se.boot)) as.vector(C.se.boot) else rep(NA_real_, length(Estimate))
  z_value <- ifelse(SE > 0, Estimate / SE, NA_real_)

  if (C.p.side == "one.sided") {
    p_value <- ifelse(is.finite(z_value), stats::pnorm(abs(z_value), lower.tail = FALSE), NA_real_)
  } else {
    p_value <- ifelse(is.finite(z_value), 2 * stats::pnorm(abs(z_value), lower.tail = FALSE), NA_real_)
  }

  rlabs <- if (!is.null(rownames(C))) rownames(C) else paste0("Resp", 1:Q)
  clabs <- if (!is.null(colnames(C))) colnames(C) else paste0("Cov", 1:R)

  coefficients <- data.frame(
    Basis     = rep(rlabs, times = R),
    Covariate = rep(clabs, each = Q),
    Estimate  = Estimate,
    SE        = SE,
    BSE       = BSE,
    z_value   = z_value,
    p_value   = p_value,
    CI_low    = if (!is.null(C.ci.lower)) as.vector(C.ci.lower) else NA_real_,
    CI_high   = if (!is.null(C.ci.upper)) as.vector(C.ci.upper) else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )

  if (print.trace) {
    message(sprintf("  Signed-Bottleneck NMF-AE inference: sandwich%s, p-value=%s",
                    if (isTRUE(wild.bootstrap)) " + wild bootstrap" else "",
                    C.p.side))
  }

  object$sigma2.used  <- sigma2.used
  object$C.se         <- C.se
  object$C.se.boot    <- C.se.boot
  object$C.ci.lower   <- C.ci.lower
  object$C.ci.upper   <- C.ci.upper
  object$coefficients <- coefficients
  object$C.p.side     <- C.p.side
  class(object) <- c("nmfae.signed.inference", "nmfae.inference",
                     "nmfae.signed", "nmfae", "nmf")
  object
}


## ==============================================================
#' @title Predict method for nmfae.signed
#' @keywords internal
#' @description
#' Computes \eqn{\hat Y_1 = X_1 (C_{+} - C_{-}) X_2 Y_2^{\mathrm{new}}}.
#' Since \eqn{\Theta = C_{+} - C_{-}} is signed, predictions may contain
#' negative entries even when \eqn{Y_1 \ge 0} in training.
#'
#' @param object A fitted \code{"nmfae.signed"} object.
#' @param newY2 New input matrix (P2 x N_new).  If \code{NULL}, returns
#'   the training fitted values.
#' @param Y1 Optional reference Y1 for scatter / confusion plot.
#' @param type Output: \code{"response"} (raw signed) or \code{"class"}.
#' @param ... Unused.
#' @return A numeric matrix (\code{"response"}) or factor (\code{"class"}).
#' @seealso \code{\link{nmfae.signed}}
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
predict.nmfae.signed <- function(object, newY2 = NULL, Y1 = NULL,
                               type = c("response", "class"), ...) {
  type <- match.arg(type)
  if (is.null(newY2)) {
    pred_mat <- object$Y1hat
  } else {
    pred_mat <- object$X1 %*% object$C %*% object$X2 %*% as.matrix(newY2)
  }

  if (type == "response") {
    result <- pred_mat
    if (!is.null(Y1)) attr(result, "Y1") <- Y1
    class(result) <- c("predict.nmfae", class(result))
  } else {
    ## For class prediction, clip negatives then take argmax (probabilities)
    clipped <- pmax(pred_mat, 0)
    labels <- rownames(pred_mat)
    if (is.null(labels)) labels <- paste0("C", seq_len(nrow(clipped)))
    pred_class <- factor(labels[apply(clipped, 2, which.max)], levels = labels)
    result <- pred_class
    if (!is.null(Y1)) {
      Y1m <- as.matrix(Y1)
      act_labels <- rownames(Y1m); if (is.null(act_labels)) act_labels <- labels
      actual <- factor(act_labels[apply(Y1m, 2, which.max)], levels = act_labels)
      attr(result, "actual") <- actual
    }
    class(result) <- c("predict.nmfae", class(result))
    attr(result, "type") <- "class"
  }
  result
}


## ==============================================================
#' @title Plot method for nmfae.signed (convergence)
#' @keywords internal
#' @description
#' Displays the convergence trajectory of the objective function.
#' @param x An \code{nmfae.signed} object.
#' @param ... Additional graphical parameters.
#' @return Invisible \code{NULL}.
#' @seealso \code{\link{nmfae.signed}}
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
plot.nmfae.signed <- function(x, ...) {
  extra_args <- list(...)
  args <- list(x = x$objfunc.iter)
  if (is.null(extra_args$main))
    args$main <- paste0("Signed-Bottleneck NMF-AE: R2 = ", round(x$r.squared, 3))
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  if (is.null(extra_args$type)) args$type <- "l"
  do.call(graphics::plot, c(args, extra_args))
  invisible(NULL)
}


## ==============================================================
#' @title Summary method for nmfae.signed
#' @keywords internal
#' @description
#' Produces a summary with dimensions, convergence, fit statistics,
#' and structure diagnostics (sparsity and negative-mass ratio).
#' @param object An \code{nmfae.signed} object.
#' @param ... Unused.
#' @return An object of class \code{"summary.nmfae.signed"}.
#' @seealso \code{\link{nmfae.signed}}
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
summary.nmfae.signed <- function(object, ...) {
  .sparsity <- function(M, thr = 1e-4) {
    if (is.null(M)) return(NA_real_); mean(abs(M) < thr)
  }
  .range <- function(M) if (is.null(M)) c(NA_real_, NA_real_) else range(M)
  .negfrac <- function(M) {
    if (is.null(M)) return(NA_real_)
    sum(pmax(-M, 0)) / (sum(abs(M)) + 1e-16)
  }
  n_params <- prod(dim(object$X1)) + prod(dim(object$Cp)) +
              prod(dim(object$Cn)) + prod(dim(object$X2))

  ans <- list(
    call        = object$call,
    dims        = object$dims,
    Q           = unname(object$rank["Q"]),
    R           = unname(object$rank["R"]),
    n.params    = n_params,
    Y.signed    = isTRUE(object$Y.signed),
    niter       = object$niter,
    runtime     = object$runtime,
    objfunc     = object$objfunc,
    r.squared          = object$r.squared,
    r.squared.uncentered     = object$r.squared.uncentered,
    r.squared.centered = object$r.squared.centered,
    sigma       = object$sigma,
    mae         = object$mae,
    X1.sparsity = .sparsity(object$X1),
    Cp.sparsity = .sparsity(object$Cp),
    Cn.sparsity = .sparsity(object$Cn),
    C.sparsity  = .sparsity(object$C),
    X2.sparsity = .sparsity(object$X2),
    X1.range    = .range(object$X1),
    Cp.range    = .range(object$Cp),
    Cn.range    = .range(object$Cn),
    C.range     = .range(object$C),
    X2.range    = .range(object$X2),
    C.negfrac   = .negfrac(object$C),
    H.negfrac   = .negfrac(object$H)
  )
  class(ans) <- "summary.nmfae.signed"
  ans
}

#' @title Print method for summary.nmfae.signed
#' @keywords internal
#' @param x A \code{"summary.nmfae.signed"} object.
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
print.summary.nmfae.signed <- function(x,
                                     digits = max(3L, getOption("digits") - 3L),
                                     ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Signed-Bottleneck NMF-AE (Direct MU; signed bottleneck Theta = Cp - Cn)\n")

  cat("\nDimensions:\n")
  cat(sprintf("  P1 (output rows): %d\n", x$dims["P1"]))
  cat(sprintf("  P2 (input rows):  %d\n", x$dims["P2"]))
  cat(sprintf("  N  (samples):     %d\n", x$dims["N"]))
  cat(sprintf("  Q (decoder rank): %d\n", x$Q))
  cat(sprintf("  R (encoder rank): %d\n", x$R))
  cat(sprintf("  Parameters:       %d\n", x$n.params))
  if (x$Y.signed) cat("  Y1 signed: TRUE\n")

  cat("\nConvergence:\n")
  cat(sprintf("  Iterations:       %d\n", x$niter))
  cat(sprintf("  Runtime (secs):   %.2f\n", x$runtime))
  cat(sprintf("  Final objfunc:    %s\n", format(x$objfunc, digits = digits)))

  cat("\nGoodness of fit:\n")
  cat(sprintf("  R-squared (cor^2):    %s\n", format(x$r.squared, digits = digits)))
  if (!is.null(x$r.squared.uncentered))
    cat(sprintf("  R-squared (uncentered):     %s\n", format(x$r.squared.uncentered, digits = digits)))
  if (!is.null(x$r.squared.centered))
    cat(sprintf("  R-squared (centered): %s\n", format(x$r.squared.centered, digits = digits)))
  cat(sprintf("  Sigma (RMSE):         %s\n", format(x$sigma,     digits = digits)))
  cat(sprintf("  MAE:                  %s\n", format(x$mae,       digits = digits)))

  cat("\nStructure (range / sparsity / negative mass):\n")
  fmt_range <- function(r)
    sprintf("[%s, %s]",
            format(r[1], digits = digits, width = 1),
            format(r[2], digits = digits, width = 1))
  cat(sprintf("  X1 (P1 x Q):     range %s   sparsity %5.1f%%\n",
              fmt_range(x$X1.range), 100 * x$X1.sparsity))
  cat(sprintf("  Cp (Q x R):      range %s   sparsity %5.1f%%\n",
              fmt_range(x$Cp.range), 100 * x$Cp.sparsity))
  cat(sprintf("  Cn (Q x R):      range %s   sparsity %5.1f%%\n",
              fmt_range(x$Cn.range), 100 * x$Cn.sparsity))
  cat(sprintf("  C=Cp-Cn (Q x R): range %s   sparsity %5.1f%%   neg-mass %5.1f%%\n",
              fmt_range(x$C.range), 100 * x$C.sparsity, 100 * x$C.negfrac))
  cat(sprintf("  X2 (R x P2):     range %s   sparsity %5.1f%%\n",
              fmt_range(x$X2.range), 100 * x$X2.sparsity))
  cat(sprintf("  Encoding H (Q x N):                                   neg-mass %5.1f%%\n",
              100 * x$H.negfrac))
  cat("\n")
  invisible(x)
}


## ==============================================================
#' @title Rename Resp/Cov labels on nmfae.signed objects
#' @keywords internal
#' @description
#' Replaces the default \code{"Resp1", "Resp2", ...} (response basis / X1 columns
#' and Cp/Cn/C rows) and \code{"Cov1", "Cov2", ...} (covariate basis / X2 rows and
#' Cp/Cn/C columns) with user-supplied labels.  Propagates to coefficients
#' tables if present (e.g., from \code{nmfae.inference}).
#'
#' @param x An \code{"nmfae.signed"} object.
#' @param X1.colnames Character vector of length Q for decoder labels.
#' @param X2.rownames Character vector of length R for encoder labels.
#' @return The renamed object (same class).
#' @seealso \code{\link{nmfae.signed}}
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
nmfae.signed.rename <- function(x, X1.colnames = NULL, X2.rownames = NULL) {
  if (!is.null(X1.colnames)) {
    Q <- ncol(x$X1)
    if (length(X1.colnames) != Q) stop("X1.colnames must have length ", Q)
    colnames(x$X1) <- X1.colnames
    rownames(x$Cp) <- X1.colnames
    rownames(x$Cn) <- X1.colnames
    rownames(x$C)  <- X1.colnames
    if (!is.null(x$coefficients)) {
      old <- paste0("Resp", seq_len(Q))
      for (k in seq_len(Q))
        x$coefficients$Basis[x$coefficients$Basis == old[k]] <- X1.colnames[k]
    }
  }
  if (!is.null(X2.rownames)) {
    R <- nrow(x$X2)
    if (length(X2.rownames) != R) stop("X2.rownames must have length ", R)
    rownames(x$X2) <- X2.rownames
    colnames(x$Cp) <- X2.rownames
    colnames(x$Cn) <- X2.rownames
    colnames(x$C)  <- X2.rownames
    if (!is.null(x$coefficients)) {
      old <- paste0("Cov", seq_len(R))
      for (k in seq_len(R))
        x$coefficients$Covariate[x$coefficients$Covariate == old[k]] <- X2.rownames[k]
    }
  }
  x
}


## ==============================================================
#' @title Element-wise Cross-Validation for Signed-Bottleneck NMF-AE
#' @keywords internal
#' @description
#' Element-wise k-fold cross-validation for \code{\link{nmfae.signed}} to
#' select the decoder / encoder ranks \eqn{(Q, R)}.  Mirrors
#' \code{\link{nmfae.ecv}} but uses the \strong{weighted} Signed-Bottleneck NMF-AE fit
#' path (\code{Y1.weights}): test-fold elements are zero-weighted during
#' fitting, and held-out MSE is computed on those elements.
#'
#' @param Y1 Output matrix (P1 x N).
#' @param Y2 Input matrix (P2 x N).  Default \code{Y1}.
#' @param rank1 Integer vector of candidate response-basis ranks. Default \code{1:2}.
#' @param rank2 Integer vector of candidate covariate-basis ranks, or \code{NULL}
#'   (default: pair rank2 = rank1, diagonal grid).
#' @param rank,rank.encoder Deprecated aliases of \code{rank1} / \code{rank2}.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{nfolds} / \code{div}}{Number of folds.  Default 5.}
#'     \item{\code{seed}}{RNG seed for fold assignment.  Default 123.}
#'     \item{\code{nstart}}{Number of random restarts per fit.  Default 1.
#'       Signed models have more local minima (the bottleneck can carry
#'       both signs), so \code{nstart >= 10} is recommended for
#'       reproducible rank selection.}
#'     \item{Other args}{\code{epsilon}, \code{maxit}, \code{warm.start},
#'       etc.\ are passed to \code{\link{nmfae.signed}}.}
#'   }
#' @return An object of class \code{c("nmfae.signed.ecv", "nmfae.ecv")} with
#'   \code{objfunc} (MSE per pair), \code{sigma} (RMSE), \code{objfunc.fold}
#'   (per-fold MSE), \code{folds}, \code{QR}, \code{paired}.
#' @seealso \code{\link{nmfae.signed}}, \code{\link{nmfae.ecv}}
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
nmfae.signed.ecv <- function(Y1, Y2 = Y1, rank1 = 1:2, rank2 = NULL, ...,
                             rank = NULL, rank.encoder = NULL) {
  extra_ecv <- list(...)
  # rank1/rank2 = response/covariate basis ranks to sweep; legacy rank/rank.encoder/Q/R
  if (is.null(rank))          rank <- rank1
  if (!is.null(extra_ecv$Q))  rank <- extra_ecv$Q
  if (is.null(rank.encoder))  rank.encoder <- rank2
  if (!is.null(extra_ecv$R))  rank.encoder <- extra_ecv$R
  nfolds <- if (!is.null(extra_ecv$nfolds)) extra_ecv$nfolds
            else if (!is.null(extra_ecv$div)) extra_ecv$div else 5
  seed   <- if (!is.null(extra_ecv$seed)) extra_ecv$seed else 123
  Q <- rank; R <- rank.encoder
  div <- nfolds

  Y1 <- as.matrix(Y1); Y2 <- as.matrix(Y2)
  P1 <- nrow(Y1); N <- ncol(Y1)

  if (is.null(R)) {
    QR <- data.frame(Q = Q, R = Q)
  } else {
    QR <- expand.grid(Q = Q, R = R)
  }
  num_pairs <- nrow(QR)

  ## Element-wise folds (shared helper)
  folds <- .ecv.make.folds(Y1, div, seed)

  pair_labels <- sprintf("Q=%d,R=%d", QR$Q, QR$R)
  has_na <- any(is.na(Y1))
  message(sprintf("Signed-Bottleneck NMF-AE ECV: %d (Q,R) pairs, %d-fold, %d tasks...",
                  num_pairs, div, num_pairs * div))

  ## Strip ecv-specific args before passing the rest to nmfae.signed
  fit_args <- extra_ecv
  fit_args$nfolds <- NULL; fit_args$div <- NULL
  fit_args$Q <- NULL; fit_args$R <- NULL

  run_one <- function(i, k) {
    test_idx <- folds[[k]]
    W_train <- matrix(1, nrow = P1, ncol = N)
    if (has_na) W_train[is.na(Y1)] <- 0
    W_train[test_idx] <- 0
    fit <- suppressMessages(
      do.call(nmfae.signed,
              c(list(Y1 = Y1, Y2 = Y2, rank = QR$Q[i], rank.encoder = QR$R[i],
                     Y1.weights = W_train), fit_args))
    )
    mean((Y1[test_idx] - fit$Y1hat[test_idx])^2)
  }

  cv <- .ecv.run(pair_labels, div, run_one,
                 progress = function(i, o, s)
                   message(sprintf("  Q=%d, R=%d: MSE=%.6f, sigma=%.4f",
                                   QR$Q[i], QR$R[i], o, s)))

  result <- list(objfunc = cv$objfunc,
                 sigma = cv$sigma,
                 objfunc.fold = cv$objfunc.fold,
                 folds = folds,
                 QR = QR,
                 paired = is.null(R))
  class(result) <- c("nmfae.signed.ecv", "nmfae.ecv")
  result
}


#' @title Rank selection for nmfae.signed (paired rank, concise diagnostics)
#' @keywords internal
#' @description
#' Fits \code{\link{nmfae.signed}} with a \strong{paired} decoder/encoder
#' rank (\eqn{Q = R}) across a range of ranks and reports
#' \code{r.squared}, the effective rank (of the latent encoding \eqn{H}),
#' and the element-wise CV error \code{sigma.ecv}, with the same concise
#' plot as \code{\link{nmfkc.rank}}.  For a full \eqn{(Q, R)} grid use
#' \code{\link{nmfae.signed.ecv}}.
#' @param Y1 Endogenous matrix (\eqn{P_1 \times N}); may be signed.
#' @param Y2 Exogenous matrix; defaults to \code{Y1}.
#' @param rank1 Integer vector of (paired) ranks to evaluate (both bases use
#'   the same value). Legacy \code{Q} accepted via \code{...}.
#' @param rank Deprecated alias of \code{rank1}.
#' @param plot Logical; draw the diagnostics plot (default \code{TRUE}).
#' @param detail \code{"full"} (default) also runs element-wise CV
#'   (\code{sigma.ecv}); \code{"fast"} skips it (plots r.squared and
#'   eff.rank only, and recommends the R-squared elbow).
#' @param ... Passed on to \code{\link{nmfae.signed}} and
#'   \code{\link{nmfae.signed.ecv}}.
#' @return A list with \code{rank.best} and \code{criteria}
#'   (\code{rank}, \code{effective.rank}, \code{effective.rank.ratio},
#'   \code{r.squared}, \code{sigma.ecv}).
#' @seealso \code{\link{nmfae.signed}}, \code{\link{nmfae.signed.ecv}},
#'   \code{\link{nmfkc.rank}}
#' @references
#' Roy, O., & Vetterli, M. (2007).  The effective rank: A measure of
#' effective dimensionality.  \emph{Proc. EUSIPCO}, 606--610.
#' (\code{effective.rank})
#' Wold, S. (1978).  Cross-validatory estimation of the number of
#' components in factor and principal components models.
#' \emph{Technometrics}, 20(4), 397--405. (\code{sigma.ecv})
#' @export
nmfae.signed.rank <- function(Y1, Y2 = Y1, rank1 = 1:5, detail = c("full", "fast"),
                              plot = TRUE, ..., rank = NULL) {
  extra <- list(...)
  if (!is.null(rank))    rank1 <- rank
  if (!is.null(extra$Q)) rank1 <- extra$Q
  extra$Q <- NULL; extra$R <- NULL; extra$rank.encoder <- NULL
  detail <- match.arg(detail)
  Y1 <- as.matrix(Y1); Y2 <- as.matrix(Y2)
  rs <- numeric(length(rank1)); er <- numeric(length(rank1))
  for (i in seq_along(rank1)) {
    f <- suppressMessages(do.call(nmfae.signed,
           c(list(Y1, Y2, rank1 = rank1[i], rank2 = rank1[i],
                  print.trace = FALSE), extra)))
    rs[i] <- f$r.squared
    er[i] <- .effective.rank(f$H)
  }
  ecv <- if (detail == "full")
    suppressMessages(do.call(nmfae.signed.ecv, c(list(Y1, Y2, rank1 = rank1), extra)))$sigma
    else rep(NA_real_, length(rank1))
  criteria <- data.frame(rank = rank1, effective.rank = er,
                         effective.rank.ratio = er / rank1,
                         r.squared = rs, sigma.ecv = as.numeric(ecv))
  .rank.finish(criteria, plot = plot,
               main = "nmfae.signed rank selection (paired Q=R)")
}


## ==============================================================
#' @title Summary method for nmfae.signed.inference objects
#' @keywords internal
#' @description
#' Produces a summary of a fitted Signed-Bottleneck NMF-AE model with
#' inference results.  Extends \code{\link{summary.nmfae.signed}} by
#' attaching the \code{coefficients} table and p-value side from
#' \code{\link{nmfae.signed.inference}}.
#'
#' @param object An object of class \code{"nmfae.signed.inference"}.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{"summary.nmfae.signed.inference"}.
#' @seealso \code{\link{nmfae.signed.inference}}, \code{\link{summary.nmfae.signed}}
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
summary.nmfae.signed.inference <- function(object, ...) {
  ans <- summary.nmfae.signed(object, ...)
  ans$coefficients <- object$coefficients
  ans$C.p.side     <- object$C.p.side
  ans$sigma2.used  <- object$sigma2.used
  class(ans) <- c("summary.nmfae.signed.inference",
                   "summary.nmfae.signed", "summary.nmfae.inference")
  ans
}

#' @title Print method for summary.nmfae.signed.inference objects
#' @keywords internal
#' @description
#' Prints the Signed-Bottleneck NMF-AE summary followed by the
#' coefficients table of Theta.
#' @param x An object of class \code{"summary.nmfae.signed.inference"}.
#' @param digits Minimum number of significant digits.
#' @param by Character; grouping order of the coefficients table.
#'   \code{"covariate"} (default) lists all bases within each covariate
#'   (1-1, 1-2, ...); \code{"basis"} lists all covariates within each basis
#'   (1-1, 2-1, ...).
#' @param ... Additional arguments (currently unused).
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @seealso \code{\link{summary.nmfae.signed.inference}}
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
print.summary.nmfae.signed.inference <- function(x,
    digits = max(3L, getOption("digits") - 3L),
    by = c("covariate", "basis"), ...) {
  by <- match.arg(by)
  print.summary.nmfae.signed(x, digits = digits, ...)

  if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
    cf <- x$coefficients
    cf <- cf[.coef.order.by(cf, by), , drop = FALSE]   # grouping order (by)
    p_side <- if (!is.null(x$C.p.side)) x$C.p.side else "two.sided"
    p_header <- if (p_side == "one.sided") "Pr(>z)" else "Pr(>|z|)"

    sig_stars <- function(p) {
      ifelse(!is.finite(p), " ",
        ifelse(p < 0.001, "***",
          ifelse(p < 0.01, "**",
            ifelse(p < 0.05, "*",
              ifelse(p < 0.1, ".", " ")))))
    }
    format_pval <- function(p) {
      ifelse(!is.finite(p), "      NA",
        ifelse(p < 2.2e-16, "  <2e-16",
          formatC(p, format = "g", digits = 4, width = 8)))
    }

    n_total <- nrow(cf)
    n_sig <- sum(cf$p_value < 0.05, na.rm = TRUE)
    cat(sprintf("\nTheta coefficients: %d total, %d significant (p < 0.05)\n",
                n_total, n_sig))
    rnames <- paste0(cf$Covariate, ":", cf$Basis)   # Covariate:Basis (matches nmfre/nmfae)
    est <- formatC(cf$Estimate, format = "f", digits = 3, width = 9)
    se  <- formatC(cf$SE,       format = "f", digits = 3, width = 10)
    bse <- formatC(cf$BSE,      format = "f", digits = 3, width = 6)
    zv  <- formatC(cf$z_value,  format = "f", digits = 2, width = 7)
    pv_str <- format_pval(cf$p_value)
    stars <- sig_stars(cf$p_value)
    max_lw <- max(nchar(rnames))
    hdr <- sprintf("%s %s %s %s %s",
                   formatC("Estimate",   width = 9),
                   formatC("Std. Error", width = 10),
                   formatC("(Boot)",     width = 6),
                   formatC("z value",    width = 7),
                   formatC(p_header,     width = 8))
    cat(sprintf("%s %s\n", formatC("Cov:Resp", width = max_lw), hdr))
    for (i in seq_len(n_total)) {
      cat(sprintf("%s %s %s %s %s %s %s\n",
                  formatC(rnames[i], width = max_lw),
                  est[i], se[i], bse[i], zv[i], pv_str[i], stars[i]))
    }
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  cat("\n")
  invisible(x)
}


## ==============================================================
#' @title Heatmap visualization of nmfae.signed factor matrices
#' @keywords internal
#' @description
#' Displays the factor blocks of a \code{\link{nmfae.signed}} fit as
#' side-by-side heatmaps.  Non-negative blocks (\eqn{X_1, C_{+}, C_{-},
#' X_2}) use the white-orange-red palette; the signed combined
#' bottleneck \eqn{C = C_{+} - C_{-}} is rendered with a diverging
#' blue-white-red palette so positive and negative weights are visually
#' distinguishable.
#'
#' @param x An object of class \code{"nmfae.signed"}.
#' @param Y1.label Character vector for rows of \eqn{X_1}.
#' @param X1.label Decoder basis labels.
#' @param X2.label Encoder basis labels.
#' @param Y2.label Input variable labels.
#' @param palette.pos Palette for non-negative blocks. Default white-orange-red.
#' @param palette.signed Palette for signed \eqn{C}. Default blue-white-red.
#' @param show.C Logical. If \code{TRUE} (default), shows the combined
#'   signed \eqn{C = C_{+} - C_{-}} as a separate panel.
#' @param ... Not used.
#'
#' @return Invisible \code{NULL}. Called for its side effect (plot).
#' @seealso \code{\link{nmfae.signed}}, \code{\link{nmfae.heatmap}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y1 <- matrix(abs(rnorm(12)), 3, 4)
#' Y2 <- matrix(abs(rnorm(20)), 5, 4)
#' res <- nmfae.signed(Y1, Y2, rank1 = 2, rank2 = 2, maxit = 200)
#' nmfae.signed.heatmap(res)
#' }
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
nmfae.signed.heatmap <- function(x,
                                  Y1.label = NULL, X1.label = NULL,
                                  X2.label = NULL, Y2.label = NULL,
                                  palette.pos = NULL,
                                  palette.signed = NULL,
                                  show.C = TRUE, ...) {
  if (!inherits(x, "nmfae.signed"))
    stop("x must be of class 'nmfae.signed'.")
  if (is.null(palette.pos))
    palette.pos <- grDevices::colorRampPalette(c("white", "orange", "red"))(64)
  if (is.null(palette.signed))
    palette.signed <- grDevices::colorRampPalette(
      c("blue", "white", "red"))(64)

  X1 <- x$X1; Cp <- x$Cp; Cn <- x$Cn; X2 <- x$X2
  C  <- if (!is.null(x$C)) x$C else (Cp - Cn)
  Q  <- ncol(X1); R  <- nrow(X2)
  P1 <- nrow(X1); P2 <- ncol(X2)

  if (is.null(Y1.label)) Y1.label <- rownames(X1)
  if (is.null(X1.label)) X1.label <- colnames(X1)
  if (is.null(X2.label)) X2.label <- rownames(X2)
  if (is.null(Y2.label)) Y2.label <- colnames(X2)
  if (is.null(Y1.label)) Y1.label <- as.character(seq_len(P1))
  if (is.null(X1.label)) X1.label <- paste0("Resp", seq_len(Q))
  if (is.null(X2.label)) X2.label <- paste0("Cov", seq_len(R))
  if (is.null(Y2.label)) Y2.label <- as.character(seq_len(P2))

  .heat <- function(mat, title, pal, zlim = NULL) {
    m <- mat[nrow(mat):1, , drop = FALSE]
    if (is.null(zlim)) zlim <- range(m, finite = TRUE)
    graphics::image(x = seq_len(ncol(m)), y = seq_len(nrow(m)),
                    z = t(m), col = pal, zlim = zlim, axes = FALSE,
                    xlab = "", ylab = "", main = title)
    graphics::axis(1, at = seq_len(ncol(m)), labels = colnames(m),
                   las = 2, cex.axis = 0.7)
    graphics::axis(2, at = seq_len(nrow(m)), labels = rev(rownames(m)),
                   las = 1, cex.axis = 0.7)
    graphics::box()
  }

  n.panels <- if (isTRUE(show.C)) 5 else 4
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par), add = TRUE)
  graphics::par(mfrow = c(1, n.panels), mar = c(4, 3, 2, 1))

  rownames(X1) <- Y1.label; colnames(X1) <- X1.label
  rownames(Cp) <- X1.label; colnames(Cp) <- X2.label
  rownames(Cn) <- X1.label; colnames(Cn) <- X2.label
  rownames(X2) <- X2.label; colnames(X2) <- Y2.label
  rownames(C)  <- X1.label; colnames(C)  <- X2.label

  .heat(X1, expression(X[1]), palette.pos)
  .heat(Cp, expression(C["+"]), palette.pos)
  .heat(Cn, expression(C["-"]), palette.pos)
  if (isTRUE(show.C)) {
    cmax <- max(abs(C), na.rm = TRUE)
    .heat(C, expression(C == C["+"] - C["-"]),
          palette.signed, zlim = c(-cmax, cmax))
  }
  .heat(X2, expression(X[2]), palette.pos)

  invisible(NULL)
}

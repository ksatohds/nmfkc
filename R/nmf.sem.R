# nmf.sem.R — NMF-FFB (formerly NMF-SEM) canonical engines + generic DOT
# Canonical:  nmf.ffb, nmf.ffb.inference, nmf.ffb.cv, nmf.ffb.split, nmf.ffb.DOT
#             (deprecated nmf.sem* aliases live in nmf.sem-deprecated.R).
# Also hosts: nmfkc.DOT / plot.nmfkc.DOT (shared DOT utilities).

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#' @title NMF-FFB Main Estimation Algorithm (formerly NMF-SEM)
#'
#' @description
#' Fits the NMF-FFB model
#' \deqn{
#'   Y_1 \approx X \bigl( \Theta_1 Y_1 + \Theta_2 Y_2 \bigr)
#' }
#' under non-negativity constraints with orthogonality and sparsity regularization.
#' The function returns the estimated latent factors, structural coefficient matrices,
#' and the implied equilibrium (input–output) mapping.
#'
#' At equilibrium, the model can be written as
#' \deqn{
#'   Y_1 \approx (I - X \Theta_1)^{-1} X \Theta_2 Y_2
#'   \equiv M_{\mathrm{model}} Y_2,
#' }
#' where \eqn{M_{\mathrm{model}} = (I - X \Theta_1)^{-1} X \Theta_2} is a
#' Leontief-type cumulative-effect operator in latent space.
#'
#' Internally, the latent feedback and exogenous loading matrices are stored as
#' \code{C1} and \code{C2}, corresponding to \eqn{\Theta_1} and \eqn{\Theta_2},
#' respectively.
#'
#' @param Y1 A non-negative numeric matrix of endogenous variables with
#'   \strong{rows = variables (P1), columns = samples (N)}.
#' @param Y2 A non-negative numeric matrix of exogenous variables with
#'   \strong{rows = variables (P2), columns = samples (N)}.
#'   Must satisfy \code{ncol(Y1) == ncol(Y2)}.
#' @param rank Integer; number of latent factors \eqn{Q}. If \code{NULL},
#'   \eqn{Q} is taken from a hidden argument in \code{...} or defaults to
#'   \code{nrow(Y2)}.
#' @param X.init Initialization strategy for the basis matrix
#'   \code{X} (\eqn{P_1 \times Q}).  One of:
#'   \itemize{
#'     \item \code{"nndsvd"} (default): Non-negative Double SVD with
#'       additive randomness (NNDSVDar; Boutsidis & Gallopoulos 2008),
#'       computed internally via \code{.nndsvdar(Y1, Q)}.  Requires
#'       \eqn{Q \le \min(P_1, N)} (over-rank case falls back to
#'       \code{"runif"}).  Uses a full SVD of \eqn{Y_1}, so for very
#'       large \eqn{Y_1} consider switching to \code{"kmeans"} to
#'       avoid SVD memory / compute cost.
#'     \item \code{"kmeans"}: k-means on the columns of \eqn{Y_1}
#'       (samples clustered into \eqn{Q} groups); the transposed
#'       cluster centers become \eqn{X}.  Scales well for large
#'       \eqn{Y_1}; this is the default of \code{\link{nmfkc}}.
#'     \item \code{"kmeansar"}: \code{"kmeans"} followed by filling
#'       zero entries of \eqn{X} with \eqn{\mathrm{Uniform}(0,
#'       \bar Y_1 / 100)} (NNDSVDar-style additive randomness to
#'       escape trivial stationary points).
#'     \item \code{"runif"}: Uniform random entries in \eqn{[0, 1]}.
#'     \item A numeric \eqn{P_1 \times Q} matrix supplied by the user;
#'       negative entries are projected to 0.
#'     \item \code{NULL}: backward-compatible alias for \code{"nndsvd"}.
#'   }
#'   In all cases the result is column-normalized to \code{colSums(X) = 1}
#'   before iteration.  The menu mirrors \code{\link{nmfkc}}'s
#'   \code{X.init} option for consistency across the package.
#' @param X.L2.ortho L2 orthogonality penalty for \code{X}. This controls
#'   the penalty term \eqn{\lambda_X \lVert X^\top X - \mathrm{diag}(X^\top X)
#'   \rVert_F^2}. Default: \code{100}.
#' @param C1.L1 L1 sparsity penalty for \code{C1} (i.e., \eqn{\Theta_1}).
#'   Default: \code{1.0}.
#' @param C2.L1 L1 sparsity penalty for \code{C2} (i.e., \eqn{\Theta_2}).
#'   Default: \code{0.1}.
#' @param epsilon Relative convergence threshold for the objective function.
#'   Iterations stop when the relative change in reconstruction loss falls
#'   below this value. Default: \code{1e-6}.
#' @param maxit Maximum number of iterations for the multiplicative updates.
#'   Default: \code{5000} (matches \code{\link{nmfkc}} and other MU
#'   functions in the package).
#' @param seed Random seed used to initialize \code{X}, \code{C1}, and \code{C2}.
#'   Default: \code{123}.
#' @param ... Additional hidden arguments controlling the optional
#'   feedforward baseline (used both as an \eqn{X} warm-start and as
#'   the reference for \code{SC.map}, the input-output structural
#'   fidelity defined in Satoh (2025) §4.SC.map):
#'   \describe{
#'     \item{\code{nmfkc.baseline}}{Controls whether a feedforward
#'       \code{\link{nmfkc}}(Y1, A = Y2) fit is used as baseline.
#'       Possible values:
#'       \itemize{
#'         \item Default (not given) — \code{nmf.sem} runs
#'           \code{\link{nmfkc}} \strong{internally} when \code{X.init}
#'           is a string method (\code{"nndsvd"}, \code{"kmeans"},
#'           \dots) or \code{NULL}, forwarding \code{X.init},
#'           \code{X.L2.ortho}, \code{epsilon}, \code{maxit},
#'           \code{seed}.  The fitted \eqn{X} of the baseline is then
#'           used as warm-start for the nmf.sem MU iterations, and
#'           \code{SC.map} is computed.  This means
#'           \code{nmf.sem(Y1, Y2, rank = Q)} runs end-to-end without
#'           a prior \code{nmfkc} call.
#'         \item \code{TRUE} — same as above, but force the internal
#'           \code{\link{nmfkc}} call even when \code{X.init} is a
#'           user-supplied matrix (the matrix is overridden).
#'         \item \code{FALSE} — opt out; no internal call,
#'           \code{SC.map = NA} (pre-v0.6.8 behavior).
#'         \item An \code{\link{nmfkc}} result (list with \code{$X}
#'           and \code{$C}) — use as the baseline directly (no
#'           internal call); also adopted as \code{X.init} when the
#'           latter is a string / NULL.
#'       }}
#'     \item{\code{M.simple}}{Optional \eqn{P_1 \times P_2} pre-computed
#'       baseline mapping.  Takes precedence over \code{nmfkc.baseline}
#'       for the SC.map calculation but does not affect warm-start.}
#'     \item{\code{Q}}{Backward-compat alias for \code{rank}.}
#'   }
#'
#' @return A list with components:
#'   \item{X}{Estimated basis matrix (\eqn{P_1 \times Q}).}
#'   \item{C1}{Estimated latent feedback matrix (\eqn{\Theta_1}, \eqn{Q \times P_1}).}
#'   \item{C2}{Estimated exogenous loading matrix (\eqn{\Theta_2}, \eqn{Q \times P_2}).}
#'   \item{XC1}{Feedback matrix \eqn{X \Theta_1}.}
#'   \item{XC2}{Direct-effect matrix \eqn{X \Theta_2}.}
#'   \item{XC1.radius}{Spectral radius \eqn{\rho(X \Theta_1)}.}
#'   \item{XC1.norm1}{Induced 1-norm \eqn{\lVert X \Theta_1 \rVert_{1,\mathrm{op}}}.}
#'   \item{Leontief.inv}{Leontief-type inverse \eqn{(I - X \Theta_1)^{-1}.}}
#'   \item{M.model}{Equilibrium mapping
#'     \eqn{M_{\mathrm{model}} = (I - X \Theta_1)^{-1} X \Theta_2}.}
#'   \item{amplification}{Latent amplification factor
#'     \eqn{\lVert M_{\mathrm{model}} \rVert_{1,\mathrm{op}} /
#'          \bigl\lVert X \Theta_2 \bigr\rVert_{1,\mathrm{op}}}.}
#'   \item{amplification.bound}{Geometric-series upper bound
#'     \eqn{1 / (1 - \lVert X \Theta_1 \rVert_{1,\mathrm{op}})} if
#'     \eqn{\lVert X \Theta_1 \rVert_{1,\mathrm{op}} < 1}, otherwise \code{Inf}.}
#'   \item{Q}{Effective latent dimension used in the fit.}
#'   \item{SC.cov}{Correlation between sample and model-implied covariance
#'     (flattened) of \eqn{Y_1}.  See \emph{second-moment fidelity} in
#'     Satoh (2025).}
#'   \item{SC.map}{Correlation between the equilibrium operator
#'     \eqn{M_{\mathrm{model}}} and a feedforward baseline mapping
#'     \eqn{M_{\mathrm{simple}} = X_0 \Theta_0}, computed only when the
#'     baseline is supplied via \code{M.simple} or \code{nmfkc.baseline}
#'     in \code{...}; otherwise \code{NA}.  See \emph{input-output
#'     structural fidelity} in Satoh (2025).}
#'   \item{MAE}{Mean absolute error between \eqn{Y_1} and its equilibrium
#'     prediction \eqn{\hat Y_1 = M_{\mathrm{model}} Y_2}.}
#'   \item{objfunc}{Vector of reconstruction losses per iteration.}
#'   \item{objfunc.full}{Vector of penalized objective values per iteration.}
#'   \item{iter}{Number of iterations actually performed.}
#'
#' @examples
#' # Simple NMF-FFB with iris data (non-negative)
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]  # Sepal
#' Y2 <- Y[3:4, ]  # Petal
#' result <- nmf.sem(Y1, Y2, rank = 2, maxit = 500)
#' result$MAE
#'
#' @seealso \code{\link{nmf.ffb.inference}}, \code{\link{nmf.ffb.cv}},
#'   \code{\link{nmf.ffb.split}}, \code{\link{nmf.ffb.DOT}},
#'   \code{\link{summary.nmf.sem}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to structural equation modeling for blind input-output analysis.
#'   arXiv:2512.18250. \url{https://arxiv.org/abs/2512.18250}
#' @export
nmf.ffb <- function(
    Y1, Y2,
    rank = NULL,
    X.init = "nndsvd",
    X.L2.ortho = 100.0,
    C1.L1 = 1.0,
    C2.L1 = 0.1,
    epsilon = 1e-6,
    maxit = 5000,
    seed  = 123,
    ...
) {
  # ------------------------------ checks ------------------------------
  if (!is.matrix(Y1)) Y1 <- as.matrix(Y1)
  if (!is.matrix(Y2)) Y2 <- as.matrix(Y2)

  if (any(!is.finite(Y1)) || any(!is.finite(Y2)))
    stop("Y1 and Y2 must not contain NA/NaN/Inf.")
  if (min(Y1) < 0 || min(Y2) < 0)
    stop("Y1 and Y2 must be non-negative.")
  if (ncol(Y1) != ncol(Y2))
    stop("ncol(Y1) must be equal to ncol(Y2).")

  extra_args <- list(...)
  Q_hidden <- if (!is.null(extra_args$Q)) extra_args$Q else NULL
  Q0 <- if (!is.null(rank)) rank else if (!is.null(Q_hidden)) Q_hidden else nrow(Y2)

  P1 <- nrow(Y1); P2 <- nrow(Y2); N <- ncol(Y1)
  Q  <- Q0
  if (Q < 1) stop("Rank Q must be >= 1.")

  # -------------------------- labels (output) -------------------------
  Y1_labels    <- if (!is.null(rownames(Y1))) rownames(Y1) else paste0("Y1_", 1:P1)
  Y2_labels    <- if (!is.null(rownames(Y2))) rownames(Y2) else paste0("Y2_", 1:P2)
  Basis_labels <- paste0("Factor", 1:Q)

  set.seed(seed)
  .eps <- 1e-10
  .xnorm  <- function(X) sweep(X, 2, pmax(colSums(X), .eps), "/")
  mat1norm <- function(A) max(colSums(abs(A)))

  # ---------- (optional) internal nmfkc warm-start + SC.map baseline --
  ## Resolution rules for `nmfkc.baseline`:
  ##   * not given (default): if X.init is a string method (or NULL),
  ##     run nmfkc(Y1, A = Y2) INTERNALLY using the same X.init,
  ##     X.L2.ortho, epsilon, maxit, seed; use its X as warm-start
  ##     and X * C as M.simple for SC.map.  This is the typical
  ##     workflow described in Satoh (2025) and means
  ##       res <- nmf.sem(Y1, Y2, rank = Q)
  ##     can run end-to-end without first calling nmfkc().
  ##   * TRUE: same as above (explicit opt-in even when X.init is a
  ##     user-supplied matrix; the matrix is overridden by nmfkc's X).
  ##   * FALSE: opt-out — no internal nmfkc, no warm-start, SC.map = NA
  ##     (equivalent to pre-v0.6.8 behavior).
  ##   * an nmfkc result (list with $X and $C): use as the baseline
  ##     for SC.map and as warm-start for X.init (no internal call).
  baseline_for_scmap <- NULL

  user_baseline    <- extra_args$nmfkc.baseline
  user_M.simple    <- extra_args$M.simple

  baseline_is_obj  <- is.list(user_baseline) &&
                      !is.null(user_baseline$X) && !is.null(user_baseline$C)
  baseline_is_TRUE <- isTRUE(user_baseline)
  baseline_is_FALSE <- isFALSE(user_baseline)
  baseline_is_default <- is.null(user_baseline)

  ## Default auto-run when baseline is unspecified AND X.init is a string
  ## method (NULL is treated as "nndsvd" string).
  auto_nmfkc <- if (baseline_is_FALSE) {
    FALSE
  } else if (baseline_is_TRUE) {
    TRUE
  } else if (baseline_is_obj) {
    FALSE
  } else if (baseline_is_default) {
    is.null(X.init) || is.character(X.init)
  } else {
    FALSE
  }

  if (auto_nmfkc) {
    ## Internal nmfkc call.  Forward only the genuinely shared options
    ## (X.init, X.L2.ortho, epsilon, maxit, seed); the nmf.sem-specific
    ## C1.L1 / C2.L1 do not apply to the feedforward baseline model.
    nmfkc_xinit <- if (is.null(X.init)) "nndsvd" else X.init
    baseline_for_scmap <- nmfkc(
      Y = Y1, A = Y2, Q = Q,
      X.init = nmfkc_xinit,
      X.L2.ortho = X.L2.ortho,
      epsilon = epsilon,
      maxit = maxit,
      seed = seed,
      verbose = FALSE,
      print.dims = FALSE
    )
    ## Override X.init with the nmfkc-fitted X for nmf.sem warm-start
    X.init <- baseline_for_scmap$X
  } else if (baseline_is_obj) {
    baseline_for_scmap <- user_baseline
    ## When the user supplies an nmfkc result AND X.init is still a
    ## string / NULL, also use baseline$X as warm-start (rmd workflow).
    if (is.null(X.init) || is.character(X.init)) {
      X.init <- baseline_for_scmap$X
    }
  }

  # ---------------------------- init X,C1,C2 --------------------------
  ## X.init dispatch: delegate string methods to the shared internal
  ## helper .init_X_method() (defined in R/nmfkc.R).  Accepts:
  ##   "nndsvd" (default), "kmeans", "kmeansar", "runif",
  ##   a numeric P1 x Q matrix, or NULL (alias for "nndsvd").
  if (is.null(X.init)) X.init <- "nndsvd"
  if (is.character(X.init)) {
    X <- .init_X_method(X.init, Y1, Q)
  } else {
    X <- as.matrix(X.init)
    if (!all(dim(X) == c(P1, Q))) {
      stop("X.init must have dimension (nrow(Y1) x rank).")
    }
    X[X < 0] <- 0
  }
  X <- .xnorm(X)

  C1 <- matrix(stats::runif(Q * P1, 0.01, 0.1), nrow = Q, ncol = P1)
  C2 <- matrix(stats::runif(Q * P2, 0.001, 0.01), nrow = Q, ncol = P2)

  min_dim <- min(Q, P2)
  for (i in 1:min_dim) {
    C2[i, i] <- stats::runif(1, 0.1, 0.2)
  }

  objfunc      <- numeric(maxit)
  objfunc.full <- numeric(maxit)

  # ----------------------------- main loop ----------------------------
  for (it in 1:maxit) {
    M  <- C1 %*% Y1 + C2 %*% Y2
    Mt <- t(M)

    # 2.1 update X
    Numerator_X       <- Y1 %*% Mt
    Denominator_X_rec <- X %*% M %*% Mt

    if (X.L2.ortho > 0) {
      XtX <- t(X) %*% X
      XtX_offdiag <- XtX
      diag(XtX_offdiag) <- 0
      Denominator_X_ortho <- X.L2.ortho * X %*% XtX_offdiag
    } else {
      Denominator_X_ortho <- 0
    }

    X <- X * (Numerator_X / (Denominator_X_rec + Denominator_X_ortho + .eps))
    X <- .xnorm(X)

    Xt  <- t(X)
    XtX <- Xt %*% X

    # 2.2 update C1
    Numerator_C1   <- Xt %*% Y1 %*% t(Y1)
    Denominator_C1 <- XtX %*% (C1 %*% Y1 + C2 %*% Y2) %*% t(Y1) + C1.L1 + .eps
    C1 <- C1 * (Numerator_C1 / Denominator_C1)

    # 2.3 update C2
    Numerator_C2   <- Xt %*% Y1 %*% t(Y2)
    Denominator_C2 <- XtX %*% (C1 %*% Y1 + C2 %*% Y2) %*% t(Y2) + C2.L1 + .eps
    C2 <- C2 * (Numerator_C2 / Denominator_C2)

    # loss + penalties
    XB <- X %*% (C1 %*% Y1 + C2 %*% Y2)
    loss_rec <- sum((Y1 - XB)^2)
    objfunc[it] <- loss_rec

    if (X.L2.ortho > 0) {
      XtX_off <- XtX
      diag(XtX_off) <- 0
      pen_X_ortho <- 0.5 * X.L2.ortho * sum(XtX_off^2)
    } else {
      pen_X_ortho <- 0
    }
    pen_C1_L1 <- C1.L1 * sum(C1)
    pen_C2_L1 <- C2.L1 * sum(C2)
    objfunc.full[it] <- loss_rec + pen_X_ortho + pen_C1_L1 + pen_C2_L1

    if (it >= 10) {
      epsilon_iter <- abs(objfunc[it] - objfunc[it - 1]) / pmax(abs(objfunc[it]), 1)
      if (epsilon_iter <= epsilon) break
    }
  }
  ## Warn when the MU loop exhausts maxit without meeting the
  ## relative-tolerance criterion (matches the nmfkc() convention).
  if (it == maxit && exists("epsilon_iter") && epsilon_iter > abs(epsilon))
    warning(paste0("maximum iterations (", maxit, ") reached..."))

  # ------------------ reorder factors (nmfkc centroid order) ----------
  centroid <- as.numeric((1:nrow(X)) / nrow(X)) %*% X
  index <- order(centroid)
  X  <- X[, index, drop = FALSE]
  C1 <- C1[index, , drop = FALSE]
  C2 <- C2[index, , drop = FALSE]

  # ------------------------------ names -------------------------------
  colnames(X)  <- Basis_labels
  rownames(C1) <- Basis_labels
  rownames(C2) <- Basis_labels
  rownames(X)  <- Y1_labels
  colnames(C1) <- Y1_labels
  colnames(C2) <- Y2_labels

  # -------------------- feedback + stability diagnostics --------------
  XC1  <- X %*% C1
  eigs <- eigen(XC1, only.values = TRUE)$values
  rho  <- max(abs(eigs))
  if (rho >= 1)
    warning("Leontief.inv may be unstable; spectral radius >= 1.")

  XC1_norm1 <- mat1norm(XC1)

  # -------------------- Leontief inverse + equilibrium mapping --------
  I_mat <- diag(nrow(XC1))
  XC2   <- X %*% C2

  Leontief.inv <- tryCatch(
    base::solve(I_mat - XC1),
    error = function(e) {
      warning("Failed to compute Leontief.inv via solve(I - XC1). Returning NA matrices.")
      matrix(NA_real_, nrow = nrow(XC1), ncol = ncol(XC1))
    }
  )
  M.model <- Leontief.inv %*% XC2

  amplification <- mat1norm(M.model) / (mat1norm(XC2) + .eps)
  amplification.bound <- if (XC1_norm1 < 1) 1 / (1 - XC1_norm1) else Inf

  # -------------------- fit indices (equilibrium prediction) ----------
  if (anyNA(Leontief.inv)) {
    Y1_hat <- matrix(NA_real_, nrow = nrow(Y1), ncol = ncol(Y1))
    SC.cov <- NA_real_
    MAE    <- NA_real_
  } else {
    Y1_hat   <- M.model %*% Y2
    S.sample <- Y1 %*% t(Y1)
    S.model  <- Y1_hat %*% t(Y1_hat)
    SC.cov   <- stats::cor(as.numeric(S.sample), as.numeric(S.model))
    MAE      <- mean(abs(Y1 - Y1_hat))
  }

  ## -------------------- input-output structural fidelity (SC.map) -----
  ## SC.map = cor(vec(M.model), vec(M.simple)), where M.simple = X0 * Theta0
  ## is the feedforward baseline mapping (Satoh 2025 §4.SC.map).
  ## `baseline_for_scmap` was set above to either:
  ##   - the result of an internal nmfkc call (default auto path), or
  ##   - the user-supplied `nmfkc.baseline` list, or
  ##   - NULL if the user opted out (nmfkc.baseline = FALSE) or
  ##     supplied X.init as a numeric matrix without nmfkc.baseline.
  ## A user-supplied `M.simple` matrix takes precedence over both.
  SC.map <- NA_real_
  M.simple <- if (!is.null(user_M.simple)) user_M.simple
              else if (!is.null(baseline_for_scmap))
                baseline_for_scmap$X %*% baseline_for_scmap$C
              else NULL
  if (!is.null(M.simple) && !anyNA(M.model)) {
    M.simple <- as.matrix(M.simple)
    if (all(dim(M.simple) == dim(M.model))) {
      SC.map <- tryCatch(
        stats::cor(as.numeric(M.simple), as.numeric(M.model)),
        error = function(e) NA_real_
      )
    } else {
      warning("M.simple has dimension ", paste(dim(M.simple), collapse = "x"),
              " but M.model is ", paste(dim(M.model), collapse = "x"),
              "; SC.map not computed.")
    }
  }

  out <- list(
    X                   = X,
    C1                  = C1,
    C2                  = C2,
    XC1                 = XC1,
    XC2                 = XC2,
    XC1.radius          = rho,
    XC1.norm1           = XC1_norm1,
    Leontief.inv        = Leontief.inv,
    M.model             = M.model,
    amplification       = amplification,
    amplification.bound = amplification.bound,
    Q                   = Q,
    SC.cov              = SC.cov,
    SC.map              = SC.map,
    MAE                 = MAE,
    ## Effective rank of the latent scores B = C1 Y1 + C2 Y2 (Q x N).
    effective.rank      = .effective.rank(C1 %*% Y1 + C2 %*% Y2),
    objfunc             = objfunc[1:it],
    objfunc.full        = objfunc.full[1:it],
    iter                = it
  )
  ## Carry both the canonical NMF-FFB class (paper-aligned, primary) and
  ## the legacy "nmf.sem" class (back-compat).  S3 methods registered on
  ## either class are dispatched correctly via inheritance.
  class(out) <- c("nmf.ffb", "nmf.sem", "nmf")
  out
}


#' @title Statistical inference for NMF-FFB via X-fixed full pair bootstrap
#' @description
#' \code{nmf.sem.inference} performs statistical inference on the structural
#' coefficient matrices \eqn{C_1} (latent feedback, \eqn{\Theta_1}) and
#' \eqn{C_2} (exogenous loading, \eqn{\Theta_2}) from a fitted
#' \code{\link{nmf.sem}} model.
#'
#' The procedure is a \strong{full pair bootstrap} that holds the basis
#' matrix \eqn{\hat X} from the original fit fixed across all replicates
#' (which avoids label switching and gives a clean conditional
#' interpretation: ``uncertainty of the structural coefficients given the
#' measurement model''):
#' \enumerate{
#'   \item For each replicate \eqn{b = 1, \dots, B}, resample column indices
#'     \eqn{(i_1, \dots, i_N)} with replacement from \eqn{\{1, \dots, N\}}
#'     and form \eqn{Y_1^{(b)} = Y_1[, i]}, \eqn{Y_2^{(b)} = Y_2[, i]}.
#'   \item Re-estimate \eqn{(C_1^{(b)}, C_2^{(b)})} by running the
#'     \code{\link{nmf.sem}} multiplicative updates \emph{with \eqn{X = \hat X}
#'     held fixed} (no \eqn{X} update; no centroid sort), using the same
#'     \code{C1.L1}, \code{C2.L1} as the original fit.
#'   \item Discard replicates that violate stationarity
#'     (\eqn{\rho(X C_1^{(b)}) \ge 1}) or have an amplification ratio
#'     exceeding the geometric-series bound by more than 1\%.
#' }
#' Because \eqn{C_1, C_2 \ge 0} are non-negative by construction, exact
#' zeros are essentially never observed in the bootstrap distribution.
#' Significance is assessed via a \strong{support rate} at a small display
#' threshold \eqn{\delta} (default \code{0.01}):
#' \deqn{
#'   \mathrm{sup}(c) \;=\; \frac{1}{|\mathrm{valid}|}
#'     \sum_{b \in \mathrm{valid}} \mathbf{1}\!\left( \hat c^{(b)} > \delta \right).
#' }
#' This is a one-sided counterpart of the classical \eqn{p}-value:
#' large \code{support_rate} indicates strong evidence that the entry is
#' meaningfully positive.  Significance markers follow the lavaan
#' convention with the natural correspondence \eqn{p = 1 - \mathrm{sup}}:
#' \code{*} (sup > 0.95), \code{**} (sup > 0.99), \code{***}
#' (sup > 0.999).  Cutoffs use strict greater-than so the rule
#' mirrors the standard R convention for p-values (p < 0.05 / 0.01 /
#' 0.001 → */**/***), translated to support_rate via
#' \eqn{\mathrm{sup} = 1 - p}.
#'
#' @param object A fitted object returned by \code{\link{nmf.sem}}.  Must
#'   contain \code{X}, \code{C1}, \code{C2}.
#' @param Y1 Endogenous variable matrix (P1 x N).  Must match the data
#'   used in \code{nmf.sem()}.
#' @param Y2 Exogenous variable matrix (P2 x N).  Same.
#' @param B Number of bootstrap replicates.  Default \code{1000}; required
#'   for the \code{***} threshold (sup > 0.999).  Reduce to 500 for
#'   exploratory speed (only \code{*} / \code{**} stay reliable).
#' @param threshold Display threshold \eqn{\delta} for the support rate
#'   \eqn{\Pr_{\mathrm{boot}}(\hat c^{(b)} > \delta)}.  Default
#'   \code{0.01}; entries below this magnitude are treated as effectively
#'   zero in the path diagram.
#' @param ci.level Confidence level for the percentile bootstrap CI.
#'   Default \code{0.95}.
#' @param C1.L1,C2.L1 L1 sparsity penalties used by the original
#'   \code{\link{nmf.sem}} fit.  These must match the fit's hyperparameters
#'   for the bootstrap to estimate the correct model.  Defaults
#'   (\code{1.0}, \code{0.1}) match \code{nmf.sem}'s defaults but you
#'   should pass the actual values used.
#' @param seed Base RNG seed for the bootstrap.  Each replicate uses
#'   \code{seed + b} (resampling) and \code{seed + 1000 + b}
#'   (\eqn{C_1, C_2} initialization).  Default \code{123}.
#' @param ... Hidden options:
#'   \describe{
#'     \item{\code{epsilon}}{Convergence tolerance for the inner fixed-X MU
#'       loop.  Default \code{1e-6}.}
#'     \item{\code{maxit}}{Maximum iterations for the inner MU loop.
#'       Default \code{5000}.}
#'     \item{\code{ncores}}{Number of parallel workers.  Default \code{1}
#'       (serial).  Cross-platform: uses \code{parallel::mclapply} on
#'       Linux/macOS and \code{parallel::parLapply} (PSOCK cluster) on
#'       Windows.}
#'     \item{\code{print.trace}}{Logical, print progress.  Default
#'       \code{FALSE}.}
#'   }
#'
#' @return The input \code{object} with additional bootstrap inference
#'   components:
#' \item{coefficients}{Data frame with rows for every entry of \eqn{C_1}
#'   and \eqn{C_2} and columns \code{Type} ("C1" / "C2"), \code{Basis},
#'   \code{Covariate}, \code{Estimate}, \code{CI_low}, \code{CI_high},
#'   \code{support_rate}, \code{p_value} (\eqn{= 1 - \mathrm{support\_rate}},
#'   for compatibility with downstream consumers such as
#'   \code{\link{nmf.sem.DOT}}), and \code{sig}.}
#' \item{C1.support.rate, C2.support.rate}{Per-element support rates
#'   (Q x P1 and Q x P2 matrices).}
#' \item{C1.ci.lower, C1.ci.upper, C2.ci.lower, C2.ci.upper}{Per-element
#'   percentile CI bounds.}
#' \item{C1.array, C2.array}{Bootstrap distributions: 3D arrays of shape
#'   B x Q x P1 (and B x Q x P2).  Invalid replicates contain \code{NA}.}
#' \item{rho.boot, AR.boot, iter.boot}{Per-replicate spectral radius,
#'   amplification ratio, and inner-loop iteration count.}
#' \item{bootstrap.B, bootstrap.threshold, bootstrap.ci.level}{Inputs
#'   recorded for reproducibility.}
#' \item{bootstrap.n.valid, bootstrap.n.invalid}{Validity counts.}
#'
#' @section Lifecycle:
#' This function's interface changed at v0.6.8: the legacy 1-step Newton
#' wild bootstrap (with sandwich SE) has been replaced by the full pair
#' bootstrap described above, following the paper revision.  The fields
#' \code{sigma2.used}, \code{C2.se}, \code{C2.se.boot}, \code{C2.p.side}
#' that the previous implementation produced are no longer present.
#'
#' @seealso \code{\link{nmf.sem}}, \code{\link{nmf.sem.DOT}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to structural equation modeling for blind input-output analysis.
#'   arXiv:2512.18250. \url{https://arxiv.org/abs/2512.18250}
#' @export
#' @examples
#' \donttest{
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
#' res  <- nmf.sem(Y1, Y2, rank = 2)
#' res2 <- nmf.sem.inference(res, Y1, Y2, B = 200)  # quick demo
#' head(res2$coefficients)
#' }
nmf.ffb.inference <- function(object, Y1, Y2,
                               B = 1000L,
                               threshold = 0.01,
                               ci.level = 0.95,
                               C1.L1 = 1.0,
                               C2.L1 = 0.1,
                               seed = 123L,
                               ...) {
  if (is.null(object$X) || is.null(object$C1) || is.null(object$C2))
    stop("object must contain X, C1, and C2 (returned by nmf.sem).")

  extra_args  <- base::list(...)
  epsilon     <- if (!is.null(extra_args$epsilon))     extra_args$epsilon     else 1e-6
  maxit       <- if (!is.null(extra_args$maxit))       extra_args$maxit       else 5000L
  ncores      <- if (!is.null(extra_args$ncores))      extra_args$ncores      else 1L
  print.trace <- if (!is.null(extra_args$print.trace)) extra_args$print.trace else FALSE

  Y1 <- base::as.matrix(Y1)
  Y2 <- base::as.matrix(Y2)
  if (ncol(Y1) != ncol(Y2)) stop("Y1 and Y2 must have the same number of columns.")

  X.hat  <- object$X
  C1.hat <- object$C1
  C2.hat <- object$C2

  Q  <- ncol(X.hat)
  P1 <- nrow(Y1)
  P2 <- nrow(Y2)
  N  <- ncol(Y1)
  if (nrow(X.hat) != P1)
    stop("nrow(X) must equal nrow(Y1); did Y1 change shape?")
  if (Q < 1L) stop("Inferred rank Q < 1.")

  ## ----------------------------------------------------------------
  ## one_boot: a single bootstrap replicate.
  ## Closure captures Y1, Y2, X.hat, C1.L1, C2.L1, Q, P1, P2, N,
  ## epsilon, maxit, seed.  Returns a list with valid flag plus
  ## (C1, C2, rho, AR, iter) when valid, or NA placeholders otherwise.
  ## ----------------------------------------------------------------
  one_boot <- function(b) {
    ## Sample column indices with replacement (pair bootstrap)
    set.seed(seed + b)
    idx  <- sample.int(N, size = N, replace = TRUE)
    Y1_b <- Y1[, idx, drop = FALSE]
    Y2_b <- Y2[, idx, drop = FALSE]

    ## Independent random init for C1, C2 (X is held at X.hat)
    set.seed(seed + 1000L + b)
    C1 <- matrix(stats::runif(Q * P1, 0.01, 0.1),  nrow = Q)
    C2 <- matrix(stats::runif(Q * P2, 0.001, 0.01), nrow = Q, ncol = P2)
    min_dim <- min(Q, P2)
    if (min_dim >= 1L) {
      for (i in 1:min_dim) C2[i, i] <- stats::runif(1, 0.1, 0.2)
    }

    .eps_local <- 1e-10
    Xt  <- t(X.hat)
    XtX <- Xt %*% X.hat
    Y1tY1_b <- tcrossprod(Y1_b)        # P1 x P1
    XtY1_b  <- Xt %*% Y1_b              # Q x N (re-used via crossprods below)
    XtY1Y1t <- XtY1_b %*% t(Y1_b)       # Q x P1   = Xt %*% Y1_b %*% t(Y1_b)
    XtY1Y2t <- XtY1_b %*% t(Y2_b)       # Q x P2
    prev_loss <- Inf
    iter_used <- 0L

    ## Fixed-X MU loop on (C1, C2)
    for (it in seq_len(maxit)) {
      M_b <- C1 %*% Y1_b + C2 %*% Y2_b   # Q x N latent representation

      ## C1 update
      Den_C1 <- XtX %*% M_b %*% t(Y1_b) + C1.L1 + .eps_local
      C1 <- C1 * (XtY1Y1t / Den_C1)

      ## C2 update (uses updated C1 via fresh M_b)
      M_b <- C1 %*% Y1_b + C2 %*% Y2_b
      Den_C2 <- XtX %*% M_b %*% t(Y2_b) + C2.L1 + .eps_local
      C2 <- C2 * (XtY1Y2t / Den_C2)

      ## Convergence check (relative change in reconstruction loss)
      XB <- X.hat %*% (C1 %*% Y1_b + C2 %*% Y2_b)
      loss <- sum((Y1_b - XB)^2)
      if (it >= 10L) {
        rel <- abs(loss - prev_loss) / max(abs(loss), 1)
        if (rel <= epsilon) { iter_used <- it; break }
      }
      prev_loss <- loss
      iter_used <- it
    }

    ## ---- Validity checks (mirrors the paper's bootstrap helper) ----
    XC1_b <- X.hat %*% C1
    rho <- tryCatch(
      max(abs(eigen(XC1_b, only.values = TRUE, symmetric = FALSE)$values)),
      error = function(e) NA_real_
    )
    norm1_XC1 <- max(colSums(abs(XC1_b)))
    AR_bound <- if (is.finite(norm1_XC1) && norm1_XC1 < 1)
                  1 / (1 - norm1_XC1)
                else Inf

    if (!is.finite(rho) || rho >= 1) {
      return(list(valid = FALSE, C1 = NULL, C2 = NULL,
                  rho = rho, AR = NA_real_, iter = iter_used))
    }
    ## XC1_b = X (P1 x Q) %*% C1 (Q x P1) is P1 x P1, so the identity
    ## here must be P1 x P1 too.  (Earlier draft used diag(Q) which is
    ## conformable only when Q == P1, otherwise solve() throws an error
    ## that the tryCatch swallowed -- causing every replicate to be
    ## marked invalid with AR = NA.)
    I_mat   <- diag(nrow(XC1_b))
    Linv    <- tryCatch(solve(I_mat - XC1_b), error = function(e) NULL)
    if (is.null(Linv)) {
      return(list(valid = FALSE, C1 = NULL, C2 = NULL,
                  rho = rho, AR = NA_real_, iter = iter_used))
    }
    XC2_b   <- X.hat %*% C2
    M.model <- Linv %*% XC2_b
    AR <- max(colSums(abs(M.model))) / (max(colSums(abs(XC2_b))) + .eps_local)
    if (!is.finite(AR) || AR > AR_bound * 1.01) {
      return(list(valid = FALSE, C1 = NULL, C2 = NULL,
                  rho = rho, AR = AR, iter = iter_used))
    }

    list(valid = TRUE, C1 = C1, C2 = C2,
         rho = rho, AR = AR, iter = iter_used)
  }

  ## ----------------------------------------------------------------
  ## Run the B replicates (parallel on Linux/macOS via mclapply,
  ## PSOCK cluster on Windows; serial when ncores == 1).
  ## ----------------------------------------------------------------
  if (print.trace)
    base::message(sprintf("  Bootstrap: B=%d, ncores=%d, threshold=%.3g, ci.level=%.2f",
                          B, ncores, threshold, ci.level))

  if (ncores > 1L) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl,
        varlist = c("Y1", "Y2", "X.hat", "C1.L1", "C2.L1",
                    "Q", "P1", "P2", "N", "epsilon", "maxit", "seed"),
        envir = environment())
      res_list <- parallel::parLapply(cl, seq_len(B), one_boot)
    } else {
      res_list <- parallel::mclapply(seq_len(B), one_boot, mc.cores = ncores)
    }
  } else {
    res_list <- lapply(seq_len(B), one_boot)
  }

  ## ----------------------------------------------------------------
  ## Aggregate replicates into 3D arrays (B x Q x P1) and (B x Q x P2)
  ## ----------------------------------------------------------------
  C1.array <- array(NA_real_, dim = c(B, Q, P1))
  C2.array <- array(NA_real_, dim = c(B, Q, P2))
  rho.vec  <- rep(NA_real_, B)
  AR.vec   <- rep(NA_real_, B)
  iter.vec <- rep(NA_integer_, B)
  valid.vec <- logical(B)

  for (b in seq_len(B)) {
    r <- res_list[[b]]
    if (is.null(r) || inherits(r, "try-error")) next
    valid.vec[b] <- isTRUE(r$valid)
    if (!is.null(r$rho))  rho.vec[b]  <- r$rho
    if (!is.null(r$AR))   AR.vec[b]   <- r$AR
    if (!is.null(r$iter)) iter.vec[b] <- as.integer(r$iter)
    if (valid.vec[b]) {
      C1.array[b, , ] <- r$C1
      C2.array[b, , ] <- r$C2
    }
  }
  n.valid <- sum(valid.vec)

  if (n.valid < 10L) {
    warning(sprintf("Only %d / %d bootstrap replicates were valid; CIs / support rates may be unreliable. Consider checking convergence (raise maxit, lower epsilon) or the stationarity of the original fit (rho < 1).",
                    n.valid, B))
  }

  ## ----------------------------------------------------------------
  ## Per-element summary: support rate at threshold, percentile CI
  ## ----------------------------------------------------------------
  alpha <- 1 - ci.level

  apply_finite <- function(arr, FUN) {
    apply(arr, c(2, 3), function(v) {
      v <- v[is.finite(v)]
      if (!length(v)) NA_real_ else FUN(v)
    })
  }

  C1.support  <- apply_finite(C1.array, function(v) mean(v > threshold))
  C2.support  <- apply_finite(C2.array, function(v) mean(v > threshold))
  C1.ci.lower <- apply_finite(C1.array, function(v) stats::quantile(v, alpha / 2,    names = FALSE))
  C1.ci.upper <- apply_finite(C1.array, function(v) stats::quantile(v, 1 - alpha / 2, names = FALSE))
  C2.ci.lower <- apply_finite(C2.array, function(v) stats::quantile(v, alpha / 2,    names = FALSE))
  C2.ci.upper <- apply_finite(C2.array, function(v) stats::quantile(v, 1 - alpha / 2, names = FALSE))

  ## Significance markers from support rate (one-sided; lavaan / stats
  ## convention).  Cutoffs use strict greater-than so the rule mirrors
  ## the R convention for p-values (p < 0.05 / 0.01 / 0.001 → */**/***),
  ## translated to support_rate via support_rate = 1 - p.  This matches
  ## the verbal description "support rates exceeding 0.95, 0.99, 0.999"
  ## in Satoh (2025), arXiv:2512.18250.
  sig.from.support <- function(s) {
    ifelse(!is.finite(s), " ",
      ifelse(s > 0.999, "***",
        ifelse(s > 0.99,  "**",
          ifelse(s > 0.95, "*", " "))))
  }

  ## ----------------------------------------------------------------
  ## Build coefficients table (rows for both C1 and C2)
  ## p_value = 1 - support_rate is provided so that downstream
  ## consumers (e.g., nmf.sem.DOT) that filter / star by p_value
  ## continue to work without modification.
  ## ----------------------------------------------------------------
  Q_lab  <- if (!is.null(rownames(C1.hat))) rownames(C1.hat) else paste0("Factor", 1:Q)
  Y1_lab <- if (!is.null(colnames(C1.hat))) colnames(C1.hat) else paste0("Y1_", 1:P1)
  Y2_lab <- if (!is.null(colnames(C2.hat))) colnames(C2.hat) else paste0("Y2_", 1:P2)

  build_block <- function(type_label, Mhat, sup, lo, hi, basislabs, varlabs) {
    Q_local <- nrow(Mhat); P_local <- ncol(Mhat)
    s_vec <- as.vector(sup)
    data.frame(
      Type         = type_label,
      Basis        = rep(basislabs, times = P_local),
      Covariate    = rep(varlabs,   each  = Q_local),
      Estimate     = as.vector(Mhat),
      CI_low       = as.vector(lo),
      CI_high      = as.vector(hi),
      support_rate = s_vec,
      p_value      = ifelse(is.finite(s_vec), 1 - s_vec, NA_real_),
      sig          = sig.from.support(s_vec),
      stringsAsFactors = FALSE
    )
  }

  C1_block <- build_block("C1", C1.hat, C1.support, C1.ci.lower, C1.ci.upper, Q_lab, Y1_lab)
  C2_block <- build_block("C2", C2.hat, C2.support, C2.ci.lower, C2.ci.upper, Q_lab, Y2_lab)
  coefficients <- rbind(C1_block, C2_block)
  rownames(coefficients) <- NULL

  if (print.trace)
    base::message(sprintf("  Bootstrap done: %d / %d valid replicates.", n.valid, B))

  ## ----------------------------------------------------------------
  ## Append to object and return
  ## ----------------------------------------------------------------
  object$bootstrap.B          <- B
  object$bootstrap.threshold  <- threshold
  object$bootstrap.ci.level   <- ci.level
  object$bootstrap.n.valid    <- n.valid
  object$bootstrap.n.invalid  <- B - n.valid
  object$rho.boot             <- rho.vec
  object$AR.boot              <- AR.vec
  object$iter.boot            <- iter.vec
  object$C1.array             <- C1.array
  object$C2.array             <- C2.array
  object$C1.support.rate      <- C1.support
  object$C2.support.rate      <- C2.support
  object$C1.ci.lower          <- C1.ci.lower
  object$C1.ci.upper          <- C1.ci.upper
  object$C2.ci.lower          <- C2.ci.lower
  object$C2.ci.upper          <- C2.ci.upper
  object$coefficients         <- coefficients

  ## Add NMF-FFB inference class on top of the legacy SEM class.
  ## Final class vector for a typical input:
  ##   c("nmf.ffb.inference", "nmf.sem.inference", "nmf.ffb", "nmf.sem")
  ## Existing S3 methods (e.g. summary.nmf.sem) still dispatch via
  ## inheritance.
  if (!inherits(object, "nmf.sem.inference"))
    class(object) <- c("nmf.sem.inference", class(object))
  if (!inherits(object, "nmf.ffb.inference"))
    class(object) <- c("nmf.ffb.inference", class(object))
  object
}


#' @title Cross-Validation for NMF-FFB
#' @description
#' Performs K-fold cross-validation to evaluate the equilibrium mapping of
#' the NMF-FFB model.
#'
#' For each fold, \code{nmf.sem} is fitted on the training samples,
#' yielding an equilibrium mapping \eqn{\hat Y_1 = M_{\mathrm{model}} Y_2}.
#' The held-out endogenous variables \eqn{Y_1} are then predicted from \eqn{Y_2}
#' using this mapping, and the mean absolute error (MAE) over all entries in the
#' test block is computed. The returned value is the average MAE across folds.
#'
#' This implements the hyperparameter selection strategy described in the paper:
#' hyperparameters are chosen by predictive cross-validation rather than direct
#' inspection of the internal structural matrices.
#'
#' @param Y1 A non-negative numeric matrix of endogenous variables with
#'   \strong{rows = variables (P1), columns = samples (N)}.
#' @param Y2 A non-negative numeric matrix of exogenous variables with
#'   \strong{rows = variables (P2), columns = samples (N)}.
#'   Must satisfy \code{ncol(Y1) == ncol(Y2)}.
#' @param rank Integer; rank (number of latent factors) passed to \code{nmf.sem}.
#'   If \code{NULL}, \code{nmf.sem} decides the effective rank (via \code{...} or \code{nrow(Y2)}).
#' @param X.init Initialization strategy for \code{X}, forwarded to
#'   \code{\link{nmf.sem}}.  One of \code{"nndsvd"} (default),
#'   \code{"kmeans"}, \code{"kmeansar"}, \code{"runif"}, a numeric
#'   \eqn{P_1 \times Q} matrix, or \code{NULL} (alias for
#'   \code{"nndsvd"}).  See \code{\link{nmf.sem}} for details.
#' @param X.L2.ortho L2 orthogonality penalty for \code{X}.
#' @param C1.L1 L1 sparsity penalty for \code{C1} (\eqn{\Theta_1}).
#' @param C2.L1 L1 sparsity penalty for \code{C2} (\eqn{\Theta_2}).
#' @param epsilon Convergence threshold for \code{nmf.sem}.
#' @param maxit Maximum number of iterations for \code{nmf.sem}.
#' @param ... Additional arguments passed to \code{nmf.sem} (except for
#'   \code{rank}, \code{seed}, \code{div}, \code{shuffle}, which are handled here).
#'   Also accepts: \code{nfolds} (number of folds, default 5; \code{div} also accepted),
#'   \code{seed} (master random seed, default \code{NULL}),
#'   \code{shuffle} (logical, default \code{TRUE}).
#'
#' @return A numeric scalar: mean MAE across CV folds.
#'
#' @examples
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]
#' Y2 <- Y[3:4, ]
#' mae <- nmf.sem.cv(Y1, Y2, rank = 2, maxit = 500, nfolds = 3)
#' mae
#'
#' @seealso \code{\link{nmf.sem}}
#' @export
nmf.ffb.cv <- function(
    Y1, Y2,
    rank = NULL,
    X.init = "nndsvd",
    X.L2.ortho = 100.0,
    C1.L1 = 1.0,        # L1 sparsity for C1 (Theta1)
    C2.L1 = 0.1,        # L1 sparsity for C2 (Theta2)
    epsilon = 1e-6,     # Convergence tolerance passed to nmf.ffb
    maxit = 5000,
    ...
){
  extra_cv <- base::list(...)
  nfolds  <- if (!is.null(extra_cv$nfolds))  extra_cv$nfolds  else if (!is.null(extra_cv$div)) extra_cv$div else 5
  seed    <- if (!is.null(extra_cv$seed))    extra_cv$seed    else NULL
  shuffle <- if (!is.null(extra_cv$shuffle)) extra_cv$shuffle else TRUE
  div <- nfolds
  # ------------------------------------------------------------------
  # 1. Basic input checks
  #
  # NMF-FFB requires non-negative matrices. We also require that Y1 and Y2
  # share the same number of samples (columns) to allow paired CV splits.
  # ------------------------------------------------------------------
  if (!is.matrix(Y1)) Y1 <- as.matrix(Y1)
  if (!is.matrix(Y2)) Y2 <- as.matrix(Y2)
  if (any(!is.finite(Y1)) || any(!is.finite(Y2))) stop("Y1 and Y2 must not contain NA/NaN/Inf.")
  div <- as.integer(div)
  if (min(Y1) < 0 || min(Y2) < 0) stop("Y1 and Y2 must be non-negative.")
  if (ncol(Y1) != ncol(Y2)) {
    stop("Y1 and Y2 must have the same number of columns (samples).")
  }

  P1 <- nrow(Y1)
  P2 <- nrow(Y2)
  N  <- ncol(Y1)

  if (div < 2L) {
    stop("div (number of CV folds) must be >= 2.")
  }
  if (div > N) {
    stop("div (number of CV folds) must be <= number of samples.")
  }

  # ------------------------------------------------------------------
  # 2. Handle extra arguments for nmf.sem
  #
  # We collect additional arguments in 'extra_args' but explicitly remove
  # those that are managed at the CV level:
  #   - div, shuffle : used only here, not passed to nmf.sem.
  #   - rank        : passed explicitly from nmf.sem.cv.
  #   - seed        : fold-specific seeds are generated here.
  # ------------------------------------------------------------------
  extra_args <- list(...)

  extra_args$div     <- NULL
  extra_args$nfolds  <- NULL
  extra_args$shuffle <- NULL
  extra_args$rank    <- NULL
  extra_args$seed    <- NULL

  # ------------------------------------------------------------------
  # 3. Set RNG for CV partition and per-fold seeds
  #
  # If a master 'seed' is given:
  #   - it is used to define the CV partition (sample permutation),
  #   - independent seeds for each nmf.sem run are drawn.
  # If 'seed' is NULL:
  #   - CV partition uses the current RNG state,
  #   - nmf.sem runs use whatever the global RNG state is at call time.
  # ------------------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Sample indices for fold division
  if (shuffle) {
    perm_index <- sample.int(N)
  } else {
    perm_index <- seq_len(N)
  }

  # Per-fold seeds for nmf.sem (optional, only if master seed specified)
  if (!is.null(seed)) {
    seeds_fold <- sample.int(.Machine$integer.max, div)
  } else {
    seeds_fold <- rep(NA_integer_, div)
  }

  # ------------------------------------------------------------------
  # 4. Create CV folds
  #
  # We assign approximately N/div samples to each fold, distributing
  # any remainder one-by-one to the earliest folds.
  # 'block[i] = k' means sample i belongs to fold k (as test set).
  # ------------------------------------------------------------------
  remainder <- N %% div
  division  <- N %/% div
  block     <- integer(N)

  processed_count <- 0L
  for (i in 1:(div - 1L)) {
    plus       <- ifelse(i <= remainder, 1L, 0L)
    chunk_size <- division + plus
    idx_range  <- (processed_count + 1L):(processed_count + chunk_size)
    target_idx <- perm_index[idx_range]
    block[target_idx] <- i
    processed_count   <- processed_count + chunk_size
  }
  # Last fold gets all remaining samples
  target_idx <- perm_index[(processed_count + 1L):N]
  block[target_idx] <- div

  # Per-fold CV loss (MAE) will be stored here
  objfunc.block <- numeric(div)

  # ------------------------------------------------------------------
  # 5. Cross-validation loop
  #
  # For each fold j:
  #   - train on all samples not in fold j,
  #   - test on samples in fold j,
  #   - fit nmf.sem on training data,
  #   - compute MAE on test block from equilibrium mapping M.model.
  # ------------------------------------------------------------------
  for (j in 1:div) {
    # Train / Test split
    train_idx <- block != j
    test_idx  <- block == j

    Y1_train <- Y1[, train_idx, drop = FALSE]
    Y1_test  <- Y1[, test_idx,  drop = FALSE]
    Y2_train <- Y2[, train_idx, drop = FALSE]
    Y2_test  <- Y2[, test_idx,  drop = FALSE]

    # Fold-specific seed for nmf.sem
    seed_j <- if (!is.null(seed)) seeds_fold[j] else NULL

    # Assemble arguments for nmf.sem
    nmf.sem.args <- c(
      extra_args,   # User-specified additional arguments (e.g., Q)
      list(
        Y1         = Y1_train,
        Y2         = Y2_train,
        rank       = rank,
        X.init     = X.init,
        X.L2.ortho = X.L2.ortho,
        C1.L1      = C1.L1,
        C2.L1      = C2.L1,
        epsilon    = epsilon,
        maxit      = maxit
      )
    )
    # Attach seed only when it is defined
    if (!is.null(seed_j)) {
      nmf.sem.args$seed <- seed_j
    }

    # Call nmf.sem on the training data (suppress messages for cleaner CV output)
    res_j <- suppressMessages(do.call("nmf.ffb", nmf.sem.args))

    # If mapping is not usable, penalize this fold (do not crash CV)
    if (is.null(res_j$M.model) || any(!is.finite(res_j$M.model))) {
      objfunc.block[j] <- Inf
      next
    }
    Pre_test <- res_j$M.model %*% Y2_test
    if (any(!is.finite(Pre_test))) {
      objfunc.block[j] <- Inf
      next
    }
    objfunc.block[j] <- mean(abs(Y1_test - Pre_test))
  }

  # ------------------------------------------------------------------
  # 6. Aggregate CV score
  #
  # The overall CV criterion is the average MAE across all folds.
  # This is typically minimized over hyperparameter grids
  # (e.g., rank, X.L2.ortho, C1.L1, C2.L1) when tuning NMF-FFB.
  # ------------------------------------------------------------------
  objfunc <- mean(objfunc.block)
  return(objfunc)
}



#' @title Heuristic Variable Splitting for NMF-FFB
#'
#' @description
#' Infers a heuristic partition of observed variables into exogenous (\eqn{Y_2})
#' and endogenous (\eqn{Y_1}) blocks for use in NMF-FFB.
#' The method is based on positive-SEM logic, causal ordering, and optional
#' sign alignment using the first principal component (PC1).
#'
#' The procedure:
#' \itemize{
#'   \item internally standardizes variables (mean 0, sd 1),
#'   \item optionally flips signs so that most variables align positively with PC1,
#'   \item infers a causal ordering by repeatedly regressing each variable on the
#'         remaining ones and selecting the variable with the largest minimum
#'         standardized coefficient,
#'   \item determines an exogenous block by scanning the ordering from upstream
#'         and stopping at the first variable whose strongest parent coefficient
#'         exceeds \code{threshold}.
#' }
#'
#' If \code{n.exogenous} is supplied, it overrides the automatic threshold rule.
#'
#' @param x A numeric matrix or data frame with
#'   \strong{rows = samples} and \strong{columns = observed variables}.
#' @param n.exogenous Optional integer specifying the number of exogenous variables
#'   (\eqn{Y_2}). If \code{NULL}, the number is inferred automatically by the
#'   coefficient cut-off rule.
#' @param threshold Standardized regression-coefficient threshold used in the
#'   automatic exogenous–endogenous split. A variable is treated as endogenous
#'   once its maximum standardized parent coefficient exceeds this value.
#'   (Default: \code{0.1})
#' @param auto.flipped Logical; if \code{TRUE}, applies PC1-based automatic
#'   sign flipping after standardization to ensure consistent orientation.
#'   (Default: \code{TRUE})
#' @param verbose Logical; if \code{TRUE}, prints progress messages and the
#'   resulting variable split. (Default: \code{FALSE})
#' @param ... Reserved for future use; currently unused (also accepted
#'   by the \code{\link{nmf.ffb.split}} alias for argument forwarding).
#'
#' @return A list with:
#'   \item{endogenous.variables}{
#'     Character vector of variables selected as endogenous (\eqn{Y_1}).}
#'   \item{exogenous.variables}{
#'     Character vector of variables selected as exogenous (\eqn{Y_2}).}
#'   \item{ordered.variables}{
#'     Variables in inferred causal order (from exogenous to endogenous).}
#'   \item{is.flipped}{
#'     Logical vector indicating which variables were sign-flipped during processing.}
#'   \item{n.exogenous}{
#'     Integer giving the number of exogenous variables.}
#'
#' @examples
#' # Infer exogenous/endogenous split from iris
#' sp <- nmf.sem.split(iris[, -5], n.exogenous = 2)
#' sp$endogenous.variables
#' sp$exogenous.variables
#'
#' @seealso \code{\link{nmf.sem}}
#' @export
nmf.ffb.split <- function(x, n.exogenous = NULL, threshold = 0.1,
                          auto.flipped = TRUE, verbose = FALSE, ...) {

  if (!is.matrix(x) && !is.data.frame(x))
    stop("x must be a numeric matrix or data frame.")

  X_raw <- as.matrix(x)
  P <- ncol(X_raw)
  col_names <- colnames(X_raw)
  if (is.null(col_names)) {
    col_names <- paste0("V", 1:P)
    colnames(X_raw) <- col_names
  }

  # --------------------------------------------------------------------
  # Preprocessing Step 1: Standardize all variables
  #
  # Variables are centered and scaled (mean 0, sd 1). NMF-FFB requires
  # non-negative matrices, but the purpose of this function is only to
  # infer variable roles (Y1/Y2), so standardized values are allowed.
  #
  # Missing or NaN values resulting from constant columns are set to 0.
  # --------------------------------------------------------------------
  X_calc <- scale(X_raw, center = TRUE, scale = TRUE)
  X_calc[is.na(X_calc)] <- 0
  X_calc[is.nan(X_calc)] <- 0

  all_indices <- 1:P

  # --------------------------------------------------------------------
  # Preprocessing Step 2: Optional sign flipping based on PC1 alignment
  #
  # In positive SEM (and NMF-FFB), variables should ideally have
  # consistent sign orientation. To enforce this heuristic, variables
  # negatively correlated with the first principal component are flipped.
  #
  # This stabilizes the causal-ordering heuristic by avoiding mixtures
  # of arbitrary sign conventions in the raw data.
  # --------------------------------------------------------------------
  is.flipped <- rep(FALSE, P)
  names(is.flipped) <- col_names

  if (auto.flipped) {
    if (verbose) message("Step 0: Checking correlations with PC1 (on standardized data)...")

    svd_res <- svd(X_calc)
    pc1 <- svd_res$u[, 1]

    cors <- stats::cor(X_calc, pc1)
    # Ensure majority alignment with PC1
    if (stats::median(cors, na.rm = TRUE) < 0) cors <- -cors

    flip_idx <- which(cors < 0)

    if (length(flip_idx) > 0) {
      is.flipped[flip_idx] <- TRUE
      X_calc[, flip_idx] <- -X_calc[, flip_idx]

      if (verbose) {
        message(sprintf("   -> Detected %d flipped variables: %s",
                    length(flip_idx), paste(col_names[flip_idx], collapse=", ")))
      }
    }
  }

  # --------------------------------------------------------------------
  # Step 1: Causal ordering heuristic
  #
  # We infer an ordering of variables consistent with positive-SEM logic:
  # repeatedly select the variable that has the *largest minimum* coefficient
  # when regressed on the remaining variables. This favors variables that
  # are least explained by others → likely exogenous.
  #
  # The resulting order approximates a causal topological order in which
  # exogenous variables appear early and endogenous variables later.
  # --------------------------------------------------------------------
  if (verbose) message("Step 1: Inferring Causal Ordering...")

  active_set <- all_indices
  ordering_reversed <- integer(P)

  for (t in 1:(P - 1)) {
    scores <- numeric(length(active_set))

    for (i in seq_along(active_set)) {
      target_col <- active_set[i]
      pred_cols <- active_set[-i]

      y_vec <- X_calc[, target_col]
      X_mat <- X_calc[, pred_cols, drop = FALSE]

      coefs <- tryCatch({
        stats::coef(stats::lm(y_vec ~ X_mat - 1))
      }, error = function(e) rep(NA_real_, length(pred_cols)))

      # Positive SEM → keep the smallest coefficient (weakest positive predictor)
      if (all(is.na(coefs))) {
        scores[i] <- -Inf
      } else {
        scores[i] <- min(coefs, na.rm = TRUE)
      }
    }

    best_idx <- which.max(scores)
    ordering_reversed[t] <- active_set[best_idx]
    active_set <- active_set[-best_idx]
  }
  ordering_reversed[P] <- active_set[1]

  # Causal order: exogenous → endogenous
  ordering_indices <- rev(ordering_reversed)

  # --------------------------------------------------------------------
  # Step 2: Automatic identification of exogenous variables
  #
  # Sweep through the causal ordering. For each variable, regress it on all
  # earlier variables. If its strongest parent coefficient exceeds the
  # threshold, the variable is considered endogenous.
  #
  # Variables before this point → exogenous (Y2)
  # Variables after this point → endogenous (Y1)
  #
  # If n.exogenous is given, it overrides this automatic rule.
  # --------------------------------------------------------------------
  if (is.null(n.exogenous)) {
    if (verbose) message("Step 2: Detecting optimal cut-off for exogenous variables...")
    cutoff <- 1

    for (k in seq_len(max(P - 2, 0)) + 1) {
      curr_idx <- ordering_indices[k]
      parent_indices <- ordering_indices[seq_len(k - 1)]

      y_vec <- X_calc[, curr_idx]
      if (length(parent_indices) == 0) next
      X_parents <- X_calc[, parent_indices, drop = FALSE]

      coefs <- stats::coef(stats::lm(y_vec ~ X_parents - 1))
      coefs <- coefs[is.finite(coefs)]
      max_influence <- if (length(coefs) > 0) max(coefs) else 0

      if (max_influence > threshold) {
        if (verbose)
          message(sprintf("   [%d] %s : Max std.coef=%.3f -> Endogenous (Stop)",
                      k, col_names[curr_idx], max_influence))
        break
      } else {
        cutoff <- k
        if (verbose)
          message(sprintf("   [%d] %s : Max std.coef=%.3f -> Exogenous (Continue)",
                      k, col_names[curr_idx], max_influence))
      }
    }
    n.exogenous <- cutoff
  }

  # --------------------------------------------------------------------
  # Step 3: Final classification into Y1 and Y2
  #
  # Variables appearing early in the ordering (determined by cut-off) are
  # treated as exogenous (Y2). The remainder are endogenous (Y1).
  #
  # Ordered list shows the full inferred causal sequence.
  # --------------------------------------------------------------------
  if(n.exogenous < 1 || n.exogenous >= P)
    stop("n.exogenous must be between 1 and P-1.")
  idx_exo <- ordering_indices[1:n.exogenous]
  idx_endo <- ordering_indices[(n.exogenous + 1):P]

  exogenous.variables <- col_names[idx_exo]
  endogenous.variables <- col_names[idx_endo]
  ordered.variables <- col_names[ordering_indices]

  if (verbose) {
    message("\n--- Auto Split Result ---")
    message(sprintf("Exogenous (Y2, n=%d): %s",
                n.exogenous, paste(exogenous.variables, collapse=", ")))
    message(sprintf("Endogenous (Y1, n=%d): %s ...",
                length(endogenous.variables),
                paste(utils::head(endogenous.variables, 3), collapse=", ")))
  }

  return(list(
    endogenous.variables = endogenous.variables,
    exogenous.variables = exogenous.variables,
    ordered.variables = ordered.variables,
    is.flipped = is.flipped,
    n.exogenous = as.integer(n.exogenous)
  ))
}




############################################################
## Common DOT Helpers
############################################################

############################################################
## 1. nmf.sem.DOT  (for NMF-FFB visualization)
############################################################

#' Generate a Graphviz DOT Diagram for an NMF-FFB Model
#'
#' @description
#' Creates a Graphviz DOT script that visualizes the structural network
#' estimated by \code{nmf.sem}.
#' The resulting diagram displays:
#' \itemize{
#'   \item endogenous observed variables (\eqn{Y_1}),
#'   \item exogenous observed variables (\eqn{Y_2}),
#'   \item latent factors (\eqn{F_1}, \dots, \eqn{F_Q}),
#' }
#' together with the non-negative path coefficients whose magnitudes
#' exceed a user-specified threshold.
#'
#' Directed edges represent estimated relationships:
#' \itemize{
#'   \item \eqn{Y_2 \rightarrow F_q}: entries of \code{C2} (exogenous loadings),
#'   \item \eqn{F_q \rightarrow Y_1}: rows of \code{X} (factor-to-endogenous mappings),
#'   \item \eqn{Y_1 \rightarrow F_q}: entries of \code{C1} (feedback paths).
#' }
#'
#' Edge widths are scaled by coefficient magnitude, and nodes are placed
#' in optional visual clusters. Only variables participating in
#' edges above the threshold are displayed, while latent factors are always shown.
#'
#' @param result A list returned by \code{nmf.sem}, containing matrices
#'   \code{X}, \code{C1}, and \code{C2}.
#' @param weight_scale Base scaling factor for edge widths.
#' @param weight_scale_c2 Scaling factor for edges
#'   \eqn{Y_2 \rightarrow F_q} (C2 matrix). Defaults to \code{weight_scale}.
#' @param weight_scale_x1 Scaling factor for edges
#'   \eqn{F_q \rightarrow Y_1} (X matrix). Defaults to \code{weight_scale}.
#' @param weight_scale_feedback Scaling factor for feedback edges
#'   \eqn{Y_1 \rightarrow F_q} (C1 matrix). Defaults to \code{weight_scale}.
#' @param threshold Minimum coefficient value needed for an edge to be drawn.
#' @param sig.level Significance level for filtering edges by p-value
#'   (requires inference results). Edges with p-value above this level are omitted.
#' @param rankdir Graphviz rank direction (e.g., \code{"LR"}, \code{"TB"}).
#' @param fill Logical; whether to use filled node shapes.
#' @param ... For backward compatibility: accepts deprecated names
#'   \code{weight_scale_y2f} (use \code{weight_scale_c2}) and
#'   \code{weight_scale_fy1} (use \code{weight_scale_x1}).
#' @param cluster.box Character string controlling the visibility and style
#'   of cluster frames around Y2, factors, and Y1 blocks.
#'   One of \code{"normal"}, \code{"faint"}, \code{"invisible"}, \code{"none"}.
#' @param cluster.labels Optional character vector of length 3 giving custom
#'   labels for the Y2, factor, and Y1 clusters.
#' @param hide.isolated Logical. If \code{TRUE} (default), Y1 and Y2 nodes
#'   that have no edges at or above \code{threshold} are excluded from the graph.
#' @param sig.level Significance level for filtering structural edges
#'   (\eqn{C_1} feedback and \eqn{C_2} exogenous loadings) when
#'   inference results are present.  If \code{result} contains a
#'   \code{coefficients} data frame from \code{\link{nmf.sem.inference}},
#'   only edges with \code{p_value < sig.level} are drawn, with
#'   significance stars (\code{*} \code{**} \code{***}) appended to
#'   the edge label.  The \eqn{X} (factor-to-\eqn{Y_1}) edges are
#'   never starred since the basis is not the inference target.
#'   Set to \code{NULL} to disable significance filtering and fall
#'   back to the \code{threshold} magnitude filter for both \eqn{C_1}
#'   and \eqn{C_2}.  Default is \code{0.1}.
#'
#' @return A character string representing a valid Graphviz DOT script.
#'
#' @examples
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]
#' Y2 <- Y[3:4, ]
#' result <- nmf.sem(Y1, Y2, rank = 2, maxit = 500)
#' dot <- nmf.sem.DOT(result)
#' cat(dot)
#'
#' @seealso \code{\link{nmf.ffb}}, \code{\link{nmf.ffb.inference}},
#'   \code{\link{plot.nmfkc.DOT}}
#' @export
nmf.ffb.DOT <- function(result,
                        weight_scale          = 5,
                        weight_scale_c2       = weight_scale,
                        weight_scale_x1       = weight_scale,
                        weight_scale_feedback = weight_scale,
                        threshold             = 0.01,
                        sig.level             = 0.1,
                        rankdir               = "LR",
                        fill                  = TRUE,
                        cluster.box           = c("normal", "faint", "invisible", "none"),
                        cluster.labels        = NULL,
                        hide.isolated         = TRUE,
                        ...) {

  ## Backward compatibility: accept deprecated names via ...
  extra_args <- base::list(...)
  if (!base::is.null(extra_args$weight_scale_y2f)) weight_scale_c2 <- extra_args$weight_scale_y2f
  if (!base::is.null(extra_args$weight_scale_fy1)) weight_scale_x1 <- extra_args$weight_scale_fy1

  ## ---------------------------------------------------------------
  ## Cluster style selection
  ## ---------------------------------------------------------------
  cluster.box <- match.arg(cluster.box)

  cluster_style <- switch(cluster.box,
                          normal    = "rounded",
                          faint     = "rounded,dashed",
                          invisible = "rounded",
                          none      = "none")

  cluster_color <- switch(cluster.box,
                          normal    = "black",
                          faint     = "gray80",
                          invisible = "none",
                          none      = "none")

  cluster_penwidth <- switch(cluster.box,
                             normal    = 1.0,
                             faint     = 0.7,
                             invisible = 0.0,
                             none      = 0.0)

  ## ---------------------------------------------------------------
  ## Cluster titles
  ## ---------------------------------------------------------------
  default.labels <- c("Exogenous (Y2)", "Latent Factors", "Endogenous (Y1)")
  if (is.null(cluster.labels)) {
    titles <- default.labels
  } else {
    if (!is.character(cluster.labels) || length(cluster.labels) != 3L) {
      stop("cluster.labels must be a character vector of length 3: c(label_Y2, label_F, label_Y1).")
    }
    titles <- cluster.labels
  }

  ## ---------------------------------------------------------------
  ## Extract matrices and labels
  ## ---------------------------------------------------------------
  X  <- result$X
  C1 <- result$C1
  C2 <- result$C2

  if (is.null(X) || is.null(C1) || is.null(C2)) {
    stop("result must contain elements X, C1, and C2.")
  }

  Y1_labels <- rownames(X)
  Y2_labels <- colnames(C2)

  if (is.null(Y1_labels)) Y1_labels <- paste0("Y1_", seq_len(nrow(X)))
  if (is.null(Y2_labels)) Y2_labels <- paste0("Y2_", seq_len(ncol(C2)))

  P1 <- length(Y1_labels)
  P2 <- length(Y2_labels)
  Q  <- ncol(X)

  ## ---------------------------------------------------------------
  ## Identify nodes involved in edges >= threshold
  ## ---------------------------------------------------------------
  if (isTRUE(hide.isolated)) {
    used_y2 <- apply(C2, 2L, function(col) any(col >= threshold, na.rm = TRUE))
    used_y1_from_X  <- apply(X,  1L, function(row) any(row >= threshold, na.rm = TRUE))
    used_y1_from_C1 <- apply(C1, 2L, function(col) any(col >= threshold, na.rm = TRUE))
    used_y1 <- used_y1_from_X | used_y1_from_C1
    idx_y1 <- which(used_y1)
    idx_y2 <- which(used_y2)
  } else {
    idx_y1 <- seq_len(P1)
    idx_y2 <- seq_len(P2)
  }

  ## ---------------------------------------------------------------
  ## Assign internal DOT node IDs
  ## ---------------------------------------------------------------
  Y1_ids <- .nmfkc_dot_sanitize_id(paste0("Y1_", seq_len(P1)))
  Y2_ids <- .nmfkc_dot_sanitize_id(paste0("Y2_", seq_len(P2)))
  F_ids  <- .nmfkc_dot_sanitize_id(paste0("F_",  seq_len(Q)))

  ## Node colors
  COLOR_Y2_NODE     <- "lightcoral"
  COLOR_Y1_NODE     <- "lightblue"
  COLOR_FACTOR_NODE <- "wheat"

  ## ---------------------------------------------------------------
  ## Header
  ## ---------------------------------------------------------------
  dot_script <- .nmfkc_dot_header(
    graph_name = "NMF_SEM_Full_Mechanism",
    rankdir    = rankdir
  )

  ## ---------------------------------------------------------------
  ## Y2 cluster
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    "\n  // Exogenous variables (Y2)\n",
    .nmfkc_dot_cluster_nodes(
      cluster_id        = "Y2",
      title             = titles[1],
      node_ids          = Y2_ids[idx_y2],
      node_labels       = Y2_labels[idx_y2],
      shape             = "box",
      fill              = fill,
      fillcolor         = COLOR_Y2_NODE,
      line_width        = 1.5,
      cluster_style     = cluster_style,
      cluster_color     = cluster_color,
      cluster_penwidth  = cluster_penwidth
    )
  )

  ## ---------------------------------------------------------------
  ## Y1 cluster
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    "\n  // Endogenous variables (Y1)\n",
    .nmfkc_dot_cluster_nodes(
      cluster_id        = "Y1",
      title             = titles[3],
      node_ids          = Y1_ids[idx_y1],
      node_labels       = Y1_labels[idx_y1],
      shape             = "box",
      fill              = fill,
      fillcolor         = COLOR_Y1_NODE,
      line_width        = 1.5,
      cluster_style     = cluster_style,
      cluster_color     = cluster_color,
      cluster_penwidth  = cluster_penwidth
    )
  )

  ## ---------------------------------------------------------------
  ## Latent factor cluster
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    "\n  // Latent Factors (F)\n",
    .nmfkc_dot_cluster_nodes(
      cluster_id        = "F",
      title             = titles[2],
      node_ids          = F_ids,
      node_labels       = paste0("Factor ", seq_len(Q)),
      shape             = "ellipse",
      fill              = fill,
      fillcolor         = COLOR_FACTOR_NODE,
      line_width        = 1.0,
      cluster_style     = cluster_style,
      cluster_color     = cluster_color,
      cluster_penwidth  = cluster_penwidth
    )
  )

  ## ---------------------------------------------------------------
  ## Edge defaults
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  edge [fontname="Arial", fontsize=8, arrowhead=open];\n'
  )

  pw     <- .nmfkc_dot_penwidth
  digits <- .nmfkc_dot_digits_from_threshold(threshold)
  fmtc   <- function(x) .nmfkc_dot_format_coef(x, digits)

  ## ---------------------------------------------------------------
  ## Significance stars for C1 (feedback) and C2 (exogenous) edges,
  ## both from nmf.sem.inference().  X (F -> Y1) edges are NOT
  ## starred even when inference results are present, since the
  ## basis is not the inference target.
  ##
  ## Filters by the optional `Type` column ("C1" / "C2") if present
  ## (newer inference output); falls back to row-name matching for
  ## back-compatibility with pre-v0.6.8 inference output that only
  ## carried C2.
  ## ---------------------------------------------------------------
  pval_to_stars <- function(p) {
    if (!is.finite(p)) ""
    else if (p < 0.001) "***"
    else if (p < 0.01)  "**"
    else if (p < 0.05)  "*"
    else                ""
  }

  C1_stars <- NULL; C1_show <- NULL
  C2_stars <- NULL; C2_show <- NULL
  if (!is.null(result$coefficients)) {
    cf        <- result$coefficients
    fac_names <- rownames(C2)
    y1_names  <- colnames(C1)
    y2_names  <- colnames(C2)

    has_type <- !is.null(cf$Type)
    cf_C1 <- if (has_type) cf[cf$Type == "C1", , drop = FALSE] else NULL
    cf_C2 <- if (has_type) cf[cf$Type == "C2", , drop = FALSE] else cf

    ## ---- C2 stars (Y2 -> F) ----
    C2_stars <- matrix("", nrow = Q, ncol = P2)
    C2_pval  <- matrix(NA_real_, nrow = Q, ncol = P2)
    for (k in seq_len(nrow(cf_C2))) {
      q  <- match(cf_C2$Basis[k], fac_names)
      p2 <- match(cf_C2$Covariate[k], y2_names)
      if (!is.na(q) && !is.na(p2) && !is.na(cf_C2$p_value[k])) {
        C2_pval[q, p2]  <- cf_C2$p_value[k]
        C2_stars[q, p2] <- pval_to_stars(cf_C2$p_value[k])
      }
    }
    if (!is.null(sig.level)) {
      C2_show <- !is.na(C2_pval) & C2_pval < sig.level
    }

    ## ---- C1 stars (Y1 -> F, feedback) ----
    if (!is.null(cf_C1) && nrow(cf_C1) > 0L) {
      C1_stars <- matrix("", nrow = Q, ncol = P1)
      C1_pval  <- matrix(NA_real_, nrow = Q, ncol = P1)
      for (k in seq_len(nrow(cf_C1))) {
        q  <- match(cf_C1$Basis[k], fac_names)
        p1 <- match(cf_C1$Covariate[k], y1_names)
        if (!is.na(q) && !is.na(p1) && !is.na(cf_C1$p_value[k])) {
          C1_pval[q, p1]  <- cf_C1$p_value[k]
          C1_stars[q, p1] <- pval_to_stars(cf_C1$p_value[k])
        }
      }
      if (!is.null(sig.level)) {
        C1_show <- !is.na(C1_pval) & C1_pval < sig.level
      }
    }
  }

  ## ---------------------------------------------------------------
  ## 1. Y2 → F edges (C2)
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  // 1. External Driving (Y2 -> Factor) [C2]\n',
    '  edge [color=black, fontcolor=black, style=solid];\n'
  )

  max_C2 <- suppressWarnings(max(C2, na.rm = TRUE))
  if (is.finite(max_C2) && max_C2 > 0) {
    for (q in seq_len(Q)) {
      for (p2 in idx_y2) {
        weight <- C2[q, p2]
        show <- if (!is.null(C2_show)) C2_show[q, p2]
                else is.finite(weight) && weight >= threshold
        if (show) {
          pen <- pw(weight, max_C2, weight_scale_c2)
          lab <- fmtc(weight)
          if (!is.null(C2_stars)) lab <- paste0(lab, C2_stars[q, p2])
          path <- sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                          Y2_ids[p2], F_ids[q], lab, pen)
          dot_script <- paste0(dot_script, path)
        }
      }
    }
  }

  ## ---------------------------------------------------------------
  ## 2. F → Y1 edges (X)
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  // 2. Generation (Factor -> Y1) [X]\n',
    '  edge [color="gray0", fontcolor="gray0", style=solid];\n'
  )

  max_X <- suppressWarnings(max(X, na.rm = TRUE))
  if (is.finite(max_X) && max_X > 0) {
    for (q in seq_len(Q)) {
      for (p1 in idx_y1) {
        weight <- X[p1, q]
        if (is.finite(weight) && weight >= threshold) {
          pen <- pw(weight, max_X, weight_scale_x1)
          lab <- fmtc(weight)
          path <- sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                          F_ids[q], Y1_ids[p1], lab, pen)
          dot_script <- paste0(dot_script, path)
        }
      }
    }
  }

  ## ---------------------------------------------------------------
  ## 3. Y1 → F edges (C1 feedback)
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  // 3. Internal Feedback (Y1 -> Factor) [C1]\n',
    '  edge [style=dashed, color="gray0", fontcolor="gray0"];\n'
  )

  max_C1 <- suppressWarnings(max(C1, na.rm = TRUE))
  if (is.finite(max_C1) && max_C1 > 0) {
    for (q in seq_len(Q)) {
      for (p1 in idx_y1) {
        weight <- C1[q, p1]
        ## When inference results provide C1 p-values, prefer the
        ## significance-based filter (parallels the C2 branch above);
        ## otherwise fall back to the magnitude threshold.
        show <- if (!is.null(C1_show)) C1_show[q, p1]
                else is.finite(weight) && weight >= threshold
        if (show) {
          pen <- pw(weight, max_C1, weight_scale_feedback)
          lab <- fmtc(weight)
          if (!is.null(C1_stars)) lab <- paste0(lab, C1_stars[q, p1])
          path <- sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                          Y1_ids[p1], F_ids[q], lab, pen)
          dot_script <- paste0(dot_script, path)
        }
      }
    }
  }

  result <- paste0(dot_script, "}\n")
  class(result) <- c("nmf.sem.DOT", "nmfkc.DOT")
  result
}





############################################################
## 2. nmfkc.DOT  (Static NMF / NMF-with-covariates visualization)
############################################################

#' Generate Graphviz DOT Scripts for NMF or NMF-with-Covariates Models
#'
#' @description
#' Produces a Graphviz DOT script visualizing the structure of an NMF model
#' (\eqn{Y \approx X C A}) or its simplified forms.
#'
#' Supported visualization types:
#' \itemize{
#'   \item \code{"YX"} — Standard NMF view: latent factors \eqn{X} map to observations \eqn{Y}.
#'   \item \code{"YA"} — Direct regression view: covariates \eqn{A} map directly to \eqn{Y}
#'         using the combined coefficient matrix \eqn{X C}.
#'   \item \code{"YXA"} — Full tri-factorization: \eqn{A \rightarrow C \rightarrow X \rightarrow Y}.
#' }
#'
#' Edge widths are scaled by coefficient magnitude, and nodes with no edges
#' above the threshold are omitted from the visualization.
#'
#' @param result The return value from \code{nmfkc}, containing matrices
#'   \code{X}, \code{B}, and optionally \code{C}.
#' @param type Character string specifying the visualization style:
#'   one of \code{"YX"}, \code{"YA"}, \code{"YXA"}.
#' @param threshold Minimum coefficient magnitude to display an edge. When
#'   \code{C.signed = TRUE} (signed \eqn{\Theta}) this is an absolute-value cut,
#'   i.e. an edge is drawn when \eqn{|coef| \ge} \code{threshold}.
#' @param C.signed Logical or \code{NULL}. Whether the coefficient matrix
#'   \eqn{\Theta} (\code{= C}) may be signed (real-valued). When \code{TRUE},
#'   the \code{threshold} is applied to \eqn{|coef|}, edge widths are scaled by
#'   \eqn{|coef|}, and negative edges are drawn as black \strong{dashed} lines
#'   (positive edges remain solid) to distinguish them; the numeric edge labels
#'   keep their sign. When \code{FALSE}, the historical non-negative behaviour is used
#'   (negative coefficients fall below the threshold and are hidden). When
#'   \code{NULL} (default), the mode is auto-detected from the fit
#'   (\code{result$C.signed}) or from the presence of negative entries in
#'   \eqn{C} / \eqn{XC}. The basis \eqn{X} is always non-negative, so
#'   \eqn{X \rightarrow Y} edges are unaffected.
#' @param sig.level Significance level for filtering C edges when inference
#'   results are available (i.e., \code{x} is of class \code{"nmfkc.inference"}).
#'   Only edges with p-value below \code{sig.level} are shown, decorated with
#'   significance stars (\code{*}, \code{**}, \code{***}). Set to \code{NULL}
#'   to disable filtering and show all edges above \code{threshold}. Default is 0.1.
#' @param rankdir Graphviz rank direction (e.g., \code{"LR"}, \code{"TB"}).
#' @param fill Logical; whether nodes should be drawn with filled shapes.
#' @param weight_scale Base scaling factor for edge widths.
#' @param weight_scale_ax Scaling factor for edges \eqn{A \rightarrow X} (type \code{"YXA"}).
#' @param weight_scale_xy Scaling factor for edges \eqn{X \rightarrow Y}.
#' @param weight_scale_ay Scaling factor for edges \eqn{A \rightarrow Y} (type \code{"YA"}).
#' @param Y.label Optional character vector for labels of Y nodes.
#' @param X.label Optional character vector for labels of X (latent factor) nodes.
#' @param A.label Optional character vector for labels of A (covariate) nodes.
#' @param Y.title Cluster title for Y nodes.
#' @param X.title Cluster title for X nodes.
#' @param A.title Cluster title for A nodes.
#' @param hide.isolated Logical. If \code{TRUE} (default), Y and A nodes that have no
#'   edges at or above \code{threshold} are excluded from the graph.
#'
#' @return A character string representing a Graphviz DOT script.
#'
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfae.DOT}}, \code{\link{nmf.sem.DOT}},
#'   \code{\link{nmfkc.ar.DOT}}, \code{\link{plot.nmfkc.DOT}}
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, rank = 1)
#' dot <- nmfkc.DOT(result)
#' cat(dot)
#'
#' @export
nmfkc.DOT <- function(
    result,
    type = c("YX","YA","YXA"),
    threshold = 0.01,
    C.signed = NULL,
    sig.level = 0.1,
    rankdir   = "LR",
    fill      = TRUE,
    weight_scale    = 5,
    weight_scale_ax = weight_scale,
    weight_scale_xy = weight_scale,
    weight_scale_ay = weight_scale,
    Y.label = NULL, X.label = NULL, A.label = NULL,
    Y.title = "Observation (Y)",
    X.title = "Basis (X)",
    A.title = "Covariates (A)",
    hide.isolated = TRUE
) {

  type <- match.arg(type)

  ## ---------------------------------------------------------
  ## Required matrices
  ## ---------------------------------------------------------
  X <- result$X
  B <- result$B
  if (is.null(X) || is.null(B)) {
    stop("result must contain X and B.")
  }

  ## If C exists and is a proper NMF-with-covariates factor:
  hasA <- !is.null(result$C) && ncol(result$C) != ncol(B)

  P <- nrow(X)
  Q <- ncol(X)

  ## ---------------------------------------------------------
  ## Labels
  ## ---------------------------------------------------------
  Y_labels <- if (is.null(Y.label)) rownames(X) else Y.label
  X_labels <- if (is.null(X.label)) colnames(X) else X.label

  if (is.null(Y_labels)) Y_labels <- paste0("Y", seq_len(P))
  if (is.null(X_labels)) X_labels <- paste0("Factor", seq_len(Q))

  Y_ids <- .nmfkc_dot_sanitize_id(paste0("Y_", seq_len(P)))
  X_ids <- .nmfkc_dot_sanitize_id(paste0("X_", seq_len(Q)))

  ## ---------------------------------------------------------
  ## Covariates and tri-factorization
  ## ---------------------------------------------------------
  if (hasA) {
    C <- as.matrix(result$C)          # Q x R
    A_cols <- ncol(C)
    A_labels <- if (is.null(A.label)) colnames(C) else A.label
    if (is.null(A_labels)) A_labels <- paste0("A", seq_len(A_cols))
    A_ids <- .nmfkc_dot_sanitize_id(paste0("A_", seq_len(A_cols)))

    ## Combined mapping for type = "YA"
    XC_mat <- X %*% C   # P x R

  } else if (type == "YX") {
    ## No A block needed
    C <- NULL
    A_cols <- 0L
    A_labels <- NULL
    A_ids <- NULL
    XC_mat <- NULL

  } else {
    stop("The model structure (matrix C) is incompatible with type 'YXA' or 'YA'.")
  }

  ## ---------------------------------------------------------
  ## Sign mode of C (= Theta). When NULL, auto-detect from the fit
  ## (result$C.signed) or from negative entries in C / XC. In signed mode the
  ## threshold is an absolute-value cut (|coef| >= threshold), edge widths use
  ## |coef|, and negative edges are coloured to distinguish them. The basis X
  ## is always non-negative, so X -> Y edges are unaffected.
  ## ---------------------------------------------------------
  if (is.null(C.signed)) {
    C.signed <- isTRUE(result$C.signed) ||
      (!is.null(C)      && any(C      < 0, na.rm = TRUE)) ||
      (!is.null(XC_mat) && any(XC_mat < 0, na.rm = TRUE))
  }
  C.signed <- isTRUE(C.signed)
  ## magnitude used for thresholding and penwidth scaling (|.| in signed mode)
  mag <- if (C.signed) function(v) abs(v) else function(v) v
  ## edge line style by sign (negative -> dashed) in signed mode; colour stays black
  estyle <- function(v) if (C.signed && is.finite(v) && v < 0) "dashed" else "solid"

  ## ---------------------------------------------------------
  ## Filter isolated nodes (hide.isolated)
  ## ---------------------------------------------------------
  idx_Y <- seq_len(P)
  idx_A <- if (hasA) seq_len(A_cols) else integer(0)
  if (isTRUE(hide.isolated)) {
    if (type == "YX" || type == "YXA") {
      # Y is connected via X: X[i, ] >= threshold
      used_Y <- apply(X, 1L, function(row) any(row >= threshold, na.rm = TRUE))
      idx_Y <- which(used_Y)
    } else if (type == "YA" && !is.null(XC_mat)) {
      # Y is connected via XC: |XC_mat[i, ]| >= threshold
      used_Y <- apply(XC_mat, 1L, function(row) any(mag(row) >= threshold, na.rm = TRUE))
      idx_Y <- which(used_Y)
    }
    if (hasA && type %in% c("YXA", "YA")) {
      if (type == "YXA" && !is.null(C)) {
        # A is connected via C: |C[, k]| >= threshold
        used_A <- apply(C, 2L, function(col) any(mag(col) >= threshold, na.rm = TRUE))
        idx_A <- which(used_A)
      } else if (type == "YA" && !is.null(XC_mat)) {
        # A is connected via XC: |XC_mat[, k]| >= threshold
        used_A <- apply(XC_mat, 2L, function(col) any(mag(col) >= threshold, na.rm = TRUE))
        idx_A <- which(used_A)
      }
    }
  }

  ## ---------------------------------------------------------
  ## DOT header
  ## ---------------------------------------------------------
  scr <- .nmfkc_dot_header(graph_name = "NMF", rankdir = rankdir)

  ## ---------------------------------------------------------
  ## Y node cluster
  ## ---------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  // Output variables (Y)\n',
    .nmfkc_dot_cluster_nodes(
      cluster_id  = "Y",
      title       = Y.title,
      node_ids    = Y_ids[idx_Y],
      node_labels = Y_labels[idx_Y],
      shape       = "box",
      fill        = fill,
      fillcolor   = "lightblue",
      line_width  = 1.5
    )
  )

  ## ---------------------------------------------------------
  ## X node cluster
  ## (Hidden when type = "YA")
  ## ---------------------------------------------------------
  if (type != "YA") {
    scr <- paste0(
      scr,
      '\n  // Latent factors (X)\n',
      .nmfkc_dot_cluster_nodes(
        cluster_id  = "X",
        title       = X.title,
        node_ids    = X_ids,
        node_labels = X_labels,
        shape       = "ellipse",
        fill        = fill,
        fillcolor   = "wheat",
        line_width  = 1.0
      )
    )
  }

  ## ---------------------------------------------------------
  ## A node cluster
  ## (Only for "YA" and "YXA")
  ## ---------------------------------------------------------
  if (type != "YX" && hasA) {
    scr <- paste0(
      scr,
      '\n  // Covariates (A)\n',
      .nmfkc_dot_cluster_nodes(
        cluster_id  = "A",
        title       = A.title,
        node_ids    = A_ids[idx_A],
        node_labels = A_labels[idx_A],
        shape       = "box",
        fill        = fill,
        fillcolor   = "lightcoral",
        line_width  = 1.5
      )
    )
  }

  ## ---------------------------------------------------------
  ## Significance stars for C edges (from nmfkc.inference)
  ## ---------------------------------------------------------
  C_stars <- NULL
  C_show  <- NULL
  if (!is.null(result$coefficients) && hasA) {
    C_stars <- matrix("", nrow = Q, ncol = A_cols)
    C_pval  <- matrix(NA_real_, nrow = Q, ncol = A_cols)
    cf <- result$coefficients
    basis_names <- rownames(C)
    cov_names   <- colnames(C)
    for (k in seq_len(nrow(cf))) {
      q <- match(cf$Basis[k], basis_names)
      r <- match(cf$Covariate[k], cov_names)
      if (!is.na(q) && !is.na(r) && !is.na(cf$p_value[k])) {
        p <- cf$p_value[k]
        C_pval[q, r] <- p
        if (p < 0.001)      C_stars[q, r] <- "***"
        else if (p < 0.01)  C_stars[q, r] <- "**"
        else if (p < 0.05)  C_stars[q, r] <- "*"
      }
    }
    if (!is.null(sig.level)) {
      C_show <- !is.na(C_pval) & C_pval < sig.level
    }
  }

  ## ---------------------------------------------------------
  ## Edge defaults
  ## ---------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  edge [fontname="Arial", fontsize=8, arrowhead=open];\n'
  )

  pw     <- .nmfkc_dot_penwidth
  digits <- .nmfkc_dot_digits_from_threshold(threshold)
  fmtc   <- function(x) .nmfkc_dot_format_coef(x, digits)


  ## =========================================================
  ## Case 1: Full tri-factorization (A → X → Y)
  ## =========================================================
  if (type == "YXA") {

    ## ---- A → X (C) edges ----
    scr <- paste0(
      scr,
      '\n  // A -> X edges (C)\n',
      '  edge [color=black, fontcolor=black, style=solid];\n'
    )

    max_C <- suppressWarnings(max(mag(C), na.rm = TRUE))
    if (is.finite(max_C) && max_C > 0) {
      for (q in seq_len(Q)) {
        for (k in idx_A) {
          val <- C[q, k]
          show <- if (!is.null(C_show)) C_show[q, k]
                  else is.finite(val) && mag(val) >= threshold
          if (show) {
            pen <- pw(mag(val), max_C, weight_scale_ax)
            lab <- fmtc(val)
            if (!is.null(C_stars)) lab <- paste0(lab, C_stars[q, k])
            sty <- estyle(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f, style=%s];\n',
                      A_ids[k], X_ids[q], lab, pen, sty)
            )
          }
        }
      }
    }

    ## ---- X → Y edges ----
    scr <- paste0(
      scr,
      '\n  // X -> Y edges (X)\n',
      '  edge [color="gray0", fontcolor="gray0", style=solid];\n'
    )

    max_X <- suppressWarnings(max(X, na.rm = TRUE))
    if (is.finite(max_X) && max_X > 0) {
      for (i in idx_Y) {
        for (j in seq_len(Q)) {
          val <- X[i, j]
          if (is.finite(val) && val >= threshold) {
            pen <- pw(val, max_X, weight_scale_xy)
            lab <- fmtc(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                      X_ids[j], Y_ids[i], lab, pen)
            )
          }
        }
      }
    }


    ## =========================================================
    ## Case 2: Direct regression view (A → Y)
    ## =========================================================
  } else if (type == "YA") {

    scr <- paste0(
      scr,
      '\n  // A -> Y edges (X %*% C)\n',
      '  edge [color=black, fontcolor=black, style=solid];\n'
    )

    max_XC <- suppressWarnings(max(mag(XC_mat), na.rm = TRUE))
    if (is.finite(max_XC) && max_XC > 0) {
      for (i in idx_Y) {
        for (k in idx_A) {
          val <- XC_mat[i, k]
          if (is.finite(val) && mag(val) >= threshold) {
            pen <- pw(mag(val), max_XC, weight_scale_ay)
            lab <- fmtc(val)
            sty <- estyle(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f, style=%s];\n',
                      A_ids[k], Y_ids[i], lab, pen, sty)
            )
          }
        }
      }
    }


    ## =========================================================
    ## Case 3: Standard NMF (X → Y)
    ## =========================================================
  } else if (type == "YX") {

    scr <- paste0(
      scr,
      '\n  // X -> Y edges (X)\n',
      '  edge [color="gray0", fontcolor="gray0", style=solid];\n'
    )

    max_X <- suppressWarnings(max(X, na.rm = TRUE))
    if (is.finite(max_X) && max_X > 0) {
      for (i in idx_Y) {
        for (j in seq_len(Q)) {
          val <- X[i, j]
          if (is.finite(val) && val >= threshold) {
            pen <- pw(val, max_X, weight_scale_xy)
            lab <- fmtc(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                      X_ids[j], Y_ids[i], lab, pen)
            )
          }
        }
      }
    }
  }

  result <- paste0(scr, "}\n")
  class(result) <- "nmfkc.DOT"
  result
}


#' @title Plot method for nmfkc.DOT objects
#' @description
#' Renders a DOT graph string using \code{DiagrammeR::grViz}.
#' If the \pkg{DiagrammeR} package is not installed, prints the DOT source
#' to the console instead.
#'
#' This method handles all DOT objects produced by the nmfkc package:
#' \code{\link{nmfkc.DOT}}, \code{\link{nmfae.DOT}}, \code{\link{nmf.sem.DOT}},
#' and \code{\link{nmfkc.ar.DOT}}.
#'
#' @param x An object of class \code{"nmfkc.DOT"} (or a subclass thereof).
#' @param ... Not used.
#' @return Called for its side effect (rendering). Returns \code{x} invisibly.
#' @seealso \code{\link{nmfkc.DOT}}, \code{\link{nmfae.DOT}},
#'   \code{\link{nmf.sem.DOT}}, \code{\link{nmfkc.ar.DOT}}
#' @export
plot.nmfkc.DOT <- function(x, ...) {
  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    print(DiagrammeR::grViz(as.character(x)))
  } else {
    message("DiagrammeR package not installed. Printing DOT source:")
    cat(as.character(x))
  }
  invisible(x)
}

# nmf.sem.R — NMF-SEM related functions
# nmf.sem, nmf.sem.cv, nmf.sem.split, nmf.sem.DOT

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#' @title NMF-SEM Main Estimation Algorithm
#'
#' @description
#' Fits the NMF-SEM model
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
#' @param X.init Optional non-negative initialization for the basis matrix
#'   \code{X} (\eqn{P_1 \times Q}). If supplied, it is projected to be
#'   non-negative and column-normalized.
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
#'   Default: \code{20000}.
#' @param seed Random seed used to initialize \code{X}, \code{C1}, and \code{C2}.
#'   Default: \code{123}.
#' @param ... Additional arguments. Currently used to pass a hidden rank
#'   \code{Q} (e.g., via \code{Q = 3}) if \code{rank} is \code{NULL}.
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
#'     (flattened) of \eqn{Y_1}.}
#'   \item{MAE}{Mean absolute error between \eqn{Y_1} and its equilibrium
#'     prediction \eqn{\hat Y_1 = M_{\mathrm{model}} Y_2}.}
#'   \item{objfunc}{Vector of reconstruction losses per iteration.}
#'   \item{objfunc.full}{Vector of penalized objective values per iteration.}
#'   \item{iter}{Number of iterations actually performed.}
#'
#' @examples
#' # Simple NMF-SEM with iris data (non-negative)
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]  # Sepal
#' Y2 <- Y[3:4, ]  # Petal
#' result <- nmf.sem(Y1, Y2, rank = 2, maxit = 500)
#' result$MAE
#'
#' @export
nmf.sem <- function(
    Y1, Y2,
    rank = NULL,
    X.init = NULL,
    X.L2.ortho = 100.0,
    C1.L1 = 1.0,
    C2.L1 = 0.1,
    epsilon = 1e-6,
    maxit = 20000,
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

  # ---------------------------- init X,C1,C2 --------------------------
  if (is.null(X.init)) {
    if (Q <= min(P1, N)) {
      X <- .nndsvdar(Y1, Q)
    } else {
      X <- matrix(stats::runif(P1 * Q), nrow = P1, ncol = Q)
    }
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

  list(
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
    MAE                 = MAE,
    objfunc             = objfunc[1:it],
    objfunc.full        = objfunc.full[1:it],
    iter                = it
  )
}



#' @title Cross-Validation for NMF-SEM
#' @description
#' Performs K-fold cross-validation to evaluate the equilibrium mapping of
#' the NMF-SEM model.
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
#' @param X.init Optional initialization for \code{X} (as in \code{nmf.sem}).
#' @param X.L2.ortho L2 orthogonality penalty for \code{X}.
#' @param C1.L1 L1 sparsity penalty for \code{C1} (\eqn{\Theta_1}).
#' @param C2.L1 L1 sparsity penalty for \code{C2} (\eqn{\Theta_2}).
#' @param epsilon Convergence threshold for \code{nmf.sem}.
#' @param maxit Maximum number of iterations for \code{nmf.sem}.
#' @param seed Master random seed for CV splitting and fold-specific calls to \code{nmf.sem}.
#'   If \code{NULL}, RNG is not controlled within folds.
#' @param div Number of CV folds. (Default: \code{5})
#' @param shuffle Logical; if \code{TRUE}, samples are randomly permuted
#'   before assigning to folds. (Default: \code{TRUE})
#' @param ... Additional arguments passed to \code{nmf.sem} (except for
#'   \code{rank}, \code{seed}, \code{div}, \code{shuffle}, which are handled here).
#'
#' @return A numeric scalar: mean MAE across CV folds.
#'
#' @examples
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]
#' Y2 <- Y[3:4, ]
#' mae <- nmf.sem.cv(Y1, Y2, rank = 2, maxit = 500, div = 3)
#' mae
#'
#' @export
nmf.sem.cv <- function(
    Y1, Y2,
    rank = NULL,
    X.init = NULL,
    X.L2.ortho = 100.0,
    C1.L1 = 1.0,        # L1 sparsity for C1 (Theta1)
    C2.L1 = 0.1,        # L1 sparsity for C2 (Theta2)
    epsilon = 1e-6,     # Convergence tolerance passed to nmf.sem
    maxit = 20000,
    seed = NULL,        # Master seed for CV (partition + fold seeds)
    div = 5,            # Number of CV folds
    shuffle = TRUE,     # Shuffle samples before assigning folds
    ...
){
  # ------------------------------------------------------------------
  # 1. Basic input checks
  #
  # NMF-SEM requires non-negative matrices. We also require that Y1 and Y2
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
    res_j <- suppressMessages(do.call("nmf.sem", nmf.sem.args))

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
  # (e.g., rank, X.L2.ortho, C1.L1, C2.L1) when tuning NMF-SEM.
  # ------------------------------------------------------------------
  objfunc <- mean(objfunc.block)
  return(objfunc)
}



#' @title Heuristic Variable Splitting for NMF-SEM
#'
#' @description
#' Infers a heuristic partition of observed variables into exogenous (\eqn{Y_2})
#' and endogenous (\eqn{Y_1}) blocks for use in NMF-SEM.
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
#'   resulting variable split. (Default: \code{TRUE})
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
#' sp$Y1.names
#' sp$Y2.names
#'
#' @export
nmf.sem.split <- function(x, n.exogenous = NULL, threshold = 0.1,
                          auto.flipped = TRUE, verbose = TRUE) {

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
  # Variables are centered and scaled (mean 0, sd 1). NMF-SEM requires
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
  # In positive SEM (and NMF-SEM), variables should ideally have
  # consistent sign orientation. To enforce this heuristic, variables
  # negatively correlated with the first principal component are flipped.
  #
  # This stabilizes the causal-ordering heuristic by avoiding mixtures
  # of arbitrary sign conventions in the raw data.
  # --------------------------------------------------------------------
  is.flipped <- rep(FALSE, P)
  names(is.flipped) <- col_names

  if (auto.flipped) {
    if (verbose) cat("Step 0: Checking correlations with PC1 (on standardized data)...\n")

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
        cat(sprintf("   -> Detected %d flipped variables: %s\n",
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
  if (verbose) cat("Step 1: Inferring Causal Ordering...\n")

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
    if (verbose) cat("Step 2: Detecting optimal cut-off for exogenous variables...\n")
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
          cat(sprintf("   [%d] %s : Max std.coef=%.3f -> Endogenous (Stop)\n",
                      k, col_names[curr_idx], max_influence))
        break
      } else {
        cutoff <- k
        if (verbose)
          cat(sprintf("   [%d] %s : Max std.coef=%.3f -> Exogenous (Continue)\n",
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
    cat("\n--- Auto Split Result ---\n")
    cat(sprintf("Exogenous (Y2, n=%d): %s\n",
                n.exogenous, paste(exogenous.variables, collapse=", ")))
    cat(sprintf("Endogenous (Y1, n=%d): %s ...\n",
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

#' Determine the decimal digits based on a threshold
#'
#' This helper computes the number of decimal places that should be
#' used for formatting coefficient labels, based on the magnitude
#' of the threshold used for edge visualization.
#'
#' @param threshold Numeric scalar (>0).
#' @return Integer specifying the number of decimal places.
#' @keywords internal
#' @noRd

############################################################
## 1. nmf.sem.DOT  (for NMF-SEM visualization)
############################################################

#' Generate a Graphviz DOT Diagram for an NMF-SEM Model
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
#' @param weight_scale_y2f Optional override for scaling edges
#'   \eqn{Y_2 \rightarrow F_q}. Defaults to \code{weight_scale}.
#' @param weight_scale_fy1 Optional override for scaling edges
#'   \eqn{F_q \rightarrow Y_1}. Defaults to \code{weight_scale}.
#' @param weight_scale_feedback Optional override for scaling feedback edges
#'   \eqn{Y_1 \rightarrow F_q}. Defaults to \code{weight_scale}.
#' @param threshold Minimum coefficient value needed for an edge to be drawn.
#' @param rankdir Graphviz rank direction (e.g., \code{"LR"}, \code{"TB"}).
#' @param fill Logical; whether to use filled node shapes.
#' @param cluster.box Character string controlling the visibility and style
#'   of cluster frames around Y2, factors, and Y1 blocks.
#'   One of \code{"normal"}, \code{"faint"}, \code{"invisible"}, \code{"none"}.
#' @param cluster.labels Optional character vector of length 3 giving custom
#'   labels for the Y2, factor, and Y1 clusters.
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
#' @export
nmf.sem.DOT <- function(result,
                        weight_scale          = 5,
                        weight_scale_y2f      = weight_scale,
                        weight_scale_fy1      = weight_scale,
                        weight_scale_feedback = weight_scale,
                        threshold             = 0.01,
                        rankdir               = "LR",
                        fill                  = TRUE,
                        cluster.box           = c("normal", "faint", "invisible", "none"),
                        cluster.labels        = NULL) {

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
  used_y2 <- apply(C2, 2L, function(col) any(col >= threshold, na.rm = TRUE))
  used_y1_from_X  <- apply(X,  1L, function(row) any(row >= threshold, na.rm = TRUE))
  used_y1_from_C1 <- apply(C1, 2L, function(col) any(col >= threshold, na.rm = TRUE))
  used_y1 <- used_y1_from_X | used_y1_from_C1

  idx_y1 <- which(used_y1)
  idx_y2 <- which(used_y2)

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
        if (is.finite(weight) && weight >= threshold) {
          pen <- pw(weight, max_C2, weight_scale_y2f)
          lab <- fmtc(weight)
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
          pen <- pw(weight, max_X, weight_scale_fy1)
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
        if (is.finite(weight) && weight >= threshold) {
          pen <- pw(weight, max_C1, weight_scale_feedback)
          lab <- fmtc(weight)
          path <- sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                          Y1_ids[p1], F_ids[q], lab, pen)
          dot_script <- paste0(dot_script, path)
        }
      }
    }
  }

  paste0(dot_script, "}\n")
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
#' @param x The return value from \code{nmfkc}, containing matrices
#'   \code{X}, \code{B}, and optionally \code{C}.
#' @param type Character string specifying the visualization style:
#'   one of \code{"YX"}, \code{"YA"}, \code{"YXA"}.
#' @param threshold Minimum coefficient magnitude to display an edge.
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
#'
#' @return A character string representing a Graphviz DOT script.
#'
#' @seealso \code{nmfkc}
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, Q = 1)
#' dot <- nmfkc.DOT(result)
#' cat(dot)
#'
#' @export
nmfkc.DOT <- function(
    x,
    type = c("YX","YA","YXA"),
    threshold = 0.01,
    rankdir   = "LR",
    fill      = TRUE,
    weight_scale    = 5,
    weight_scale_ax = weight_scale,
    weight_scale_xy = weight_scale,
    weight_scale_ay = weight_scale,
    Y.label = NULL, X.label = NULL, A.label = NULL,
    Y.title = "Observation (Y)",
    X.title = "Basis (X)",
    A.title = "Covariates (A)"
) {

  type <- match.arg(type)

  ## ---------------------------------------------------------
  ## Required matrices
  ## ---------------------------------------------------------
  X <- x$X
  B <- x$B
  if (is.null(X) || is.null(B)) {
    stop("x must contain X and B.")
  }

  ## If C exists and is a proper NMF-with-covariates factor:
  hasA <- !is.null(x$C) && ncol(x$C) != ncol(B)

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
    C <- as.matrix(x$C)          # Q x R
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
      node_ids    = Y_ids,
      node_labels = Y_labels,
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
        node_ids    = A_ids,
        node_labels = A_labels,
        shape       = "box",
        fill        = fill,
        fillcolor   = "lightcoral",
        line_width  = 1.5
      )
    )
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

    max_C <- suppressWarnings(max(C, na.rm = TRUE))
    if (is.finite(max_C) && max_C > 0) {
      for (q in seq_len(Q)) {
        for (k in seq_len(A_cols)) {
          val <- C[q, k]
          if (is.finite(val) && val >= threshold) {
            pen <- pw(val, max_C, weight_scale_ax)
            lab <- fmtc(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                      A_ids[k], X_ids[q], lab, pen)
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
      for (i in seq_len(P)) {
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

    max_XC <- suppressWarnings(max(XC_mat, na.rm = TRUE))
    if (is.finite(max_XC) && max_XC > 0) {
      for (i in seq_len(P)) {
        for (k in seq_len(A_cols)) {
          val <- XC_mat[i, k]
          if (is.finite(val) && val >= threshold) {
            pen <- pw(val, max_XC, weight_scale_ay)
            lab <- fmtc(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                      A_ids[k], Y_ids[i], lab, pen)
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
      for (i in seq_len(P)) {
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

  paste0(scr, "}\n")
}

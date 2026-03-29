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
#' @param ... Additional arguments (reserved for future use).
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
#' @seealso \code{\link{nmf.sem.inference}}, \code{\link{nmf.sem.cv}},
#'   \code{\link{nmf.sem.split}}, \code{\link{nmf.sem.DOT}},
#'   \code{\link{summary.nmf.sem}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to structural equation modeling for blind input-output analysis.
#'   arXiv:2512.18250. \url{https://arxiv.org/abs/2512.18250}
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
    MAE                 = MAE,
    objfunc             = objfunc[1:it],
    objfunc.full        = objfunc.full[1:it],
    iter                = it
  )
  class(out) <- "nmf.sem"
  out
}


#' @title Statistical inference for the exogenous parameter matrix C2
#' @description
#' \code{nmf.sem.inference} performs statistical inference on the exogenous
#' parameter matrix \eqn{C_2} from a fitted \code{nmf.sem} model, conditional
#' on the estimated basis matrix \eqn{\hat{X}} and the endogenous parameter
#' matrix \eqn{\hat{C}_1}.
#'
#' Under the working model \eqn{R = Y_1 - X C_1 Y_1 \approx X C_2 Y_2 + \varepsilon},
#' inference on \eqn{C_2} is conducted via sandwich covariance estimation and
#' one-step wild bootstrap with non-negative projection.
#'
#' @param object A list returned by \code{\link{nmf.sem}}, containing at least
#'   \code{X}, \code{C1}, and \code{C2}.
#' @param Y1 Endogenous variable matrix (P1 x N). Must match the data used in
#'   \code{nmf.sem()}.
#' @param Y2 Exogenous variable matrix (P2 x N). Must match the data used in
#'   \code{nmf.sem()}.
#' @param wild.bootstrap Logical. If \code{TRUE} (default), performs wild bootstrap
#'   for confidence intervals and bootstrap standard errors.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{wild.B}}{Number of bootstrap replicates. Default is 1000.}
#'     \item{\code{wild.seed}}{Seed for bootstrap. Default is 42.}
#'     \item{\code{wild.level}}{Confidence level for bootstrap CI. Default is 0.95.}
#'     \item{\code{sandwich}}{Logical. Use sandwich covariance. Default is \code{TRUE}.}
#'     \item{\code{C.p.side}}{P-value type: \code{"one.sided"} (default) or \code{"two.sided"}.}
#'     \item{\code{cov.ridge}}{Ridge stabilization for information matrix inversion. Default is 1e-8.}
#'     \item{\code{print.trace}}{Logical. If \code{TRUE}, prints progress. Default is \code{FALSE}.}
#'   }
#'
#' @return The input \code{object} with additional inference components:
#' \item{sigma2.used}{Estimated \eqn{\sigma^2} used for inference.}
#' \item{C2.se}{Sandwich standard errors for \eqn{C_2} (Q x P2 matrix).}
#' \item{C2.se.boot}{Bootstrap standard errors for \eqn{C_2} (Q x P2 matrix).}
#' \item{C2.ci.lower}{Lower CI bounds for \eqn{C_2} (Q x P2 matrix).}
#' \item{C2.ci.upper}{Upper CI bounds for \eqn{C_2} (Q x P2 matrix).}
#' \item{coefficients}{Data frame with Estimate, SE, BSE, z, p-value for each element of \eqn{C_2}.}
#' \item{C2.p.side}{P-value type used.}
#'
#' @seealso \code{\link{nmf.sem}}, \code{\link{nmf.sem.DOT}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to structural equation modeling for blind input-output analysis.
#'   arXiv:2512.18250. \url{https://arxiv.org/abs/2512.18250}
#' @export
#' @examples
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
#' res <- nmf.sem(Y1, Y2, rank = 2)
#' res2 <- nmf.sem.inference(res, Y1, Y2)
#' res2$coefficients
#'
nmf.sem.inference <- function(object, Y1, Y2, wild.bootstrap = TRUE, ...) {
  if (is.null(object$X) || is.null(object$C1) || is.null(object$C2))
    stop("object must contain X, C1, and C2 (returned by nmf.sem).")

  extra_args <- base::list(...)
  wild.B      <- if (!is.null(extra_args$wild.B))      extra_args$wild.B      else 500
  wild.seed   <- if (!is.null(extra_args$wild.seed))   extra_args$wild.seed   else 123
  wild.level  <- if (!is.null(extra_args$wild.level))  extra_args$wild.level  else 0.95
  sandwich    <- if (!is.null(extra_args$sandwich))     extra_args$sandwich    else TRUE
  C.p.side    <- if (!is.null(extra_args$C.p.side))    extra_args$C.p.side    else "one.sided"
  cov.ridge   <- if (!is.null(extra_args$cov.ridge))   extra_args$cov.ridge   else 1e-8
  print.trace <- if (!is.null(extra_args$print.trace)) extra_args$print.trace else FALSE

  Y1 <- base::as.matrix(Y1)
  Y2 <- base::as.matrix(Y2)

  X     <- object$X    # P1 x Q
  C1    <- object$C1   # Q x P1
  C2    <- object$C2   # Q x P2

  Q  <- base::ncol(X)
  P1 <- base::nrow(Y1)
  P2 <- base::nrow(Y2)
  N  <- base::ncol(Y1)

  # Residual after removing endogenous feedback: R = Y1 - X*C1*Y1
  R_endo <- Y1 - X %*% C1 %*% Y1          # P1 x N
  R_C2   <- R_endo - X %*% C2 %*% Y2      # P1 x N  (full residual)

  # sigma2 estimate
  denom <- base::max(P1 * N - Q * P2, 1)
  sigma2.used <- base::sum(R_C2^2) / denom

  # Information matrix: I = sigma^{-2} (Y2*Y2' x X'X)
  XtX  <- base::crossprod(X)                # Q x Q
  Y2Y2 <- base::tcrossprod(Y2)              # P2 x P2
  Info_core <- base::kronecker(Y2Y2, XtX)   # QP2 x QP2
  Info <- Info_core / base::max(sigma2.used, 1e-12)
  Info <- Info + base::diag(cov.ridge, base::nrow(Info))

  Hinv <- base::tryCatch(base::solve(Info), error = function(e) {
    if (base::requireNamespace("MASS", quietly = TRUE)) MASS::ginv(Info)
    else base::stop("Information matrix singular; install MASS package.")
  })

  # Sandwich covariance: V = Hinv J Hinv
  V_sand <- NULL
  if (base::isTRUE(sandwich)) {
    Xt <- base::t(X)
    J <- base::matrix(0, Q * P2, Q * P2)
    for (n in 1:N) {
      y2_n <- Y2[, n, drop = FALSE]       # P2 x 1
      r_n  <- R_C2[, n, drop = FALSE]     # P1 x 1
      g_n  <- Xt %*% r_n                  # Q x 1
      S_n  <- -(g_n %*% base::t(y2_n)) / base::max(sigma2.used, 1e-12)  # Q x P2
      s_n  <- base::as.vector(S_n)
      J    <- J + base::tcrossprod(s_n)
    }
    if (N > 1) J <- (N / (N - 1)) * J    # CR1 correction
    V_sand <- Hinv %*% J %*% Hinv
  }

  C.vec.cov <- if (!is.null(V_sand)) V_sand else Hinv

  # Sandwich SE
  se_vec <- base::sqrt(base::pmax(base::diag(C.vec.cov), 0))
  C2.se <- base::matrix(se_vec, nrow = Q, ncol = P2, byrow = FALSE)

  # ---- Wild bootstrap (one-step Newton) ----
  C2.se.boot  <- NULL
  C2.ci.lower <- NULL
  C2.ci.upper <- NULL

  if (base::isTRUE(wild.bootstrap)) {
    base::set.seed(wild.seed)
    Xt <- base::t(X)
    score_mat <- base::matrix(0, Q * P2, N)
    for (n in 1:N) {
      y2_n <- Y2[, n, drop = FALSE]
      r_n  <- R_C2[, n, drop = FALSE]
      g_n  <- Xt %*% r_n
      G_n  <- -(g_n %*% base::t(y2_n)) / base::max(sigma2.used, 1e-12)
      score_mat[, n] <- base::as.vector(G_n)
    }

    C2_hat_vec <- base::as.vector(C2)
    C2_boot <- base::matrix(NA_real_, nrow = Q * P2, ncol = wild.B)
    for (b in 1:wild.B) {
      w <- stats::rexp(N, rate = 1) - 1       # Exp(1)-centered multiplier
      grad_b <- base::as.vector(score_mat %*% w)
      c_b <- C2_hat_vec - base::as.vector(Hinv %*% grad_b)
      c_b <- base::pmax(c_b, 0)               # project onto C2 >= 0
      C2_boot[, b] <- c_b
    }

    # Bootstrap SE
    sd_vec <- base::apply(C2_boot, 1, stats::sd, na.rm = TRUE)
    C2.se.boot <- base::matrix(sd_vec, nrow = Q, ncol = P2, byrow = FALSE)

    # Bootstrap CI
    alpha <- 1 - wild.level
    lo <- base::apply(C2_boot, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
    hi <- base::apply(C2_boot, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
    C2.ci.lower <- base::matrix(lo, nrow = Q, ncol = P2, byrow = FALSE)
    C2.ci.upper <- base::matrix(hi, nrow = Q, ncol = P2, byrow = FALSE)
  }

  # ---- Coefficients table ----
  Estimate <- base::as.vector(C2)
  SE  <- base::as.vector(C2.se)
  BSE <- if (!is.null(C2.se.boot)) base::as.vector(C2.se.boot) else base::rep(NA_real_, base::length(Estimate))
  z_value <- base::ifelse(SE > 0, Estimate / SE, NA_real_)

  if (C.p.side == "one.sided") {
    p_value <- base::ifelse(base::is.finite(z_value), stats::pnorm(z_value, lower.tail = FALSE), NA_real_)
  } else {
    p_value <- base::ifelse(base::is.finite(z_value), 1 - stats::pchisq(z_value^2, df = 1), NA_real_)
  }

  rlabs <- if (!is.null(base::rownames(C2))) base::rownames(C2) else base::paste0("Factor", 1:Q)
  clabs <- if (!is.null(base::colnames(C2))) base::colnames(C2) else base::paste0("Y2_", 1:P2)

  coefficients <- base::data.frame(
    Basis     = base::rep(rlabs, times = P2),
    Covariate = base::rep(clabs, each = Q),
    Estimate  = Estimate,
    SE        = SE,
    BSE       = BSE,
    z_value   = z_value,
    p_value   = p_value,
    CI_low    = if (!is.null(C2.ci.lower)) base::as.vector(C2.ci.lower) else NA_real_,
    CI_high   = if (!is.null(C2.ci.upper)) base::as.vector(C2.ci.upper) else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )

  if (print.trace) {
    if (base::isTRUE(wild.bootstrap)) {
      base::message("  Inference: sandwich SE + wild bootstrap done.")
    } else {
      base::message("  Inference: sandwich SE done (wild bootstrap skipped).")
    }
  }

  object$sigma2.used  <- sigma2.used
  object$C2.se        <- C2.se
  object$C2.se.boot   <- C2.se.boot
  object$C2.ci.lower  <- C2.ci.lower
  object$C2.ci.upper  <- C2.ci.upper
  object$coefficients <- coefficients
  object$C2.p.side    <- C.p.side
  return(object)
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
nmf.sem.cv <- function(
    Y1, Y2,
    rank = NULL,
    X.init = NULL,
    X.L2.ortho = 100.0,
    C1.L1 = 1.0,        # L1 sparsity for C1 (Theta1)
    C2.L1 = 0.1,        # L1 sparsity for C2 (Theta2)
    epsilon = 1e-6,     # Convergence tolerance passed to nmf.sem
    maxit = 20000,
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
#'   resulting variable split. (Default: \code{FALSE})
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
nmf.sem.split <- function(x, n.exogenous = NULL, threshold = 0.1,
                          auto.flipped = TRUE, verbose = FALSE) {

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
#' @param hide.isolated Logical. If \code{TRUE} (default), Y1 and Y2 nodes
#'   that have no edges at or above \code{threshold} are excluded from the graph.
#' @param sig.level Significance level for filtering C2 edges when inference
#'   results are present. If \code{result} contains a \code{coefficients} data
#'   frame (from \code{\link{nmf.sem.inference}}), only edges with
#'   \code{p_value < sig.level} are drawn, with significance stars appended.
#'   Set to \code{NULL} to disable filtering. Default is \code{0.1}.
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
#' @seealso \code{\link{nmf.sem}}, \code{\link{nmf.sem.inference}},
#'   \code{\link{plot.nmfkc.DOT}}
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
                        cluster.labels        = NULL,
                        hide.isolated         = TRUE,
                        sig.level             = 0.1) {

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
  ## Significance stars for C2 edges (from nmf.sem.inference)
  ## ---------------------------------------------------------------
  C2_stars <- NULL
  C2_show  <- NULL
  if (!is.null(result$coefficients)) {
    C2_stars <- matrix("", nrow = Q, ncol = P2)
    C2_pval  <- matrix(NA_real_, nrow = Q, ncol = P2)
    cf <- result$coefficients
    fac_names <- rownames(C2)
    exo_names <- colnames(C2)
    for (k in seq_len(nrow(cf))) {
      q  <- match(cf$Basis[k], fac_names)
      p2 <- match(cf$Covariate[k], exo_names)
      if (!is.na(q) && !is.na(p2) && !is.na(cf$p_value[k])) {
        p <- cf$p_value[k]
        C2_pval[q, p2] <- p
        if (p < 0.001)      C2_stars[q, p2] <- "***"
        else if (p < 0.01)  C2_stars[q, p2] <- "**"
        else if (p < 0.05)  C2_stars[q, p2] <- "*"
      }
    }
    if (!is.null(sig.level)) {
      C2_show <- !is.na(C2_pval) & C2_pval < sig.level
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
          pen <- pw(weight, max_C2, weight_scale_y2f)
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
#' @param threshold Minimum coefficient magnitude to display an edge.
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
      # Y is connected via XC: XC_mat[i, ] >= threshold
      used_Y <- apply(XC_mat, 1L, function(row) any(row >= threshold, na.rm = TRUE))
      idx_Y <- which(used_Y)
    }
    if (hasA && type %in% c("YXA", "YA")) {
      if (type == "YXA" && !is.null(C)) {
        # A is connected via C: C[, k] >= threshold
        used_A <- apply(C, 2L, function(col) any(col >= threshold, na.rm = TRUE))
        idx_A <- which(used_A)
      } else if (type == "YA" && !is.null(XC_mat)) {
        # A is connected via XC: XC_mat[, k] >= threshold
        used_A <- apply(XC_mat, 2L, function(col) any(col >= threshold, na.rm = TRUE))
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

    max_C <- suppressWarnings(max(C, na.rm = TRUE))
    if (is.finite(max_C) && max_C > 0) {
      for (q in seq_len(Q)) {
        for (k in idx_A) {
          val <- C[q, k]
          show <- if (!is.null(C_show)) C_show[q, k]
                  else is.finite(val) && val >= threshold
          if (show) {
            pen <- pw(val, max_C, weight_scale_ax)
            lab <- fmtc(val)
            if (!is.null(C_stars)) lab <- paste0(lab, C_stars[q, k])
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

    max_XC <- suppressWarnings(max(XC_mat, na.rm = TRUE))
    if (is.finite(max_XC) && max_XC > 0) {
      for (i in idx_Y) {
        for (k in idx_A) {
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

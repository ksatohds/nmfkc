# nmfae.R
# 3-layer NMF (NMF-AE): Y1 ≈ X1 C X2 Y2
# Author: Kenichi Satoh
# Date: 2026-03-07
#
# Depends: nmfkc
# Suggests: DiagrammeR, DiagrammeRsvg, rsvg (for DOT graph rendering)

#' @title Three-Layer Non-negative Matrix Factorization (NMF-AE)
#' @description
#' \code{nmfae} fits a three-layer nonnegative matrix factorization model
#' \eqn{Y_1 \approx X_1 \Theta X_2 Y_2}, where \eqn{X_1} is a decoder basis
#' (column sum 1), \eqn{\Theta} is a bottleneck parameter matrix,
#' \eqn{X_2} is an encoder basis (row sum 1), and \eqn{Y_2} is the input matrix.
#'
#' When \code{Y2 = Y1}, the model acts as a non-negative autoencoder.
#' When \code{Y1 != Y2}, it acts as a heteroencoder.
#'
#' Initialization uses a three-step NMF procedure via \code{\link{nmfkc}}:
#' (1) \code{nmfkc(Y1, rank=Q)} to obtain \eqn{X_1},
#' (2) \code{nmfkc(Y1, A=Y2, rank=Q)} with fixed \eqn{X_1} to obtain \eqn{C = \Theta X_2},
#' (3) \code{nmfkc(Y2, rank=R)} to factor \eqn{C} into \eqn{\Theta} and \eqn{X_2}.
#'
#' @param Y1 Output matrix \eqn{Y_1} (P1 x N). Non-negative. May contain \code{NA}s
#'   (handled via \code{Y1.weights}).
#' @param Y2 Input matrix \eqn{Y_2} (P2 x N). Non-negative. Default is \code{Y1} (autoencoder).
#' @param rank Integer. Rank of the decoder basis \eqn{X_1} (P1 x Q). Default is 2.
#'   For backward compatibility, \code{Q} is accepted via \code{...}.
#' @param rank.encoder Integer. Rank of the encoder basis \eqn{X_2} (R x P2).
#'   Default is \code{rank}. For backward compatibility, \code{R} is accepted via \code{...}.
#' @param epsilon Positive convergence tolerance. Default is \code{1e-4}.
#' @param maxit Maximum number of multiplicative update iterations. Default is 5000.
#' @param verbose Logical. If \code{TRUE}, prints progress messages during fitting. Default is \code{FALSE}.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{Y1.weights}}{Optional non-negative weight matrix
#'       (P1 x N) or vector for \eqn{Y_1}, analogous to the
#'       \code{weights} argument of \code{\link[stats]{lm}}.  Loss becomes
#'       \eqn{\sum W_{ij} \, (Y_{1,ij} - \hat Y_{1,ij})^2}
#'       (\code{lm()}-style, \strong{linear} in \eqn{W}).  Logical
#'       matrices (\code{TRUE} / \code{FALSE}) are also accepted.
#'       Typical ECV / CV usage passes a binary mask
#'       \eqn{W \in \{0,1\}} for held-out elements; real-valued weights
#'       for importance weighting are also supported.  Default: if
#'       \code{Y1} has \code{NA}, a binary mask is auto-generated
#'       (0 for \code{NA}, 1 elsewhere).}
#'     \item{\code{C.L1}}{L1 regularization parameter for \eqn{C}. Default is 0.}
#'     \item{\code{X1.L2.ortho}}{L2 orthogonality regularization for \eqn{X_1} columns. Default is 0.}
#'     \item{\code{X2.L2.ortho}}{L2 orthogonality regularization for \eqn{X_2} rows. Default is 0.}
#'     \item{\code{seed}}{Integer seed for reproducibility. Default is 123.}
#'     \item{\code{print.trace}}{Logical. If \code{TRUE}, prints progress. Default is \code{FALSE}.}
#'   }
#'
#' @return An object of class \code{"nmfae"}, a list with components:
#' \item{X1}{Decoder basis matrix (P1 x Q), column sum 1.}
#' \item{C}{Parameter matrix (Q x R).}
#' \item{X2}{Encoder basis matrix (R x P2), row sum 1.}
#' \item{Y1hat}{Fitted values \eqn{X_1 \Theta X_2 Y_2} (P1 x N).}
#' \item{rank}{Named integer vector \code{c(Q, R)}.}
#' \item{objfunc}{Final objective value.}
#' \item{objfunc.iter}{Objective values by iteration.}
#' \item{r.squared}{\eqn{\mathrm{cor}(Y, \widehat Y)^2} (Pearson; in \eqn{[0,1]}).}
#' \item{r.squared.uncentered}{Uncentered \eqn{R^2 = 1 - \|Y - \widehat Y\|_F^2 / \|Y\|_F^2} (baseline = zero matrix).}
#' \item{r.squared.centered}{Row-mean centered \eqn{1 - \|Y - \widehat Y\|_F^2 / \|Y - \bar Y_{p\cdot}\|_F^2}.}
#' \item{niter}{Number of iterations performed.}
#' \item{runtime}{Elapsed time as a \code{difftime} object.}
#' \item{n.missing}{Number of missing (or zero-weighted) elements in \eqn{Y_1}.}
#' \item{n.total}{Total number of elements in \eqn{Y_1} (P1 x N).}
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in future versions.
#'
#' @seealso \code{\link{nmfae.inference}}, \code{\link{predict.nmfae}}, \code{\link{nmfae.ecv}}, \code{\link{nmfae.DOT}}, \code{\link{nmfkc}}
#' @export
#' @source Satoh, K. (2025). Applying Non-negative Matrix Factorization with Covariates
#'   to Multivariate Time Series. \emph{Japanese Journal of Statistics and Data Science}.
#' @references
#' Lee, D. D. and Seung, H. S. (2001). Algorithms for Non-negative Matrix Factorization.
#'   \emph{Advances in Neural Information Processing Systems}, 13.
#'
#' Saha, S. et al. (2022). Hierarchical Deep Learning Neural Network (HiDeNN):
#'   An Artificial Intelligence (AI) Framework for Computational Science and Engineering.
#'   \emph{Computer Methods in Applied Mechanics and Engineering}, 399.
#' @examples
#' # Autoencoder example
#' Y <- matrix(c(1,0,1,0, 0,1,0,1, 1,1,0,0), nrow=3, byrow=TRUE)
#' res <- nmfae(Y, rank=2, rank.encoder=2)
#' res$r.squared
#'
#' # Heteroencoder example
#' Y1 <- matrix(c(1,0,0,1), nrow=2)
#' Y2 <- matrix(runif(8), nrow=4)
#' res2 <- nmfae(Y1, Y2, rank=2, rank.encoder=2)
#'
nmfae <- function(Y1, Y2 = Y1, rank = 2, rank.encoder = rank,
                  epsilon = 1e-4, maxit = 5000, verbose = FALSE, ...) {

  cl <- match.call()

  extra_args <- list(...)
  # backward compatibility: Q -> rank, R -> rank.encoder
  if (!is.null(extra_args$Q)) rank <- extra_args$Q
  if (!is.null(extra_args$R)) rank.encoder <- extra_args$R
  Q <- rank
  R <- rank.encoder
  Y1.weights  <- if (!is.null(extra_args$Y1.weights))  extra_args$Y1.weights  else NULL
  C.L1        <- if (!is.null(extra_args$C.L1))        extra_args$C.L1        else 0
  X1.L2.ortho <- if (!is.null(extra_args$X1.L2.ortho)) extra_args$X1.L2.ortho else 0
  X2.L2.ortho <- if (!is.null(extra_args$X2.L2.ortho)) extra_args$X2.L2.ortho else 0
  seed        <- if (!is.null(extra_args$seed))        extra_args$seed        else 123
  print.trace <- verbose
  if (!is.null(extra_args$print.trace)) print.trace <- extra_args$print.trace  # backward compat

  # --- Input validation ---
  Y1 <- as.matrix(Y1); storage.mode(Y1) <- "double"
  Y2 <- as.matrix(Y2); storage.mode(Y2) <- "double"
  if (min(Y2, na.rm = TRUE) < 0) stop("Y2 must be non-negative.")
  if (ncol(Y1) != ncol(Y2)) stop("Y1 and Y2 must have the same number of columns (N).")

  P1 <- nrow(Y1); P2 <- nrow(Y2); N <- ncol(Y1)
  eps <- 1e-10
  start.time <- Sys.time()

  # --- Y1.weights handling (cf. nmfkc Y.weights) ---
  # Vector -> matrix expansion
  if (!is.null(Y1.weights) && is.vector(Y1.weights)) {
    if (length(Y1.weights) == N) {
      Y1.weights <- matrix(Y1.weights, nrow = P1, ncol = N, byrow = TRUE)
    } else if (length(Y1.weights) == 1) {
      Y1.weights <- matrix(Y1.weights, nrow = P1, ncol = N)
    } else {
      stop("Length of Y1.weights vector must match ncol(Y1) or be 1.")
    }
  }
  # NULL: auto-detect NAs
  if (is.null(Y1.weights)) {
    if (any(is.na(Y1))) {
      Y1.weights <- matrix(1, nrow = P1, ncol = N)
      Y1.weights[is.na(Y1)] <- 0
      Y1[is.na(Y1)] <- 0
      if (print.trace) message("  Notice: NAs in Y1 detected, Y1.weights set automatically.")
    }
  } else {
    if (!is.matrix(Y1.weights)) Y1.weights <- as.matrix(Y1.weights)
    if (!all(dim(Y1.weights) == dim(Y1))) stop("Dimension mismatch between Y1 and Y1.weights.")
    Y1.weights[is.na(Y1.weights)] <- 0
    Y1[is.na(Y1)] <- 0
    Y1[Y1.weights == 0] <- 0
  }
  # Determine code path
  has.weights <- !is.null(Y1.weights) && any(Y1.weights != 1)
  if (!has.weights) Y1.weights <- NULL

  if (min(Y1, na.rm = TRUE) < 0) stop("Y1 must be non-negative.")

  if (print.trace) {
    message(sprintf("Y1(%d,%d) ~ X1(%d,%d) C(%d,%d) X2(%d,%d) Y2(%d,%d)",
                    P1, N, P1, Q, Q, R, R, P2, P2, N))
  }

  # === Initialization using nmfkc ===
  # Step 1: X1 from nmfkc(Y1, rank=Q)
  if (print.trace) message("  Init step 1: nmfkc(Y1, rank=Q)...")
  res1 <- nmfkc(Y1, rank = Q, seed = seed, print.dims = FALSE,
                Y.weights = Y1.weights)
  X1 <- res1$X  # P1 x Q, column sum 1

  # Step 2: CX2 = C X2 with X1 fixed
  if (print.trace) message("  Init step 2: nmfkc(Y1, A=Y2, X.restriction='fixed')...")
  res2 <- nmfkc(Y1, A = Y2, rank = Q,
                 X.init = X1, X.restriction = "fixed",
                 seed = seed, print.dims = FALSE,
                 Y.weights = Y1.weights)
  CX2_init <- res2$C  # Q x P2

  # Step 3: X2 from nmfkc(Y2, rank=R), then C from CX2 ≈ C X2
  if (print.trace) message("  Init step 3: nmfkc(Y2, rank=R)...")
  res3 <- nmfkc(Y2, rank = R, seed = seed, print.dims = FALSE)
  X2 <- t(res3$X)             # R x P2, row sum 1
  C <- CX2_init %*% t(X2)     # Q x R, non-negative

  # Precompute (only needed for unweighted path)
  if (!has.weights) Y2Y2t <- tcrossprod(Y2)  # P2 x P2

  # === Multiplicative updates ===
  objfunc.iter <- numeric(maxit)

  W <- Y1.weights  # alias for readability (NULL when unweighted)

  for (iter in 1:maxit) {
    # 1. Update X1: Y1 ≈ X1 F, F = C X2 Y2
    F_mat <- C %*% X2 %*% Y2                           # Q x N
    if (has.weights) {
      num_X1 <- (W * Y1) %*% t(F_mat)
      den_X1 <- (W * (X1 %*% F_mat)) %*% t(F_mat) + eps
    } else {
      num_X1 <- tcrossprod(Y1, F_mat)
      den_X1 <- X1 %*% tcrossprod(F_mat) + eps
    }
    if (X1.L2.ortho > 0) {
      X1tX1 <- crossprod(X1); diag(X1tX1) <- 0
      den_X1 <- den_X1 + X1.L2.ortho * (X1 %*% X1tX1)
    }
    X1 <- X1 * (num_X1 / den_X1)

    # 2. Normalize X1 columns -> absorb scale into C (left)
    cs <- colSums(X1)
    X1 <- sweep(X1, 2, cs, "/")
    C <- sweep(C, 1, cs, "*")

    # 3. Update C: Y1 ≈ X1 C G, G = X2 Y2
    G <- X2 %*% Y2                                     # R x N
    if (has.weights) {
      num_C <- crossprod(X1, W * Y1) %*% t(G)
      den_C <- crossprod(X1, W * (X1 %*% C %*% G)) %*% t(G) + eps
    } else {
      num_C <- crossprod(X1, Y1) %*% t(G)
      den_C <- crossprod(X1) %*% C %*% tcrossprod(G) + eps
    }
    if (C.L1 > 0) den_C <- den_C + (C.L1 / 2)
    C <- C * (num_C / den_C)

    # 4. Update X2: Y1 ≈ H X2 Y2, H = X1 C
    H <- X1 %*% C                                      # P1 x R
    if (has.weights) {
      num_X2 <- crossprod(H, W * Y1) %*% t(Y2)
      den_X2 <- crossprod(H, W * (H %*% X2 %*% Y2)) %*% t(Y2) + eps
    } else {
      num_X2 <- crossprod(H, Y1) %*% t(Y2)
      den_X2 <- crossprod(H) %*% X2 %*% Y2Y2t + eps
    }
    if (X2.L2.ortho > 0) {
      X2X2t <- tcrossprod(X2); diag(X2X2t) <- 0
      den_X2 <- den_X2 + X2.L2.ortho * (X2X2t %*% X2)
    }
    X2 <- X2 * (num_X2 / den_X2)

    # 5. Normalize X2 rows -> absorb scale into C (right)
    rs <- rowSums(X2)
    X2 <- sweep(X2, 1, rs, "/")
    C <- sweep(C, 2, rs, "*")

    # Objective function (with regularization penalties)
    # lm()-style weighted least squares: L = sum(W * (Y1 - Y1hat)^2).
    # The MU (num_X1, num_C, etc.) carries W linearly, so reporting the
    # linear-W objective here keeps MU target and reported loss consistent.
    # For binary W in {0,1} (standard ECV / NA-mask case) this is identical
    # to sum((W*(Y1-Y1hat))^2) since W == W^2.
    Y1hat <- X1 %*% C %*% X2 %*% Y2
    if (has.weights) {
      obj <- sum(W * (Y1 - Y1hat)^2)
    } else {
      obj <- sum((Y1 - Y1hat)^2)
    }
    if (C.L1 > 0) obj <- obj + C.L1 * sum(C)
    if (X1.L2.ortho > 0) {
      X1tX1 <- crossprod(X1); diag(X1tX1) <- 0
      obj <- obj + (X1.L2.ortho / 2) * sum(X1tX1^2)
    }
    if (X2.L2.ortho > 0) {
      X2X2t <- tcrossprod(X2); diag(X2X2t) <- 0
      obj <- obj + (X2.L2.ortho / 2) * sum(X2X2t^2)
    }
    objfunc.iter[iter] <- obj

    if (print.trace && (iter %% 100 == 0 || iter == 1)) {
      message(sprintf("  iter %5d: objfunc = %.6f", iter, objfunc.iter[iter]))
    }

    # Convergence check
    if (iter > 1) {
      rel_change <- abs(objfunc.iter[iter] - objfunc.iter[iter-1]) /
                       (objfunc.iter[iter-1] + eps)
      if (rel_change < epsilon) {
        if (print.trace) message(sprintf("  Converged at iter %d", iter))
        break
      }
    }
  }
  ## Warn when the MU loop exhausts maxit without meeting the
  ## relative-tolerance criterion (matches nmfkc() / nmf.sem() convention).
  if (iter == maxit && exists("rel_change") && rel_change >= epsilon)
    warning(paste0("maximum iterations (", maxit, ") reached..."))

  niter <- iter
  objfunc.iter <- objfunc.iter[1:niter]
  diff.time <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))

  # --- Reorder bases by centroid position (cf. nmfkc) ---
  # X1 (P1 x Q): sort columns by weighted centroid of row indices
  if (Q > 1) {
    idx1 <- order(matrix(seq_len(P1) / P1, nrow = 1) %*% X1)
    X1 <- X1[, idx1, drop = FALSE]
    C <- C[idx1, , drop = FALSE]
  }
  # X2 (R x P2): sort rows by weighted centroid of column indices
  if (R > 1) {
    idx2 <- order(X2 %*% matrix(seq_len(P2) / P2, ncol = 1))
    X2 <- X2[idx2, , drop = FALSE]
    C <- C[, idx2, drop = FALSE]
  }

  # --- Assign Dec/Enc names ---
  colnames(X1) <- paste0("Dec", 1:Q)
  rownames(C)  <- paste0("Dec", 1:Q)
  colnames(C)  <- paste0("Enc", 1:R)
  rownames(X2) <- paste0("Enc", 1:R)

  Y1hat <- X1 %*% C %*% X2 %*% Y2
  objfunc <- utils::tail(objfunc.iter, 1)

  # R-squared, sigma, mae, and missing count
  if (has.weights) {
    valid <- (W > 0)
    n.missing <- sum(!valid)
    n.valid <- sum(valid)
    r2_all <- .r.squared.all(Y1, Y1hat, Y.weights = W)
    sigma <- sqrt(objfunc / n.valid)
    mae <- mean(abs(Y1[valid] - Y1hat[valid]))
  } else {
    n.missing <- 0L
    n.valid <- P1 * N
    r2_all <- .r.squared.all(Y1, Y1hat)
    sigma <- sqrt(objfunc / n.valid)
    mae <- mean(abs(Y1 - Y1hat))
  }
  r.squared          <- r2_all$r.squared
  r.squared.uncentered     <- r2_all$r.squared.uncentered
  r.squared.centered <- r2_all$r.squared.centered

  if (print.trace) {
    msg <- sprintf("  Done: %d iterations, %.1f sec, R2 = %.4f", niter, diff.time, r.squared)
    if (n.missing > 0) msg <- paste0(msg, sprintf(", missing=%d", n.missing))
    message(msg)
  }

  # --- Soft/hard clustering (cf. nmfkc B.prob / B.cluster) ---
  H <- C %*% X2 %*% Y2   # Q x N encoding
  eps_bp <- 1e-16
  B.prob <- t( t(H) / (colSums(H) + eps_bp) )   # column-normalized
  B.cluster <- apply(B.prob, 2, which.max)
  B.cluster[colSums(H) == 0] <- NA

  result <- list(
    call = cl,
    X1 = X1,
    C = C,
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
    runtime = diff.time,
    n.missing = n.missing,
    n.total = P1 * N
  )
  class(result) <- c("nmfae", "nmf")
  return(result)
}

#' @title Statistical Inference for NMF-AE Parameter Matrix
#' @description
#' Performs post-estimation inference for \eqn{\Theta} in the three-layer NMF model
#' \eqn{Y_1 \approx X_1 \Theta X_2 Y_2}, conditional on \eqn{(\hat{X}_1, \hat{X}_2)}.
#' Uses sandwich covariance estimation and one-step wild bootstrap with
#' non-negative projection.
#'
#' @param object An object of class \code{"nmfae"} returned by \code{\link{nmfae}}.
#' @param Y1 Output matrix \eqn{Y_1} (P1 x N). Must match the data used in \code{nmfae()}.
#' @param Y2 Input matrix \eqn{Y_2} (P2 x N). Default is \code{Y1} (autoencoder).
#' @param wild.bootstrap Logical. If \code{TRUE} (default), performs wild bootstrap
#'   for bootstrap SE and confidence intervals. If \code{FALSE}, only sandwich SE
#'   and z-test p-values are computed (faster).
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
#' @return An object of class \code{c("nmfae.inference", "nmfae")}, inheriting all
#' components from the input \code{object}, with additional inference components:
#' \item{sigma2.used}{Estimated \eqn{\sigma^2} used for inference.}
#' \item{C.se}{Sandwich standard errors for \eqn{\Theta} (Q x R matrix).}
#' \item{C.se.boot}{Bootstrap standard errors for \eqn{\Theta} (Q x R matrix).}
#' \item{C.ci.lower}{Lower CI bounds for \eqn{\Theta} (Q x R matrix).}
#' \item{C.ci.upper}{Upper CI bounds for \eqn{\Theta} (Q x R matrix).}
#' \item{coefficients}{Data frame with Estimate, SE, BSE, z, p-value for each element of \eqn{\Theta}.}
#' \item{C.p.side}{P-value type used.}
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @seealso \code{\link{nmfae}}, \code{\link{summary.nmfae.inference}}
#' @export
#' @examples
#' Y <- matrix(c(1,0,1,0, 0,1,0,1, 1,1,0,0), nrow=3, byrow=TRUE)
#' res <- nmfae(Y, rank=2, rank.encoder=2)
#' res2 <- nmfae.inference(res, Y)
#' summary(res2)
#'
nmfae.inference <- function(object, Y1, Y2 = Y1,
                            wild.bootstrap = TRUE, ...) {
  if (!inherits(object, "nmfae")) stop("object must be of class 'nmfae'")

  extra_args <- list(...)
  wild.B     <- if (!is.null(extra_args$wild.B))     extra_args$wild.B     else 500
  wild.seed  <- if (!is.null(extra_args$wild.seed))  extra_args$wild.seed  else 123
  wild.level <- if (!is.null(extra_args$wild.level)) extra_args$wild.level else 0.95
  sandwich   <- if (!is.null(extra_args$sandwich))   extra_args$sandwich   else TRUE
  C.p.side   <- if (!is.null(extra_args$C.p.side))   extra_args$C.p.side   else "one.sided"
  cov.ridge  <- if (!is.null(extra_args$cov.ridge))  extra_args$cov.ridge  else 1e-8
  print.trace <- if (!is.null(extra_args$print.trace)) extra_args$print.trace else FALSE

  X1 <- object$X1
  C  <- object$C
  X2 <- object$X2
  Y1hat <- object$Y1hat
  Q  <- object$rank["Q"]
  R  <- object$rank["R"]
  P1 <- nrow(Y1)
  N  <- ncol(Y1)

  Z   <- X2 %*% Y2                        # R x N
  R_C <- Y1 - Y1hat                       # P1 x N  residuals

  # sigma2 estimate (denominator: PN - QR)
  denom <- max(P1 * N - Q * R, 1)
  sigma2.used <- sum(R_C^2) / denom

  # Information matrix: I = sigma^{-2} (ZZ' x X1'X1)
  X1tX1 <- crossprod(X1)                  # Q x Q
  ZZt   <- tcrossprod(Z)                  # R x R
  Info_core <- kronecker(ZZt, X1tX1)      # QR x QR
  Info <- Info_core / max(sigma2.used, 1e-12)
  Info <- Info + diag(cov.ridge, nrow(Info))

  Hinv <- tryCatch(solve(Info), error = function(e) {
    if (requireNamespace("MASS", quietly = TRUE)) MASS::ginv(Info)
    else stop("Information matrix singular; install MASS package.")
  })

  # Sandwich covariance: V = Hinv J Hinv
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
    if (N > 1) J <- (N / (N - 1)) * J     # CR1 correction
    V_sand <- Hinv %*% J %*% Hinv
  }

  C.vec.cov <- if (!is.null(V_sand)) V_sand else Hinv

  # Sandwich SE
  se_vec <- sqrt(pmax(diag(C.vec.cov), 0))
  C.se <- matrix(se_vec, nrow = Q, ncol = R, byrow = FALSE)

  # ---- Wild bootstrap (one-step Newton) ----
  C.se.boot <- NULL
  C.ci.lower <- NULL
  C.ci.upper <- NULL

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

    C_hat_vec <- as.vector(C)
    C_boot <- matrix(NA_real_, nrow = Q * R, ncol = wild.B)
    for (b in 1:wild.B) {
      w <- stats::rexp(N, rate = 1) - 1     # Exp(1)-centered multiplier
      grad_b <- as.vector(score_mat %*% w)
      c_b <- C_hat_vec - as.vector(Hinv %*% grad_b)
      c_b <- pmax(c_b, 0)                   # project onto Theta >= 0
      C_boot[, b] <- c_b
    }

    # Bootstrap SE
    sd_vec <- apply(C_boot, 1, stats::sd, na.rm = TRUE)
    C.se.boot <- matrix(sd_vec, nrow = Q, ncol = R, byrow = FALSE)

    # Bootstrap CI
    alpha <- 1 - wild.level
    lo <- apply(C_boot, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
    hi <- apply(C_boot, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
    C.ci.lower <- matrix(lo, nrow = Q, ncol = R, byrow = FALSE)
    C.ci.upper <- matrix(hi, nrow = Q, ncol = R, byrow = FALSE)
  }

  # ---- Coefficients table ----
  Estimate <- as.vector(C)
  SE <- as.vector(C.se)
  BSE <- if (!is.null(C.se.boot)) as.vector(C.se.boot) else rep(NA_real_, length(Estimate))
  z_value <- ifelse(SE > 0, Estimate / SE, NA_real_)

  if (C.p.side == "one.sided") {
    p_value <- ifelse(is.finite(z_value), stats::pnorm(z_value, lower.tail = FALSE), NA_real_)
  } else {
    p_value <- ifelse(is.finite(z_value), 1 - stats::pchisq(z_value^2, df = 1), NA_real_)
  }

  # Row/column labels for C
  rlabs <- if (!is.null(rownames(C))) rownames(C) else paste0("Dec", 1:Q)
  clabs <- if (!is.null(colnames(C))) colnames(C) else paste0("Enc", 1:R)

  coefficients <- data.frame(
    Basis    = rep(rlabs, times = R),
    Covariate = rep(clabs, each = Q),
    Estimate = Estimate,
    SE       = SE,
    BSE      = BSE,
    z_value  = z_value,
    p_value  = p_value,
    CI_low   = if (!is.null(C.ci.lower)) as.vector(C.ci.lower) else NA_real_,
    CI_high  = if (!is.null(C.ci.upper)) as.vector(C.ci.upper) else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )

  if (print.trace) {
    if (isTRUE(wild.bootstrap)) {
      message("  Inference: sandwich SE + wild bootstrap done.")
    } else {
      message("  Inference: sandwich SE done (wild bootstrap skipped).")
    }
  }

  # Add inference fields to the object
  object$sigma2.used  <- sigma2.used
  object$C.se         <- C.se
  object$C.se.boot    <- C.se.boot
  object$C.ci.lower   <- C.ci.lower
  object$C.ci.upper   <- C.ci.upper
  object$coefficients <- coefficients
  object$C.p.side     <- C.p.side
  class(object) <- c("nmfae.inference", "nmf.inference", "nmfae", "nmf")
  return(object)
}

#' @title Rename decoder and encoder bases
#' @description
#' Assigns user-specified names to the decoder (X1 columns) and encoder
#' (X2 rows) bases of an \code{nmfae} object.  The names propagate to
#' \eqn{\Theta}, the coefficients table, and all downstream displays
#' such as \code{summary}, \code{nmfae.DOT}, and \code{nmfae.heatmap}.
#'
#' @param x An object of class \code{"nmfae"} returned by \code{\link{nmfae}}.
#' @param X1.colnames Character vector of length \eqn{Q} for decoder bases
#'   (columns of \eqn{X_1} / rows of \eqn{\Theta}).  If \code{NULL}
#'   (default), the decoder names are left unchanged.
#' @param X2.rownames Character vector of length \eqn{R} for encoder bases
#'   (rows of \eqn{X_2} / columns of \eqn{\Theta}).  If \code{NULL}
#'   (default), the encoder names are left unchanged.
#'
#' @return A modified copy of \code{x} with updated names.
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(runif(15), nrow = 3)
#' res <- nmfae(Y, rank = 2, rank.encoder = 2)
#' res <- nmfae.rename(res,
#'   X1.colnames = c("Basis1", "Basis2"),
#'   X2.rownames = c("Enc1", "Enc2"))
#' summary(res)
#' }
#' @seealso \code{\link{nmfae}}
#' @export
nmfae.rename <- function(x, X1.colnames = NULL, X2.rownames = NULL) {
  if (!is.null(X1.colnames)) {
    Q <- ncol(x$X1)
    if (length(X1.colnames) != Q) stop("X1.colnames must have length ", Q)
    colnames(x$X1) <- X1.colnames
    rownames(x$C)  <- X1.colnames
    if (!is.null(x$coefficients)) {
      old <- paste0("Dec", 1:Q)
      for (k in seq_len(Q))
        x$coefficients$Basis[x$coefficients$Basis == old[k]] <- X1.colnames[k]
    }
  }
  if (!is.null(X2.rownames)) {
    R <- nrow(x$X2)
    if (length(X2.rownames) != R) stop("X2.rownames must have length ", R)
    rownames(x$X2) <- X2.rownames
    colnames(x$C)  <- X2.rownames
    if (!is.null(x$coefficients)) {
      old <- paste0("Enc", 1:R)
      for (k in seq_len(R))
        x$coefficients$Covariate[x$coefficients$Covariate == old[k]] <- X2.rownames[k]
    }
  }
  x
}

#' \code{plot.nmfae} displays the convergence trajectory of the objective function
#' across iterations. The title shows the achieved \eqn{R^2}.
#'
#' @param x An object of class \code{"nmfae"} returned by \code{\link{nmfae}}.
#' @param ... Additional graphical parameters passed to \code{plot}.
#'
#' @return Invisible \code{NULL}. Called for its side effect (plot).
#' @seealso \code{\link{nmfae}}, \code{\link{nmfae.heatmap}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(runif(20), nrow = 4)
#' res <- nmfae(Y, rank = 2)
#' plot(res)
#' }
#' @export
plot.nmfae <- function(x, ...) {
  extra_args <- list(...)
  args <- list(x = x$objfunc.iter)
  if (is.null(extra_args$main))
    args$main <- paste0("R2 = ", round(x$r.squared, 3))
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  if (is.null(extra_args$type)) args$type <- "l"
  all_args <- c(args, extra_args)
  do.call("plot", all_args)
  invisible(NULL)
}

#' @title Summary method for nmfae objects
#' @description
#' \code{summary.nmfae} produces a summary of a fitted NMF-AE model,
#' including dimensions, convergence status, goodness-of-fit statistics,
#' and structure diagnostics (sparsity of factor matrices).
#'
#' @param object An object of class \code{"nmfae"} returned by \code{\link{nmfae}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.nmfae"}, a list with components:
#' \item{call}{The matched call.}
#' \item{dims}{Named vector \code{c(P1, P2, N)}.}
#' \item{Q}{Decoder rank.}
#' \item{R}{Encoder rank.}
#' \item{n.params}{Total number of parameters (P1*Q + Q*R + R*P2).}
#' \item{autoencoder}{Logical; TRUE if P1 == P2 and Y1 was used as Y2.}
#' \item{niter}{Number of iterations.}
#' \item{runtime}{Elapsed time.}
#' \item{objfunc}{Final objective value.}
#' \item{r.squared}{R-squared.}
#' \item{sigma}{Residual standard error (RMSE).}
#' \item{mae}{Mean absolute error.}
#' \item{n.missing}{Number of missing elements.}
#' \item{prop.missing}{Percentage of missing elements.}
#' \item{X1.sparsity}{Proportion of near-zero elements in X1.}
#' \item{C.sparsity}{Proportion of near-zero elements in C.}
#' \item{X2.sparsity}{Proportion of near-zero elements in X2.}
#'
#' @seealso \code{\link{nmfae}}, \code{\link{print.summary.nmfae}}
#' @export
summary.nmfae <- function(object, ...) {
  ans <- list()
  ans$call <- object$call
  ans$dims <- object$dims
  Q <- object$rank["Q"]
  R <- object$rank["R"]
  ans$Q <- Q
  ans$R <- R

  P1 <- object$dims["P1"]
  P2 <- object$dims["P2"]
  N  <- object$dims["N"]
  ans$n.params <- P1 * Q + Q * R + R * P2

  # Autoencoder detection
  ans$autoencoder <- (P1 == P2) &&
    !is.null(object$call$Y2) == FALSE
  # More robust: check if Y2 argument was not explicitly provided
  if (is.null(object$call$Y2)) {
    ans$autoencoder <- TRUE
  } else {
    ans$autoencoder <- FALSE
  }

  ans$niter <- object$niter
  ans$runtime <- object$runtime

  ans$objfunc <- object$objfunc
  ans$r.squared          <- object$r.squared
  ans$r.squared.uncentered     <- object$r.squared.uncentered
  ans$r.squared.centered <- object$r.squared.centered
  ans$sigma <- object$sigma
  ans$mae <- object$mae
  ## Effective rank of the latent encoding H (Q x N).
  ans$effective.rank <- if (!is.null(object$H)) .effective.rank(object$H) else NA_real_
  ans$rank <- if (!is.null(object$rank)) object$rank[1] else
              if (!is.null(object$H)) nrow(object$H) else NA

  # Missing values
  ans$n.missing <- object$n.missing
  if (!is.null(object$n.missing) && !is.null(object$n.total)) {
    ans$prop.missing <- object$n.missing / object$n.total * 100
  } else {
    ans$prop.missing <- 0
  }

  # Sparsity diagnostics
  if (!is.null(object$X1) && is.matrix(object$X1)) {
    ans$X1.sparsity <- mean(object$X1 < 1e-4)
  }
  if (!is.null(object$C) && is.matrix(object$C)) {
    ans$C.sparsity <- mean(object$C < 1e-4)
  }
  if (!is.null(object$X2) && is.matrix(object$X2)) {
    ans$X2.sparsity <- mean(object$X2 < 1e-4)
  }

  # Inference
  ans$coefficients <- object$coefficients
  ans$C.p.side <- object$C.p.side

  class(ans) <- "summary.nmfae"
  return(ans)
}

#' @title Print method for summary.nmfae objects
#' @description
#' Prints a formatted summary of an NMF-AE model fit.
#'
#' @param x An object of class \code{"summary.nmfae"}.
#' @param digits Minimum number of significant digits to be used.
#' @param max.coef Maximum number of coefficient rows to display. If the table
#'   has more rows, only significant rows (p < 0.05) are shown. Default is 20.
#' @param ... Additional arguments (currently unused).
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @seealso \code{\link{summary.nmfae}}
#' @export
print.summary.nmfae <- function(x, digits = max(3L, getOption("digits") - 3L),
                                max.coef = 20, ...) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  P1 <- x$dims["P1"]; P2 <- x$dims["P2"]; N <- x$dims["N"]
  Q <- x$Q; R <- x$R

  cat("Dimensions:\n")
  type_str <- if (x$autoencoder) "  (autoencoder)" else "  (heteroencoder)"
  cat("  Y1:             ", sprintf("%d x %d", P1, N), "\n")
  cat("  Y2:             ", sprintf("%d x %d", P2, N), type_str, "\n")
  cat("  Decoder rank Q: ", Q, "\n")
  cat("  Encoder rank R: ", R, "\n")
  cat("  Parameters:     ",
      sprintf("X1(%dx%d) + C(%dx%d) + X2(%dx%d) = %d",
              P1, Q, Q, R, R, P2, x$n.params), "\n")

  cat("\nConvergence:\n")
  cat("  Iterations:     ", x$niter, "\n")
  cat("  Runtime:        ", sprintf("%.1f secs", x$runtime), "\n")

  if (x$n.missing > 0) {
    cat("  Missing:        ", x$n.missing,
        sprintf("(%.1f%%)", x$prop.missing), "\n")
  }

  .print.fit.statistics(x, header = "Goodness of fit:", digits = digits)

  .print.structure.diagnostics(
    sparsity = c("Decoder (X1)" = x$X1.sparsity,
                 "Bottleneck (C)" = x$C.sparsity,
                 "Encoder (X2)" = x$X2.sparsity))

  # Coefficients table (inference)
  if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
    cf <- x$coefficients
    n_total <- nrow(cf)
    rnames <- paste0(cf$Covariate, ":", cf$Basis)

    # Determine which rows to display
    if (n_total <= max.coef) {
      show_idx <- seq_len(n_total)
      truncated <- FALSE
    } else {
      sig_idx <- which(cf$p_value < 0.05)
      if (length(sig_idx) == 0) {
        show_idx <- seq_len(min(max.coef, n_total))
        truncated <- TRUE
      } else if (length(sig_idx) <= max.coef) {
        show_idx <- sig_idx
        truncated <- FALSE
      } else {
        show_idx <- sig_idx[seq_len(max.coef)]
        truncated <- TRUE
      }
    }

    n_sig <- sum(cf$p_value < 0.05, na.rm = TRUE)
    cat(sprintf("\nCoefficients (conditional on X1, X2): %d total, %d significant\n",
                n_total, n_sig))
    if (n_total > max.coef) {
      cat(sprintf("  (showing %d significant rows; use res$coefficients for full table)\n",
                  length(show_idx)))
    }

    p_side <- if (!is.null(x$C.p.side)) x$C.p.side else "one.sided"
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

    est <- formatC(cf$Estimate[show_idx], format = "f", digits = 3, width = 9)
    se  <- formatC(cf$SE[show_idx], format = "f", digits = 3, width = 10)
    bse <- formatC(cf$BSE[show_idx], format = "f", digits = 3, width = 6)
    zv  <- formatC(cf$z_value[show_idx], format = "f", digits = 2, width = 7)
    pv_str <- format_pval(cf$p_value[show_idx])
    stars <- sig_stars(cf$p_value[show_idx])
    show_names <- rnames[show_idx]

    max_lw <- max(nchar(show_names))
    hdr <- sprintf("%s %s %s %s %s %s",
                   formatC("Estimate", width = 9),
                   formatC("Std. Error", width = 10),
                   formatC("(Boot)", width = 6),
                   formatC("z value", width = 7),
                   formatC(p_header, width = 8), "")
    cat(sprintf("%s %s\n", formatC("Enc:Dec", width = max_lw), hdr))
    for (i in seq_along(show_names)) {
      cat(sprintf("%s %s %s %s %s %s %s\n",
                  formatC(show_names[i], width = max_lw),
                  est[i], se[i], bse[i], zv[i], pv_str[i], stars[i]))
    }
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  cat("\n")
  invisible(x)
}

#' @title Summary method for nmfae.inference objects
#' @description
#' Produces a summary of a fitted NMF-AE model with inference results,
#' including the coefficients table for \eqn{\Theta}.
#'
#' @param object An object of class \code{"nmfae.inference"}.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{"summary.nmfae.inference"}.
#' @seealso \code{\link{nmfae.inference}}, \code{\link{summary.nmfae}}
#' @export
summary.nmfae.inference <- function(object, ...) {
  ans <- summary.nmfae(object, ...)
  class(ans) <- "summary.nmfae.inference"
  return(ans)
}

#' @title Print method for summary.nmfae.inference objects
#' @description
#' Prints a formatted summary including the coefficients table.
#' @param x An object of class \code{"summary.nmfae.inference"}.
#' @param digits Minimum number of significant digits.
#' @param max.coef Maximum coefficient rows to display. Default is 20.
#' @param ... Additional arguments (currently unused).
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @seealso \code{\link{summary.nmfae.inference}}
#' @export
print.summary.nmfae.inference <- function(x, digits = max(3L, getOption("digits") - 3L),
                                          max.coef = 20, ...) {
  print.summary.nmfae(x, digits = digits, max.coef = max.coef, ...)
}



#' @title Heatmap visualization of nmfae factor matrices
#' @description
#' \code{nmfae.heatmap} displays the three factor matrices \eqn{X_1}, \eqn{\Theta},
#' and \eqn{X_2} as side-by-side heatmaps. This provides an alternative to DOT graph
#' visualization, especially when \eqn{Y_2} has many variables (e.g., kernel matrix).
#'
#' @param x An object of class \code{"nmfae"} returned by \code{\link{nmfae}}.
#' @param Y1.label Character vector of output variable names (rows of \eqn{X_1}).
#' @param X1.label Character vector of decoder basis labels (columns of \eqn{X_1}).
#' @param X2.label Character vector of encoder basis labels (rows of \eqn{X_2}).
#' @param Y2.label Character vector of input variable names (columns of \eqn{X_2}).
#' @param palette Color palette vector. Default is white-orange-red (64 colors).
#' @param ... Not used.
#'
#' @return Invisible \code{NULL}. Called for its side effect (plot).
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @seealso \code{\link{nmfae}}, \code{\link{plot.nmfae}}, \code{\link{nmfae.DOT}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(runif(20), nrow = 4)
#' res <- nmfae(Y, rank = 2)
#' nmfae.heatmap(res)
#' }
#' @export
nmfae.heatmap <- function(x,
                          Y1.label = NULL, X1.label = NULL,
                          X2.label = NULL, Y2.label = NULL,
                          palette = NULL, ...) {
  X1 <- x$X1
  C_mat <- x$C
  X2 <- x$X2
  P1 <- nrow(X1); Q <- ncol(X1)
  R <- nrow(X2); P2 <- ncol(X2)

  if (is.null(palette))
    palette <- grDevices::colorRampPalette(c("white", "orange", "red"))(64)

  # Default labels: use dimnames from the object, then fallback to prefix
  resolve_lab <- function(lab, dimnames_default, expected_len, prefix) {
    if (!is.null(lab) && length(lab) != expected_len) {
      warning(sprintf("Label length (%d) does not match dimension (%d). Labels ignored.",
                      length(lab), expected_len))
      lab <- NULL
    }
    if (is.null(lab)) lab <- dimnames_default
    if (is.null(lab)) lab <- paste0(prefix, 1:expected_len)
    lab
  }
  Y1.label <- resolve_lab(Y1.label, rownames(X1), P1, "Y1.")
  X1.label <- resolve_lab(X1.label, colnames(X1), Q,  "Q")
  X2.label <- resolve_lab(X2.label, rownames(X2), R,  "R")
  Y2.label <- resolve_lab(Y2.label, colnames(X2), P2, "Y2.")

  # Helper: plot a matrix as heatmap (row 1 at top)
  plot_mat <- function(mat, main, row.lab, col.lab, mar) {
    opar <- graphics::par(mar = mar)
    on.exit(graphics::par(opar))
    nr <- nrow(mat); nc <- ncol(mat)
    graphics::image(1:nc, 1:nr, t(mat[nr:1, , drop = FALSE]),
                    col = palette, axes = FALSE,
                    xlab = "", ylab = "", main = main)
    graphics::box()
    if (!is.null(col.lab)) {
      cex_col <- min(0.8, 20 / length(col.lab))
      graphics::axis(1, at = 1:nc, labels = col.lab, las = 2, cex.axis = cex_col)
    }
    if (!is.null(row.lab)) {
      cex_row <- min(0.8, 20 / length(row.lab))
      graphics::axis(2, at = 1:nr, labels = rev(row.lab), las = 1, cex.axis = cex_row)
    }
  }

  opar <- graphics::par(mfrow = c(1, 3))
  on.exit(graphics::par(opar))

  plot_mat(X1, sprintf("X1 (%d x %d)", P1, Q), Y1.label, X1.label,
           mar = c(5, 6, 3, 1))
  plot_mat(C_mat, sprintf("C (%d x %d)", Q, R), X1.label, X2.label,
           mar = c(5, 4, 3, 1))
  plot_mat(X2, sprintf("X2 (%d x %d)", R, P2), X2.label, Y2.label,
           mar = c(5, 4, 3, 1))

  invisible(NULL)
}

#' @title Predict method for nmfae objects
#' @description
#' \code{predict.nmfae} computes fitted or predicted values from a three-layer NMF model.
#' Without \code{newY2}, returns the in-sample fitted values \eqn{X_1 \Theta X_2 Y_2}.
#' With \code{newY2}, computes out-of-sample predictions \eqn{X_1 \Theta X_2 \cdot \mathrm{newY2}}.
#'
#' When \code{type = "class"}, each column is classified to the row with the
#' maximum predicted value (useful when \eqn{Y_1} is a one-hot class matrix
#' from \code{\link{nmfkc.class}}).
#'
#' If \code{Y1} (actual values) is provided, it is stored as an attribute so that
#' \code{plot.predict.nmfae} can produce an observed-vs-predicted scatter plot
#' (for \code{type = "response"}) or a confusion matrix heatmap
#' (for \code{type = "class"}).
#'
#' @param object An object of class \code{"nmfae"} returned by \code{\link{nmfae}}.
#' @param newY2 Optional new input matrix (P2 x M) for prediction.
#'   If \code{NULL}, returns in-sample fitted values.
#' @param Y1 Optional actual output matrix for comparison plotting.
#' @param type Character. \code{"response"} (default) returns the predicted matrix.
#'   \code{"class"} returns a factor of predicted class labels (row with max value).
#' @param ... Not used.
#'
#' @return For \code{type = "response"}: a matrix of class \code{"predict.nmfae"}.
#'   For \code{type = "class"}: a factor of class \code{"predict.nmfae"} with
#'   predicted class labels. If \code{Y1} was provided, actual classes are stored
#'   in \code{attr(result, "actual")}.
#' @seealso \code{\link{nmfae}}, \code{\link{plot.predict.nmfae}},
#'   \code{\link{nmfkc.class}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(runif(20), nrow = 4)
#' res <- nmfae(Y, rank = 2)
#' pred <- predict(res)
#' }
#' @export
predict.nmfae <- function(object, newY2 = NULL, Y1 = NULL,
                          type = c("response", "class"), ...) {
  type <- match.arg(type)
  if (is.null(newY2)) pred_mat <- object$Y1hat
  else pred_mat <- object$X1 %*% object$C %*% object$X2 %*% as.matrix(newY2)

  if (type == "response") {
    result <- pred_mat
    if (!is.null(Y1)) attr(result, "Y1") <- Y1
    class(result) <- c("predict.nmfae", class(result))
  } else {
    # Class prediction: row with max value per column
    labels <- rownames(pred_mat)
    if (is.null(labels)) labels <- paste0("C", seq_len(nrow(pred_mat)))
    pred_class <- factor(labels[apply(pred_mat, 2, which.max)], levels = labels)
    result <- pred_class
    if (!is.null(Y1)) {
      Y1 <- as.matrix(Y1)
      act_labels <- rownames(Y1)
      if (is.null(act_labels)) act_labels <- labels
      actual <- factor(act_labels[apply(Y1, 2, which.max)], levels = act_labels)
      attr(result, "actual") <- actual
    }
    class(result) <- c("predict.nmfae", class(result))
    attr(result, "type") <- "class"
  }
  result
}

#' @title Plot method for predict.nmfae objects
#' @description
#' For \code{type = "response"}: if actual values \eqn{Y_1} were stored,
#' displays an observed-vs-predicted scatter plot with \eqn{R^2} in the title.
#' Otherwise, displays the predicted matrix as a heatmap.
#'
#' For \code{type = "class"}: if actual classes were stored, displays a
#' confusion matrix heatmap with accuracy (ACC) in the title.
#'
#' @param x An object of class \code{"predict.nmfae"} returned by \code{\link{predict.nmfae}}.
#' @param ... Additional graphical parameters passed to \code{plot} or \code{image}.
#'
#' @return Invisible \code{NULL}. Called for its side effect (plot).
#' @seealso \code{\link{predict.nmfae}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(runif(20), nrow = 4)
#' res <- nmfae(Y, rank = 2)
#' pred <- predict(res)
#' plot(pred)
#' }
#' @export
plot.predict.nmfae <- function(x, ...) {
  extra_args <- list(...)
  pred_type <- attr(x, "type")

  if (!is.null(pred_type) && pred_type == "class") {
    # --- Confusion matrix heatmap ---
    actual <- attr(x, "actual")
    if (is.null(actual)) {
      message("No actual classes stored. Use Y1 argument in predict().")
      return(invisible(NULL))
    }
    pred_class <- factor(x, levels = levels(actual))
    cm <- table(Actual = actual, Predicted = pred_class)
    acc <- sum(diag(cm)) / sum(cm)

    # Heatmap of confusion matrix
    mat <- as.matrix(cm)
    nr <- nrow(mat); nc <- ncol(mat)
    pal <- grDevices::colorRampPalette(c("white", "orange", "red"))(64)
    opar <- graphics::par(mar = c(5, 5, 3, 1))
    on.exit(graphics::par(opar))
    graphics::image(1:nc, 1:nr, t(mat[nr:1, , drop = FALSE]),
                    col = pal, axes = FALSE, xlab = "Predicted", ylab = "Actual")
    graphics::box()
    labs <- rownames(mat)
    graphics::axis(1, at = 1:nc, labels = colnames(mat), las = 2, cex.axis = 0.8)
    graphics::axis(2, at = 1:nr, labels = rev(labs), las = 1, cex.axis = 0.8)
    # Cell labels
    for (i in 1:nr) {
      for (j in 1:nc) {
        val <- mat[i, j]
        if (val > 0) {
          yi <- nr - i + 1
          graphics::text(j, yi, val, cex = 0.7,
                         col = ifelse(val > stats::median(mat), "white", "black"))
        }
      }
    }
    if (is.null(extra_args$main))
      graphics::title(main = sprintf("ACC = %.1f%%", acc * 100))
    else
      graphics::title(main = extra_args$main)

  } else {
    # --- Response type ---
    Y1 <- attr(x, "Y1")

    if (!is.null(Y1)) {
      # Scatter plot: observed vs predicted
      obs <- as.numeric(as.matrix(Y1))
      pred <- as.numeric(unclass(x))
      r2 <- stats::cor(obs, pred)^2
      args <- list(x = obs, y = pred)
      if (is.null(extra_args$main))
        args$main <- sprintf("R2 = %.3f", r2)
      if (is.null(extra_args$xlab)) args$xlab <- "Observed"
      if (is.null(extra_args$ylab)) args$ylab <- "Predicted"
      if (is.null(extra_args$pch))  args$pch <- 16
      if (is.null(extra_args$cex))  args$cex <- 0.8
      if (is.null(extra_args$col))  args$col <- grDevices::rgb(0, 0, 0, 0.4)
      all_args <- c(args, extra_args)
      do.call("plot", all_args)
      graphics::abline(0, 1, col = "red")
    } else {
      # Heatmap of predicted matrix
      mat <- unclass(x)
      nr <- nrow(mat); nc <- ncol(mat)
      pal <- grDevices::colorRampPalette(c("white", "orange", "red"))(64)
      args <- list(x = 1:nc, y = 1:nr, z = t(mat[nr:1, , drop = FALSE]),
                   col = pal, axes = FALSE, xlab = "Sample", ylab = "Variable")
      if (is.null(extra_args$main)) args$main <- "Predicted Y1"
      all_args <- c(args, extra_args)
      do.call("image", all_args)
      graphics::box()
    }
  }
  invisible(NULL)
}

#' @title Element-wise Cross-Validation for nmfae (Wold's CV)
#' @description
#' \code{nmfae.ecv} performs k-fold element-wise cross-validation by randomly
#' holding out individual elements of \eqn{Y_1}, assigning them a weight of 0
#' via \code{Y1.weights}, and evaluating the reconstruction error on those
#' held-out elements.
#'
#' This method (also known as Wold's CV) is suitable for determining the optimal
#' rank pair \eqn{(Q, R)} in three-layer NMF. Both \code{rank} and \code{rank.encoder} accept
#' vector inputs. When \code{rank.encoder = NULL} (default), \code{rank.encoder} is set equal to \code{rank}
#' and pairs are evaluated element-wise (i.e., \eqn{(Q_1, R_1), (Q_2, R_2), \dots}).
#' When \code{rank.encoder} is explicitly specified, all combinations of \code{rank} and \code{rank.encoder}
#' are evaluated via \code{expand.grid}.
#'
#' @param Y1 Output matrix \eqn{Y_1} (P1 x N).
#' @param Y2 Input matrix \eqn{Y_2} (P2 x N). Default is \code{Y1}.
#' @param rank Integer vector of decoder ranks to evaluate. Default is \code{1:2}.
#' @param rank.encoder Integer vector of encoder ranks to evaluate. Default is \code{NULL},
#'   which sets \code{rank.encoder = rank} and evaluates element-wise pairs.
#'   When explicitly specified, all combinations with \code{rank} are evaluated.
#' @param ... Additional arguments passed to \code{\link{nmfae}} (e.g., \code{epsilon}, \code{maxit}).
#'   Also accepts: \code{nfolds} (number of folds, default 5; \code{div} also accepted),
#'   \code{seed} (integer seed, default 123).
#'   For backward compatibility, \code{Q} and \code{R} are accepted as aliases for
#'   \code{rank} and \code{rank.encoder}.
#'
#' @return A list with components:
#' \item{objfunc}{Named numeric vector of mean MSE for each (Q, R) pair.}
#' \item{sigma}{Named numeric vector of RMSE (square root of MSE) for each pair.}
#' \item{objfunc.fold}{Named list of per-fold MSE vectors for each pair.}
#' \item{folds}{List of length \code{div} containing the held-out element indices for each fold.}
#' \item{QR}{Data frame with columns \code{Q} and \code{R} listing the evaluated pairs.}
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @seealso \code{\link{nmfae}}, \code{\link{nmfkc.ecv}}
#' @export
#' @examples
#' Y <- t(iris[1:30, 1:4])
#' # Default: rank.encoder=NULL -> paired rank=rank.encoder
#' res <- nmfae.ecv(Y, rank = 1:3, nfolds = 3, maxit = 500)
#' res$sigma
#' # Explicit rank.encoder: full grid
#' res2 <- nmfae.ecv(Y, rank = 1:3, rank.encoder = 1:3, nfolds = 3, maxit = 500)
#' res2$sigma
#'
nmfae.ecv <- function(Y1, Y2 = Y1, rank = 1:2, rank.encoder = NULL, ...) {
  extra_ecv <- list(...)
  if (!is.null(extra_ecv$Q)) rank <- extra_ecv$Q
  if (!is.null(extra_ecv$R)) rank.encoder <- extra_ecv$R
  nfolds <- if (!is.null(extra_ecv$nfolds)) extra_ecv$nfolds else if (!is.null(extra_ecv$div)) extra_ecv$div else 5
  seed   <- if (!is.null(extra_ecv$seed))   extra_ecv$seed   else 123
  Q <- rank; R <- rank.encoder
  div <- nfolds

  Y1 <- as.matrix(Y1); Y2 <- as.matrix(Y2)
  P1 <- nrow(Y1); N <- ncol(Y1)

  # R=NULL -> paired with Q; R specified -> full grid
  if (is.null(R)) {
    QR <- data.frame(Q = Q, R = Q)
  } else {
    QR <- expand.grid(Q = Q, R = R)
  }
  num_pairs <- nrow(QR)

  # Create folds (element-wise on Y1; shared helper)
  folds <- .ecv.make.folds(Y1, div, seed)

  # Prepare result storage
  pair_labels <- sprintf("Q=%d,R=%d", QR$Q, QR$R)
  has_na <- any(is.na(Y1))

  message(sprintf("Element-wise CV: %d (Q,R) pairs, %d-fold, %d tasks...",
                  num_pairs, div, num_pairs * div))

  extra_args <- list(...)

  # Model-specific worker: mask fold k, refit at pair i, held-out loss
  run_one <- function(i, k) {
    test_idx <- folds[[k]]
    weights_train <- matrix(1, nrow = P1, ncol = N)
    if (has_na) weights_train[is.na(Y1)] <- 0
    weights_train[test_idx] <- 0
    fit <- suppressMessages(
      do.call(nmfae, c(list(Y1 = Y1, Y2 = Y2, Q = QR$Q[i], R = QR$R[i],
                            Y1.weights = weights_train), extra_args))
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
  class(result) <- "nmfae.ecv"
  return(result)
}


#' @title Rank selection for nmfae (paired rank, concise diagnostics)
#' @description
#' Fits \code{\link{nmfae}} with a \strong{paired} decoder/encoder rank
#' (\eqn{Q = R}) across a range of ranks and reports \code{r.squared},
#' the effective rank (of the latent encoding \eqn{H}), and the
#' element-wise CV error \code{sigma.ecv}, with the same concise plot as
#' \code{\link{nmfkc.rank}}.  For a full \eqn{(Q, R)} grid use
#' \code{\link{nmfae.ecv}} with \code{rank.encoder} and its heatmap.
#' @param Y1 Endogenous matrix (\eqn{P_1 \times N}).
#' @param Y2 Exogenous matrix; defaults to \code{Y1} (autoencoder).
#' @param rank Integer vector of (paired) ranks to evaluate.
#' @param detail \code{"full"} (default) also runs element-wise CV
#'   (\code{sigma.ecv}); \code{"fast"} skips it (plots r.squared and
#'   eff.rank only, and recommends the R-squared elbow).
#' @param plot Logical; draw the diagnostics plot (default \code{TRUE}).
#' @param ... Passed on to \code{\link{nmfae}} and \code{\link{nmfae.ecv}}
#'   (e.g.\ \code{maxit}, \code{nfolds}, \code{seed}).
#' @return A list with \code{rank.best} and \code{criteria}
#'   (\code{rank}, \code{effective.rank}, \code{effective.rank.ratio},
#'   \code{r.squared}, \code{sigma.ecv}).
#' @seealso \code{\link{nmfae}}, \code{\link{nmfae.ecv}},
#'   \code{\link{nmfkc.rank}}
#' @export
nmfae.rank <- function(Y1, Y2 = Y1, rank = 1:5, detail = c("full", "fast"),
                       plot = TRUE, ...) {
  detail <- match.arg(detail)
  Y1 <- as.matrix(Y1); Y2 <- as.matrix(Y2)
  rs <- numeric(length(rank)); er <- numeric(length(rank))
  for (i in seq_along(rank)) {
    f <- suppressMessages(nmfae(Y1, Y2, rank = rank[i],
                                rank.encoder = rank[i],
                                print.trace = FALSE, ...))
    rs[i] <- f$r.squared
    er[i] <- .effective.rank(f$H)
  }
  ecv <- if (detail == "full")
    suppressMessages(nmfae.ecv(Y1, Y2, rank = rank, ...))$sigma
    else rep(NA_real_, length(rank))
  criteria <- data.frame(rank = rank, effective.rank = er,
                         effective.rank.ratio = er / rank,
                         r.squared = rs, sigma.ecv = as.numeric(ecv))
  .rank.finish(criteria, plot = plot, main = "nmfae rank selection (paired Q=R)")
}

#' @title Plot method for nmfae.ecv objects
#' @description
#' Visualizes element-wise cross-validation results.
#' When \code{rank.encoder} was \code{NULL} (paired), a line plot of sigma vs rank is drawn.
#' When \code{rank.encoder} was explicitly specified (grid), a heatmap of sigma over the (rank, rank.encoder) grid is drawn.
#'
#' @param x An object of class \code{"nmfae.ecv"} returned by \code{\link{nmfae.ecv}}.
#' @param ... Additional graphical parameters (currently unused).
#'
#' @return Invisible \code{NULL}. Called for its side effect of producing a plot.
#' @seealso \code{\link{nmfae.ecv}}
#' @export
plot.nmfae.ecv <- function(x, ...) {
  Qs <- sort(unique(x$QR$Q))
  Rs <- sort(unique(x$QR$R))
  use_line <- x$paired || length(Qs) == 1 || length(Rs) == 1

  if (use_line) {
    # --- Line plot ---
    if (x$paired) {
      xvals <- x$QR$Q; xlab_str <- "Q = R"
    } else if (length(Qs) == 1) {
      xvals <- x$QR$R; xlab_str <- paste0("R (Q=", Qs, ")")
    } else {
      xvals <- x$QR$Q; xlab_str <- paste0("Q (R=", Rs, ")")
    }
    plot(xvals, x$sigma, type = "b", pch = 16,
         xlab = xlab_str, ylab = "sigma",
         main = "Element-wise CV",
         xaxt = "n")
    graphics::axis(1, at = xvals)
    best <- which.min(x$sigma)
    graphics::points(xvals[best], x$sigma[best], pch = 16, col = "red", cex = 1.5)
  } else {
    # --- Heatmap: sigma over (Q, R) grid ---
    Qs <- sort(unique(x$QR$Q))
    Rs <- sort(unique(x$QR$R))
    sigma_mat <- matrix(NA, length(Qs), length(Rs))
    for (i in 1:nrow(x$QR)) {
      qi <- which(Qs == x$QR$Q[i])
      ri <- which(Rs == x$QR$R[i])
      sigma_mat[qi, ri] <- x$sigma[i]
    }

    ncolors <- 64
    col_palette <- grDevices::hcl.colors(ncolors, "YlOrRd", rev = TRUE)

    opar <- graphics::par(mar = c(5, 5, 3, 8))
    on.exit(graphics::par(opar))

    graphics::image(Qs, Rs, sigma_mat, col = col_palette,
          xlab = "Q (decoder rank)", ylab = "R (encoder rank)",
          main = "Element-wise CV", axes = FALSE)
    graphics::axis(1, at = Qs)
    graphics::axis(2, at = Rs)
    graphics::box()

    # Cell values
    for (qi in seq_along(Qs)) {
      for (ri in seq_along(Rs)) {
        val <- sigma_mat[qi, ri]
        txt_col <- "black"
        graphics::text(Qs[qi], Rs[ri], sprintf("%.2f", val), cex = 0.65, col = txt_col)
      }
    }

    # Color bar — use cell step as unit for positioning
    old_par <- graphics::par(xpd = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    q_step <- if (length(Qs) > 1) diff(Qs)[1] else 1
    r_step <- if (length(Rs) > 1) diff(Rs)[1] else 1
    bar_x <- max(Qs) + q_step * 0.7
    bar_w <- q_step * 0.3
    bar_y <- seq(min(Rs), max(Rs), length.out = ncolors + 1)
    for (i in 1:ncolors) {
      graphics::rect(bar_x, bar_y[i], bar_x + bar_w, bar_y[i + 1],
           col = col_palette[i], border = NA)
    }
    graphics::rect(bar_x, min(Rs), bar_x + bar_w, max(Rs))
    sigma_range <- range(sigma_mat)
    ticks <- pretty(sigma_range, 5)
    ticks <- ticks[ticks >= sigma_range[1] & ticks <= sigma_range[2]]
    tick_y <- min(Rs) + (ticks - sigma_range[1]) / diff(sigma_range) * diff(range(Rs))
    graphics::text(bar_x + bar_w + q_step * 0.15, tick_y,
         sprintf("%.1f", ticks), cex = 0.7, adj = 0)
    graphics::text(bar_x + bar_w / 2, max(Rs) + r_step * 0.5,
         "sigma", cex = 0.8)
  }
}

#' @title Sample-wise k-fold Cross-Validation for nmfae
#' @description
#' \code{nmfae.cv} performs k-fold cross-validation by splitting columns (samples)
#' of \eqn{Y_1} and \eqn{Y_2} into \code{div} folds. For each fold, the model
#' \eqn{Y_1 \approx X_1 \Theta X_2 Y_2} is fitted on the training samples and
#' predictive performance is evaluated on the held-out samples.
#'
#' When \code{Y2} is a kernel matrix created by \code{\link{nmfkc.kernel}}
#' (detected via attributes), the symmetric kernel splitting convention is used:
#' \code{Y2[train, train]} for training and \code{Y2[train, test]} for prediction.
#'
#' @param Y1 Output matrix \eqn{Y_1} (P1 x N). Non-negative.
#' @param Y2 Input matrix \eqn{Y_2} (P2 x N), or a kernel matrix (N x N).
#'   Default is \code{Y1} (autoencoder).
#' @param rank Integer. Rank of the decoder basis. Default is 2.
#' @param rank.encoder Integer. Rank of the encoder basis. Default is \code{rank}.
#' @param ... Additional arguments passed to \code{\link{nmfae}}
#'   (e.g., \code{epsilon}, \code{maxit}, \code{Y1.weights}).
#'   Also accepts: \code{nfolds} (number of folds, default 5; \code{div} also accepted),
#'   \code{seed} (integer seed, default 123), \code{shuffle} (logical, default \code{TRUE}).
#'   For backward compatibility, \code{Q}, \code{R} are accepted as aliases for
#'   \code{rank}, \code{rank.encoder}.
#'
#' @return A list with components:
#' \item{objfunc}{Mean squared error per valid element over all folds.}
#' \item{sigma}{Residual standard error (RMSE), same scale as \eqn{Y_1}.}
#' \item{objfunc.block}{Per-fold squared error totals.}
#' \item{block}{Integer vector of fold assignments (1, ..., \code{div}) for each column.}
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @seealso \code{\link{nmfae}}, \code{\link{nmfae.ecv}}, \code{\link{nmfae.kernel.beta.cv}},
#'   \code{\link{nmfkc.cv}}
#' @export
#' @examples
#' Y <- t(iris[1:30, 1:4])
#' res <- nmfae.cv(Y, rank = 2, rank.encoder = 2, nfolds = 5, maxit = 500)
#' res$sigma
#'
nmfae.cv <- function(Y1, Y2 = Y1, rank = 2, rank.encoder = rank, ...) {
  extra_cv <- list(...)
  if (!is.null(extra_cv$Q)) rank <- extra_cv$Q
  if (!is.null(extra_cv$R)) rank.encoder <- extra_cv$R
  nfolds  <- if (!is.null(extra_cv$nfolds))  extra_cv$nfolds  else if (!is.null(extra_cv$div)) extra_cv$div else 5
  seed    <- if (!is.null(extra_cv$seed))    extra_cv$seed    else 123
  shuffle <- if (!is.null(extra_cv$shuffle)) extra_cv$shuffle else TRUE
  Q <- rank; R <- rank.encoder
  div <- nfolds

  extra_args <- list(...)

  Y1 <- as.matrix(Y1); storage.mode(Y1) <- "double"
  Y2 <- as.matrix(Y2); storage.mode(Y2) <- "double"
  P1 <- nrow(Y1); N <- ncol(Y1)

  # --- Y1.weights preparation (cf. nmfkc.cv) ---
  Y1.weights <- extra_args$Y1.weights
  if (!is.null(Y1.weights) && is.vector(Y1.weights)) {
    if (length(Y1.weights) == N) {
      Y1.weights <- matrix(Y1.weights, nrow = P1, ncol = N, byrow = TRUE)
    } else if (length(Y1.weights) == 1) {
      Y1.weights <- matrix(Y1.weights, nrow = P1, ncol = N)
    } else {
      stop("Length of Y1.weights vector must match ncol(Y1) (or be 1).")
    }
  }
  if (is.null(Y1.weights)) {
    if (any(is.na(Y1))) {
      Y1.weights <- matrix(1, nrow = P1, ncol = N)
      Y1.weights[is.na(Y1)] <- 0
      Y1[is.na(Y1)] <- 0
    } else {
      Y1.weights <- matrix(1, nrow = P1, ncol = N)
    }
  } else {
    if (!is.matrix(Y1.weights)) Y1.weights <- as.matrix(Y1.weights)
    Y1.weights[is.na(Y1.weights)] <- 0
    Y1[is.na(Y1)] <- 0
  }

  # --- Detect kernel matrix (cf. nmfkc.cv) ---
  A.function <- attr(Y2, "function.name")
  is_kernel_matrix <- !is.null(A.function) && A.function == "nmfkc.kernel"
  is_symmetric <- isSymmetric(Y2, tol = .Machine$double.eps)

  # --- Create folds (sample-wise, same logic as nmfkc.cv) ---
  remainder <- N %% div
  division  <- N %/% div
  block <- integer(N)

  if (shuffle) {
    set.seed(seed)
    perm_index <- sample(1:N, N, replace = FALSE)
  } else {
    perm_index <- 1:N
  }

  processed <- 0
  for (i in 1:(div - 1)) {
    chunk <- division + ifelse(i <= remainder, 1, 0)
    idx <- perm_index[(processed + 1):(processed + chunk)]
    block[idx] <- i
    processed <- processed + chunk
  }
  block[perm_index[(processed + 1):N]] <- div

  # --- Prepare nmfae pass-through args ---
  nmfae_extra <- extra_args
  nmfae_extra$Y1.weights <- NULL
  nmfae_extra$print.trace <- NULL
  nmfae_extra$seed <- NULL

  objfunc.block <- numeric(div)
  total_valid <- 0

  for (j in 1:div) {
    train <- (block != j)
    test  <- (block == j)

    # Y1 split
    Y1_train <- Y1[, train, drop = FALSE]
    Y1_test  <- Y1[, test,  drop = FALSE]
    W_train  <- Y1.weights[, train, drop = FALSE]
    W_test   <- Y1.weights[, test,  drop = FALSE]

    # Y2 split
    if (is_symmetric && is_kernel_matrix) {
      Y2_train <- Y2[train, train, drop = FALSE]
      Y2_test  <- Y2[train, test,  drop = FALSE]
    } else {
      Y2_train <- Y2[, train, drop = FALSE]
      Y2_test  <- Y2[, test,  drop = FALSE]
    }

    # Fit on training set
    nmfae_args <- c(
      list(Y1 = Y1_train, Y2 = Y2_train, Q = Q, R = R,
           Y1.weights = W_train, print.trace = FALSE),
      nmfae_extra
    )
    res_j <- suppressMessages(do.call("nmfae", nmfae_args))

    # Predict on test set
    Y1hat_test <- res_j$X1 %*% res_j$C %*% res_j$X2 %*% Y2_test

    # Evaluate weighted error (lm-style: sum(W * resid^2))
    objfunc.block[j] <- sum(W_test * (Y1_test - Y1hat_test)^2)
    total_valid <- total_valid + sum(W_test > 0)
  }

  objfunc <- sum(objfunc.block) / max(total_valid, 1)
  sigma <- sqrt(objfunc)

  result <- list(objfunc = objfunc, sigma = sigma,
                 objfunc.block = objfunc.block, block = block)
  class(result) <- "nmfae.cv"
  result
}

#' @title Plot method for nmfae.cv objects
#' @description
#' Displays a bar chart of per-fold cross-validation errors from
#' \code{\link{nmfae.cv}}. The overall RMSE (sigma) is shown in the title.
#'
#' @param x An object of class \code{"nmfae.cv"} returned by \code{\link{nmfae.cv}}.
#' @param ... Additional graphical parameters passed to \code{barplot}.
#'
#' @return Invisible \code{NULL}. Called for its side effect (plot).
#' @seealso \code{\link{nmfae.cv}}
#' @export
plot.nmfae.cv <- function(x, ...) {
  extra_args <- list(...)
  args <- list(height = x$objfunc.block)
  if (is.null(extra_args$main))
    args$main <- sprintf("sigma = %.4f", x$sigma)
  if (is.null(extra_args$xlab)) args$xlab <- "Fold"
  if (is.null(extra_args$ylab)) args$ylab <- "Error"
  if (is.null(extra_args$col))  args$col <- "steelblue"
  if (is.null(extra_args$names.arg))
    args$names.arg <- seq_along(x$objfunc.block)
  all_args <- c(args, extra_args)
  do.call("barplot", all_args)
  invisible(NULL)
}

#' @title Optimize kernel beta for nmfae by cross-validation
#' @description
#' \code{nmfae.kernel.beta.cv} selects the optimal \code{beta} parameter of the
#' kernel function by evaluating \code{\link{nmfae.cv}} for each candidate value.
#' The kernel matrix \eqn{A = K(U, V; \beta)} replaces \eqn{Y_2} in the three-layer
#' NMF model.
#'
#' When \code{beta = NULL}, candidate values are automatically generated via
#' \code{\link{nmfkc.kernel.beta.nearest.med}}.
#'
#' @param Y1 Output matrix \eqn{Y_1} (P1 x N). Non-negative.
#' @param rank Integer. Rank of the decoder basis. Default is 2.
#' @param rank.encoder Integer. Rank of the encoder basis. Default is \code{rank}.
#' @param U Covariate matrix \eqn{U} (K x M). Rows are features, columns are samples
#'   (or knot points for non-symmetric kernels).
#' @param V Covariate matrix \eqn{V} (K x N). If \code{NULL} (default), \code{V = U}
#'   and a symmetric kernel is used.
#' @param beta Numeric vector of candidate beta values. If \code{NULL}, automatically
#'   determined via \code{\link{nmfkc.kernel.beta.nearest.med}}.
#' @param plot Logical. If \code{TRUE} (default), plots the objective function curve.
#' @param ... Additional arguments. Kernel-specific args (\code{kernel}, \code{degree})
#'   are passed to \code{\link{nmfkc.kernel}}; all others
#'   (\code{div}, \code{seed}, \code{shuffle}, \code{epsilon}, \code{maxit}, etc.)
#'   are passed to \code{\link{nmfae.cv}}.
#'   For backward compatibility, \code{Q} and \code{R} are accepted as aliases for
#'   \code{rank} and \code{rank.encoder}.
#'
#' @return A list with components:
#' \item{beta}{The beta value that minimizes the cross-validation objective.}
#' \item{objfunc}{Named numeric vector of objective function values for each candidate beta.}
#'
#' @seealso \code{\link{nmfae.cv}}, \code{\link{nmfkc.kernel}},
#'   \code{\link{nmfkc.kernel.beta.cv}}
#' @export
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' U <- matrix(cars$speed, nrow = 1)
#' res <- nmfae.kernel.beta.cv(Y, rank = 1, rank.encoder = 1, U = U,
#'                              beta = c(0.01, 0.02, 0.05), nfolds = 5)
#' res$beta
#'
nmfae.kernel.beta.cv <- function(Y1, rank = 2, rank.encoder = rank, U, V = NULL,
                                  beta = NULL, plot = TRUE, ...) {

  extra_args <- list(...)
  # backward compat: Q -> rank, R -> rank.encoder
  if (!is.null(extra_args$Q)) rank <- extra_args$Q
  if (!is.null(extra_args$R)) rank.encoder <- extra_args$R
  extra_args <- extra_args[!names(extra_args) %in% c("Q", "R")]

  # Separate kernel-specific args from cv/nmfae args
  kernel_only <- c("kernel", "degree")
  kernel_args <- extra_args[names(extra_args) %in% kernel_only]
  cv_args     <- extra_args[!names(extra_args) %in% kernel_only]

  # Auto-generate beta candidates if not provided
  if (is.null(beta)) {
    if (is.null(V)) V <- U
    result.beta <- nmfkc.kernel.beta.nearest.med(V)
    beta <- result.beta$beta_candidates
    if (is.null(beta) || length(beta) == 0)
      stop("Failed to determine beta candidates from nearest-neighbor median.")
  }

  objfuncs <- numeric(length(beta))
  for (i in seq_along(beta)) {
    start.time <- Sys.time()
    message(paste0("beta=", beta[i], "..."), appendLF = FALSE)

    kernel_call <- c(list(U = U, V = V, beta = beta[i]), kernel_args)
    A <- do.call("nmfkc.kernel", kernel_call)

    cv_call <- c(list(Y1 = Y1, Y2 = A, rank = rank, rank.encoder = rank.encoder), cv_args)
    result <- do.call("nmfae.cv", cv_call)

    objfuncs[i] <- result$objfunc

    end.time <- Sys.time()
    diff.time <- difftime(end.time, start.time, units = "sec")
    diff.time.st <- ifelse(diff.time <= 180,
                           paste0(round(diff.time, 1), "sec"),
                           paste0(round(diff.time / 60, 1), "min"))
    message(diff.time.st)
  }

  i0 <- which.min(objfuncs)
  beta.best <- beta[i0]

  names(objfuncs) <- beta
  result <- list(beta = beta.best, objfunc = objfuncs)
  class(result) <- "nmfae.kernel.beta.cv"
  if (plot) plot(result)
  result
}

#' @title Plot method for nmfae.kernel.beta.cv objects
#' @description
#' Displays the cross-validation objective function across candidate
#' \code{beta} values (log scale). The optimal beta is highlighted in red.
#'
#' @param x An object of class \code{"nmfae.kernel.beta.cv"} returned by
#'   \code{\link{nmfae.kernel.beta.cv}}.
#' @param ... Additional graphical parameters passed to \code{plot}.
#'
#' @return Invisible \code{NULL}. Called for its side effect (plot).
#' @seealso \code{\link{nmfae.kernel.beta.cv}}
#' @export
plot.nmfae.kernel.beta.cv <- function(x, ...) {
  beta <- as.numeric(names(x$objfunc))
  objfuncs <- x$objfunc
  i0 <- which.min(objfuncs)

  extra_args <- list(...)
  args <- list(x = beta, y = objfuncs, type = "l", col = 2, log = "x")
  if (is.null(extra_args$xlab)) args$xlab <- "beta"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  if (is.null(extra_args$main))
    args$main <- sprintf("Best beta = %g", x$beta)
  all_args <- c(args, extra_args)
  do.call("plot", all_args)
  graphics::points(beta[i0], objfuncs[i0], cex = 3, col = 2)
  graphics::text(beta, objfuncs,
                 format(beta, scientific = TRUE, digits = 5))
  invisible(NULL)
}

#' @title DOT graph visualization for nmfae objects
#' @description
#' \code{nmfae.DOT} generates a DOT language string for visualizing the structure
#' of a three-layer NMF model. Two graph types are supported:
#' \code{"XCX"} shows encoder factors, \eqn{\Theta}, and decoder factors;
#' \code{"YXCXY"} shows the full structure from \eqn{Y_2} through \eqn{X_2}, \eqn{\Theta},
#' \eqn{X_1} to \eqn{Y_1}.
#'
#' Edge widths are proportional to matrix element values, and edges below
#' \code{threshold} are omitted for clarity.
#'
#' @param result An object of class \code{"nmfae"} returned by \code{\link{nmfae}}.
#' @param type Character. Graph type: \code{"XCX"} (default) or \code{"YXCXY"}.
#' @param threshold Numeric. Edges with values below this are omitted. Default is 0.01.
#' @param sig.level Numeric or \code{NULL}. Significance level for filtering C edges
#'   when inference results are available. Only edges with p-value below \code{sig.level}
#'   are shown, with significance stars. Set to \code{NULL} to disable. Default is 0.1.
#' @param rankdir Character. Graph direction for DOT layout. Default is \code{"LR"} (left to right).
#' @param fill Logical. If \code{TRUE}, nodes are filled with color. Default is \code{TRUE}.
#' @param weight_scale Numeric. Base scale factor for edge widths. Default is 5.
#' @param weight_scale_x1 Numeric. Scale factor for \eqn{X_1} edges.
#' @param weight_scale_theta Numeric. Scale factor for \eqn{\Theta} edges.
#' @param weight_scale_x2 Numeric. Scale factor for \eqn{X_2} edges.
#' @param Y1.label Character vector of output variable labels.
#' @param X1.label Character vector of decoder basis labels.
#' @param X2.label Character vector of encoder basis labels.
#' @param Y2.label Character vector of input variable labels.
#' @param Y1.title Character. Title for output node group. Default is \code{"Output (Y1)"}.
#' @param X1.title Character. Title for decoder node group. Default is \code{"Decoder (X1)"}.
#' @param X2.title Character. Title for encoder node group. Default is \code{"Encoder (X2)"}.
#' @param Y2.title Character. Title for input node group. Default is \code{"Input (Y2)"}.
#' @param hide.isolated Logical. If \code{TRUE} (default), Y1 and Y2 nodes that have no
#'   edges at or above \code{threshold} are excluded from the graph. Only
#'   applies when \code{type = "YXCXY"}.
#'
#' @return A character string containing the DOT graph specification.
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @seealso \code{\link{nmfae}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(runif(20), nrow = 4)
#' res <- nmfae(Y, rank = 2)
#' dot <- nmfae.DOT(res)
#' }
#' @export
nmfae.DOT <- function(result,
                      type = c("XCX", "YXCXY"),
                      threshold = 0.01,
                      sig.level = 0.1,
                      rankdir = "LR",
                      fill = TRUE,
                      weight_scale = 5,
                      weight_scale_x1 = weight_scale,
                      weight_scale_theta = weight_scale,
                      weight_scale_x2 = weight_scale,
                      Y1.label = NULL, X1.label = NULL,
                      X2.label = NULL, Y2.label = NULL,
                      Y1.title = "Output (Y1)",
                      X1.title = "Decoder (X1)",
                      X2.title = "Encoder (X2)",
                      Y2.title = "Input (Y2)",
                      hide.isolated = TRUE) {

  type <- match.arg(type)

  X1 <- result$X1        # P1 x Q
  C_mat <- result$C   # Q x R
  X2 <- result$X2        # R x P2

  P1 <- nrow(X1); Q <- ncol(X1)
  R <- nrow(X2); P2 <- ncol(X2)

  # --- Labels ---
  Y1_labels <- if (!is.null(Y1.label)) Y1.label else rownames(X1)
  if (is.null(Y1_labels)) Y1_labels <- paste0("Y1_", seq_len(P1))
  X1_labels <- if (!is.null(X1.label)) X1.label else colnames(X1)
  if (is.null(X1_labels)) X1_labels <- paste0("D", seq_len(Q))
  X2_labels <- if (!is.null(X2.label)) X2.label else rownames(X2)
  if (is.null(X2_labels)) X2_labels <- paste0("E", seq_len(R))
  Y2_labels <- if (!is.null(Y2.label)) Y2.label else colnames(X2)
  if (is.null(Y2_labels)) Y2_labels <- paste0("Y2_", seq_len(P2))

  # --- Node IDs ---
  sanitize <- function(s) gsub("[^[:alnum:]_.]", "_", s, perl = TRUE)
  Y1_ids <- sanitize(paste0("Y1_", seq_len(P1)))
  X1_ids <- sanitize(paste0("X1_", seq_len(Q)))
  X2_ids <- sanitize(paste0("X2_", seq_len(R)))
  Y2_ids <- sanitize(paste0("Y2_", seq_len(P2)))

  # --- Filter isolated Y1/Y2 nodes (hide.isolated) ---
  idx_Y1 <- seq_len(P1)
  idx_Y2 <- seq_len(P2)
  if (isTRUE(hide.isolated) && type == "YXCXY") {
    # Y1 is connected if any X1[i, ] >= threshold
    used_Y1 <- apply(X1, 1L, function(row) any(row >= threshold, na.rm = TRUE))
    idx_Y1 <- which(used_Y1)
    # Y2 is connected if any X2[, j] >= threshold
    used_Y2 <- apply(X2, 2L, function(col) any(col >= threshold, na.rm = TRUE))
    idx_Y2 <- which(used_Y2)
  }

  # --- Helpers ---
  pw <- function(value, max_value, ws) {
    if (!is.finite(max_value) || max_value <= 0) return(0.5)
    if (!is.finite(value) || value <= 0) return(0.5)
    max(0.5, value * ws / max_value)
  }
  dg <- max(0L, nchar(sub("^0\\.", "", format(threshold, scientific = FALSE, trim = TRUE))))
  fmtc <- function(v) sprintf(paste0("%.", dg, "f"), v)

  # --- Cluster builder ---
  make_cluster <- function(id, title, ids, labels, shape, color) {
    style_str <- if (fill) {
      sprintf('    node [shape=%s, style="filled,rounded", fillcolor="%s", color=black, penwidth=1.5];\n', shape, color)
    } else {
      sprintf('    node [shape=%s, style="rounded", color=black, penwidth=1.5];\n', shape)
    }
    s <- paste0('  subgraph cluster_', id, ' {label="', title,
                '" style="rounded" color="black" penwidth=1.0;\n', style_str)
    for (i in seq_along(ids)) {
      s <- paste0(s, sprintf('    %s [label="%s"];\n', ids[i], labels[i]))
    }
    paste0(s, "  }\n")
  }

  # --- Edge builder ---
  make_edges <- function(from_ids, to_ids, mat, ws, comment,
                         dir = "row_to_col", stars_mat = NULL,
                         show_mat = NULL) {
    # dir = "row_to_col": mat[i,j] -> from_ids[j] -> to_ids[i]
    # dir = "col_to_row": mat[i,j] -> from_ids[i] -> to_ids[j]
    # stars_mat: optional matrix of "" or "*"/"**"/"***"
    # show_mat: optional logical matrix; if provided, overrides threshold filter
    s <- paste0('\n  // ', comment, '\n',
                '  edge [color="gray0", fontcolor="gray0", style=solid];\n')
    max_val <- suppressWarnings(max(mat, na.rm = TRUE))
    if (!is.finite(max_val) || max_val <= 0) return(s)
    nr <- nrow(mat); nc <- ncol(mat)
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        val <- mat[i, j]
        show <- if (!is.null(show_mat)) show_mat[i, j]
                else is.finite(val) && val >= threshold
        if (show) {
          pen <- pw(val, max_val, ws)
          lab <- fmtc(val)
          if (!is.null(stars_mat)) lab <- paste0(lab, stars_mat[i, j])
          if (dir == "row_to_col") {
            s <- paste0(s, sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                                   from_ids[j], to_ids[i], lab, pen))
          } else {
            s <- paste0(s, sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                                   from_ids[i], to_ids[j], lab, pen))
          }
        }
      }
    }
    s
  }

  # --- Significance stars and p-value filter for Theta edges ---
  C_stars <- NULL
  C_show  <- NULL  # logical matrix for sig.level filter
  if (!is.null(result$coefficients)) {
    C_stars <- matrix("", nrow = Q, ncol = R)
    C_pval  <- matrix(NA_real_, nrow = Q, ncol = R)
    cf <- result$coefficients
    dec_names <- rownames(C_mat)  # current basis names (may be renamed)
    enc_names <- colnames(C_mat)  # current covariate names (may be renamed)
    for (k in seq_len(nrow(cf))) {
      q <- match(cf$Basis[k], dec_names)
      r <- match(cf$Covariate[k], enc_names)
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

  # === Build DOT ===
  scr <- paste0(
    "digraph NMF_AE {\n",
    "  graph [rankdir=", rankdir, " compound=true];\n",
    "  splines=true; nodesep=0.4; ranksep=0.7; fontname=\"Arial\";\n",
    '  edge [fontname="Arial", fontsize=8, arrowhead=open];\n'
  )

  if (type == "YXCXY") {
    # Clusters: Y2 -> X2 -> X1 -> Y1
    scr <- paste0(scr,
      make_cluster("Y2", Y2.title, Y2_ids[idx_Y2], Y2_labels[idx_Y2], "box", "lightcoral"),
      make_cluster("X2", X2.title, X2_ids, X2_labels, "ellipse", "wheat"),
      make_cluster("X1", X1.title, X1_ids, X1_labels, "ellipse", "wheat"),
      make_cluster("Y1", Y1.title, Y1_ids[idx_Y1], Y1_labels[idx_Y1], "box", "lightblue")
    )
    # Edges: Y2 -> X2 (X2 matrix: R x P2, X2[r,j] = weight from Y2_j to X2_r)
    scr <- paste0(scr,
      make_edges(Y2_ids[idx_Y2], X2_ids, X2[, idx_Y2, drop = FALSE], weight_scale_x2, "Y2 -> X2 (encoder)", "row_to_col"))
    # Edges: X2 -> X1 (C: Q x R, C[q,r] = weight from X2_r to X1_q)
    scr <- paste0(scr,
      make_edges(X2_ids, X1_ids, C_mat, weight_scale_theta, "X2 -> X1 (C)",
                 "row_to_col", stars_mat = C_stars, show_mat = C_show))
    # Edges: X1 -> Y1 (X1: P1 x Q, X1[i,q] = weight from X1_q to Y1_i)
    scr <- paste0(scr,
      make_edges(X1_ids, Y1_ids[idx_Y1], X1[idx_Y1, , drop = FALSE], weight_scale_x1, "X1 -> Y1 (decoder)", "row_to_col"))

  } else {
    # XCX: X2 -> X1 only
    scr <- paste0(scr,
      make_cluster("X2", X2.title, X2_ids, X2_labels, "ellipse", "wheat"),
      make_cluster("X1", X1.title, X1_ids, X1_labels, "ellipse", "wheat")
    )
    scr <- paste0(scr,
      make_edges(X2_ids, X1_ids, C_mat, weight_scale_theta, "X2 -> X1 (C)",
                 "row_to_col", stars_mat = C_stars, show_mat = C_show))
  }

  result <- paste0(scr, "}\n")
  class(result) <- c("nmfae.DOT", "nmfkc.DOT")
  result
}

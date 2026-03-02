# ============================================================
# NMF-RE: Non-negative Matrix Factorization with Random Effects
# ============================================================

# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

#' Safe element-wise division with fallback
#' @param num Numerator matrix.
#' @param den Denominator matrix.
#' @param eps Small constant to avoid division by zero.
#' @return Matrix of element-wise quotients.
#' @keywords internal
#' @noRd
.nmfre.safe.div <- function(num, den, eps = 1e-10) {
  res <- num / pmax(den, eps)
  res[is.na(res) | is.infinite(res)] <- 0
  res
}

#' Normalize columns of X to sum 1, rescaling C and U
#'
#' Applies X <- X D^{-1}, C <- D C, U <- D U where D = diag(colSums(X)).
#' @param X Basis matrix (P x Q).
#' @param C Coefficient matrix (Q x K).
#' @param U Random effects matrix (Q x N), or NULL.
#' @param eps Small constant to avoid division by zero.
#' @return A list with components \code{X}, \code{C}, \code{U} (if supplied), and \code{scale}.
#' @keywords internal
#' @noRd
.nmfre.normalize.X <- function(X, C, U = NULL, eps = 1e-12) {
  col_sums <- pmax(colSums(X), eps)
  Dinv <- diag(1 / col_sums, nrow = length(col_sums))
  D    <- diag(col_sums,     nrow = length(col_sums))

  Xn <- X %*% Dinv
  Cn <- D %*% C

  if (!is.null(U)) {
    Un <- D %*% U
    return(list(X = Xn, C = Cn, U = Un, scale = col_sums))
  }
  list(X = Xn, C = Cn, scale = col_sums)
}

#' Compute effective degrees of freedom for U from eigenvalues and lambda
#' @param d Eigenvalues of X'X.
#' @param N Number of observations.
#' @param lambda Ridge penalty parameter.
#' @return Scalar dfU value.
#' @keywords internal
#' @noRd
.nmfre.dfU.from.lambda <- function(d, N, lambda) {
  N * sum(d / (d + lambda))
}

#' Find lambda that achieves a target dfU cap via bisection
#' @param d Eigenvalues of X'X.
#' @param N Number of observations.
#' @param target Target dfU value.
#' @param lambda_min Lower bound for search.
#' @param lambda_max Upper bound for search.
#' @param tol Convergence tolerance.
#' @return Scalar lambda value.
#' @keywords internal
#' @noRd
.nmfre.lambda.for.dfU.cap <- function(d, N, target,
                                       lambda_min = 0,
                                       lambda_max = 1e12,
                                       tol = 1e-8) {
  if (.nmfre.dfU.from.lambda(d, N, lambda_min) <= target) return(lambda_min)

  hi <- lambda_min
  df_hi <- Inf
  while (hi < lambda_max) {
    hi <- max(1e-12, if (hi == 0) 1e-6 else hi * 10)
    df_hi <- .nmfre.dfU.from.lambda(d, N, hi)
    if (df_hi <= target) break
  }
  if (df_hi > target) return(lambda_max)

  f <- function(lam) .nmfre.dfU.from.lambda(d, N, lam) - target
  out <- stats::uniroot(f, lower = lambda_min, upper = hi, tol = tol)
  out$root
}


# ============================================================
# nmfre() - Main estimation function
# ============================================================

#' @title Non-negative Matrix Factorization with Random Effects
#' @description
#' Estimates the NMF-RE model
#' \deqn{Y = X(\Theta A + U) + \mathcal{E}}
#' where \eqn{Y} (\eqn{P \times N}) is a non-negative observation matrix,
#' \eqn{X} (\eqn{P \times Q}) is a non-negative basis matrix learned from the data,
#' \eqn{\Theta} (\eqn{Q \times K}) is a non-negative coefficient matrix capturing
#' systematic covariate effects on latent scores,
#' \eqn{A} (\eqn{K \times N}) is a covariate matrix, and
#' \eqn{U} (\eqn{Q \times N}) is a random effects matrix capturing
#' unit-specific deviations in the latent score space.
#'
#' NMF-RE can be viewed as a mixed-effects latent-variable model defined on a
#' reconstruction (mean) structure. The non-negativity constraint on \eqn{X}
#' induces sparse, parts-based loadings, achieving measurement-side variable
#' selection without an explicit sparsity penalty. Inference on \eqn{\Theta}
#' provides covariate-side variable selection by identifying which covariates
#' significantly affect which components.
#'
#' Estimation alternates ridge-type BLUP-like closed-form updates for \eqn{U}
#' with multiplicative non-negative updates for \eqn{X} and \eqn{\Theta}.
#' The effective degrees of freedom consumed by \eqn{U} are monitored and a
#' df-based cap can be enforced to prevent near-saturated fits.
#'
#' When \code{wild.bootstrap = TRUE}, inference on \eqn{\Theta} is performed
#' conditional on \eqn{(\hat{X}, \hat{U})} via asymptotic linearization,
#' a one-step Newton update, and a multiplier (wild) bootstrap, yielding
#' standard errors, z-values, p-values, and confidence intervals without
#' repeated constrained re-optimization.
#'
#' @param Y Observation matrix (P x N), non-negative.
#' @param A Covariate matrix (K x N). Default is a row of ones (intercept only).
#' @param Q Integer. Rank of the basis matrix \eqn{X}.
#' @param dfU.cap.rate Rate for computing the dfU cap (\code{cap = rate * N * Q}).
#'   If \code{NULL} (default), runs \code{\link{nmfre.dfU.scan}} internally and
#'   selects the minimum rate where the cap is not binding.
#'   Use \code{\link{nmfre.dfU.scan}} beforehand to examine dfU behavior across rates
#'   and choose an appropriate value.
#' @param wild.bootstrap Logical. If \code{TRUE} (default), perform wild bootstrap inference
#'   on \eqn{\Theta}.
#' @param epsilon Convergence tolerance for relative change in objective (default 1e-5).
#' @param maxit Maximum number of iterations (default 50000).
#' @param ... Additional arguments for initialization, variance control, dfU control,
#'   optimization, and inference settings.
#'   \itemize{
#'     \item \code{X.init}: Initial basis matrix (P x Q), or \code{NULL}.
#'       When \code{NULL}, \code{\link{nmfkc}} is called internally to generate initial values.
#'     \item \code{C.init}: Initial coefficient matrix (Q x K), or \code{NULL}.
#'       When \code{NULL}, \code{\link{nmfkc}} is called internally to generate initial values.
#'     \item \code{U.init}: Initial random effects matrix (Q x N), or \code{NULL} (all zeros).
#'     \item \code{prefix}: Prefix for basis names (default \code{"Basis"}).
#'     \item \code{sigma2}: Initial residual variance (default 1).
#'     \item \code{sigma2.update}: Logical. Update \eqn{\sigma^2} during iterations (default \code{TRUE}).
#'     \item \code{tau2}: Initial random effect variance (default 1).
#'     \item \code{tau2.update}: Logical. Update \eqn{\tau^2} by moment matching (default \code{TRUE}).
#'       Disabled when dfU cap is active.
#'     \item \code{dfU.control}: Either \code{"cap"} (default) to enforce a cap on dfU,
#'       or \code{"off"} for no cap.
#'     \item \code{print.trace}: Logical. If \code{TRUE}, print progress every 100 iterations (default \code{FALSE}).
#'     \item \code{seed}: Integer seed for reproducibility (default 1).
#'     \item \code{C.p.side}: P-value sidedness: \code{"one.sided"} (default, for boundary null
#'       H0: C=0 vs H1: C>0) or \code{"two.sided"}.
#'     \item \code{wild.B}: Number of wild bootstrap replicates (default 500).
#'     \item \code{wild.seed}: Seed for wild bootstrap (default 123).
#'   }
#' @return A list of class \code{"nmfre"} with components.
#'   The model is \eqn{Y = X(\Theta A + U) + \mathcal{E}}.
#'
#'   \strong{Core matrices}
#'   \describe{
#'     \item{\code{X}}{Basis matrix \eqn{X} (\eqn{P \times Q}), columns normalized to sum to 1.}
#'     \item{\code{C}}{Coefficient matrix \eqn{\Theta} (\eqn{Q \times K}).}
#'     \item{\code{U}}{Random effects matrix \eqn{U} (\eqn{Q \times N}).}
#'   }
#'
#'   \strong{Variance components}
#'   \describe{
#'     \item{\code{sigma2}}{Residual variance \eqn{\hat{\sigma}^2}.}
#'     \item{\code{tau2}}{Random effect variance \eqn{\hat{\tau}^2}.}
#'     \item{\code{lambda}}{Ridge penalty \eqn{\lambda = \sigma^2 / \tau^2}.}
#'   }
#'
#'   \strong{Convergence diagnostics}
#'   \describe{
#'     \item{\code{converged}}{Logical. Whether the algorithm converged.}
#'     \item{\code{stop.reason}}{Character string describing why iteration stopped.}
#'     \item{\code{iter}}{Number of iterations performed.}
#'     \item{\code{maxit}}{Maximum iterations setting used.}
#'     \item{\code{epsilon}}{Convergence tolerance used.}
#'     \item{\code{objfunc}}{Final objective function value
#'       \eqn{\|Y - X(\Theta A + U)\|^2 + \lambda \|U\|^2}.}
#'     \item{\code{rel.change.final}}{Final relative change in objective.}
#'     \item{\code{objfunc.iter}}{Numeric vector of objective values per iteration.}
#'     \item{\code{rss.trace}}{Numeric vector of \eqn{\|Y - X(\Theta A + U)\|^2} per iteration.}
#'   }
#'
#'   \strong{Effective degrees of freedom (dfU) diagnostics}
#'   \describe{
#'     \item{\code{dfU}}{Final effective degrees of freedom
#'       \eqn{\mathrm{df}_U = N \sum_q d_q / (d_q + \lambda)},
#'       where \eqn{d_q} are eigenvalues of \eqn{X'X}.}
#'     \item{\code{dfU.cap}}{Upper bound imposed on \eqn{\mathrm{df}_U}.}
#'     \item{\code{dfU.cap.rate}}{Rate used to compute the cap.}
#'     \item{\code{dfU.cap.scan}}{Result of \code{\link{nmfre.dfU.scan}}, or \code{NULL}.}
#'     \item{\code{lambda.enforced}}{Final \eqn{\lambda} enforced to satisfy the cap.}
#'     \item{\code{dfU.hit.cap}}{Logical. Whether the cap was binding.}
#'     \item{\code{dfU.hit.iter}}{Iteration at which the cap first bound.}
#'     \item{\code{dfU.frac}}{\eqn{\mathrm{df}_U / (NQ)}, fraction of maximum df.}
#'     \item{\code{dfU.cap.frac}}{\eqn{\mathrm{df}_U^{\mathrm{cap}} / (NQ)}.}
#'   }
#'
#'   \strong{Fitted matrices}
#'   \describe{
#'     \item{\code{B}}{Fixed-effect scores \eqn{\Theta A} (\eqn{Q \times N}).}
#'     \item{\code{B.prob}}{Column-normalized probabilities from
#'       \eqn{\max(\Theta A, 0)}.}
#'     \item{\code{B.blup}}{BLUP scores \eqn{\Theta A + U} (\eqn{Q \times N}).}
#'     \item{\code{B.blup.pos}}{Non-negative BLUP scores
#'       \eqn{\max(\Theta A + U, 0)} (\eqn{Q \times N}).}
#'     \item{\code{B.blup.prob}}{Column-normalized probabilities from
#'       \eqn{\max(\Theta A + U, 0)}.}
#'     \item{\code{XB}}{Fitted values from fixed effects
#'       \eqn{X \Theta A} (\eqn{P \times N}).}
#'     \item{\code{XB.blup}}{Fitted values including random effects
#'       \eqn{X(\Theta A + U)} (\eqn{P \times N}).}
#'   }
#'
#'   \strong{Fit statistics}
#'   \describe{
#'     \item{\code{r.squared}}{Coefficient of determination \eqn{R^2} for
#'       \eqn{Y} vs \eqn{X(\Theta A + U)} (BLUP).}
#'     \item{\code{r.squared.fixed}}{Coefficient of determination \eqn{R^2} for
#'       \eqn{Y} vs \eqn{X \Theta A} (fixed effects only).}
#'     \item{\code{ICC}}{Trace-based Intraclass Correlation Coefficient.
#'       In the NMF-RE model, the conditional covariance of the \eqn{n}-th
#'       observation column is
#'       \eqn{\mathrm{Var}(Y_n) = \tau^2 X X^\top + \sigma^2 I_P},
#'       a \eqn{P \times P} matrix. Unlike a standard random intercept model
#'       where the design matrix \eqn{Z} is a simple indicator (so the ICC
#'       reduces to \eqn{\tau^2 / (\sigma^2 + \tau^2)}), the basis matrix
#'       \eqn{X} plays the role of \eqn{Z} in a random slopes model,
#'       making the variance contribution of \eqn{U} depend on \eqn{X}.
#'       To obtain a scalar summary, we take the trace of each component:
#'       \deqn{\mathrm{ICC} = \frac{\tau^2 \, \mathrm{tr}(X^\top X)}
#'       {\tau^2 \, \mathrm{tr}(X^\top X) + \sigma^2 P}.}
#'       This equals the average (over \eqn{P} dimensions) proportion of
#'       per-column variance attributable to the random effects.}
#'   }
#'
#'   \strong{Inference on \eqn{\Theta} (wild bootstrap)}
#'   \describe{
#'     \item{\code{sigma2.used}}{\eqn{\hat{\sigma}^2} used for inference.}
#'     \item{\code{C.vec.cov}}{Variance-covariance matrix for
#'       \eqn{\mathrm{vec}(\Theta)} (\eqn{QK \times QK}).}
#'     \item{\code{C.se}}{Standard error matrix for \eqn{\Theta} (\eqn{Q \times K}).}
#'     \item{\code{C.se.hess}}{Sandwich (Hessian-based) SE matrix for \eqn{\Theta}.}
#'     \item{\code{C.se.boot}}{Bootstrap SE matrix for \eqn{\Theta}.}
#'     \item{\code{coefficients}}{Data frame with columns Estimate, Std. Error, z value,
#'       Pr(>|z|), and confidence interval bounds for each element of \eqn{\Theta}.}
#'     \item{\code{C.ci.lower}}{Lower confidence interval matrix for \eqn{\Theta} (\eqn{Q \times K}).}
#'     \item{\code{C.ci.upper}}{Upper confidence interval matrix for \eqn{\Theta} (\eqn{Q \times K}).}
#'     \item{\code{C.boot.sd}}{Bootstrap standard deviation matrix for \eqn{\Theta} (\eqn{Q \times K}).}
#'     \item{\code{C.p.side}}{P-value sidedness used: \code{"one.sided"} or \code{"two.sided"}.}
#'   }
#' @export
#' @examples
#' # Example 1. cars data
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(intercept = 1, speed = cars$speed)
#' res <- nmfre(Y, A, Q = 1, maxit = 5000)
#' summary(res)
#'
#' \donttest{
#' # Example 2. Orthodont data (nlme)
#' library(nlme)
#' Y <- matrix(Orthodont$distance, 4, 27)
#' male <- ifelse(Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
#' A <- rbind(intercept = 1, male = male)
#'
#' # Scan dfU cap rates to choose an appropriate value
#' nmfre.dfU.scan(1:10/10, Y, A, Q = 1)
#'
#' # Fit with chosen rate
#' res <- nmfre(Y, A, Q = 1, dfU.cap.rate = 0.2)
#' summary(res)
#' }
#'
nmfre <- function(Y, A = NULL, Q = 2, dfU.cap.rate = NULL,
                  wild.bootstrap = TRUE, epsilon = 1e-5,
                  maxit = 50000, ...) {

  extra_args <- base::list(...)

  # --- Parameter Extraction from ... ---
  # initialization
  X.init    <- if (!is.null(extra_args$X.init)) extra_args$X.init else NULL
  C.init    <- if (!is.null(extra_args$C.init)) extra_args$C.init else NULL
  U.init    <- if (!is.null(extra_args$U.init)) extra_args$U.init else NULL
  prefix    <- if (!is.null(extra_args$prefix)) extra_args$prefix else "Basis"

  # variance handling
  sigma2        <- if (!is.null(extra_args$sigma2)) extra_args$sigma2 else 1
  sigma2.update <- if (!is.null(extra_args$sigma2.update)) extra_args$sigma2.update else TRUE
  tau2          <- if (!is.null(extra_args$tau2)) extra_args$tau2 else 1
  tau2.update   <- if (!is.null(extra_args$tau2.update)) extra_args$tau2.update else TRUE

  # variance update internals (fixed at defaults)
  sigma2.update.start <- 50
  sigma2.update.every <- 10
  sigma2.update.rate  <- 0.05
  sigma2.min <- 1e-12; sigma2.max <- 1e12
  tau2.update.start <- 1
  tau2.update.every <- 1
  tau2.update.rate  <- 0.2
  tau2.min <- 1e-12; tau2.max <- 1e12

  # dfU control
  dfU.control  <- if (!is.null(extra_args$dfU.control)) extra_args$dfU.control else "cap"

  # dfU internals (fixed at defaults)
  dfU.cap <- NULL
  dfU.enforce.every <- 1
  dfU.lambda.max <- 1e12

  # optimization extras
  print.trace <- if (!is.null(extra_args$print.trace)) extra_args$print.trace else FALSE
  seed        <- if (!is.null(extra_args$seed)) extra_args$seed else 1

  # inference settings
  C.p.side <- if (!is.null(extra_args$C.p.side)) extra_args$C.p.side else "one.sided"

  # inference internals (fixed at defaults)
  sigma2.hat  <- NULL
  df.sigma    <- "PN-df"
  C.info.mode <- "IminusH"
  sandwich    <- TRUE
  sandwich.cr1 <- TRUE
  se.rule     <- "sandwich"
  cov.ridge   <- 1e-8

  # wild bootstrap
  wild.B    <- if (!is.null(extra_args$wild.B)) extra_args$wild.B else 500
  wild.seed <- if (!is.null(extra_args$wild.seed)) extra_args$wild.seed else 123

  # wild bootstrap internals (fixed at defaults)
  wild.level <- 0.95

  set.seed(seed)
  .eps <- 1e-10

  dfU.control <- match.arg(dfU.control, choices = c("cap", "off"))
  C.p.side    <- match.arg(C.p.side, choices = c("one.sided", "two.sided"))

  # ---- dimensions ----
  P <- nrow(Y); N <- ncol(Y)
  if (is.null(A)) A <- matrix(1, 1, N)
  if (!is.matrix(A)) A <- as.matrix(A)
  stopifnot(ncol(A) == N)
  K <- nrow(A)

  # ---- init via nmfkc when X.init/C.init are NULL ----
  if (is.null(X.init) || is.null(C.init)) {
    res0 <- nmfkc(Y, A, Q = Q, epsilon = epsilon, seed = seed,
                          print.trace = print.trace, print.dims = FALSE)
    if (is.null(X.init)) X.init <- res0$X
    if (is.null(C.init)) C.init <- res0$C
  }

  X <- X.init
  stopifnot(nrow(X) == P, ncol(X) == Q)

  C_mat <- C.init
  stopifnot(nrow(C_mat) == Q, ncol(C_mat) == K)

  if (is.null(U.init)) {
    U <- matrix(0, Q, N)
  } else {
    U <- U.init
    stopifnot(nrow(U) == Q, ncol(U) == N)
  }

  # ---- dfU.cap.rate auto-selection via scan ----
  dfU.cap.scan <- NULL
  if (is.null(dfU.cap.rate) && dfU.control == "cap") {
    dfU.cap.scan <- nmfre.dfU.scan(
      Y = Y, A = A, Q = Q,
      X.init = X, C.init = C_mat, U.init = U,
      print.trace = FALSE, wild.bootstrap = FALSE,
      maxit = maxit, epsilon = epsilon, seed = seed
    )
    selected_rate <- dfU.cap.scan$cap.rate
    if (is.null(selected_rate) || is.na(selected_rate)) {
      print(dfU.cap.scan$table)
      stop("No safeguard rate found in auto scan. Try specifying dfU.cap.rate manually.")
    }
    dfU.cap.rate <- selected_rate
  }
  if (is.null(dfU.cap.rate)) dfU.cap.rate <- 0.10

  if (!is.numeric(dfU.cap.rate) || length(dfU.cap.rate) != 1 ||
      !is.finite(dfU.cap.rate) || dfU.cap.rate <= 0) {
    stop("'dfU.cap.rate' must be a positive finite number (e.g., 0.05).")
  }

  # ---- normalize X columns to sum 1 (rescale C AND U) ----
  normed <- .nmfre.normalize.X(pmax(X, .eps), pmax(C_mat, .eps), U)
  X <- normed$X
  C_mat <- normed$C
  U <- normed$U

  clip_val <- function(x, xmin, xmax) min(max(x, xmin), xmax)

  # ---- convergence bookkeeping ----
  obj <- NA_real_
  obj_prev <- Inf
  rel_change <- NA_real_

  converged <- FALSE
  stop_reason <- "maxit"
  iter_done <- 0

  obj_trace <- rep(NA_real_, maxit)
  rss_trace <- rep(NA_real_, maxit)

  # ---- dfU diagnostics ----
  dfU_last <- NA_real_
  dfU_cap_last <- NA_real_
  lambda_last <- NA_real_

  dfU_hit_cap <- FALSE
  dfU_hit_iter <- integer(0)

  # =========================================================
  # main loop
  # =========================================================
  for (iter in 1:maxit) {
    iter_done <- iter

    # (0) current components
    CA <- C_mat %*% A  # Q x N

    # (1) U-step: ridge BLUP with lambda = sigma2/tau2
    tau2 <- clip_val(tau2, tau2.min, tau2.max)
    sigma2 <- clip_val(sigma2, sigma2.min, sigma2.max)
    lambda <- sigma2 / tau2

    XtX <- crossprod(X)  # Q x Q

    # ---- dfU cap enforcement ----
    if (dfU.control == "cap" && (iter %% dfU.enforce.every == 0)) {

      dfU_cap_now <- if (!is.null(dfU.cap)) dfU.cap else dfU.cap.rate * (N * Q)
      dfU_cap_now <- max(dfU_cap_now, 1)

      d <- eigen(XtX, symmetric = TRUE, only.values = TRUE)$values
      d <- pmax(d, 0)

      dfU_now <- .nmfre.dfU.from.lambda(d, N, lambda)

      if (dfU_now > dfU_cap_now) {
        lambda_new <- .nmfre.lambda.for.dfU.cap(
          d, N, target = dfU_cap_now,
          lambda_min = pmax(lambda, 0),
          lambda_max = dfU.lambda.max
        )
        lambda <- lambda_new
        tau2 <- clip_val(sigma2 / pmax(lambda, 1e-12), tau2.min, tau2.max)
        dfU_now <- .nmfre.dfU.from.lambda(d, N, lambda)

        dfU_hit_cap <- TRUE
        dfU_hit_iter <- c(dfU_hit_iter, iter)
      }

      dfU_last <- dfU_now
      dfU_cap_last <- dfU_cap_now
    }

    lambda_last <- lambda

    M <- XtX + diag(pmax(lambda, 1e-12), Q)
    cholM <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(cholM)) stop("Cholesky failed: XtX + lambda I not SPD. Increase lambda (or decrease tau2).")

    # residual without U
    fit_fixed <- X %*% CA
    R0 <- Y - fit_fixed

    # solve for each column u_n
    for (n in 1:N) {
      rhs <- crossprod(X, R0[, n, drop = FALSE])  # Q x 1
      u <- backsolve(cholM, forwardsolve(t(cholM), rhs))
      U[, n] <- as.numeric(u)
    }

    # center U (identifiability)
    U <- sweep(U, 1, rowMeans(U), "-")

    # (1.5) tau2 update by moment matching (disabled when dfU cap is on)
    if (isTRUE(tau2.update) && dfU.control != "cap" &&
        iter >= tau2.update.start && (iter %% tau2.update.every == 0)) {

      term1 <- sum(U^2) / (N * Q)

      Minv <- tryCatch(chol2inv(cholM), error = function(e) NULL)
      if (!is.null(Minv)) {
        term2 <- sigma2 * sum(diag(Minv)) / Q
        tau2_new <- clip_val(term1 + term2, tau2.min, tau2.max)
        tau2 <- (1 - tau2.update.rate) * tau2 + tau2.update.rate * tau2_new
      }
    }

    # (2) X-step (EU MU) with positive-part stabilization
    Y_tilde <- pmax(Y, 0)
    B_sem <- CA + U
    B_pos <- pmax(B_sem, 0)

    numX <- Y_tilde %*% t(B_pos)
    denX <- X %*% (B_pos %*% t(B_pos))
    X <- X * .nmfre.safe.div(numX, denX, eps = .eps)
    X <- pmax(X, .eps)

    # normalize columns sum 1; absorb scale into C AND U
    normed <- .nmfre.normalize.X(X, C_mat, U)
    X <- pmax(normed$X, .eps)
    C_mat <- pmax(normed$C, .eps)
    U <- normed$U

    # (3) C-step (EU MU) with positive-part stabilization
    Y_star <- Y_tilde - X %*% U
    Y_star_pos <- pmax(Y_star, 0)

    numC <- crossprod(X, Y_star_pos) %*% t(A)
    denC <- (crossprod(X, X) %*% C_mat) %*% (A %*% t(A))
    C_mat <- C_mat * .nmfre.safe.div(numC, denC, eps = .eps)
    C_mat <- pmax(C_mat, .eps)

    # (3.5) sigma2 update from residual
    if (isTRUE(sigma2.update) && iter >= sigma2.update.start &&
        (iter %% sigma2.update.every == 0)) {
      CA <- C_mat %*% A
      fit <- X %*% (CA + U)
      R <- Y - fit
      sigma2_new <- clip_val(mean(R^2), sigma2.min, sigma2.max)
      sigma2 <- (1 - sigma2.update.rate) * sigma2 + sigma2.update.rate * sigma2_new
    }

    # (4) objective & convergence
    CA <- C_mat %*% A
    fit <- X %*% (CA + U)
    R <- Y - fit

    lambda_obj <- sigma2 / clip_val(tau2, tau2.min, tau2.max)
    obj <- sum(R^2) + lambda_obj * sum(U^2)
    rss <- sum(R^2)

    obj_trace[iter] <- obj
    rss_trace[iter] <- rss

    if (is.finite(obj_prev)) {
      rel_change <- abs(obj_prev - obj) / (abs(obj_prev) + .eps)
    } else {
      rel_change <- NA_real_
    }

    if (!is.finite(obj)) {
      stop_reason <- "nonfinite_obj"
      converged <- FALSE
      break
    }

    if (is.finite(obj_prev) && rel_change < epsilon) {
      stop_reason <- "epsilon"
      converged <- TRUE
      break
    }

    obj_prev <- obj

    if (isTRUE(print.trace) && iter %% 100 == 0) {
      message(sprintf("iter=%d  obj=%.6g  rel=%.3g  sigma2=%.4g  tau2=%.4g  lambda=%.4g",
                      iter, obj, rel_change, sigma2, tau2, sigma2 / tau2))
    }
  }

  # ---- post loop: finalize ----
  obj_final <- obj
  rel_change_final <- rel_change

  if (iter_done < maxit && identical(stop_reason, "maxit")) {
    stop_reason <- "unknown_break"
    converged <- FALSE
  }

  # ---- reorder basis ----
  if (ncol(X) > 1) {
    w_ord <- matrix((1:P) / P, nrow = 1)
    score <- as.numeric(w_ord %*% X)
    index <- order(score)
    X <- X[, index, drop = FALSE]
    C_mat <- C_mat[index, , drop = FALSE]
    U <- U[index, , drop = FALSE]
  }

  # ---- final fitted matrices ----
  CA <- C_mat %*% A
  B_fixed <- CA
  B_blup_raw <- CA + U

  XB <- X %*% B_fixed
  XB.blup <- X %*% B_blup_raw

  # ---- final dfU recomputation ----
  lambda_final <- sigma2 / clip_val(tau2, tau2.min, tau2.max)
  XtX_final <- crossprod(X)
  d_final <- eigen(XtX_final, symmetric = TRUE, only.values = TRUE)$values
  d_final <- pmax(d_final, 0)
  dfU_final <- .nmfre.dfU.from.lambda(d_final, N, lambda_final)

  if (dfU.control == "cap") {
    dfU_cap_last <- if (!is.null(dfU.cap)) max(dfU.cap, 1) else max(dfU.cap.rate * (N * Q), 1)
  } else {
    dfU_cap_last <- NA_real_
  }
  dfU_last <- dfU_final
  lambda_last <- lambda_final

  # ---- names ----
  colnames(X) <- paste0(prefix, 1:ncol(X))
  rownames(C_mat) <- colnames(X)
  colnames(C_mat) <- rownames(A)
  rownames(U) <- rownames(C_mat)
  colnames(U) <- colnames(Y)

  # ---- r.squared ----
  ss_total <- sum((Y - mean(Y))^2)
  R_blup <- Y - XB.blup
  R_fixed <- Y - XB
  r.squared <- if (ss_total > 0) 1 - sum(R_blup^2) / ss_total else NA_real_
  r.squared.fixed <- if (ss_total > 0) 1 - sum(R_fixed^2) / ss_total else NA_real_

  # ---- ICC (trace-based) ----
  trXtX <- sum(d_final)  # tr(X'X), using eigenvalues already computed
  ICC <- tau2 * trXtX / (tau2 * trXtX + sigma2 * P)

  # =========================================================
  # Inference on C
  # =========================================================
  C.vec.cov <- NULL
  C.se <- NULL
  C.se.hess <- NULL
  C.se.boot <- NULL
  coefficients <- NULL
  sigma2_used <- sigma2.hat

  C.ci.lower <- NULL
  C.ci.upper <- NULL
  C.boot.sd <- NULL

  if (isTRUE(wild.bootstrap)) {

    # Y_star for inference (conditioning on X,U)
    Y_star_inf <- Y - X %*% U
    R_C <- Y_star_inf - X %*% (C_mat %*% A)
    RSS_inf <- sum(R_C^2)

    lambda_inf <- sigma2 / pmax(tau2, 1e-12)

    XtX_now <- crossprod(X)
    AAt <- A %*% t(A)

    M_inf <- XtX_now + diag(pmax(lambda_inf, 1e-12), Q)
    cholM_inf <- tryCatch(chol(M_inf), error = function(e) NULL)
    if (!is.null(cholM_inf)) {
      Minv <- chol2inv(cholM_inf)
    } else {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for generalized inverse when Cholesky fails.")
      }
      Minv <- tryCatch(solve(M_inf), error = function(e) MASS::ginv(M_inf))
    }

    trH <- sum(diag(Minv %*% XtX_now))
    dfU_inf <- N * trH

    if (is.null(sigma2_used)) {
      if (df.sigma == "PN") {
        denom <- P * N
      } else if (df.sigma == "PN-QR") {
        denom <- max(P * N - Q * K, 1)
      } else { # "PN-df"
        dfC <- Q * K
        denom <- max(P * N - dfU_inf - dfC, 1)
      }
      sigma2_used <- RSS_inf / denom
    }

    if (C.info.mode == "plain") {
      Info_core <- kronecker(AAt, XtX_now)
    } else {
      Xt_IH_X <- XtX_now - XtX_now %*% Minv %*% XtX_now
      Info_core <- kronecker(AAt, Xt_IH_X)
    }

    Info <- Info_core / max(sigma2_used, 1e-12)
    Info <- Info + diag(cov.ridge, nrow(Info))

    if (!requireNamespace("MASS", quietly = TRUE)) {
      Hinv <- tryCatch(solve(Info), error = function(e) {
        stop("Package 'MASS' is required for generalized inverse.")
      })
    } else {
      Hinv <- tryCatch(solve(Info), error = function(e) MASS::ginv(Info))
    }
    V_naive <- Hinv

    V_sand <- NULL
    if (isTRUE(sandwich)) {
      Xt <- t(X)
      J <- matrix(0, Q * K, Q * K)

      for (n in 1:N) {
        a_n <- A[, n, drop = FALSE]
        g_n <- Xt %*% R_C[, n, drop = FALSE]
        S_n <- -(g_n %*% t(a_n)) / max(sigma2_used, 1e-12)
        s_n <- as.vector(S_n)
        J <- J + tcrossprod(s_n)
      }

      if (isTRUE(sandwich.cr1) && N > 1) {
        J <- (N / (N - 1)) * J
      }
      V_sand <- Hinv %*% J %*% Hinv
    }

    if (se.rule == "naive") {
      C.vec.cov <- V_naive
    } else if (se.rule == "sandwich") {
      C.vec.cov <- if (is.null(V_sand)) V_naive else V_sand
    } else { # "max"
      if (is.null(V_sand)) {
        C.vec.cov <- V_naive
      } else {
        d1 <- pmax(diag(V_naive), 0)
        d2 <- pmax(diag(V_sand), 0)
        dmax <- pmax(d1, d2)
        C.vec.cov <- V_naive
        diag(C.vec.cov) <- dmax
      }
    }

    se_vec <- sqrt(pmax(diag(C.vec.cov), 0))
    C.se.hess <- matrix(se_vec, nrow = Q, ncol = K, byrow = FALSE)
    dimnames(C.se.hess) <- dimnames(C_mat)

    C.se <- C.se.hess

    # ---- Multiplier (wild) bootstrap ----
    set.seed(wild.seed)

    Xt <- t(X)
    score_mat <- matrix(0, Q * K, N)

    for (n in 1:N) {
      a_n <- A[, n, drop = FALSE]
      r_n <- R_C[, n, drop = FALSE]
      g_n <- Xt %*% r_n
      G_n <- -(g_n %*% t(a_n)) / max(sigma2_used, 1e-12)
      score_mat[, n] <- as.vector(G_n)
    }

    C_hat_vec <- as.vector(C_mat)
    C_boot <- matrix(NA_real_, nrow = Q * K, ncol = wild.B)

    for (b in 1:wild.B) {
      w <- stats::rexp(N, rate = 1) - 1
      grad_b <- as.vector(score_mat %*% w)
      c_b <- C_hat_vec - as.vector(Hinv %*% grad_b)
      c_b <- pmax(c_b, 0)
      C_boot[, b] <- c_b
    }

    alpha <- 1 - wild.level
    lo_q <- alpha / 2
    hi_q <- 1 - alpha / 2

    lo <- apply(C_boot, 1, stats::quantile, probs = lo_q, na.rm = TRUE, names = FALSE)
    hi <- apply(C_boot, 1, stats::quantile, probs = hi_q, na.rm = TRUE, names = FALSE)

    C.ci.lower <- matrix(lo, nrow = Q, ncol = K, byrow = FALSE)
    C.ci.upper <- matrix(hi, nrow = Q, ncol = K, byrow = FALSE)
    dimnames(C.ci.lower) <- dimnames(C_mat)
    dimnames(C.ci.upper) <- dimnames(C_mat)

    sd_vec <- apply(C_boot, 1, stats::sd, na.rm = TRUE)
    C.boot.sd <- matrix(sd_vec, nrow = Q, ncol = K, byrow = FALSE)
    dimnames(C.boot.sd) <- dimnames(C_mat)

    C.se.boot <- C.boot.sd

    # ---- coefficients table ----
    Estimate <- as.vector(C_mat)
    SE_hess_vec <- as.vector(C.se.hess)
    SE <- SE_hess_vec
    BSE_vec <- as.vector(C.boot.sd)

    z_value <- rep(NA_real_, length(Estimate))
    ok <- is.finite(Estimate) & is.finite(SE) & (SE > 0)
    z_value[ok] <- Estimate[ok] / SE[ok]

    Wald <- rep(NA_real_, length(Estimate))
    Wald[ok] <- z_value[ok]^2

    p_two <- rep(NA_real_, length(Estimate))
    p_two[ok] <- 1 - stats::pchisq(Wald[ok], df = 1)

    p_one <- rep(NA_real_, length(Estimate))
    p_one[ok] <- stats::pnorm(z_value[ok], lower.tail = FALSE)

    p_value <- if (C.p.side == "one.sided") p_one else p_two

    coefficients <- data.frame(
      Variable = rep(colnames(C_mat), each = Q),
      Basis    = rep(rownames(C_mat), times = K),
      Estimate = Estimate,
      SE       = SE,
      BSE      = BSE_vec,
      z_value  = z_value,
      Wald     = Wald,
      p_value  = p_value,
      p_two_sided = p_two,
      p_one_sided = p_one,
      p_side_used = C.p.side,
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    if (!is.null(C.ci.lower) && !is.null(C.ci.upper)) {
      coefficients$CI_low  <- as.vector(C.ci.lower)
      coefficients$CI_high <- as.vector(C.ci.upper)
    }
  }

  # ---- convenience probabilities ----
  colnorm_prob <- function(M, eps = 1e-12) {
    cs <- colSums(M)
    sweep(M, 2, pmax(cs, eps), "/")
  }

  out <- list(
    X = X,
    C = C_mat,
    U = U,

    sigma2 = sigma2,
    tau2 = tau2,
    lambda = sigma2 / tau2,

    # convergence diagnostics
    converged = converged,
    stop.reason = stop_reason,
    iter = iter_done,
    maxit = maxit,
    epsilon = epsilon,
    objfunc = obj_final,
    rel.change.final = rel_change_final,
    objfunc.iter = obj_trace[seq_len(iter_done)],
    rss.trace = rss_trace[seq_len(iter_done)],

    # dfU diagnostics
    dfU = dfU_last,
    dfU.cap = dfU_cap_last,
    dfU.cap.rate = dfU.cap.rate,
    dfU.cap.scan = dfU.cap.scan,
    lambda.enforced = lambda_last,
    dfU.hit.cap  = dfU_hit_cap,
    dfU.hit.iter = dfU_hit_iter,
    dfU.frac     = if (is.finite(dfU_last)) dfU_last / (N * Q) else NA_real_,
    dfU.cap.frac = if (is.finite(dfU_cap_last)) dfU_cap_last / (N * Q) else NA_real_,

    # fitted matrices
    B = B_fixed,
    B.prob = colnorm_prob(pmax(B_fixed, 0)),
    B.blup = B_blup_raw,
    B.blup.pos = pmax(B_blup_raw, 0),
    B.blup.prob = colnorm_prob(pmax(B_blup_raw, 0)),
    XB = XB,
    XB.blup = XB.blup,

    # fit statistics
    r.squared = r.squared,
    r.squared.fixed = r.squared.fixed,
    ICC = ICC,

    # inference outputs
    sigma2.used = sigma2_used,
    C.vec.cov = C.vec.cov,
    C.se = C.se,
    C.se.hess = C.se.hess,
    C.se.boot = C.se.boot,
    coefficients = coefficients,
    C.ci.lower = C.ci.lower,
    C.ci.upper = C.ci.upper,
    C.boot.sd  = C.boot.sd,

    C.p.side = C.p.side
  )

  class(out) <- "nmfre"
  out
}


# ============================================================
# summary.nmfre() - S3 summary method
# ============================================================

#' @title Summary method for objects of class \code{nmfre}
#' @description
#' Displays a concise summary of an NMF-RE model fit, including dimensions,
#' convergence, variance components, and a coefficient table following
#' standard R regression output conventions.
#'
#' @param object An object of class \code{nmfre}, returned by \code{\link{nmfre}}.
#' @param show_ci Logical. If \code{TRUE}, show confidence interval columns (default \code{FALSE}).
#' @param ... Additional arguments (currently unused).
#' @return The input object, invisibly.
#' @export
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(intercept = 1, speed = cars$speed)
#' res <- nmfre(Y, A, Q = 1, maxit = 5000)
#' summary(res)
#'
summary.nmfre <- function(object, show_ci = FALSE, ...) {

  x <- object
  P <- nrow(x$X); Q <- ncol(x$X)
  K <- ncol(x$C); N <- ncol(x$U)

  # ---- header ----
  cat(sprintf("NMF-RE: Y(%d,%d) = X(%d,%d) [C(%d,%d) A + U(%d,%d)]\n",
              P, N, P, Q, Q, K, Q, N))

  # ---- convergence ----
  conv_str <- if (isTRUE(x$converged)) "converged" else "NOT converged"
  cat(sprintf("Iterations: %d (%s, epsilon = %s)\n",
              x$iter, conv_str, format(x$epsilon, digits = 1)))

  # ---- R-squared ----
  if (!is.null(x$r.squared) && is.finite(x$r.squared)) {
    r2_fixed <- if (!is.null(x$r.squared.fixed) && is.finite(x$r.squared.fixed)) {
      sprintf(", %.4f (XB)", x$r.squared.fixed)
    } else ""
    cat(sprintf("R-squared: %.4f (XB+blup)%s\n", x$r.squared, r2_fixed))
  }

  # ---- variance components ----
  cat(sprintf("\nVariance components:\n"))
  cat(sprintf("  sigma2 = %.4g  (residual)\n", x$sigma2))
  cat(sprintf("  tau2   = %.4g  (random effect)\n", x$tau2))
  cat(sprintf("  lambda = %.4g  (sigma2 / tau2)\n", x$lambda))
  cat(sprintf("  ICC    = %.4f  (tau2*tr(X'X) / (tau2*tr(X'X) + sigma2*P))\n", x$ICC))
  if (!is.null(x$dfU.cap) && is.finite(x$dfU.cap)) {
    cat(sprintf("  dfU    = %.2f <= %.2f (cap rate = %.2f)\n",
                x$dfU, x$dfU.cap, x$dfU.cap.rate))
  } else if (is.finite(x$dfU)) {
    cat(sprintf("  dfU    = %.2f\n", x$dfU))
  }

  # ---- coefficients ----
  if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
    cat("\nCoefficients:\n")

    cf <- x$coefficients

    # row names: Variable:Basis
    rnames <- paste0(cf$Variable, ":", cf$Basis)

    # p-value formatting (lm-style)
    p_side <- if (!is.null(x$C.p.side)) x$C.p.side else "one.sided"
    pv <- cf$p_value

    format_pval <- function(p) {
      ifelse(!is.finite(p), "      NA",
        ifelse(p < 2.2e-16, "  <2e-16",
          formatC(p, format = "g", digits = 4, width = 8)))
    }

    # significance stars
    sig_stars <- function(p) {
      ifelse(!is.finite(p), " ",
        ifelse(p < 0.001, "***",
          ifelse(p < 0.01, "**",
            ifelse(p < 0.05, "*",
              ifelse(p < 0.1, ".", " ")))))
    }

    # build display matrix
    p_header <- if (p_side == "one.sided") "Pr(>z)" else "Pr(>|z|)"

    # columns
    est <- formatC(cf$Estimate, format = "f", digits = 3, width = 9)
    se  <- formatC(cf$SE, format = "f", digits = 3, width = 10)
    zv  <- formatC(cf$z_value, format = "f", digits = 2, width = 7)
    pv_str <- format_pval(pv)
    stars <- sig_stars(pv)

    has_bse <- "BSE" %in% names(cf) && any(is.finite(cf$BSE))
    if (has_bse) {
      bse <- formatC(cf$BSE, format = "f", digits = 3, width = 6)
    }

    # header line
    if (has_bse) {
      hdr <- sprintf("%s %s %s %s %s %s",
                     formatC("Estimate", width = 9),
                     formatC("Std. Error", width = 10),
                     formatC("(Boot)", width = 6),
                     formatC("z value", width = 7),
                     formatC(p_header, width = 8),
                     "")
    } else {
      hdr <- sprintf("%s %s %s %s %s",
                     formatC("Estimate", width = 9),
                     formatC("Std. Error", width = 10),
                     formatC("z value", width = 7),
                     formatC(p_header, width = 8),
                     "")
    }

    # max label width for alignment
    max_lw <- max(nchar(rnames))
    cat(sprintf("%s %s\n", formatC("", width = max_lw), hdr))

    for (i in seq_along(rnames)) {
      lab <- formatC(rnames[i], width = max_lw)
      if (has_bse) {
        cat(sprintf("%s %s %s %s %s %s %s\n",
                    lab, est[i], se[i], bse[i], zv[i], pv_str[i], stars[i]))
      } else {
        cat(sprintf("%s %s %s %s %s %s\n",
                    lab, est[i], se[i], zv[i], pv_str[i], stars[i]))
      }
    }

    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

    # CI
    if (isTRUE(show_ci) && all(c("CI_low", "CI_high") %in% names(cf))) {
      level <- if (!is.null(x$wild.level)) x$wild.level else 0.95
      cat(sprintf("\n%.0f%% Bootstrap CI:\n", level * 100))
      ci_lo <- formatC(cf$CI_low, format = "f", digits = 3, width = 9)
      ci_hi <- formatC(cf$CI_high, format = "f", digits = 3, width = 9)
      for (i in seq_along(rnames)) {
        cat(sprintf("  %s: [%s, %s]\n", rnames[i], ci_lo[i], ci_hi[i]))
      }
    }
  }

  invisible(x)
}


# ============================================================
# nmfre.dfU.scan() - dfU scan utility
# ============================================================

#' @title Scan dfU cap rates for NMF-RE
#' @description
#' Fits the NMF-RE model across a range of \code{dfU.cap.rate} values and
#' returns a diagnostic table showing the resulting effective degrees of freedom,
#' variance components, and convergence diagnostics for each rate.
#'
#' The dfU cap limits the effective degrees of freedom consumed by the random
#' effects \eqn{U}. The cap is computed as \code{rate * N * Q}, where \eqn{N}
#' is the number of observations and \eqn{Q} is the rank. A suitable rate is
#' one where the final \eqn{\mathrm{df}_U} is below the cap
#' (\code{safeguard = TRUE}) and the model has converged
#' (\code{converged = TRUE}).
#'
#' When called automatically by \code{\link{nmfre}} (i.e., \code{dfU.cap.rate = NULL}),
#' the minimum rate satisfying both \code{safeguard = TRUE} and
#' \code{converged = TRUE} is selected.
#'
#' @param rates Numeric vector of cap rates to scan (default \code{(1:10)/10}).
#' @param Y Observation matrix (P x N).
#' @param A Covariate matrix (K x N).
#' @param Q Integer. Rank of the basis matrix.
#' @param X.init Initial basis matrix, or \code{NULL}.
#' @param C.init Initial coefficient matrix, or \code{NULL}.
#' @param U.init Initial random effects matrix, or \code{NULL}.
#' @param print.trace Logical. Print progress for each fit (default \code{FALSE}).
#' @param ... Additional arguments passed to \code{\link{nmfre}}.
#' @return An object of class \code{"nmfre.dfU.scan"} with two components:
#' \describe{
#'   \item{\code{table}}{A data frame with the following columns:
#'   \describe{
#'     \item{\code{rate}}{Cap rate used. The dfU cap is \code{rate * N * Q}.}
#'     \item{\code{dfU.cap}}{The dfU cap value (upper bound on effective degrees of freedom).}
#'     \item{\code{dfU}}{Final effective degrees of freedom for \eqn{U} at convergence.}
#'     \item{\code{safeguard}}{Logical. \code{TRUE} if the dfU cap is functioning
#'       as a safeguard (\code{dfU / dfU.cap < 0.99}): the cap prevents
#'       random-effects saturation without over-constraining \eqn{U}.
#'       \code{FALSE} if dfU is at or near the cap, indicating the cap is
#'       binding and the rate may be too small.}
#'     \item{\code{hit}}{Logical. \code{TRUE} if the cap was reached at least once
#'       during iteration, even if dfU later decreased below the cap.}
#'     \item{\code{converged}}{Logical. \code{TRUE} if the algorithm converged
#'       within the maximum number of iterations.}
#'     \item{\code{tau2}}{Final random effect variance \eqn{\hat{\tau}^2}.}
#'     \item{\code{sigma2}}{Final residual variance \eqn{\hat{\sigma}^2}.}
#'     \item{\code{ICC}}{Trace-based Intraclass Correlation Coefficient
#'       \eqn{\tau^2 \, \mathrm{tr}(X^\top X) /
#'       (\tau^2 \, \mathrm{tr}(X^\top X) + \sigma^2 P)}.
#'       See \code{\link{nmfre}} for details.}
#'   }}
#'   \item{\code{cap.rate}}{Optimal cap rate selected automatically.
#'     If rows with \code{safeguard = TRUE} and \code{hit = TRUE} exist,
#'     the maximum rate among them is chosen (safeguard activated but
#'     giving \eqn{U} the most freedom). Otherwise, the minimum rate
#'     with \code{safeguard = TRUE} and \code{hit = FALSE} is chosen.
#'     \code{NA} if no suitable rate is found.}
#' }
#'
#' When printed, only the \code{table} is displayed. Access \code{cap.rate}
#' directly from the returned object.
#' @export
#' @examples
#' # Example 1. cars data (small maxit for speed)
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(intercept = 1, speed = cars$speed)
#' tab <- nmfre.dfU.scan(rates = c(0.1, 0.2), Y = Y, A = A, Q = 1, maxit = 1000)
#' print(tab)
#'
#' \donttest{
#' # Example 2. Orthodont data (nlme)
#' library(nlme)
#' Y <- matrix(Orthodont$distance, 4, 27)
#' male <- ifelse(Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
#' A <- rbind(intercept = 1, male = male)
#'
#' nmfre.dfU.scan(1:10/10, Y, A, Q = 1)
#' }
#'
nmfre.dfU.scan <- function(
  rates = (1:10) / 10,
  Y, A, Q,
  X.init = NULL, C.init = NULL, U.init = NULL,
  print.trace = FALSE,
  ...
) {
  N <- ncol(Y)
  NQ <- N * Q

  one_run <- function(rate) {

    r <- nmfre(
      Y = Y, A = A, Q = Q,
      X.init = X.init, C.init = C.init, U.init = U.init,
      print.trace = print.trace,
      dfU.control = "cap",
      dfU.cap.rate = rate,
      ...
    )

    safeguard <- if (is.finite(r$dfU) && is.finite(r$dfU.cap) && r$dfU.cap > 0) {
      r$dfU / r$dfU.cap < 0.99
    } else {
      NA
    }
    data.frame(
      rate = rate,
      dfU.cap = r$dfU.cap,
      dfU = round(r$dfU, 2),
      safeguard = safeguard,
      hit = isTRUE(r$dfU.hit.cap),
      converged = isTRUE(r$converged),
      tau2 = round(r$tau2, 3),
      sigma2 = round(r$sigma2, 3),
      ICC = round(r$ICC, 4),
      stringsAsFactors = FALSE
    )
  }
  tbl <- do.call(rbind, lapply(rates, one_run))
  rownames(tbl) <- NULL

  # optimal rate selection
  sg_hit <- tbl[tbl$safeguard == TRUE & tbl$hit == TRUE & tbl$converged, ]
  sg_free <- tbl[tbl$safeguard == TRUE & tbl$hit == FALSE & tbl$converged, ]
  cap.rate <- if (nrow(sg_hit) > 0) {
    max(sg_hit$rate)
  } else if (nrow(sg_free) > 0) {
    min(sg_free$rate)
  } else {
    NA_real_
  }

  result <- list(table = tbl, cap.rate = cap.rate)
  class(result) <- "nmfre.dfU.scan"
  result
}

#' @export
print.nmfre.dfU.scan <- function(x, ...) {
  print(x$table, ...)
  invisible(x)
}

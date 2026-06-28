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

#' Semi-NMF multiplicative update for X >= 0 with a sign-free score matrix
#'
#' Solves (monotonically) \eqn{\min_{X\ge0}\|Y-XB\|_F^2} for a real-valued
#' \eqn{B} (Ding, Li & Jordan, 2010): with \eqn{G=BB^\top},
#' \eqn{X\leftarrow X\odot\sqrt{([YB^\top]^+ + X[G]^-)/([YB^\top]^- + X[G]^+)}}.
#' When a symmetric PD matrix \code{S} is supplied (the EM posterior-variance
#' term), it is added to \eqn{G}, solving \eqn{\min_{X\ge0}\|Y-XB\|_F^2 +
#' \mathrm{tr}(XSX^\top)} (still monotone).
#' @keywords internal
#' @noRd
.nmfre.seminmf.X <- function(X, Y, B, S = NULL, eps = 1e-10) {
  YBt <- Y %*% t(B)
  G   <- B %*% t(B)
  if (!is.null(S)) G <- G + S
  pos <- function(M) pmax(M, 0)
  neg <- function(M) pmax(-M, 0)
  numer <- pos(YBt) + X %*% neg(G)
  denom <- neg(YBt) + X %*% pos(G)
  X * sqrt(.nmfre.safe.div(numer, denom, eps = eps))
}

#' Marginal negative log-likelihood of the NMF-RE model
#'
#' \eqn{\ell = \tfrac12\big[N\,(P\log 2\pi + \log|\Sigma|) + \sum_n
#' (y_n - X\Theta a_n)^\top \Sigma^{-1}(y_n - X\Theta a_n)\big]}, with
#' \eqn{\Sigma = \sigma^2 I_P + \tau^2 X X^\top}. The random effects \eqn{U}
#' are integrated out, so this is the quantity the ECM decreases monotonically
#' (unlike the fixed-\eqn{\lambda} penalized objective, which jumps when
#' \eqn{\lambda=\sigma^2/\tau^2} is updated).
#'
#' Evaluated in the Woodbury / \eqn{Q\times Q} form (efficient and numerically
#' stable for \eqn{Q \ll P}); identical to the standalone NMF-RE core:
#' \eqn{\log|\Sigma| = P\log\sigma^2 + \log|I_Q + (\tau^2/\sigma^2)X^\top X|} and
#' \eqn{\mathrm{quad} = \sigma^{-2}\big[\|R\|^2 -
#' \mathrm{tr}\{R^\top X(X^\top X + \lambda I_Q)^{-1}X^\top R\}\big]},
#' \eqn{R = Y - X\Theta A}, \eqn{\lambda = \sigma^2/\tau^2}.
#' @keywords internal
#' @noRd
.nmfre.marginal.nll <- function(Y, X, C, A, sigma2, tau2, eps = 1e-12) {
  P <- nrow(Y); N <- ncol(Y); Q <- ncol(X)
  s2 <- max(sigma2, eps); t2 <- max(tau2, eps)
  lambda <- s2 / t2
  Rm   <- Y - X %*% (C %*% A)                       # P x N fixed-effects residual
  XtRm <- crossprod(X, Rm)                          # Q x N
  Mll  <- crossprod(X) + diag(pmax(lambda, 1e-12), Q)
  quad <- tryCatch((sum(Rm^2) - sum(XtRm * solve(Mll, XtRm))) / s2,
                   error = function(e) NA_real_)
  logdetQ <- tryCatch(
    as.numeric(determinant(diag(Q) + (t2 / s2) * crossprod(X), logarithm = TRUE)$modulus),
    error = function(e) NA_real_)
  if (!is.finite(quad) || !is.finite(logdetQ)) return(NA_real_)
  0.5 * (N * (P * log(2 * pi) + P * log(s2) + logdetQ) + quad)
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
#' Estimation is an outer--inner ECM: an inner loop (fixed \eqn{\lambda}) runs a
#' ridge/BLUP update for \eqn{U}, a complete-EM semi-NMF update for \eqn{X}, and
#' a fixed-effect update for \eqn{\Theta}; an outer loop runs the EM M-steps for
#' \eqn{(\sigma^2, \tau^2)}. The variance components are estimated from the data;
#' the effective degrees of freedom \eqn{df_U} are reported only as a diagnostic.
#'
#' \code{nmfre()} performs \strong{optimization only}. Hypothesis tests and
#' standard errors for \eqn{\Theta} are obtained separately with
#' \code{\link{nmfre.inference}} (sandwich SE + wild bootstrap), mirroring the
#' \code{\link{nmfkc}} / \code{nmfkc.inference} split.
#'
#' @param Y Observation matrix (P x N), non-negative.
#' @param A Covariate matrix (K x N). Default is a row of ones (intercept only).
#' @param rank Integer. Rank of the basis matrix \eqn{X}. Default is 2.
#'   For backward compatibility, \code{Q} is accepted via \code{...}.
#' @param C.signed Logical. Whether the fixed-effect coefficients \eqn{C}
#'   (\eqn{= \Theta} in the paper) are sign-free; the single switch that selects
#'   the whole estimation scheme. \code{TRUE} (default, recommended, matches the
#'   paper) treats \eqn{C} as real-valued and updates it by exact least squares,
#'   with the basis \eqn{X} estimated by the complete-EM semi-NMF step and a
#'   two-sided test (interior null). \code{FALSE} constrains \eqn{C \ge 0} via a
#'   multiplicative update, with \eqn{X} estimated by the positive-part
#'   multiplicative update and a one-sided test (boundary null). The basis
#'   \eqn{X} is always non-negative. A character value (\code{"signed"} /
#'   \code{"nonneg"}) is also accepted for backward compatibility.
#' @param epsilon Convergence tolerance for relative change in objective (default 1e-5).
#' @param maxit Maximum number of iterations.  Default \code{5000}
#'   (matches \code{\link{nmfkc}} and the other MU functions in the
#'   package).  When the cap is hit without meeting the relative-
#'   tolerance criterion, a \code{"maximum iterations (...) reached..."}
#'   warning is emitted so users notice unconverged fits.
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
#'     \item \code{x.postvar}: Logical. Include the posterior-variance term in
#'       the semi-NMF \eqn{X}-step (default \code{TRUE}; advanced). Applies only
#'       to the sign-free / semi-NMF path (\code{C.signed = TRUE}).
#'     \item \code{dfU.control}: Deprecated and inert. The algorithm imposes no
#'       cap on \eqn{df_U} (\code{"off"}, the only behaviour); \eqn{df_U} is
#'       reported as a diagnostic only.
#'     \item \code{print.trace}: Logical. If \code{TRUE}, print progress every 100 iterations (default \code{FALSE}).
#'     \item \code{seed}: Integer seed for reproducibility (default 1).
#'     \item \code{nstart}: Number of random restarts for the \code{nmfkc()}
#'       initialisation step (passed to the k-means initialiser). Default
#'       \code{1} (single start; historical behaviour). A larger value
#'       (e.g.\ 10-20) gives a more stable initialisation.
#'     \item \code{inner.maxit}, \code{outer.maxit}: Maximum inner (fixed-\eqn{\lambda}
#'       block-coordinate) and outer (EM variance) iterations (defaults
#'       \code{10000} and \code{500}).
#'     \item \code{epsilon.outer}: Convergence tolerance for the outer EM loop on
#'       \eqn{\lambda} (default \code{1e-6}).
#'   }
#' @return A list of class \code{"nmfre"} with components.
#'   The model is \eqn{Y = X(\Theta A + U) + \mathcal{E}}.
#'
#'   \strong{Core matrices}
#'   \describe{
#'     \item{\code{X}}{Basis matrix \eqn{X} (\eqn{P \times Q}), columns normalized to sum to 1.}
#'     \item{\code{X.prob}}{Row-wise soft-clustering probabilities from the
#'       non-negative \eqn{X} (each row normalized to sum to 1), as in \code{\link{nmfkc}}.}
#'     \item{\code{X.cluster}}{Hard-clustering label for each row of \eqn{X}
#'       (argmax over \code{X.prob}).}
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
#'     \item{\code{objfunc.iter}}{Numeric vector of the fixed-\eqn{\lambda}
#'       penalized objective \eqn{\|Y - X(\Theta A + U)\|^2 + \lambda\|U\|^2}
#'       per iteration. Monotone within an inner loop but \emph{not} across
#'       outer iterations (the penalty jumps when \eqn{\lambda} is updated).}
#'     \item{\code{rss.trace}}{Numeric vector of \eqn{\|Y - X(\Theta A + U)\|^2} per iteration.}
#'     \item{\code{nll.trace}}{Numeric vector of the marginal negative
#'       log-likelihood \eqn{\ell(X,\Theta,\sigma^2,\tau^2)} per iteration
#'       (random effects integrated out). This is the ECM-monotone quantity and
#'       is what \code{\link{plot.nmfre}} displays.}
#'   }
#'
#'   \strong{Effective degrees of freedom (dfU) diagnostics}
#'   \describe{
#'     \item{\code{dfU}}{Final effective degrees of freedom
#'       \eqn{\mathrm{df}_U = N \sum_q d_q / (d_q + \lambda)},
#'       where \eqn{d_q} are eigenvalues of \eqn{X'X}.}
#'     \item{\code{dfU.cap}}{Upper bound imposed on \eqn{\mathrm{df}_U}.}
#'     \item{\code{dfU.cap.rate}}{Rate used to compute the cap.}
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
#'     \item{\code{r.squared}}{Pearson \eqn{\mathrm{cor}(Y, X(\Theta A + U))^2}
#'       (BLUP prediction).}
#'     \item{\code{r.squared.uncentered}}{Uncentered
#'       \eqn{1 - \|Y - X(\Theta A + U)\|_F^2 / \|Y\|_F^2} (BLUP;
#'       baseline = zero matrix).}
#'     \item{\code{r.squared.centered}}{Row-mean centered
#'       \eqn{1 - \|Y - X(\Theta A + U)\|_F^2 / \|Y - \bar Y_{p\cdot}\|_F^2}
#'       (BLUP; baseline = per-row mean).}
#'     \item{\code{r.squared.fixed}}{Pearson \eqn{\mathrm{cor}(Y, X\Theta A)^2}
#'       (fixed-only prediction).}
#'     \item{\code{r.squared.fixed.uncentered}, \code{r.squared.fixed.centered}}{Uncentered
#'       and centered \eqn{R^2} for the fixed-only prediction.}
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
#'   \strong{Sign convention}
#'   \describe{
#'     \item{\code{C.signed}}{Logical. Whether \eqn{C} was sign-free (\code{TRUE}) or non-negative (\code{FALSE}).}
#'   }
#'
#'   Standard errors, z-values, p-values, and confidence intervals for
#'   \eqn{\Theta} are \strong{not} computed here; obtain them by passing the fit
#'   to \code{\link{nmfre.inference}}.
#' @seealso \code{\link{nmfre.inference}},
#'   \code{\link{nmfkc.DOT}}, \code{\link{summary.nmfre}}
#' @export
#' @references
#' Satoh, K. (2026). Wild Bootstrap Inference for Non-Negative Matrix
#'   Factorization with Random Effects. arXiv:2603.01468.
#'   \url{https://arxiv.org/abs/2603.01468}
#' @examples
#' # Example 1. cars data
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(intercept = 1, speed = cars$speed)
#' res <- nmfre(Y, A, rank = 1, maxit = 5000)
#' summary(res)
#'
#' \donttest{
#' # Example 2. Orthodont data (nlme)
#' if (requireNamespace("nlme", quietly = TRUE)) {
#'   Y <- matrix(nlme::Orthodont$distance, 4, 27)
#'   male <- ifelse(nlme::Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
#'   A <- rbind(intercept = 1, male = male)
#'
#'   # Fit (sign-free Theta by default; variances estimated by EM/ECM)
#'   res <- nmfre(Y, A, rank = 1)
#'   summary(res)
#' }}
#'
nmfre <- function(Y, A = NULL, rank = 2, C.signed = TRUE,
                  epsilon = 1e-5, maxit = 5000, ...) {

  extra_args <- base::list(...)
  # backward compatibility
  if (!is.null(extra_args$Q)) rank <- extra_args$Q
  Q <- rank
  ## legacy dfU.cap.rate via ... is tolerated but inert (no cap is applied)
  dfU.cap.rate <- if (!is.null(extra_args$dfU.cap.rate)) extra_args$dfU.cap.rate else NULL

  # --- Parameter Extraction from ... ---
  # initialization
  X.init    <- if (!is.null(extra_args$X.init)) extra_args$X.init else NULL
  C.init    <- if (!is.null(extra_args$C.init)) extra_args$C.init else NULL
  U.init    <- if (!is.null(extra_args$U.init)) extra_args$U.init else NULL
  ## Multi-start for the nmfkc() initialisation (k-means nstart).
  ## Default 1 keeps the historical single-start behaviour.
  nstart    <- if (!is.null(extra_args$nstart)) extra_args$nstart else 1
  prefix    <- if (!is.null(extra_args$prefix)) extra_args$prefix else "Basis"
  ## Sign constraint on C (= Theta in the paper). C.signed = TRUE (default)
  ## gives real-valued, sign-free covariate effects (the semi-NMF reading);
  ## C.signed = FALSE constrains C >= 0 (compositional/intensity scores).
  ## A character C.signed ("signed"/"nonneg") is also accepted (back-compat).
  ## C.mode is the internal string used throughout the optimizer.
  C.mode <- if (is.character(C.signed)) base::match.arg(C.signed, c("signed", "nonneg"))
            else if (isTRUE(C.signed)) "signed" else "nonneg"
  ## The X-step rule follows the sign convention of C directly: sign-free (LS) C
  ## uses the complete-EM semi-NMF step (Ding-Li-Jordan 2010); non-negative (MU)
  ## C uses the legacy positive-part multiplicative update.
  ## Complete-EM X-step: add the posterior-variance term tr(X S X'),
  ## S = N*sigma2*(X'X+lambda I)^{-1}; FALSE recovers the conditional X-step.
  x.postvar <- if (!is.null(extra_args$x.postvar)) isTRUE(extra_args$x.postvar) else TRUE

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

  # dfU control (DEPRECATED: no cap is applied -- the variance components are
  # estimated by ECM and df_U/(NQ) is reported only as a complexity diagnostic).
  # The argument is accepted (inert) for backward compatibility; "cap" behaves
  # identically to "off".
  dfU.control  <- if (!is.null(extra_args$dfU.control)) extra_args$dfU.control else "off"

  # dfU internals (fixed at defaults)
  dfU.cap <- NULL
  dfU.enforce.every <- 1
  dfU.lambda.max <- 1e12

  # optimization extras
  print.trace <- if (!is.null(extra_args$print.trace)) extra_args$print.trace else FALSE
  seed        <- if (!is.null(extra_args$seed)) extra_args$seed else 1

  ## NOTE: nmfre() performs optimization only. Hypothesis tests / standard errors
  ## for Theta (= C) are obtained separately via nmfre.inference(), mirroring the
  ## nmfkc() / nmfkc.inference() split.

  set.seed(seed)
  .eps <- 1e-10

  dfU.control <- match.arg(dfU.control, choices = c("cap", "off"))

  # ---- dimensions ----
  P <- nrow(Y); N <- ncol(Y)
  if (is.null(A)) A <- matrix(1, 1, N)
  if (!is.matrix(A)) A <- as.matrix(A)
  stopifnot(ncol(A) == N)
  K <- nrow(A)

  # ---- init via nmfkc when X.init/C.init are NULL ----
  ## X.init may be NULL (default), a character init-method name forwarded to
  ## nmfkc() (e.g. "kmeans", "runif", "nndsvd"), or a numeric basis matrix
  ## used as-is (with C then estimated given that fixed X).  A character
  ## X.init previously fell through unresolved and crashed downstream.
  if (is.null(X.init) || is.null(C.init)) {
    init_args <- list(Y, A, Q = Q, epsilon = epsilon, seed = seed,
                      nstart = nstart, print.trace = print.trace,
                      print.dims = FALSE)
    if (is.character(X.init)) {
      init_args$X.init <- X.init                  # forward the named method
    } else if (is.matrix(X.init) || is.data.frame(X.init)) {
      init_args$X.init <- as.matrix(X.init)       # fix the supplied basis
      init_args$X.restriction <- "fixed"
    }
    res0 <- do.call(nmfkc, init_args)
    if (is.null(X.init) || is.character(X.init)) X.init <- res0$X
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

  # (dfU cap removed; dfU.cap.rate is inert, kept only for the diagnostic frac)
  if (is.null(dfU.cap.rate)) dfU.cap.rate <- 0.10

  ## outer-inner ECM internals (no user knobs; df_U is uncapped)
  inner.maxit   <- if (!is.null(extra_args$inner.maxit)) extra_args$inner.maxit else 10000L
  outer.maxit   <- if (!is.null(extra_args$outer.maxit)) extra_args$outer.maxit else 500L
  epsilon.outer <- if (!is.null(extra_args$epsilon.outer)) extra_args$epsilon.outer else 1e-6

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
  nll_trace <- rep(NA_real_, maxit)   # marginal NLL (ECM-monotone; see plot.nmfre)

  # ---- dfU diagnostics ----
  dfU_last <- NA_real_
  dfU_cap_last <- NA_real_
  lambda_last <- NA_real_

  dfU_hit_cap <- FALSE
  dfU_hit_iter <- integer(0)

  # =========================================================
  # main loop: outer-inner ECM
  #   inner = fixed-lambda block-coordinate descent (U, X, C); monotone
  #   outer = EM M-steps for (sigma2, tau2); repeat until lambda stabilizes
  # =========================================================
  any_var_update <- isTRUE(sigma2.update) || isTRUE(tau2.update)
  total_iter <- 0L
  lambda     <- sigma2 / tau2

  for (outer in 1:outer.maxit) {

    tau2   <- clip_val(tau2,   tau2.min,   tau2.max)
    sigma2 <- clip_val(sigma2, sigma2.min, sigma2.max)
    lambda <- sigma2 / tau2
    lambda_last <- lambda

    # ----- INNER: fixed-lambda block-coordinate descent -----
    obj_prev_inner <- Inf
    inner_done <- 0L
    for (inner in 1:inner.maxit) {
      total_iter <- total_iter + 1L
      iter_done  <- total_iter
      inner_done <- inner

      CA <- C_mat %*% A

      # (1) U-step: ridge / BLUP at the fixed lambda (E-step posterior mean)
      XtX   <- crossprod(X)
      M     <- XtX + diag(pmax(lambda, 1e-12), Q)
      cholM <- tryCatch(chol(M), error = function(e) NULL)
      if (is.null(cholM)) stop("Cholesky failed: XtX + lambda I not SPD. Increase lambda (or decrease tau2).")
      R0 <- Y - X %*% CA
      for (n in 1:N) {
        rhs <- crossprod(X, R0[, n, drop = FALSE])
        U[, n] <- as.numeric(backsolve(cholM, forwardsolve(t(cholM), rhs)))
      }
      U <- sweep(U, 1, rowMeans(U), "-")          # center for identifiability

      # (2) X-step: X >= 0 with sign-free score matrix B = CA + U
      B_sem <- CA + U
      if (C.mode == "signed") {
        ## complete-EM X-step: add posterior-variance S = N*sigma2*(X'X+lambda I)^{-1}
        S_pv <- if (isTRUE(x.postvar))
          N * clip_val(sigma2, sigma2.min, sigma2.max) * chol2inv(cholM) else NULL
        X <- .nmfre.seminmf.X(X, Y, B_sem, S = S_pv, eps = .eps)
      } else {
        Y_tilde <- pmax(Y, 0); B_pos <- pmax(B_sem, 0)
        numX <- Y_tilde %*% t(B_pos)
        denX <- X %*% (B_pos %*% t(B_pos))
        X <- X * .nmfre.safe.div(numX, denX, eps = .eps)
      }
      X <- pmax(X, .eps)
      normed <- .nmfre.normalize.X(X, C_mat, U)
      X     <- pmax(normed$X, .eps)
      C_mat <- if (C.mode == "nonneg") pmax(normed$C, .eps) else normed$C
      U     <- normed$U

      # (3) C-step
      if (C.mode == "signed") {
        ## exact least squares (sign-free): C = (X'X)^{-1} X'(Y - X U) A'(A A')^{-1}
        Y_star <- Y - X %*% U
        XtX_t  <- crossprod(X) + diag(1e-10, Q)
        AAt_t  <- tcrossprod(A) + diag(1e-10, K)
        rhs    <- crossprod(X, Y_star) %*% t(A)
        C_mat  <- solve(XtX_t, rhs)
        C_mat  <- t(solve(AAt_t, t(C_mat)))
      } else {
        ## non-negative MU (Ding 2006), positive-part stabilization
        Y_tilde <- pmax(Y, 0); Y_star <- Y_tilde - X %*% U; Y_star_pos <- pmax(Y_star, 0)
        numC <- crossprod(X, Y_star_pos) %*% t(A)
        denC <- (crossprod(X, X) %*% C_mat) %*% (A %*% t(A))
        C_mat <- C_mat * .nmfre.safe.div(numC, denC, eps = .eps)
        C_mat <- pmax(C_mat, .eps)
      }

      # (4) objective at the FIXED lambda (what the inner loop minimizes)
      CA  <- C_mat %*% A
      R   <- Y - X %*% (CA + U)
      obj <- sum(R^2) + lambda * sum(U^2)
      if (total_iter <= maxit) {
        obj_trace[total_iter] <- obj
        rss_trace[total_iter] <- sum(R^2)
        ## marginal NLL (U integrated out): the ECM-monotone diagnostic
        nll_trace[total_iter] <- .nmfre.marginal.nll(Y, X, C_mat, A, sigma2, tau2)
      }

      if (!is.finite(obj)) { stop_reason <- "nonfinite_obj"; converged <- FALSE; break }
      rel_change <- if (is.finite(obj_prev_inner))
        abs(obj_prev_inner - obj) / (abs(obj_prev_inner) + .eps) else NA_real_
      if (is.finite(obj_prev_inner) && rel_change < epsilon) break   # inner converged
      obj_prev_inner <- obj
      if (total_iter >= maxit) break
    }

    if (identical(stop_reason, "nonfinite_obj")) break

    # ----- OUTER: EM M-steps for (sigma2, tau2) from the converged inner fit -----
    CA <- C_mat %*% A
    R  <- Y - X %*% (CA + U)
    sigma2_old <- sigma2
    Minv <- if (isTRUE(sigma2.update) || isTRUE(tau2.update))
      tryCatch(chol2inv(chol(crossprod(X) + diag(pmax(lambda, 1e-12), Q))),
               error = function(e) NULL) else NULL
    if (isTRUE(sigma2.update)) {
      if (!is.null(Minv)) {
        trH   <- Q - pmax(lambda, 1e-12) * sum(diag(Minv))   # tr(H_lambda) in [0,Q]
        sigma2 <- clip_val(mean(R^2) + sigma2_old * (N * trH) / (P * N), sigma2.min, sigma2.max)
      } else sigma2 <- clip_val(mean(R^2), sigma2.min, sigma2.max)
    }
    if (isTRUE(tau2.update) && !is.null(Minv)) {
      tau2 <- clip_val(sum(U^2) / (N * Q) + sigma2_old * sum(diag(Minv)) / Q, tau2.min, tau2.max)
    }

    lambda_new <- clip_val(sigma2, sigma2.min, sigma2.max) / clip_val(tau2, tau2.min, tau2.max)

    if (isTRUE(print.trace)) {
      message(sprintf("[outer %d] inner=%d obj=%.6g sigma2=%.4g tau2=%.4g lambda %.4g -> %.4g",
                      outer, inner_done, obj, sigma2, tau2, lambda, lambda_new))
    }

    # ----- outer convergence -----
    if (!any_var_update) { converged <- TRUE; stop_reason <- "fixed_lambda"; break }
    if (abs(lambda_new - lambda) / (abs(lambda) + .eps) < epsilon.outer) {
      converged <- TRUE; stop_reason <- "outer_lambda"; break
    }
    if (total_iter >= maxit) { stop_reason <- "maxit"; converged <- FALSE; break }
  }

  # ---- post loop: finalize ----
  obj_final <- obj
  rel_change_final <- rel_change

  if (iter_done < maxit && identical(stop_reason, "maxit")) {
    stop_reason <- "unknown_break"
    converged <- FALSE
  }

  ## Warn when the MU loop exhausted maxit without meeting the
  ## relative-tolerance criterion (matches nmfkc() / nmf.sem()
  ## convention).  The non-finite-objective break path issues its
  ## own diagnostic via stop_reason; here we only fire on maxit.
  if (identical(stop_reason, "maxit") && iter_done == maxit)
    warning(paste0("maximum iterations (", maxit, ") reached..."))

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

  dfU_cap_last <- NA_real_   # no cap is applied (df_U is a diagnostic only)
  dfU_last <- dfU_final
  lambda_last <- lambda_final

  # ---- names ----
  colnames(X) <- paste0(prefix, 1:ncol(X))
  rownames(C_mat) <- colnames(X)
  colnames(C_mat) <- rownames(A)
  rownames(U) <- rownames(C_mat)
  colnames(U) <- colnames(Y)

  # ---- Three-variant R^2 (cor^2, uncentered, row-mean
  # centered), separately for BLUP-prediction (XB.blup = XB + U) and
  # fixed-only prediction (XB).  Consistent with nmfkc / nmfae /
  # nmfae.signed / nmfkc.net / nmfkc.signed naming convention.
  r2_blup  <- .r.squared.all(Y, XB.blup)
  r2_fixed <- .r.squared.all(Y, XB)
  r.squared                <- r2_blup$r.squared
  r.squared.uncentered           <- r2_blup$r.squared.uncentered
  r.squared.centered       <- r2_blup$r.squared.centered
  r.squared.fixed          <- r2_fixed$r.squared
  r.squared.fixed.uncentered     <- r2_fixed$r.squared.uncentered
  r.squared.fixed.centered <- r2_fixed$r.squared.centered

  # ---- ICC (trace-based) ----
  trXtX <- sum(d_final)  # tr(X'X), using eigenvalues already computed
  ICC <- tau2 * trXtX / (tau2 * trXtX + sigma2 * P)

  # ---- convenience probabilities ----
  colnorm_prob <- function(M, eps = 1e-12) {
    cs <- colSums(M)
    sweep(M, 2, pmax(cs, eps), "/")
  }

  ## row-wise soft clustering of the (non-negative) basis X, mirroring nmfkc():
  ## each row of X is normalized to sum to 1 and assigned to its argmax factor.
  X.prob <- X / (base::rowSums(X) + .eps)
  X.cluster <- base::apply(X.prob, 1, base::which.max)
  X.cluster[base::rowSums(X) == 0] <- NA

  out <- list(
    X = X,
    X.prob = X.prob,
    X.cluster = X.cluster,
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
    nll.trace = nll_trace[seq_len(iter_done)],

    # dfU diagnostics
    dfU = dfU_last,
    dfU.cap = dfU_cap_last,
    dfU.cap.rate = dfU.cap.rate,
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
    r.squared                = r.squared,
    r.squared.uncentered           = r.squared.uncentered,
    r.squared.centered       = r.squared.centered,
    r.squared.fixed          = r.squared.fixed,
    r.squared.fixed.uncentered     = r.squared.fixed.uncentered,
    r.squared.fixed.centered = r.squared.fixed.centered,
    ICC = ICC,

    ## sign convention of C (used by nmfre.inference() to pick the test side)
    C.signed = (C.mode == "signed")
  )

  class(out) <- c("nmfre", "nmf")
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
#' @seealso \code{\link{nmfre}}, \code{\link{nmfre.inference}}
#' @export
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(intercept = 1, speed = cars$speed)
#' res <- nmfre(Y, A, rank = 1, maxit = 5000)
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
    cat(sprintf("R-squared (cor^2):    %.4f (XB+blup)%s\n", x$r.squared, r2_fixed))
    if (!is.null(x$r.squared.uncentered) && is.finite(x$r.squared.uncentered)) {
      r2f_fixed <- if (!is.null(x$r.squared.fixed.uncentered) && is.finite(x$r.squared.fixed.uncentered))
        sprintf(", %.4f (XB)", x$r.squared.fixed.uncentered) else ""
      cat(sprintf("R-squared (uncentered):     %.4f (XB+blup)%s\n", x$r.squared.uncentered, r2f_fixed))
    }
    if (!is.null(x$r.squared.centered) && is.finite(x$r.squared.centered)) {
      r2c_fixed <- if (!is.null(x$r.squared.fixed.centered) && is.finite(x$r.squared.fixed.centered))
        sprintf(", %.4f (XB)", x$r.squared.fixed.centered) else ""
      cat(sprintf("R-squared (centered): %.4f (XB+blup)%s\n", x$r.squared.centered, r2c_fixed))
    }
  }

  # ---- effective rank (BLUP scores B = Theta A + U, Q x N) ----
  eff <- if (!is.null(x$B.blup)) .effective.rank(x$B.blup) else NA_real_
  if (is.finite(eff)) cat(sprintf("Effective Rank:       %.2f / %d  (%.1f%%)\n",
                                  eff, Q, 100 * eff / Q))

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
  if (is.null(x$coefficients) || !is.data.frame(x$coefficients)) {
    cat("\nCoefficients (Theta): run nmfre.inference(fit, Y, A) for SE / p-values.\n")
  }
  if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
    cat("\nCoefficients:\n")

    cf <- x$coefficients

    # row names: Covariate:Basis
    rnames <- paste0(cf$Covariate, ":", cf$Basis)

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
    ## sign convention from the fit: logical C.signed preferred; tolerate the
    ## legacy character C.signed and the short-lived logical C.nonnegative.
    C.nn <- if (!is.null(x$C.signed)) {
              if (is.logical(x$C.signed)) !isTRUE(x$C.signed) else identical(x$C.signed, "nonneg")
            } else if (!is.null(x$C.nonnegative)) isTRUE(x$C.nonnegative) else NA
    if (!is.na(C.nn)) {
      conv <- if (!C.nn) "sign-free (real-valued)" else "non-negative"
      cat(sprintf("C (= Theta) update: %s; p-values %s\n",
                  conv, if (!is.null(x$C.p.side)) x$C.p.side else "one.sided"))
    }

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
# nmfre.inference() - statistical inference on C (separated)
# ============================================================

#' @title Statistical inference for the coefficient matrix C from NMF-RE
#'
#' @description
#' \code{nmfre.inference} performs statistical inference on the coefficient
#' matrix \eqn{C} (\eqn{\Theta}) from a fitted \code{nmfre} model,
#' conditional on the estimated basis matrix \eqn{\hat{X}} and random
#' effects \eqn{\hat{U}}.
#'
#' Under the working model \eqn{Y^* = Y - X\hat{U} \approx X C A + \varepsilon},
#' inference is conducted via sandwich covariance estimation and
#' one-step wild bootstrap with non-negative projection.
#'
#' The result is compatible with \code{\link{nmfkc.DOT}} for visualization
#' (pass the result directly as \code{x} with \code{type = "YXA"}).
#'
#' @param object An object of class \code{"nmfre"} returned by
#'   \code{\link{nmfre}}.
#' @param Y Observation matrix (P x N). Must match the data used in
#'   \code{nmfre()}.
#' @param A Covariate matrix (K x N). Default is \code{NULL} (intercept only).
#' @param wild.bootstrap Logical. If \code{TRUE} (default), performs wild
#'   bootstrap for confidence intervals and bootstrap standard errors.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{wild.B}}{Number of bootstrap replicates. Default is 500.}
#'     \item{\code{wild.seed}}{Seed for bootstrap. Default is 123.}
#'     \item{\code{wild.level}}{Confidence level for bootstrap CI. Default is 0.95.}
#'     \item{\code{C.p.side}}{P-value type: \code{"one.sided"} (default) or \code{"two.sided"}.}
#'     \item{\code{cov.ridge}}{Ridge stabilization. Default is 1e-8.}
#'     \item{\code{print.trace}}{Logical. Default is \code{FALSE}.}
#'   }
#'
#' @return The input \code{object} with additional inference components:
#' \item{sigma2.used}{Estimated \eqn{\sigma^2} used for inference.}
#' \item{C.vec.cov}{Full covariance matrix for \eqn{vec(C)}.}
#' \item{C.se}{Sandwich standard errors for \eqn{C}.}
#' \item{C.se.boot}{Bootstrap standard errors for \eqn{C}.}
#' \item{C.ci.lower}{Lower CI bounds for \eqn{C}.}
#' \item{C.ci.upper}{Upper CI bounds for \eqn{C}.}
#' \item{coefficients}{Data frame with Basis, Covariate, Estimate, SE, BSE,
#'   z_value, p_value, CI_low, CI_high.}
#' \item{C.p.side}{P-value type used.}
#'
#' @seealso \code{\link{nmfre}}, \code{\link{nmfkc.DOT}},
#'   \code{\link{summary.nmfre}}
#' @references
#' Satoh, K. (2026). Wild Bootstrap Inference for Non-Negative Matrix
#'   Factorization with Random Effects. arXiv:2603.01468.
#'   \url{https://arxiv.org/abs/2603.01468}
#' @export
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(intercept = 1, speed = cars$speed)
#' res <- nmfre(Y, A, rank = 1, wild.bootstrap = FALSE)
#' res2 <- nmfre.inference(res, Y, A)
#' res2$coefficients
#'
nmfre.inference <- function(object, Y, A = NULL, wild.bootstrap = TRUE, ...) {
  if (!inherits(object, "nmfre"))
    stop("object must be of class 'nmfre' (returned by nmfre).")

  extra_args <- base::list(...)
  wild.B      <- if (!is.null(extra_args$wild.B))      extra_args$wild.B      else 500
  wild.seed   <- if (!is.null(extra_args$wild.seed))   extra_args$wild.seed   else 123
  wild.level  <- if (!is.null(extra_args$wild.level))  extra_args$wild.level  else 0.95
  ## sign convention of C (= Theta) from the fit; resolves the default p-side
  ## (two-sided for sign-free C, one-sided for the non-negative variant) and
  ## whether the bootstrap replicates are projected onto C >= 0. Logical
  ## C.signed preferred; tolerate the legacy character C.signed and the
  ## short-lived logical C.nonnegative. C.mode is the internal string.
  C.mode      <- if (!is.null(object$C.signed)) {
                   if (is.logical(object$C.signed)) (if (isTRUE(object$C.signed)) "signed" else "nonneg")
                   else base::match.arg(object$C.signed, c("signed", "nonneg"))
                 } else if (!is.null(object$C.nonnegative)) {
                   if (isTRUE(object$C.nonnegative)) "nonneg" else "signed"
                 } else "nonneg"
  C.p.side    <- if (!is.null(extra_args$C.p.side))    extra_args$C.p.side
                 else if (C.mode == "nonneg") "one.sided" else "two.sided"
  cov.ridge   <- if (!is.null(extra_args$cov.ridge))   extra_args$cov.ridge   else 1e-8
  print.trace <- if (!is.null(extra_args$print.trace)) extra_args$print.trace else FALSE

  Y <- base::as.matrix(Y)
  if (is.null(A)) A <- base::matrix(1, nrow = 1, ncol = base::ncol(Y))
  A <- base::as.matrix(A)

  X     <- object$X       # P x Q
  C_mat <- object$C       # Q x K
  U     <- object$U       # Q x N
  sigma2 <- object$sigma2
  tau2   <- object$tau2

  P <- base::nrow(Y)
  N <- base::ncol(Y)
  Q <- base::ncol(X)
  K <- base::nrow(A)

  lambda_inf <- sigma2 / base::pmax(tau2, 1e-12)

  # Y_star = Y - XU (remove random effects)
  Y_star_inf <- Y - X %*% U
  R_C <- Y_star_inf - X %*% (C_mat %*% A)
  RSS_inf <- base::sum(R_C^2)

  XtX_now <- base::crossprod(X)       # Q x Q
  AAt <- A %*% base::t(A)             # K x K

  M_inf <- XtX_now + base::diag(base::pmax(lambda_inf, 1e-12), Q)
  cholM_inf <- base::tryCatch(base::chol(M_inf), error = function(e) NULL)
  if (!is.null(cholM_inf)) {
    Minv <- base::chol2inv(cholM_inf)
  } else {
    Minv <- base::tryCatch(base::solve(M_inf), error = function(e) {
      if (base::requireNamespace("MASS", quietly = TRUE)) MASS::ginv(M_inf)
      else base::stop("Information matrix singular; install MASS package.")
    })
  }

  trH <- base::sum(base::diag(Minv %*% XtX_now))
  dfU_inf <- N * trH
  dfC <- Q * K
  denom <- base::max(P * N - dfU_inf - dfC, 1)
  sigma2_used <- RSS_inf / denom

  # Information matrix with (I-H) correction
  Xt_IH_X <- XtX_now - XtX_now %*% Minv %*% XtX_now
  Info_core <- base::kronecker(AAt, Xt_IH_X)
  Info <- Info_core / base::max(sigma2_used, 1e-12)
  Info <- Info + base::diag(cov.ridge, base::nrow(Info))

  Hinv <- base::tryCatch(base::solve(Info), error = function(e) {
    if (base::requireNamespace("MASS", quietly = TRUE)) MASS::ginv(Info)
    else base::stop("Information matrix singular; install MASS package.")
  })

  # Sandwich covariance
  V_sand <- NULL
  Xt <- base::t(X)
  J <- base::matrix(0, Q * K, Q * K)
  for (n in 1:N) {
    a_n <- A[, n, drop = FALSE]
    r_n <- R_C[, n, drop = FALSE]
    g_n <- Xt %*% r_n
    S_n <- -(g_n %*% base::t(a_n)) / base::max(sigma2_used, 1e-12)
    s_n <- base::as.vector(S_n)
    J <- J + base::tcrossprod(s_n)
  }
  if (N > 1) J <- (N / (N - 1)) * J    # CR1 correction
  V_sand <- Hinv %*% J %*% Hinv

  C.vec.cov <- V_sand
  se_vec <- base::sqrt(base::pmax(base::diag(C.vec.cov), 0))
  C.se <- base::matrix(se_vec, nrow = Q, ncol = K, byrow = FALSE)
  base::dimnames(C.se) <- base::dimnames(C_mat)

  # ---- Wild bootstrap (one-step Newton) ----
  C.se.boot  <- NULL
  C.ci.lower <- NULL
  C.ci.upper <- NULL

  if (base::isTRUE(wild.bootstrap)) {
    base::set.seed(wild.seed)
    score_mat <- base::matrix(0, Q * K, N)
    for (n in 1:N) {
      a_n <- A[, n, drop = FALSE]
      r_n <- R_C[, n, drop = FALSE]
      g_n <- Xt %*% r_n
      G_n <- -(g_n %*% base::t(a_n)) / base::max(sigma2_used, 1e-12)
      score_mat[, n] <- base::as.vector(G_n)
    }

    ## project to C >= 0 only for the non-negative variant; sign-free C is interior.
    C_boot <- .boot.onestep(base::as.vector(C_mat), score_mat, Hinv, wild.B,
                            dist = "exp", seed = wild.seed,
                            project = (C.mode == "nonneg"))

    alpha <- 1 - wild.level
    lo <- base::apply(C_boot, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE)
    hi <- base::apply(C_boot, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)

    C.ci.lower <- base::matrix(lo, nrow = Q, ncol = K, byrow = FALSE)
    C.ci.upper <- base::matrix(hi, nrow = Q, ncol = K, byrow = FALSE)
    base::dimnames(C.ci.lower) <- base::dimnames(C_mat)
    base::dimnames(C.ci.upper) <- base::dimnames(C_mat)

    sd_vec <- base::apply(C_boot, 1, stats::sd, na.rm = TRUE)
    C.se.boot <- base::matrix(sd_vec, nrow = Q, ncol = K, byrow = FALSE)
    base::dimnames(C.se.boot) <- base::dimnames(C_mat)
  }

  # ---- Coefficients table (nmfkc.DOT compatible: Basis/Covariate) ----
  Estimate <- base::as.vector(C_mat)
  SE <- base::as.vector(C.se)
  BSE <- if (!is.null(C.se.boot)) base::as.vector(C.se.boot) else base::rep(NA_real_, base::length(Estimate))
  z_value <- base::ifelse(SE > 0, Estimate / SE, NA_real_)

  if (C.p.side == "one.sided") {
    p_value <- base::ifelse(base::is.finite(z_value), stats::pnorm(z_value, lower.tail = FALSE), NA_real_)
  } else {
    p_value <- base::ifelse(base::is.finite(z_value), 1 - stats::pchisq(z_value^2, df = 1), NA_real_)
  }

  rlabs <- if (!is.null(base::rownames(C_mat))) base::rownames(C_mat) else base::paste0("Basis", 1:Q)
  clabs <- if (!is.null(base::colnames(C_mat))) base::colnames(C_mat) else base::paste0("Cov", 1:K)

  coefficients <- base::data.frame(
    Basis     = base::rep(rlabs, times = K),
    Covariate = base::rep(clabs, each = Q),
    Estimate  = Estimate,
    SE        = SE,
    BSE       = BSE,
    z_value   = z_value,
    p_value   = p_value,
    CI_low    = if (!is.null(C.ci.lower)) base::as.vector(C.ci.lower) else NA_real_,
    CI_high   = if (!is.null(C.ci.upper)) base::as.vector(C.ci.upper) else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )

  if (print.trace) {
    base::message("  nmfre.inference: sandwich SE + wild bootstrap done.")
  }

  object$sigma2.used  <- sigma2_used
  object$C.vec.cov    <- C.vec.cov
  object$C.se         <- C.se
  object$C.se.hess    <- C.se
  object$C.se.boot    <- C.se.boot
  object$C.ci.lower   <- C.ci.lower
  object$C.ci.upper   <- C.ci.upper
  object$C.boot.sd    <- C.se.boot
  object$coefficients <- coefficients
  object$C.p.side     <- C.p.side
  return(object)
}

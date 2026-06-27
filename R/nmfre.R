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
#' @param df.rate Rate for computing the dfU cap (\code{cap = rate * N * Q}).
#'   For backward compatibility, \code{dfU.cap.rate} is accepted via \code{...}.
#'   If \code{NULL} (default), runs \code{\link{nmfre.dfU.scan}} internally and
#'   selects the minimum rate where the cap is not binding.
#'   Use \code{\link{nmfre.dfU.scan}} beforehand to examine dfU behavior across rates
#'   and choose an appropriate value.
#' @param wild.bootstrap Logical. If \code{TRUE} (default), perform wild bootstrap inference
#'   on \eqn{\Theta}.
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
#'     \item \code{C.p.side}: P-value sidedness. Defaults to \code{"two.sided"}
#'       when \code{C.signed = TRUE} (interior null H0: C=0) and
#'       \code{"one.sided"} when \code{C.signed = FALSE} (boundary null
#'       H0: C=0 vs H1: C>0). Override explicitly if desired.
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
#'     \item{\code{C.signed}}{Logical. Whether \eqn{C} was sign-free (\code{TRUE}) or non-negative (\code{FALSE}).}
#'     \item{\code{x.update}}{Basis update rule used: \code{"seminmf"} or \code{"mu"}.}
#'   }
#' @seealso \code{\link{nmfre.inference}}, \code{\link{nmfre.dfU.scan}},
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
#'   # Scan dfU cap rates to choose an appropriate value
#'   nmfre.dfU.scan(1:10/10, Y, A, rank = 1)
#'
#'   # Fit with chosen rate
#'   res <- nmfre(Y, A, rank = 1, df.rate = 0.2)
#'   summary(res)
#' }}
#'
nmfre <- function(Y, A = NULL, rank = 2, C.signed = TRUE, df.rate = NULL,
                  wild.bootstrap = TRUE, epsilon = 1e-5,
                  maxit = 5000, ...) {

  extra_args <- base::list(...)
  # backward compatibility
  if (!is.null(extra_args$Q)) rank <- extra_args$Q
  if (!is.null(extra_args$dfU.cap.rate)) df.rate <- extra_args$dfU.cap.rate
  Q <- rank
  dfU.cap.rate <- df.rate

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
  ## X-step rule is selected by the sign convention of C: the sign-free (LS) C
  ## pairs with the complete-EM semi-NMF step (Ding-Li-Jordan 2010), and the
  ## non-negative (MU) C pairs with the legacy positive-part multiplicative
  ## update.  An explicit x.update in ... overrides this default (advanced use).
  x.update  <- if (!is.null(extra_args$x.update)) extra_args$x.update
               else if (C.mode == "nonneg") "mu" else "seminmf"
  x.update  <- base::match.arg(x.update, c("seminmf", "mu"))
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

  # inference settings
  ## C.p.side default resolved from C.mode: two-sided for the sign-free
  ## default (C is an interior parameter), one-sided for the non-negative
  ## variant (boundary null H0: C = 0 vs H1: C > 0).
  C.p.side <- if (!is.null(extra_args$C.p.side)) extra_args$C.p.side
              else if (C.mode == "nonneg") "one.sided" else "two.sided"

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
  dfU.cap.scan <- NULL
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
      if (x.update == "seminmf") {
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
      if (total_iter <= maxit) { obj_trace[total_iter] <- obj; rss_trace[total_iter] <- sum(R^2) }

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

    ## project to C >= 0 only for the non-negative variant; sign-free C is interior.
    C_boot <- .boot.onestep(as.vector(C_mat), score_mat, Hinv, wild.B,
                            dist = "exp", seed = wild.seed,
                            project = (C.mode == "nonneg"))

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
      Basis     = rep(rownames(C_mat), times = K),
      Covariate = rep(colnames(C_mat), each = Q),
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
    r.squared                = r.squared,
    r.squared.uncentered           = r.squared.uncentered,
    r.squared.centered       = r.squared.centered,
    r.squared.fixed          = r.squared.fixed,
    r.squared.fixed.uncentered     = r.squared.fixed.uncentered,
    r.squared.fixed.centered = r.squared.fixed.centered,
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

    C.p.side = C.p.side,
    C.signed = (C.mode == "signed"),
    x.update = x.update
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
#' @param rank Integer. Rank of the basis matrix.
#'   For backward compatibility, \code{Q} is accepted via \code{...}.
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
#' @seealso \code{\link{nmfre}}
#' @export
#' @examples
#' # Example 1. cars data (small maxit for speed)
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(intercept = 1, speed = cars$speed)
#' tab <- nmfre.dfU.scan(rates = c(0.1, 0.2), Y = Y, A = A, rank = 1, maxit = 1000)
#' print(tab)
#'
#' \donttest{
#' # Example 2. Orthodont data (nlme)
#' if (requireNamespace("nlme", quietly = TRUE)) {
#'   Y <- matrix(nlme::Orthodont$distance, 4, 27)
#'   male <- ifelse(nlme::Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
#'   A <- rbind(intercept = 1, male = male)
#'   nmfre.dfU.scan(1:10/10, Y, A, rank = 1)
#' }}
#'
nmfre.dfU.scan <- function(
  rates = (1:10) / 10,
  Y, A, rank = NULL,
  X.init = NULL, C.init = NULL, U.init = NULL,
  print.trace = FALSE,
  ...
) {
  extra_scan <- base::list(...)
  if (!is.null(extra_scan$Q)) rank <- extra_scan$Q
  Q <- rank
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

#' @seealso \code{\link{nmfre.dfU.scan}}
#' @export
print.nmfre.dfU.scan <- function(x, ...) {
  print(x$table, ...)
  invisible(x)
}

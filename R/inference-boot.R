# ============================================================
# Shared wild (multiplier) bootstrap engine for *.inference()
# ------------------------------------------------------------
# Two flavours, both wild multiplier bootstraps, distinguished by the
# UPDATE scheme (not by "wild"-ness, which they share):
#   method = "onestep" : one-step Newton linearisation around the fit,
#                        c_b = C_hat - Hinv (score %*% w).  Needs Hinv.
#   method = "refit"   : perturb the response and RE-ESTIMATE C to
#                        convergence with the basis X held FIXED.  Needs
#                        no information matrix, so it stays valid when the
#                        Fisher information AA' is singular (over-
#                        parameterised covariates) or C sits on the >= 0
#                        boundary.
# The multiplier distribution (wild.dist) is orthogonal to the method.
# ============================================================

#' @title Wild bootstrap multipliers (Internal)
#' @description Draw \code{n} i.i.d. mean-0, variance-1 multipliers for a
#'   wild (multiplier) bootstrap.
#' @param n Number of multipliers.
#' @param dist One of \code{"rademacher"} (\eqn{\pm 1}; keeps the residual
#'   magnitude, flips sign), \code{"mammen"} (Mammen's 1993 two-point,
#'   matches the third moment), or \code{"exp"} (\eqn{\mathrm{Exp}(1)-1};
#'   the historical default of the one-step engine).
#' @return Numeric vector of length \code{n}.
#' @keywords internal
#' @noRd
.wild.multipliers <- function(n, dist = c("rademacher", "mammen", "exp")) {
  dist <- base::match.arg(dist)
  if (dist == "rademacher") {
    base::sample(c(-1, 1), n, replace = TRUE)
  } else if (dist == "exp") {
    stats::rexp(n, rate = 1) - 1
  } else {                                  # Mammen (1993) two-point
    a <- -(base::sqrt(5) - 1) / 2           # ~ -0.618  (low value)
    b <-  (base::sqrt(5) + 1) / 2           # ~  1.618  (high value)
    pa <- (base::sqrt(5) + 1) / (2 * base::sqrt(5))  # P(w = a)
    base::ifelse(stats::runif(n) < pa, a, b)
  }
}


#' @title One-step Newton wild bootstrap draws (Internal)
#' @description Linearised (one-step) wild bootstrap of \eqn{vec(C)} around
#'   the fitted value, using the score matrix and the inverse information
#'   \code{Hinv}.  Reproduces the historical loop
#'   \code{c_b = C_hat - Hinv (score \%*\% w)} with optional non-negative
#'   projection.
#' @param C.hat.vec Numeric vector \eqn{vec(\hat C)} (length \eqn{QK}).
#' @param score.mat \eqn{QK \times N} per-observation score matrix.
#' @param Hinv \eqn{QK \times QK} inverse (or pseudo-inverse) information.
#' @param B Number of replicates.
#' @param dist Multiplier distribution (see \code{.wild.multipliers}).
#'   Default \code{"exp"} for backward compatibility.
#' @param seed RNG seed (set once before the loop).
#' @param project Logical; if \code{TRUE}, clip each draw at 0 (for
#'   non-negative \eqn{C}); \code{FALSE} for signed \eqn{C}.
#' @return \eqn{QK \times B} matrix of bootstrap draws.
#' @keywords internal
#' @noRd
.boot.onestep <- function(C.hat.vec, score.mat, Hinv, B,
                          dist = "exp", seed = 123, project = TRUE) {
  base::set.seed(seed)
  N  <- base::ncol(score.mat)
  QK <- base::length(C.hat.vec)
  C.boot <- base::matrix(NA_real_, nrow = QK, ncol = B)
  for (b in base::seq_len(B)) {
    w      <- .wild.multipliers(N, dist)
    grad_b <- base::as.vector(score.mat %*% w)
    c_b    <- C.hat.vec - base::as.vector(Hinv %*% grad_b)
    if (project) c_b <- base::pmax(c_b, 0)
    C.boot[, b] <- c_b
  }
  C.boot
}


#' @title X-fixed re-estimation of C for the working model (Internal)
#' @description Solve \eqn{\min_{C \ge 0} \lVert Y_p - X C A \rVert_F^2} with
#'   the basis \eqn{X} (and design \eqn{A}) held FIXED, by multiplicative
#'   updates.  With \eqn{X} fixed the problem is convex in \eqn{C}; the
#'   sign-split numerator/denominator keep \eqn{C \ge 0} even when \eqn{A}
#'   has negative entries.  Warm-started from \code{C.init} it converges in
#'   a few iterations.
#' @param Yp \eqn{P \times N} (perturbed) response.
#' @param X \eqn{P \times Q} fixed basis.
#' @param A \eqn{K \times N} fixed design.
#' @param C.init \eqn{Q \times K} warm start.
#' @param maxit,eps Iteration cap and relative-objective tolerance.
#' @return \eqn{Q \times K} re-estimated \eqn{C}.
#' @keywords internal
#' @noRd
.refit.C.MU <- function(Yp, X, A, C.init, maxit = 300, eps = 1e-8) {
  Xt  <- base::t(X)
  XtX <- Xt %*% X                       # Q x Q (>= 0 since X >= 0)
  XtY <- Xt %*% Yp %*% base::t(A)        # Q x K
  AAt <- A %*% base::t(A)                # K x K
  Gp  <- base::pmax(XtY, 0); Gn <- base::pmax(-XtY, 0)
  Mp  <- base::pmax(AAt, 0); Mn <- base::pmax(-AAt, 0)
  small <- 1e-12
  C <- base::pmax(C.init, small)
  obj_prev <- Inf
  for (it in base::seq_len(maxit)) {
    XtXC <- XtX %*% C                    # Q x K
    num  <- Gp + XtXC %*% Mn
    den  <- Gn + XtXC %*% Mp + small
    C <- C * num / den
    if (it %% 10 == 0 || it == maxit) {
      resid <- Yp - X %*% C %*% A
      obj   <- base::sum(resid * resid)
      if (base::is.finite(obj_prev) &&
          base::abs(obj_prev - obj) < eps * base::max(obj_prev, small)) break
      obj_prev <- obj
    }
  }
  C
}


#' @title Re-fit wild bootstrap draws (Internal)
#' @description Residual wild bootstrap that RE-FITS the coefficient on each
#'   replicate.  The response is perturbed as
#'   \eqn{Y^* = \mathrm{clip}_{\ge 0}(\hat Y + W \odot R)} with wild
#'   multipliers \eqn{W}, then \code{refit.fun(Ys)} re-estimates and returns
#'   \eqn{vec(C_b)}.  No information matrix is used.
#' @param Yhat \eqn{P \times N} fitted values.
#' @param R \eqn{P \times N} residuals \eqn{Y - \hat Y}.
#' @param refit.fun Function of one argument (the perturbed \eqn{P \times N}
#'   response) returning a numeric vector \eqn{vec(C_b)}.
#' @param B Number of replicates.
#' @param dist Multiplier distribution. Default \code{"rademacher"}.
#' @param seed RNG seed (set once before the loop).
#' @param unit \code{"element"} (i.i.d. multiplier per matrix cell) or
#'   \code{"column"} (one multiplier per sample column, shared over rows).
#' @param clipY Logical; clip the perturbed response at 0 (NMF needs
#'   \eqn{Y \ge 0}). Default \code{TRUE}.
#' @return \eqn{QK \times B} matrix of bootstrap draws (\code{NA} column if a
#'   replicate fails).
#' @keywords internal
#' @noRd
.boot.refit <- function(Yhat, R, refit.fun, B,
                        dist = "rademacher", seed = 123,
                        unit = c("element", "column"), clipY = TRUE) {
  unit <- base::match.arg(unit)
  base::set.seed(seed)
  P <- base::nrow(R); N <- base::ncol(R)
  cols <- base::vector("list", B)
  len  <- NA_integer_
  for (b in base::seq_len(B)) {
    if (unit == "column") {
      w <- .wild.multipliers(N, dist)
      W <- base::matrix(w, nrow = P, ncol = N, byrow = TRUE)
    } else {
      W <- base::matrix(.wild.multipliers(P * N, dist), nrow = P, ncol = N)
    }
    Ys <- Yhat + W * R
    if (clipY) Ys <- base::pmax(Ys, 0)
    cb <- base::tryCatch(base::as.vector(refit.fun(Ys)),
                         error = function(e) NULL)
    if (!base::is.null(cb)) { cols[[b]] <- cb; len <- base::length(cb) }
  }
  if (base::is.na(len)) base::stop("all re-fit bootstrap replicates failed")
  out <- base::matrix(NA_real_, nrow = len, ncol = B)
  for (b in base::seq_len(B)) if (!base::is.null(cols[[b]])) out[, b] <- cols[[b]]
  out
}


#' @title Summarise bootstrap draws into SE / CI / bootstrap p-values (Internal)
#' @description From a \eqn{QK \times B} draw matrix, compute the bootstrap
#'   SE (row sd), percentile CI, and a two-sided bootstrap p-value
#'   \eqn{2 \min(\widehat P[c_b > 0], \widehat P[c_b < 0])} (useful when the
#'   analytical z is unavailable, e.g. method = "refit").
#' @param C.boot \eqn{QK \times B} draws.
#' @param level Confidence level for the percentile CI.
#' @return List with \code{se} (length \eqn{QK}), \code{ci.lower},
#'   \code{ci.upper}, and \code{p.boot} vectors.
#' @keywords internal
#' @noRd
.boot.summarize <- function(C.boot, level = 0.95) {
  alpha <- 1 - level
  se <- base::apply(C.boot, 1, stats::sd, na.rm = TRUE)
  lo <- base::apply(C.boot, 1, stats::quantile, probs = alpha / 2,
                    na.rm = TRUE, names = FALSE)
  hi <- base::apply(C.boot, 1, stats::quantile, probs = 1 - alpha / 2,
                    na.rm = TRUE, names = FALSE)
  p.boot <- base::apply(C.boot, 1, function(v) {
    v <- v[base::is.finite(v)]
    if (!base::length(v)) return(NA_real_)
    pg <- base::mean(v > 0); pl <- base::mean(v < 0)
    base::min(1, 2 * base::min(pg, pl))
  })
  base::list(se = se, ci.lower = lo, ci.upper = hi, p.boot = p.boot)
}

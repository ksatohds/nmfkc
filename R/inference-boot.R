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


#' @title Block wild bootstrap multipliers (Internal)
#' @description Draw \code{n} mean-0, variance-1 multipliers that are constant
#'   within blocks of \code{L} consecutive columns, for a block wild bootstrap
#'   of serially dependent residuals (columns = time points).  Residuals at
#'   nearby time points then share a multiplier, so their serial dependence
#'   within a block is kept.  The block grid starts at a random offset in
#'   \code{0, ..., L-1}, so no column is always at a block edge.  With
#'   \code{L = 1} this is exactly \code{.wild.multipliers(n, dist)} and draws
#'   the same random numbers.
#' @param n Number of multipliers (columns).
#' @param L Block length, a positive integer.
#' @param dist Multiplier distribution (see \code{.wild.multipliers}).
#' @return Numeric vector of length \code{n}.
#' @keywords internal
#' @noRd
.wild.block.multipliers <- function(n, L, dist = c("rademacher", "mammen", "exp")) {
  dist <- base::match.arg(dist)
  L <- base::suppressWarnings(base::as.integer(L))
  if (base::length(L) != 1L || base::is.na(L) || L < 1L)
    base::stop("The block length must be a single positive integer.")
  if (L == 1L) return(.wild.multipliers(n, dist))
  s  <- base::sample.int(L, 1L) - 1L                    # offset in 0..(L-1)
  id <- (base::seq_len(n) - 1L + s) %/% L + 1L           # block of each column
  v  <- .wild.multipliers(base::max(id), dist)
  v[id]
}


#' @title Default block length for the block wild bootstrap (Internal)
#' @description \code{max(1, round(n^(1/3)))}, the usual rate for block
#'   resampling of a stationary series of length \code{n}.
#' @param n Number of columns (time points).
#' @return Positive integer.
#' @keywords internal
#' @noRd
.wild.block.default <- function(n) base::max(1L, base::as.integer(base::round(n^(1 / 3))))


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
#' @param unit \code{"element"} (i.i.d. multiplier per matrix cell),
#'   \code{"column"} (one multiplier per sample column, shared over rows) or
#'   \code{"block"} (one multiplier per block of \code{block} consecutive
#'   columns, shared over rows; see \code{.wild.block.multipliers}).
#' @param block Block length for \code{unit = "block"}.
#' @param clipY Logical; clip the perturbed response at 0 (NMF needs
#'   \eqn{Y \ge 0}). Default \code{TRUE}.
#' @return \eqn{QK \times B} matrix of bootstrap draws (\code{NA} column if a
#'   replicate fails).
#' @keywords internal
#' @noRd
.boot.refit <- function(Yhat, R, refit.fun, B,
                        dist = "rademacher", seed = 123,
                        unit = c("element", "column", "block"), block = 1L,
                        clipY = TRUE) {
  unit <- base::match.arg(unit)
  base::set.seed(seed)
  P <- base::nrow(R); N <- base::ncol(R)
  cols <- base::vector("list", B)
  len  <- NA_integer_
  for (b in base::seq_len(B)) {
    if (unit == "column") {
      w <- .wild.multipliers(N, dist)
      W <- base::matrix(w, nrow = P, ncol = N, byrow = TRUE)
    } else if (unit == "block") {
      w <- .wild.block.multipliers(N, block, dist)
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
#'   SE (row sd), percentile CI, and a bootstrap p-value for
#'   \eqn{H_0: c_{qj} = 0} obtained by basic-bootstrap inversion,
#'   \eqn{\widehat P[c^*_b \ge 2\hat c]} (useful when the analytical z is
#'   unavailable, e.g. method = "refit").  The replicates are centred at the
#'   point estimate: comparing the raw replicates with zero is degenerate under
#'   a non-negativity constraint, where every replicate is \eqn{\ge 0}.
#' @param C.boot \eqn{QK \times B} draws.
#' @param est Length-\eqn{QK} vector of point estimates the draws are centred on.
#' @param level Confidence level for the percentile CI.
#' @param p.side \code{"one.sided"} (default, natural for a non-negative
#'   coefficient) or \code{"two.sided"}.
#' @return List with \code{se} (length \eqn{QK}), \code{ci.lower},
#'   \code{ci.upper}, and \code{p.boot} vectors.
#' @keywords internal
#' @noRd
.boot.summarize <- function(C.boot, est, level = 0.95,
                            p.side = c("one.sided", "two.sided")) {
  p.side <- base::match.arg(p.side)
  alpha <- 1 - level
  se <- base::apply(C.boot, 1, stats::sd, na.rm = TRUE)
  qs <- base::apply(C.boot, 1, stats::quantile,
                    probs = c(alpha / 2, 1 - alpha / 2),
                    na.rm = TRUE, names = FALSE)
  lo <- qs[1L, ]
  hi <- qs[2L, ]
  ## Basic-bootstrap inversion of H0: C_qj = 0.  The replicates must be
  ## CENTRED at the point estimate: (Chat* - Chat) approximates (Chat - C), so
  ## under H0 the observed value has upper-tail probability
  ## P*(Chat* - Chat >= Chat) = P*(Chat* >= 2 Chat).
  ## Testing the raw replicates against zero instead is degenerate here, since
  ## the non-negativity constraint makes every replicate >= 0 and would return
  ## p = 0 for every coefficient.
  p.boot <- base::vapply(base::seq_len(base::nrow(C.boot)), function(i) {
    v <- C.boot[i, ]; v <- v[base::is.finite(v)]
    if (!base::length(v)) return(NA_real_)
    p1 <- base::mean(v >= 2 * est[i])
    if (p.side == "one.sided") base::min(1, p1)
    else base::min(1, 2 * base::min(p1, 1 - p1))
  }, base::numeric(1))
  base::list(se = se, ci.lower = lo, ci.upper = hi, p.boot = p.boot)
}

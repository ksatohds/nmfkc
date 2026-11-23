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
#' @param cores Number of workers for the re-fits (\code{.nmfkc.parlapply}).
#'   With more than one, the multipliers are drawn up front in the order the
#'   sequential loop draws them; \code{refit.fun} must leave the random stream
#'   as it found it (\code{nmfkc()} does), so the draws do not depend on
#'   \code{cores}.
#' @return \eqn{QK \times B} matrix of bootstrap draws (\code{NA} column if a
#'   replicate fails).
#' @keywords internal
#' @noRd
.boot.refit <- function(Yhat, R, refit.fun, B,
                        dist = "rademacher", seed = 123,
                        unit = c("element", "column", "block"), block = 1L,
                        clipY = TRUE, cores = 1L) {
  unit <- base::match.arg(unit)
  base::set.seed(seed)
  P <- base::nrow(R); N <- base::ncol(R)
  cols <- base::vector("list", B)
  len  <- NA_integer_
  cores <- base::suppressWarnings(base::as.integer(cores))
  if (base::length(cores) != 1L || base::is.na(cores) || cores <= 1L) {
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
  } else {
    ## the multipliers of all replicates, drawn as the loop above draws them
    ws <- base::lapply(base::seq_len(B), function(b) {
      if (unit == "column") .wild.multipliers(N, dist)
      else if (unit == "block") .wild.block.multipliers(N, block, dist)
      else .wild.multipliers(P * N, dist)
    })
    one <- function(w) {
      W <- if (unit == "element") base::matrix(w, nrow = P, ncol = N)
           else base::matrix(w, nrow = P, ncol = N, byrow = TRUE)
      Ys <- Yhat + W * R
      if (clipY) Ys <- base::pmax(Ys, 0)
      base::tryCatch(base::as.vector(refit.fun(Ys)), error = function(e) NULL)
    }
    cols <- .nmfkc.parlapply(ws, one, cores = cores, envir = base::environment())
    ## a worker that died returns an error object, not a vector: count it as failed
    cols <- base::lapply(cols, function(cb) if (base::is.numeric(cb)) cb else NULL)
    for (cb in cols) if (!base::is.null(cb)) len <- base::length(cb)
  }
  if (base::is.na(len)) base::stop("all re-fit bootstrap replicates failed")
  out <- base::matrix(NA_real_, nrow = len, ncol = B)
  for (b in base::seq_len(B)) if (!base::is.null(cols[[b]])) out[, b] <- cols[[b]]
  out
}


#' @title All permutations of the bases (Internal)
#' @description The \code{factorial(Q)} permutations of \code{1:Q}, one per
#'   row, for \code{Q <= 6}; \code{NULL} beyond (\code{.align.bases} then
#'   matches the bases greedily).
#' @param Q Number of bases.
#' @keywords internal
#' @noRd
.all.perms <- function(Q) {
  if (Q <= 6)
    do.call(rbind, lapply(seq_len(factorial(Q)), function(i) .nmfvar.perm(Q, i)))
  else NULL
}


#' @title Align the bases of a re-fit to the original basis (Internal)
#' @description The permutation \code{p} that makes \code{Xs[, p]} closest to
#'   \code{X} in squared distance: over all permutations when \code{perms} is
#'   given, otherwise greedily, basis by basis.  Shared by the bootstraps that
#'   re-estimate the basis (\code{nmfkc.ar.latent.inference},
#'   \code{nmfkc.inference(method = "refit.basis")}).
#' @param Xs Basis of a re-fit (\eqn{P \times Q}).
#' @param X Original basis (\eqn{P \times Q}).
#' @param perms From \code{.all.perms(Q)}.
#' @return Integer permutation of \code{1:Q}: basis \code{q} of the original
#'   is basis \code{p[q]} of the re-fit.
#' @keywords internal
#' @noRd
.align.bases <- function(Xs, X, perms) {
  Q <- ncol(X)
  if (!is.null(perms)) {
    d <- apply(perms, 1, function(p) sum((Xs[, p, drop = FALSE] - X)^2))
    perms[which.min(d), ]
  } else {
    p <- integer(Q); avail <- seq_len(Q)
    for (q in seq_len(Q)) {
      d <- vapply(avail, function(j) sum((Xs[, j] - X[, q])^2), numeric(1))
      p[q] <- avail[which.min(d)]; avail <- setdiff(avail, p[q])
    }
    p
  }
}


#' @title Arguments for bootstrap re-fits of the whole model (Internal)
#' @description A re-fit must use the same estimator as the original fit, or
#'   the bootstrap describes a different one.  Starts from what the fit records
#'   (\code{fit.args}), overrides with the caller's \code{refit.args}, inherits
#'   a literal numeric \code{seed} from the recorded call, and lists the other
#'   arguments of that call that cannot be reproduced (the caller warns).
#' @param fit.args Named list recorded from the fit (\code{NULL} entries dropped).
#' @param call The matched call of the fit, or \code{NULL}.
#' @param user The caller's \code{refit.args}, or \code{NULL}.
#' @return List with \code{args} (the arguments for the re-fits) and
#'   \code{unseen} (names of arguments in the call that the re-fits cannot see).
#' @keywords internal
#' @noRd
.refit.args.collect <- function(fit.args, call, user) {
  refit.args <- fit.args[!vapply(fit.args, is.null, logical(1))]
  if (!is.null(user)) refit.args[names(user)] <- user
  unseen <- character(0)
  if (!is.null(call)) {
    ## Arguments that cannot change the estimator (labels and reporting) plus the
    ## ones set or inherited here.
    seen <- c("Y", "A", "rank", "data", "epsilon", "maxit", "verbose",
              "prefix", "print.trace", "detail", "", names(refit.args))
    cl <- as.list(call)
    ## A literal seed can be inherited outright; anything else has to be passed.
    if (!is.null(cl$seed) && is.numeric(cl$seed) && is.null(refit.args$seed)) {
      refit.args$seed <- cl$seed
      seen <- c(seen, "seed")
    }
    unseen <- setdiff(names(cl)[-1], seen)
  }
  list(args = refit.args, unseen = unseen)
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

# nmf.ffb.fiml.R -- NMF-FFB, likelihood-based (FIML) estimator
#
# Internal engine behind nmf.ffb(method = "fiml") and the fiml branch of
# nmf.ffb.inference().  Nothing in this file is exported; the user-facing
# wrappers live in R/nmf.sem.R.
#
# Model (conditional on a fixed non-negative, column-stochastic basis X,
# P1 x Q, estimated in a first stage by nmfkc(Y1, A = Y2)):
#
#   Y1 = X B + E,   B = Theta1 Y1 + Theta2 Y2 + U,
#   U ~ N(0, Phi)  (Q x Q, full via Cholesky or diagonal),
#   E ~ N(0, diag(psi))  (P1 x P1).
#
# Reduced form:  Y1 | Y2 ~ N(M Y2, Sigma),  A = X Theta1,  L = (I - A)^{-1},
#   M = L X Theta2,  Sigma = L (X Phi X' + diag(psi)) L'.
#
# The feed-forward (FF) null is Theta1 = 0, a non-negative MIMIC factor model
# with correlated factors; the feedback (FFB) alternative frees the entries of
# Theta1 admitted by an exclusion mask.  Theta1 is not recoverable from the
# reduced form alone (the structural and reduced forms have the same
# second moments once X is free), which is why X is fixed first and the
# self-loop / dominant-factor exclusion restriction is imposed; the L1 path
# and BIC then select a sparse Theta1.  The likelihood-ratio statistic against
# the FF null is returned WITHOUT a p-value: Theta1 >= 0 puts the null on the
# boundary and the BIC refit is a post-selection statistic, so the reference
# distribution is obtained by parametric bootstrap in nmf.ffb.inference().
#
# Objective (negative log-likelihood up to the constant (N P1 / 2) log 2 pi):
#   nll = (N/2) [log det Sigma + tr(Sigma^{-1} S)] + lambda1 * sum(Theta1),
#   S = R R' / N,  R = Y1 - M Y2.
# Gradient (Omega = Sigma^{-1}, G = Omega - Omega S Omega, Yhat = M Y2):
#   dTheta1 = N X' L' (I - Omega S) - X' L' Omega R Yhat' + lambda1   (free entries)
#   dTheta2 = - X' L' Omega R Y2'
#   Gphi    = (N/2) X' L' G L X       (w.r.t. symmetric Phi)
#     Cholesky Phi = Lc Lc':  dLc = 2 Gphi Lc (lower triangle); the diagonal
#     is parameterized as log Lc_ii, so its entry is multiplied by Lc_ii
#     diagonal Phi = diag(exp(v)):  dv = diag(Gphi) exp(v)
#   dw (w = log psi) = (N/2) diag(L' G L) psi
# The term -X' L' Omega R Yhat' carries no trailing L' (Yhat' already contains
# it); this was verified against central finite differences.
#
# S, R Y2' and R Yhat' are formed from the sufficient statistics Y1 Y1',
# Y1 Y2', Y2 Y2', so one objective / gradient evaluation is O(P^3), not
# O(P^2 N).

## ---------------------------------------------------------------------------
## small helpers
## ---------------------------------------------------------------------------

#' Spectral radius of a square matrix (Internal)
#' @param A A square numeric matrix.
#' @return \code{max(Mod(eigen(A)$values))}.
#' @keywords internal
#' @noRd
.ffb.spectral.radius <- function(A) {
  base::max(base::Mod(base::eigen(A, only.values = TRUE)$values))
}

#' Exclusion mask for Theta1 (Internal)
#'
#' @param X Basis matrix (P1 x Q), column-stochastic.
#' @param mask \code{"block"}: entry (q, i) is free unless q is the dominant
#'   factor of outcome i (\code{which.max(X[i, ])}), i.e. an outcome may not
#'   feed back into the factor that generates it (no self-loops);
#'   \code{"cross"}: free where \code{X[i, q] < cross.threshold};
#'   \code{"none"}: all entries free; or a user Q x P1 0/1 matrix.
#' @param cross.threshold Loading threshold for \code{mask = "cross"}.
#' @return A Q x P1 matrix of 0 / 1 (1 = free).
#' @keywords internal
#' @noRd
.ffb.fiml.mask <- function(X, mask = "block", cross.threshold = 0.05) {
  P1 <- base::nrow(X); Q <- base::ncol(X)
  if (base::is.matrix(mask) || (base::is.numeric(mask) && base::length(mask) == Q * P1)) {
    M <- base::matrix(base::as.numeric(mask), Q, P1)
    if (base::any(!(M %in% c(0, 1))))
      base::stop("a user-supplied `mask` must be a Q x P1 matrix of 0 / 1.")
  } else {
    mask <- base::match.arg(mask, c("block", "cross", "none"))
    M <- base::matrix(1, Q, P1)
    if (mask == "block") {
      dom <- base::apply(X, 1L, base::which.max)
      for (i in base::seq_len(P1)) M[dom[i], i] <- 0
    } else if (mask == "cross") {
      M <- (base::t(X) < cross.threshold) * 1
    }
  }
  base::dimnames(M) <- base::list(base::colnames(X), base::rownames(X))
  M
}

#' Objective and analytic gradient of the working-model FIML (Internal)
#'
#' Builds closures for one design (data, basis, mask, penalty).  The
#' parameter vector is \code{c(free entries of Theta1, Theta2, phi-parameters,
#' log psi)}; masked-out entries of Theta1 are not parameters at all.
#'
#' @param Y1,Y2 Data blocks (P1 x N, P2 x N).
#' @param X Basis (P1 x Q).
#' @param feedback Logical; \code{FALSE} fits the FF null (Theta1 = 0).
#' @param T1.mask Q x P1 0/1 matrix of free entries (used when
#'   \code{feedback = TRUE}).
#' @param lambda1 L1 penalty on Theta1 (its entries are non-negative, so the
#'   penalty is \code{lambda1 * sum(Theta1)}).
#' @param phi.full Logical; full (Cholesky) or diagonal Phi.
#' @return A list of closures and index bookkeeping (\code{fn}, \code{gr},
#'   \code{unpack}, \code{phi.par}, \code{npar}, ...).
#' @keywords internal
#' @noRd
.ffb.fiml.objective <- function(Y1, Y2, X, feedback = TRUE, T1.mask = NULL,
                                lambda1 = 0, phi.full = TRUE) {
  N <- base::ncol(Y1); P1 <- base::nrow(Y1); P2 <- base::nrow(Y2); Q <- base::ncol(X)
  if (base::is.null(T1.mask)) T1.mask <- base::matrix(1, Q, P1)
  free1 <- if (feedback) base::which(T1.mask != 0) else base::integer(0)
  n1 <- base::length(free1); n2 <- Q * P2
  nphi <- if (phi.full) Q * (Q + 1) / 2 else Q
  idx1 <- base::seq_len(n1)
  idx2 <- n1 + base::seq_len(n2)
  idxphi <- n1 + n2 + base::seq_len(nphi)
  idxpsi <- n1 + n2 + nphi + base::seq_len(P1)
  npar <- n1 + n2 + nphi + P1
  lowtri <- base::which(base::lower.tri(base::diag(Q), diag = TRUE))
  mkLc <- function(v) {
    Lc <- base::matrix(0, Q, Q); Lc[lowtri] <- v
    base::diag(Lc) <- base::exp(base::diag(Lc)); Lc
  }
  mkPhi <- function(v) {
    if (!phi.full) return(base::diag(base::exp(v), Q))
    Lc <- mkLc(v); Lc %*% base::t(Lc)
  }
  unpack <- function(th) {
    T1 <- base::matrix(0, Q, P1)
    if (n1 > 0) T1[free1] <- th[idx1]
    base::list(T1 = T1, T2 = base::matrix(th[idx2], Q), Phi = mkPhi(th[idxphi]),
               psi = base::exp(th[idxpsi]))
  }
  I_P <- base::diag(P1)
  ## sufficient statistics
  Y1Y1t <- base::tcrossprod(Y1); Y1Y2t <- base::tcrossprod(Y1, Y2); Y2Y2t <- base::tcrossprod(Y2)
  bad <- base::list(value = 1e10, grad = base::rep(0, npar))
  fg <- function(th) {
    p <- unpack(th); A <- X %*% p$T1
    if (n1 > 0 && .ffb.spectral.radius(A) >= 0.999) return(bad)
    L <- base::tryCatch(base::solve(I_P - A), error = function(e) NULL)
    if (base::is.null(L)) return(bad)
    LX <- L %*% X                                     # P1 x Q
    M <- LX %*% p$T2                                  # P1 x P2
    RY2t <- Y1Y2t - M %*% Y2Y2t                       # R Y2'
    S <- (Y1Y1t - Y1Y2t %*% base::t(M) - M %*% base::t(Y1Y2t) + M %*% Y2Y2t %*% base::t(M)) / N
    S <- (S + base::t(S)) / 2
    Sig <- LX %*% p$Phi %*% base::t(LX) + L %*% (p$psi * base::t(L))
    Sig <- (Sig + base::t(Sig)) / 2
    ch <- base::tryCatch(base::chol(Sig), error = function(e) NULL)
    if (base::is.null(ch)) return(bad)
    Om <- base::chol2inv(ch); OmS <- Om %*% S
    val <- (N / 2) * (2 * base::sum(base::log(base::diag(ch))) + base::sum(base::diag(OmS))) +
      lambda1 * base::sum(p$T1)
    G <- Om - OmS %*% Om; G <- (G + base::t(G)) / 2   # Omega - Omega S Omega
    XtLtOm <- base::t(LX) %*% Om                      # X' L' Omega  (Q x P1)
    g <- base::numeric(npar)
    if (n1 > 0) {
      dT1 <- N * base::t(LX) %*% (I_P - OmS) - XtLtOm %*% RY2t %*% base::t(M) + lambda1
      g[idx1] <- dT1[free1]
    }
    g[idx2] <- base::as.numeric(-XtLtOm %*% RY2t)
    Gphi <- (N / 2) * base::t(LX) %*% G %*% LX
    if (phi.full) {
      Lc <- mkLc(th[idxphi]); dLc <- 2 * Gphi %*% Lc
      base::diag(dLc) <- base::diag(dLc) * base::diag(Lc)
      g[idxphi] <- dLc[lowtri]
    } else {
      g[idxphi] <- base::diag(Gphi) * base::exp(th[idxphi])
    }
    g[idxpsi] <- (N / 2) * base::diag(base::t(L) %*% G %*% L) * p$psi
    base::list(value = val, grad = g)
  }
  ## optim() calls fn() and gr() separately; memoise the last evaluation so
  ## the two share one computation.
  cache <- base::new.env(); cache$th <- NULL; cache$res <- NULL
  getfg <- function(th) {
    if (!base::is.null(cache$th) && base::identical(cache$th, th)) return(cache$res)
    r <- fg(th); cache$th <- th; cache$res <- r; r
  }
  fn <- function(th) getfg(th)$value
  gr <- function(th) getfg(th)$grad
  phi.par <- function(Phi) {
    if (!phi.full) return(base::log(base::pmax(base::diag(Phi), 1e-6)))
    Lc <- base::t(base::chol(Phi + base::diag(1e-8, Q))); d <- base::diag(Lc)
    Lc2 <- Lc; base::diag(Lc2) <- base::log(base::pmax(d, 1e-6)); Lc2[lowtri]
  }
  base::list(fn = fn, gr = gr, unpack = unpack, phi.par = phi.par, npar = npar,
             n1 = n1, n2 = n2, nphi = nphi, free1 = free1, I_P = I_P)
}

#' One FIML fit of the working model by L-BFGS-B (Internal)
#'
#' @inheritParams .ffb.fiml.objective
#' @param init Optional list with \code{T1}, \code{T2}, \code{Phi} (or
#'   \code{phi}), \code{psi} used as starting values.
#' @param maxit Maximum L-BFGS-B iterations.
#' @param factr \code{optim} tolerance (\code{control$factr}).
#' @return A list: \code{T1, T2, Phi, psi, loglik, conv, iter, rho, M, MAE,
#'   npar, nnz}.  \code{loglik} is the full Gaussian log-likelihood of the
#'   reduced form (penalty removed, constant included).
#' @keywords internal
#' @noRd
.ffb.fiml.fit <- function(Y1, Y2, X, feedback = TRUE, T1.mask = NULL, lambda1 = 0,
                          init = NULL, maxit = 3000, phi.full = TRUE, factr = 1e3) {
  N <- base::ncol(Y1); P1 <- base::nrow(Y1); P2 <- base::nrow(Y2); Q <- base::ncol(X)
  if (base::is.null(T1.mask)) T1.mask <- base::matrix(1, Q, P1)
  ob <- .ffb.fiml.objective(Y1, Y2, X, feedback = feedback, T1.mask = T1.mask,
                            lambda1 = lambda1, phi.full = phi.full)
  n1 <- ob$n1; nphi <- ob$nphi
  if (base::is.null(init)) {
    T2_0 <- base::matrix(0.1, Q, P2); Phi0 <- base::diag(0.01, Q)
    psi0 <- base::rep(0.01, P1); T1_0 <- base::matrix(0.01, Q, P1)
  } else {
    T2_0 <- base::pmax(init$T2, 1e-4); psi0 <- base::pmax(init$psi, 1e-6)
    Phi0 <- if (!base::is.null(init$Phi)) init$Phi else base::diag(base::pmax(init$phi, 1e-6), Q)
    T1_0 <- if (!base::is.null(init$T1)) base::pmax(init$T1, 1e-4) else base::matrix(0.01, Q, P1)
  }
  th0 <- c(if (n1 > 0) base::as.numeric(T1_0)[ob$free1], base::as.numeric(T2_0),
           ob$phi.par(Phi0), base::log(psi0))
  lower <- c(base::rep(0, n1), base::rep(0, ob$n2), base::rep(-20, nphi), base::rep(-20, P1))
  upper <- c(base::rep(5, n1), base::rep(50, ob$n2), base::rep(5, nphi), base::rep(5, P1))
  o <- stats::optim(th0, ob$fn, ob$gr, method = "L-BFGS-B", lower = lower, upper = upper,
                    control = base::list(maxit = maxit, factr = factr))
  p <- ob$unpack(o$par); A <- X %*% p$T1
  L <- base::solve(ob$I_P - A); M <- L %*% X %*% p$T2
  ll <- -(o$value - lambda1 * base::sum(p$T1)) - (N * P1 / 2) * base::log(2 * base::pi)
  base::dimnames(p$T1) <- base::list(base::colnames(X), base::rownames(Y1))
  base::dimnames(p$T2) <- base::list(base::colnames(X), base::rownames(Y2))
  base::dimnames(p$Phi) <- base::list(base::colnames(X), base::colnames(X))
  base::names(p$psi) <- base::rownames(Y1)
  base::list(T1 = p$T1, T2 = p$T2, Phi = p$Phi, psi = p$psi, loglik = ll,
             conv = o$convergence, iter = base::as.integer(o$counts[["function"]]),
             rho = .ffb.spectral.radius(A), M = M,
             MAE = base::mean(base::abs(Y1 - M %*% Y2)), npar = ob$npar,
             nnz = base::sum(p$T1 > 1e-3))
}

#' Starting values for the FF null from a feed-forward fit (Internal)
#'
#' \code{Theta2} is taken from the nmfkc coefficient matrix; \code{phi} from
#' the variance of the free scores around \code{C Y2}; \code{psi} from the
#' residual variance of \code{Y1 - X C Y2}.
#' @param Y1,Y2 Data blocks.
#' @param X Basis (P1 x Q).
#' @param C Coefficient matrix (Q x P2) of the feed-forward fit.
#' @return A list \code{T2, phi, psi}.
#' @keywords internal
#' @noRd
.ffb.ff.init <- function(Y1, Y2, X, C) {
  B <- C %*% Y2
  Bfree <- base::tryCatch(base::solve(base::crossprod(X), base::t(X) %*% Y1),
                          error = function(e) B)
  phi <- base::pmax(base::apply(Bfree - B, 1, stats::var), 1e-4)
  psi <- base::pmax(base::apply(Y1 - X %*% C %*% Y2, 1, stats::var), 1e-4)
  base::list(T2 = C, phi = phi, psi = psi)
}

#' FF null -> unpenalized FFB -> multi-start L1 path -> BIC selection (Internal)
#'
#' The penalized problem is non-convex: L-BFGS-B started from different
#' points lands in different local optima and proposes different supports.
#' Every point of the path is therefore fitted from several starts
#' (\code{starts}), every distinct support proposed by any (lambda1, start)
#' is re-estimated without penalty, and BIC is minimized over ALL distinct
#' candidates (plus the null and the unpenalized full model).  The refit on
#' a support is itself run from three starts (the penalized solution, the
#' unpenalized fit restricted to the support, a small constant) and the best
#' log-likelihood is kept.
#'
#' @param Y1,Y2 Data blocks.
#' @param X Fixed basis (P1 x Q).
#' @param C.init Q x P2 starting value for Theta2 in the null fit (the nmfkc
#'   coefficient matrix, or the null Theta2 of an observed fit in the
#'   bootstrap).
#' @param mask Q x P1 0/1 matrix of admitted Theta1 entries.
#' @param phi.full Logical; full or diagonal Phi.
#' @param lambda1 Numeric vector of L1 penalties (the path; sorted
#'   increasingly, non-positive / infinite values dropped).
#' @param select \code{"BIC"} or \code{"none"} (the unpenalized fit is the
#'   selected model).
#' @param maxit,factr Passed to \code{.ffb.fiml.fit}.
#' @param starts Character vector, subset of \code{c("full", "path", "null",
#'   "soft")}: starting points of the penalized fits.  \code{"full"}: the
#'   unpenalized feedback fit; \code{"path"}: continuation from the best
#'   penalized solution at the previous (smaller) lambda1; \code{"null"}:
#'   the FF null with Theta1 = 0.05 on the free entries; \code{"soft"}: the
#'   unpenalized Theta1 soft-thresholded at lambda1 / N.
#' @return A list \code{f0} (null), \code{f1} (unpenalized), \code{fsel}
#'   (selected), \code{path} (data.frame, one row per (lambda1, start)),
#'   \code{candidates} (data.frame, one row per distinct support),
#'   \code{supports} (list of the distinct supports, Q x P1 logical),
#'   \code{sel} (selected path row), \code{lambda1.selected},
#'   \code{support.selected} (id into \code{supports}).
#' @keywords internal
#' @noRd
.ffb.fiml.pipeline <- function(Y1, Y2, X, C.init, mask, phi.full = TRUE, lambda1,
                               select = "BIC", maxit = 3000, factr = 1e3,
                               starts = c("full", "path", "null", "soft")) {
  N <- base::ncol(Y1); P1 <- base::nrow(Y1); Q <- base::ncol(X)
  starts <- base::match.arg(starts, c("full", "path", "null", "soft"), several.ok = TRUE)
  mask <- (mask != 0) * 1
  f0 <- .ffb.fiml.fit(Y1, Y2, X, feedback = FALSE, init = .ffb.ff.init(Y1, Y2, X, C.init),
                      maxit = maxit, phi.full = phi.full, factr = factr)
  k0 <- f0$npar
  no_feedback <- base::sum(mask) == 0
  f1 <- if (no_feedback) f0 else
    .ffb.fiml.fit(Y1, Y2, X, feedback = TRUE, T1.mask = mask,
                  init = base::list(T2 = f0$T2, Phi = f0$Phi, psi = f0$psi),
                  maxit = maxit, phi.full = phi.full, factr = factr)
  fit1 <- function(T1.mask, l, init)
    .ffb.fiml.fit(Y1, Y2, X, feedback = TRUE, T1.mask = T1.mask, lambda1 = l, init = init,
                  maxit = maxit, phi.full = phi.full, factr = factr)
  as_init <- function(f, T1 = f$T1) base::list(T1 = T1, T2 = f$T2, Phi = f$Phi, psi = f$psi)
  bic_of <- function(f, nnz) -2 * f$loglik + base::log(N) * (nnz + k0)

  ## ---- registry of distinct supports: key -> id; one refit per id ----
  key_of <- function(supp) base::paste(base::as.integer(supp != 0), collapse = "")
  keys <- base::character(0); supports <- base::list(); cand_fits <- base::list()
  cand_first <- base::list()      # (lambda1, start) that first proposed the support
  register <- function(supp, fit, l, start) {
    k <- key_of(supp); id <- base::match(k, keys)
    if (base::is.na(id)) {
      id <- base::length(keys) + 1L
      keys[id] <<- k; cand_fits[[id]] <<- fit
      sm <- (supp != 0); base::dimnames(sm) <- base::dimnames(mask); supports[[id]] <<- sm
      cand_first[[id]] <<- base::list(lambda1 = l, start = start)
    }
    id
  }
  id_null <- register(0 * mask, f0, Inf, "null")
  id_full <- if (no_feedback) id_null else register(mask, f1, 0, "full")
  ## unpenalized refit of a proposed support, best of three starts
  refit_support <- function(supp, fl) {
    nnz <- base::sum(supp)
    if (nnz == 0) return(f0)
    if (nnz == base::sum(mask)) return(f1)
    inits <- base::list(as_init(fl), as_init(f1),
                        as_init(f0, T1 = base::matrix(0.05, Q, P1)))
    best <- NULL
    for (ini in inits) {
      fr <- fit1(supp, 0, ini)
      if (base::is.null(best) || fr$loglik > best$loglik) best <- fr
    }
    best
  }
  row <- function(l, start, id, f, nnz, pen = NA_real_) base::data.frame(
    lambda1 = l, start = start, support_id = id, nnz = nnz, rho = f$rho, loglik = f$loglik,
    BIC = bic_of(f, nnz), MAE = f$MAE, pen.value = pen, stringsAsFactors = FALSE)
  ## lambda1 = 0: the unpenalized fit re-estimated on its own non-zero
  ## entries (the model with every admitted entry free is candidate id_full)
  supp1 <- (f1$T1 > 1e-3) * mask
  id_thr <- if (no_feedback) id_null else {
    id <- base::match(key_of(supp1), keys)
    if (base::is.na(id)) register(supp1, refit_support(supp1, f1), 0, "full") else id
  }
  path <- base::list(row(0, "full", id_thr, cand_fits[[id_thr]], base::sum(supp1)))

  ## ---- the multi-start L1 path ----
  if (select == "BIC" && !no_feedback) {
    lambda1 <- base::sort(base::unique(lambda1[base::is.finite(lambda1) & lambda1 > 0]))
    prev <- NULL                  # best penalized solution at the previous lambda1
    for (l in lambda1) {
      inits <- base::list()
      if ("full" %in% starts) inits$full <- as_init(f1)
      if ("path" %in% starts && !base::is.null(prev)) inits$path <- as_init(prev)
      if ("null" %in% starts) inits$null <- as_init(f0, T1 = base::matrix(0.05, Q, P1))
      if ("soft" %in% starts) inits$soft <- as_init(f1, T1 = base::pmax(f1$T1 - l / N, 0))
      best_pen <- NULL
      for (st in base::names(inits)) {
        fl <- fit1(mask, l, inits[[st]])
        pen <- -fl$loglik + l * base::sum(fl$T1)      # penalized objective (up to a constant)
        supp <- (fl$T1 > 1e-3) * mask
        id <- base::match(key_of(supp), keys)
        if (base::is.na(id)) id <- register(supp, refit_support(supp, fl), l, st)
        path[[base::length(path) + 1]] <- row(l, st, id, cand_fits[[id]], base::sum(supp), pen)
        if (base::is.null(best_pen) || pen < best_pen$pen) best_pen <- base::list(pen = pen, fit = fl)
      }
      prev <- best_pen$fit
    }
  }
  path[[base::length(path) + 1]] <- row(Inf, "null", id_null, f0, 0)
  path <- base::do.call(base::rbind, path)
  base::rownames(path) <- NULL
  path$duplicate <- base::duplicated(path$support_id)

  ## ---- candidates (one row per distinct support) and selection ----
  candidates <- base::do.call(base::rbind, base::lapply(base::seq_along(keys), function(id) {
    f <- cand_fits[[id]]; nnz <- base::sum(supports[[id]])
    base::data.frame(support_id = id, nnz = nnz, rho = f$rho, loglik = f$loglik,
                     BIC = bic_of(f, nnz), MAE = f$MAE,
                     lambda1.first = cand_first[[id]]$lambda1, start.first = cand_first[[id]]$start,
                     stringsAsFactors = FALSE)
  }))
  base::rownames(candidates) <- NULL
  if (select == "BIC") {
    ## smallest BIC over every distinct support; ties go to the sparser model
    ord <- base::order(candidates$BIC, candidates$nnz)
    id_sel <- candidates$support_id[ord[1]]
  } else {
    id_sel <- id_full             # select = "none": the unpenalized fit
  }
  fsel <- cand_fits[[id_sel]]
  candidates$selected <- candidates$support_id == id_sel
  b <- candidates[id_sel, ]
  base::list(f0 = f0, f1 = f1, fsel = fsel, path = path, candidates = candidates,
             supports = supports, sel = b, lambda1.selected = b$lambda1.first,
             support.selected = id_sel)
}

#' Equilibrium / Leontief quantities shared with the MU estimator (Internal)
#' @param X Basis (P1 x Q).
#' @param C1,C2 Coefficient matrices.
#' @param Y1,Y2 Data blocks (for MAE).
#' @return A list of the legacy fields (\code{XC1, XC2, XC1.radius, ...}).
#' @keywords internal
#' @noRd
.ffb.equilibrium <- function(X, C1, C2, Y1, Y2) {
  .eps <- 1e-10
  mat1norm <- function(A) base::max(base::colSums(base::abs(A)))
  XC1 <- X %*% C1; XC2 <- X %*% C2
  rho <- .ffb.spectral.radius(XC1)
  XC1_norm1 <- mat1norm(XC1)
  Leontief.inv <- base::tryCatch(base::solve(base::diag(base::nrow(XC1)) - XC1),
                                 error = function(e) base::matrix(NA_real_, base::nrow(XC1), base::ncol(XC1)))
  M.model <- Leontief.inv %*% XC2
  Y1_hat <- M.model %*% Y2
  base::list(XC1 = XC1, XC2 = XC2, XC1.radius = rho, XC1.norm1 = XC1_norm1,
             Leontief.inv = Leontief.inv, M.model = M.model,
             amplification = mat1norm(M.model) / (mat1norm(XC2) + .eps),
             amplification.bound = if (XC1_norm1 < 1) 1 / (1 - XC1_norm1) else Inf,
             MAE = if (base::anyNA(M.model)) NA_real_ else base::mean(base::abs(Y1 - Y1_hat)),
             effective.rank = .effective.rank(C1 %*% Y1 + C2 %*% Y2))
}

## ---------------------------------------------------------------------------
## nmf.ffb(method = "fiml")
## ---------------------------------------------------------------------------

#' Two-stage likelihood-based NMF-FFB fit (Internal)
#'
#' Stage 1 estimates the basis with \code{\link{nmfkc}} (unless \code{X} is
#' supplied); stage 2 runs \code{.ffb.fiml.pipeline} conditional on it and
#' assembles an object of class \code{c("nmf.ffb", "nmf.sem", "nmf")} with the
#' legacy fields plus the likelihood fields.  See \code{\link{nmf.ffb}}.
#' @keywords internal
#' @noRd
.nmf.ffb.fiml <- function(Y1, Y2, rank, X.init, X.L2.ortho, epsilon, maxit, seed,
                          X = NULL, mask = "block", cross.threshold = 0.05,
                          phi = "full", lambda1 = NULL, select = "BIC",
                          starts = c("full", "path", "null", "soft"), cl = NULL, ...) {
  extra_args <- base::list(...)
  starts <- base::match.arg(starts, c("full", "path", "null", "soft"), several.ok = TRUE)
  if (!base::is.matrix(Y1)) Y1 <- base::as.matrix(Y1)
  if (!base::is.matrix(Y2)) Y2 <- base::as.matrix(Y2)
  if (base::any(!base::is.finite(Y1)) || base::any(!base::is.finite(Y2)))
    base::stop("Y1 and Y2 must not contain NA/NaN/Inf.")
  if (base::min(Y1) < 0 || base::min(Y2) < 0) base::stop("Y1 and Y2 must be non-negative.")
  if (base::ncol(Y1) != base::ncol(Y2)) base::stop("ncol(Y1) must be equal to ncol(Y2).")
  P1 <- base::nrow(Y1); P2 <- base::nrow(Y2); N <- base::ncol(Y1)
  Y1_labels <- if (!base::is.null(base::rownames(Y1))) base::rownames(Y1) else base::paste0("Y1_", 1:P1)
  Y2_labels <- if (!base::is.null(base::rownames(Y2))) base::rownames(Y2) else base::paste0("Y2_", 1:P2)
  base::rownames(Y1) <- Y1_labels; base::rownames(Y2) <- Y2_labels

  fiml.maxit <- if (!base::is.null(extra_args$fiml.maxit)) extra_args$fiml.maxit else 3000L
  factr      <- if (!base::is.null(extra_args$factr))      extra_args$factr      else 1e3
  phi.full   <- base::identical(phi, "full")

  ## Keep the self-seeding of stage 1 out of the caller's random stream.
  .rng <- .nmfkc.rng.save(seed)
  base::on.exit(.nmfkc.rng.restore(.rng), add = TRUE)

  ## ---- stage 1: basis ----
  stage1 <- NULL; C.init <- NULL
  if (base::is.null(X)) {
    Q_hidden <- if (!base::is.null(extra_args$Q)) extra_args$Q else NULL
    Q <- if (!base::is.null(rank)) rank else if (!base::is.null(Q_hidden)) Q_hidden else P2
    if (Q < 1) base::stop("Rank Q must be >= 1.")
    if (base::is.null(X.init)) X.init <- "nndsvd"
    stage1 <- nmfkc(Y = Y1, A = Y2, Q = Q, X.init = X.init, X.L2.ortho = X.L2.ortho,
                    epsilon = epsilon, maxit = maxit, seed = seed,
                    verbose = FALSE, print.dims = FALSE)
    Xb <- stage1$X; C.init <- stage1$C
  } else {
    Xb <- if (base::is.list(X) && !base::is.null(X$X)) X$X else base::as.matrix(X)
    if (base::nrow(Xb) != P1) base::stop("X must have nrow(Y1) rows.")
    Q <- base::ncol(Xb)
    ## starting Theta2: the nmfkc coefficient matrix, or the null Theta2 of a
    ## previous fiml fit (an nmf.ffb object carries no $C)
    if (base::is.list(X) && !base::is.null(X$C) && base::all(base::dim(X$C) == c(Q, P2)))
      C.init <- X$C
    else if (base::is.list(X) && !base::is.null(X$null$C2) && base::all(base::dim(X$null$C2) == c(Q, P2)))
      C.init <- X$null$C2
  }
  Xb[Xb < 0] <- 0
  Xb <- base::sweep(Xb, 2, base::pmax(base::colSums(Xb), 1e-10), "/")
  Basis_labels <- base::paste0("Factor", 1:Q)
  base::dimnames(Xb) <- base::list(Y1_labels, Basis_labels)
  if (base::is.null(C.init)) {
    ## least-squares start for Theta2 given X: B = (X'X)^-1 X'Y1, C = B Y2^+
    B0 <- base::tryCatch(base::solve(base::crossprod(Xb), base::crossprod(Xb, Y1)),
                         error = function(e) base::matrix(0.1, Q, N))
    C.init <- base::tryCatch(B0 %*% base::t(Y2) %*% base::solve(base::tcrossprod(Y2)),
                             error = function(e) base::matrix(0.1, Q, P2))
    C.init <- base::pmax(C.init, 1e-4)
  }
  M1 <- .ffb.fiml.mask(Xb, mask, cross.threshold)
  if (base::is.null(lambda1)) lambda1 <- N * c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  if (base::sum(M1) == 0)
    base::warning("the exclusion mask admits no feedback entry (all Theta1 entries excluded); ",
                  "the fit reduces to the FF null.")

  ## ---- stage 2: FIML pipeline ----
  pp <- .ffb.fiml.pipeline(Y1, Y2, Xb, C.init, M1, phi.full = phi.full, lambda1 = lambda1,
                           select = select, maxit = fiml.maxit, factr = factr, starts = starts)
  f0 <- pp$f0; f1 <- pp$f1; fs <- pp$fsel
  eq <- .ffb.equilibrium(Xb, fs$T1, fs$T2, Y1, Y2)
  if (eq$XC1.radius >= 1) base::warning("Leontief.inv may be unstable; spectral radius >= 1.")
  eq0 <- .ffb.equilibrium(Xb, 0 * fs$T1, f0$T2, Y1, Y2)
  nnz_sel <- base::sum(fs$T1 > 1e-3)
  LR <- c(full = 2 * (f1$loglik - f0$loglik), selected = 2 * (fs$loglik - f0$loglik))
  LR.df <- c(full = base::sum(M1), selected = nnz_sel)
  base::attr(LR, "df") <- LR.df
  k0 <- f0$npar
  ll <- c(null = f0$loglik, full = f1$loglik, selected = fs$loglik)
  kk <- c(null = k0, full = k0 + base::sum(M1), selected = k0 + nnz_sel)

  out <- base::list(
    call = cl,
    dims = base::sprintf("Y1(%d,%d)~X(%d,%d)[C1(%d,%d)Y1+C2(%d,%d)Y2]",
                         P1, N, P1, Q, Q, P1, Q, P2),
    X = Xb, C1 = fs$T1, C2 = fs$T2,
    XC1 = eq$XC1, XC2 = eq$XC2, XC1.radius = eq$XC1.radius, XC1.norm1 = eq$XC1.norm1,
    Leontief.inv = eq$Leontief.inv, M.model = eq$M.model,
    amplification = eq$amplification, amplification.bound = eq$amplification.bound,
    Q = Q, SC.cov = NULL, SC.map = NULL, MAE = eq$MAE,
    effective.rank = eq$effective.rank,
    objfunc = NULL, objfunc.full = NULL,
    iter = fs$iter, maxit = fiml.maxit, epsilon = epsilon, converged = fs$conv == 0,
    method = "fiml",
    Phi = fs$Phi, psi = fs$psi, loglik = fs$loglik, npar = kk[["selected"]],
    null = base::list(C2 = f0$T2, Phi = f0$Phi, psi = f0$psi, loglik = f0$loglik,
                      npar = k0, M.model = eq0$M.model, conv = f0$conv),
    full = base::list(C1 = f1$T1, C2 = f1$T2, Phi = f1$Phi, psi = f1$psi, loglik = f1$loglik,
                      npar = kk[["full"]], XC1.radius = f1$rho, conv = f1$conv),
    path = pp$path, candidates = pp$candidates, supports = pp$supports,
    support.selected = pp$support.selected,
    mask = M1, lambda1 = lambda1, lambda1.selected = pp$lambda1.selected,
    support = fs$T1 > 1e-3,
    LR = LR, LR.df = LR.df,
    BIC = -2 * ll + base::log(N) * kk,
    AIC = -2 * ll + 2 * kk,
    phi = phi, select = select, starts = starts, cross.threshold = cross.threshold,
    factr = factr, X.L2.ortho = X.L2.ortho,
    stage1 = if (base::is.null(stage1)) NULL else
      base::list(iter = stage1$iter, converged = stage1$converged, objfunc = stage1$objfunc)
  )
  base::class(out) <- c("nmf.ffb", "nmf.sem", "nmf")
  out
}

## ---------------------------------------------------------------------------
## nmf.ffb.inference(method = "fiml"): two parametric bootstraps
## ---------------------------------------------------------------------------

#' Draw Y1* from the working model (Internal)
#'
#' \eqn{Y_1^* = L X (\Theta_2 Y_2 + U^*) + L E^*}, \eqn{L = (I - X\Theta_1)^{-1}}
#' (\eqn{L = I} under the FF null), negatives clipped at 0.
#' @param X Basis; \code{T1} Theta1 (or \code{NULL} for the null);
#'   \code{T2, Phi, psi} the generating parameters; \code{Y2} the fixed
#'   exogenous block.
#' @return A P1 x N matrix.
#' @keywords internal
#' @noRd
.ffb.fiml.sim <- function(X, T1, T2, Phi, psi, Y2) {
  Q <- base::ncol(X); N <- base::ncol(Y2); P1 <- base::nrow(X)
  Lphi <- base::t(base::chol(Phi + base::diag(1e-10, Q)))
  U <- Lphi %*% base::matrix(stats::rnorm(Q * N), Q, N)
  E <- base::matrix(stats::rnorm(P1 * N), P1, N) * base::sqrt(psi)
  Y1s <- X %*% (T2 %*% Y2 + U) + E
  if (!base::is.null(T1) && base::any(T1 != 0)) {
    L <- base::solve(base::diag(P1) - X %*% T1)
    Y1s <- L %*% Y1s
  }
  Y1s <- base::pmax(Y1s, 0)
  base::rownames(Y1s) <- base::rownames(X)
  Y1s
}

#' Parametric-bootstrap inference for a fiml NMF-FFB fit (Internal)
#'
#' Implements the \code{method = "fiml"} branch of
#' \code{\link{nmf.ffb.inference}}; see that help page for the returned
#' fields.
#' @keywords internal
#' @noRd
.nmf.ffb.inference.fiml <- function(object, Y1, Y2, B = 1000L, threshold = 0.01,
                                    ci.level = 0.95, seed = 123L, ...) {
  extra_args <- base::list(...)
  cores <- if (!base::is.null(extra_args$cores)) extra_args$cores
           else if (!base::is.null(extra_args$ncores)) extra_args$ncores
           else base::getOption("mc.cores", 1L)
  print.trace <- if (!base::is.null(extra_args$print.trace)) extra_args$print.trace else FALSE
  factr <- if (!base::is.null(extra_args$factr)) extra_args$factr
           else if (!base::is.null(object$factr)) object$factr else 1e3
  fiml.maxit <- if (!base::is.null(extra_args$fiml.maxit)) extra_args$fiml.maxit
                else if (!base::is.null(object$maxit)) object$maxit else 3000L
  boot.null <- if (!base::is.null(extra_args$boot.null)) base::isTRUE(extra_args$boot.null) else TRUE
  B <- base::as.integer(B)

  .rng <- .nmfkc.rng.save(seed)
  base::on.exit(.nmfkc.rng.restore(.rng), add = TRUE)

  Y1 <- base::as.matrix(Y1); Y2 <- base::as.matrix(Y2)
  if (base::ncol(Y1) != base::ncol(Y2)) base::stop("Y1 and Y2 must have the same number of columns.")
  X <- object$X; Q <- base::ncol(X); P1 <- base::nrow(Y1); P2 <- base::nrow(Y2); N <- base::ncol(Y1)
  if (base::nrow(X) != P1) base::stop("nrow(X) must equal nrow(Y1); did Y1 change shape?")
  base::rownames(Y1) <- base::rownames(X)
  base::rownames(Y2) <- base::colnames(object$C2)
  mask <- object$mask
  phi.full <- base::identical(object$phi, "full")
  lambda1 <- object$lambda1
  select <- if (!base::is.null(object$select)) object$select else "BIC"
  starts <- if (!base::is.null(object$starts)) object$starts else c("full", "path", "null", "soft")
  support <- (object$C1 > 1e-3) * 1
  nnz_obs <- base::sum(support)

  ## ---- (i) null bootstrap: LR calibration ----
  T2_0 <- object$null$C2; Phi_0 <- object$null$Phi; psi_0 <- object$null$psi
  LR_obs <- object$LR
  boot_null_one <- function(b) {
    base::set.seed(seed + b)
    Y1s <- .ffb.fiml.sim(X, NULL, T2_0, Phi_0, psi_0, Y2)
    r <- base::tryCatch(
      .ffb.fiml.pipeline(Y1s, Y2, X, T2_0, mask, phi.full = phi.full, lambda1 = lambda1,
                         select = select, maxit = fiml.maxit, factr = factr, starts = starts),
      error = function(e) NULL)
    if (base::is.null(r))
      return(c(full = NA_real_, selected = NA_real_, nnz = NA_real_,
               conv.null = NA_real_, conv.full = NA_real_))
    ## conv.* are the L-BFGS-B convergence codes of the two fits (0 = converged).
    ## They are carried out of the worker so that the caller can report how many
    ## null replicates hit `maxit`: on a flat likelihood (small N, full Phi) that
    ## can be a large share, and it is not visible from the LR values alone.
    c(full = 2 * (r$f1$loglik - r$f0$loglik), selected = 2 * (r$fsel$loglik - r$f0$loglik),
      nnz = base::sum(r$fsel$T1 > 1e-3),
      conv.null = r$f0$conv, conv.full = r$f1$conv)
  }
  ## ---- (ii) selected-model bootstrap: coefficient uncertainty ----
  T1_s <- object$C1; T2_s <- object$C2; Phi_s <- object$Phi; psi_s <- object$psi
  boot_sel_one <- function(b) {
    base::set.seed(seed + B + b)
    Y1s <- .ffb.fiml.sim(X, T1_s, T2_s, Phi_s, psi_s, Y2)
    r <- base::tryCatch(
      .ffb.fiml.fit(Y1s, Y2, X, feedback = nnz_obs > 0, T1.mask = support,
                    init = base::list(T1 = T1_s, T2 = T2_s, Phi = Phi_s, psi = psi_s),
                    maxit = fiml.maxit, phi.full = phi.full, factr = factr),
      error = function(e) NULL)
    if (base::is.null(r) || !base::is.finite(r$rho) || r$rho >= 1)
      return(base::list(valid = FALSE, C1 = NULL, C2 = NULL, rho = if (base::is.null(r)) NA_real_ else r$rho))
    base::list(valid = TRUE, C1 = r$T1, C2 = r$T2, rho = r$rho)
  }

  if (print.trace)
    base::message(base::sprintf("  Parametric bootstrap (fiml): B=%d, cores=%d, threshold=%.3g, ci.level=%.2f",
                                B, base::as.integer(cores), threshold, ci.level))
  if (boot.null) {
    res_null <- .nmfkc.parlapply(base::seq_len(B), boot_null_one, cores = cores, envir = base::environment())
    LR.boot <- base::do.call(base::rbind, res_null)
    nnz.boot <- LR.boot[, "nnz"]
    conv.null.boot <- LR.boot[, "conv.null"]; conv.full.boot <- LR.boot[, "conv.full"]
    LR.boot <- LR.boot[, c("full", "selected"), drop = FALSE]
    ok <- base::is.finite(LR.boot[, "full"]) & base::is.finite(LR.boot[, "selected"])
    n.ok <- base::sum(ok)
    ## (1 + #)/(1 + B_ok), not the raw proportion: with a strongly significant
    ## statistic no replicate exceeds LR_obs and the raw proportion is exactly 0,
    ## which is not a valid bootstrap p-value (Davison & Hinkley 1997, sec. 4.2).
    ## CONVENTIONS.md 6 already forbids the degenerate "returns p = 0" form for
    ## coefficients; this is the same rule for the LR statistic.  The floor is
    ## 1/(1 + B_ok), so B controls the smallest reportable p-value.
    LR.p.boot <- c(
      full     = (1 + base::sum(LR.boot[ok, "full"]     >= LR_obs[["full"]]))     / (1 + n.ok),
      selected = (1 + base::sum(LR.boot[ok, "selected"] >= LR_obs[["selected"]])) / (1 + n.ok))
    LR.null.quantile <- c(full = stats::quantile(LR.boot[ok, "full"], 0.95, type = 8, names = FALSE),
                          selected = stats::quantile(LR.boot[ok, "selected"], 0.95, type = 8, names = FALSE))
    prob.select.null <- base::mean(nnz.boot[ok] > 0)
    ## How many null replicates failed to meet the optimizer tolerance.  These are
    ## kept in the calibration (an early stop is still a draw from the procedure),
    ## but the user is told, because a large share means the null likelihood is
    ## flat and the calibration deserves a sensitivity check on `fiml.maxit`/`factr`.
    n.nonconv <- c(null = base::sum(conv.null.boot[ok] != 0, na.rm = TRUE),
                   full = base::sum(conv.full.boot[ok] != 0, na.rm = TRUE))
    if (base::max(n.nonconv) > 0.1 * base::max(n.ok, 1L))
      base::warning(base::sprintf(
        "%d / %d null-bootstrap replicates did not meet the optimizer tolerance (null fit) and %d / %d (feedback fit); the null likelihood may be flat. Consider raising `fiml.maxit` or lowering `factr` and checking that LR.p.boot is stable.",
        n.nonconv[["null"]], n.ok, n.nonconv[["full"]], n.ok))
  } else {
    LR.boot <- NULL; nnz.boot <- NULL; LR.p.boot <- NULL; LR.null.quantile <- NULL
    prob.select.null <- NULL; n.nonconv <- NULL
    conv.null.boot <- NULL; conv.full.boot <- NULL
  }
  res_sel <- .nmfkc.parlapply(base::seq_len(B), boot_sel_one, cores = cores, envir = base::environment())

  C1.array <- base::array(NA_real_, dim = c(B, Q, P1))
  C2.array <- base::array(NA_real_, dim = c(B, Q, P2))
  rho.vec <- base::rep(NA_real_, B); valid.vec <- base::logical(B)
  for (b in base::seq_len(B)) {
    r <- res_sel[[b]]
    valid.vec[b] <- base::isTRUE(r$valid)
    if (!base::is.null(r$rho)) rho.vec[b] <- r$rho
    if (valid.vec[b]) { C1.array[b, , ] <- r$C1; C2.array[b, , ] <- r$C2 }
  }
  n.valid <- base::sum(valid.vec)
  if (n.valid < 10L && n.valid < B)
    base::warning(base::sprintf("Only %d / %d bootstrap replicates were valid; CIs / support rates may be unreliable.", n.valid, B))

  alpha <- 1 - ci.level
  apply_finite <- function(arr, FUN) {
    base::apply(arr, c(2, 3), function(v) { v <- v[base::is.finite(v)]; if (!base::length(v)) NA_real_ else FUN(v) })
  }
  ## Centred (basic) percentile interval: [2 hat - q_{1-a/2}, 2 hat - q_{a/2}]
  ## (CONVENTIONS.md 6: invert the centred replicate distribution).
  q_lo <- function(arr) apply_finite(arr, function(v) stats::quantile(v, alpha / 2, names = FALSE))
  q_hi <- function(arr) apply_finite(arr, function(v) stats::quantile(v, 1 - alpha / 2, names = FALSE))
  C1.ci.lower <- 2 * T1_s - q_hi(C1.array); C1.ci.upper <- 2 * T1_s - q_lo(C1.array)
  C2.ci.lower <- 2 * T2_s - q_hi(C2.array); C2.ci.upper <- 2 * T2_s - q_lo(C2.array)
  C1.support <- apply_finite(C1.array, function(v) base::mean(base::abs(v) > threshold))
  C1.support[support == 0] <- 0        # excluded / unselected entries are fixed at 0
  C2.support <- apply_finite(C2.array, function(v) base::mean(base::abs(v) > threshold))
  sig.from.support <- function(s) {
    base::ifelse(!base::is.finite(s), " ",
      base::ifelse(s > 0.999, "***", base::ifelse(s > 0.99, "**", base::ifelse(s > 0.95, "*", " "))))
  }
  Q_lab <- base::rownames(object$C1); Y1_lab <- base::colnames(object$C1); Y2_lab <- base::colnames(object$C2)
  build_block <- function(type_label, Mhat, sup, lo, hi, basislabs, varlabs) {
    s_vec <- base::as.vector(sup)
    base::data.frame(
      Type = type_label,
      Basis = base::rep(basislabs, times = base::ncol(Mhat)),
      Covariate = base::rep(varlabs, each = base::nrow(Mhat)),
      Estimate = base::as.vector(Mhat),
      CI_low = base::as.vector(lo), CI_high = base::as.vector(hi),
      support_rate = s_vec,
      prob.unsupported = base::ifelse(base::is.finite(s_vec), 1 - s_vec, NA_real_),
      p_value = base::ifelse(base::is.finite(s_vec), 1 - s_vec, NA_real_),
      sig = sig.from.support(s_vec),
      stringsAsFactors = FALSE)
  }
  coefficients <- base::rbind(
    build_block("C1", T1_s, C1.support, C1.ci.lower, C1.ci.upper, Q_lab, Y1_lab),
    build_block("C2", T2_s, C2.support, C2.ci.lower, C2.ci.upper, Q_lab, Y2_lab))
  base::rownames(coefficients) <- NULL
  if (print.trace)
    base::message(base::sprintf("  Bootstrap done: %d / %d valid replicates.", n.valid, B))

  object$bootstrap.B <- B
  object$bootstrap.threshold <- threshold
  object$bootstrap.ci.level <- ci.level
  object$bootstrap.n.valid <- n.valid
  object$bootstrap.n.invalid <- B - n.valid
  object$bootstrap.type <- "parametric (fiml)"
  object$LR.boot <- LR.boot
  object$nnz.boot <- nnz.boot
  object$LR.p.boot <- LR.p.boot
  object$LR.null.quantile <- LR.null.quantile
  object$prob.select.null <- prob.select.null
  object$LR.boot.n.nonconv <- n.nonconv
  object$LR.boot.n.ok <- if (boot.null) n.ok else NULL
  object$rho.boot <- rho.vec
  object$C1.array <- C1.array
  object$C2.array <- C2.array
  object$C1.support.rate <- C1.support
  object$C2.support.rate <- C2.support
  object$C1.ci.lower <- C1.ci.lower; object$C1.ci.upper <- C1.ci.upper
  object$C2.ci.lower <- C2.ci.lower; object$C2.ci.upper <- C2.ci.upper
  object$coefficients <- coefficients
  if (!base::inherits(object, "nmf.sem.inference"))
    base::class(object) <- c("nmf.sem.inference", base::class(object))
  if (!base::inherits(object, "nmf.ffb.inference"))
    base::class(object) <- c("nmf.ffb.inference", "nmf.inference", base::class(object))
  object
}

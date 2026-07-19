# =====================================================================
# nmf.cox.R -- NMF-COX: separating a proportional-hazards effect from a
#              low-rank, non-negative time-varying structure in the Cox model.
#
# Model
#   lambda_n(t) = h0(t) * exp( z_n' gamma + a_n' beta(t) ),  beta(t) = C_mat' x(t)
#   with a non-negative time basis X >= 0 (P x Q, columns summing to one) and a
#   sign-free loading C_mat (Q x R): a semi-NMF representation estimated directly
#   on the Cox partial likelihood (Breslow ties).
#
#   The covariates of interest z keep a time-constant effect gamma (the estimand);
#   the non-proportional covariates a receive the time-varying offset a' beta(t).
#   The offset is linear in a at each event time, which is what makes the
#   cross-fit (Neyman-orthogonal) inference for gamma in nmf.cox.cf applicable.
#
# Entry points (nmfkc house style: optimization / inference split)
#   nmf.cox           full-sample fit, optimization only (alternating block
#                     updates; projected block-Newton for the X >= 0 block)
#   nmf.cox.inference given-basis Wald inference on C and beta(t) = X C
#   nmf.cox.cv        selection of rank Q and smoothing X.L2.smooth (CVPL/AIC/BIC)
#   nmf.cox.cf        two-stage cross-fit estimator for gamma (cluster-robust SE)
#   nmf.cox.phtest    cross-fitted Wald test of H0: beta_r(t) = 0
#
#   The parameter matrix is called C (Q x R), matching the nmfkc family
#   (Theta in the paper). Object classes are single ("nmf.cox", ...).
#   Do NOT re-add an "nmf" class: nmfkc's coef.nmf/fitted.nmf look for
#   $coefficients/$XB, which these objects do not have.
#
# Reference: Satoh K. NMF-COX (2026).
# =====================================================================

# --------------------------------------------------------------------
# nmf.cox : fit
# --------------------------------------------------------------------

#' @title Fit NMF-COX: Cox model with a low-rank non-negative time-varying offset
#' @description
#' \code{nmf.cox} fits the model
#' \deqn{\lambda_n(t) = h_0(t)\exp(\bm{z}_n^\top\gamma + \bm{a}_n^\top\beta(t)),
#'       \qquad \beta(t) = \Theta^\top \bm{x}(t),\quad \bm{x}(t)\ge 0,}
#' on the Cox partial likelihood (Breslow ties). The covariates of interest
#' \code{z} (right-hand side of \code{formula}) keep a time-constant effect
#' \eqn{\gamma}; the non-proportional covariates \code{a} (given by \code{A})
#' receive a low-rank time-varying offset built from a non-negative time basis
#' \eqn{X} (event times by rank) and a sign-free loading \eqn{\Theta}
#' (a semi-NMF representation).
#'
#' @details
#' The penalised partial log-likelihood
#' \deqn{\ell(\gamma,X,\Theta) - \tfrac{\lambda_X}{2}\sum_{p\ge2} w_p
#'       \|(\bm{x}_p-\bm{x}_{p-1})^\top\Theta\|^2}
#' is optimised by \emph{alternating block updates}:
#' \itemize{
#'   \item the \eqn{(\gamma,\Theta)} block, a penalised Cox (Newton) update with
#'         the interactions \eqn{x_{pq}a_{rn}} as covariates (concave);
#'   \item the \eqn{X\ge0} block, solved by a projected block-Newton method with an
#'         active set, Armijo backtracking and a projected-gradient KKT stopping
#'         rule (\code{X.update = "newton"}, the default), or by L-BFGS-B
#'         (\code{"optim"}).
#' }
#' The \eqn{(\gamma,\Theta)} Newton update is not line-searched, so monotone
#' ascent is not guaranteed in principle; the objective is monitored at every
#' iteration, the best iterate is kept, and iteration stops at a plateau.
#'
#' Only \code{ties = "breslow"} is supported (Efron is not implemented in the
#' \eqn{X} block). Identification: each row of \code{A} is standardised across
#' subjects and each column of \eqn{X} is normalised (\code{X.restriction}),
#' moving scale into \eqn{\Theta}; \eqn{\beta(t)=X\Theta} is invariant.
#'
#' \code{nmf.cox} performs \strong{optimization only}, mirroring the
#' \code{\link{nmfkc}} / \code{\link{nmfkc.inference}} split: run
#' \code{\link{nmf.cox.inference}} on the fitted object for standard errors of
#' \eqn{\beta(t)}, the coefficients table of \eqn{C} and the Wald tests.
#'
#' @param formula A model formula whose left-hand side is a \code{\link[survival]{Surv}}
#'   object and whose right-hand side gives the proportional-hazards covariates \code{z}.
#' @param data A data frame containing the variables in \code{formula} and \code{A}.
#' @param A The non-proportional covariates \code{a}: either a one-sided formula
#'   (e.g. \code{~ age + grade}) or a numeric matrix with \code{N} columns (R x N).
#' @param rank Rank \eqn{Q} of the time basis. Default \code{2}.
#' @param X.L2.smooth Non-negative roughness penalty \eqn{\lambda_X} on
#'   \eqn{\beta(t) = XC} (adjacent-event-time differences). Default \code{0}.
#' @param epsilon Convergence tolerance on the relative change of the objective.
#' @param maxit Maximum number of outer iterations. Default \code{50}.
#' @param verbose Logical; print per-iteration progress.
#' @param ... Additional arguments:
#'   \itemize{
#'     \item \code{X.init}: initialisation of \eqn{X}: \code{"kmeans++"}
#'       (default), \code{"nndsvd"}, a numeric matrix, or \code{NULL} for a
#'       deterministic start.
#'     \item \code{X.restriction}: column normalisation of \eqn{X}:
#'       \code{"colSums"} (default, composition probabilities),
#'       \code{"colSqSums"}, \code{"totalSum"}, \code{"none"} or \code{"fixed"}.
#'     \item \code{X.update}: solver for the \eqn{X} block: \code{"newton"}
#'       (default), \code{"optim"}, \code{"poisson"} or \code{"MU"}.
#'     \item \code{X.maxit}: maximum inner iterations of the \eqn{X}-block
#'       solver (default 100).
#'     \item \code{nstart}: number of restarts passed to the initialiser (default 1).
#'     \item \code{seed}: random seed for the initialiser (default 123).
#'     \item \code{prefix}: prefix for basis names (default \code{"Basis"}).
#'     \item \code{smooth.time}: spacing used by the penalty weights:
#'       \code{"gap"} (default), \code{"log"} or \code{"uniform"}.
#'     \item \code{A.normalize}: scaling of \code{A}: \code{"center"} (default,
#'       per-SD), \code{"minmax"} or \code{"none"}.
#'     \item \code{relax}: relaxation factor for the \eqn{\gamma} update (default 1).
#'     \item \code{patience}: iterations without improvement before stopping (default 5).
#'     \item \code{nonneg}: logical; \code{FALSE} drops the sign constraint on
#'       \eqn{X} (a sign-free reduced-rank variant; default \code{TRUE}).
#'     \item \code{aic}: logical; also return an effective-df partial-likelihood
#'       AIC (default \code{FALSE}).
#'     \item \code{ties}: tie handling; only \code{"breslow"} is supported.
#'   }
#'
#' @return An object of class \code{"nmf.cox"}: a list with components including
#'   \code{gamma} (the proportional-hazards estimand), \code{C} (the sign-free
#'   parameter matrix \eqn{\Theta} of the paper, Q x R), \code{X} (non-negative
#'   time basis), \code{X.prob} (columns summing to one), \code{X.cluster},
#'   \code{beta.t} (\eqn{\beta(t)=XC}, per-SD log hazard ratio),
#'   \code{cox.fit}, \code{event.times}, \code{objfunc}, \code{objfunc.iter},
#'   \code{iter}, \code{converged}, \code{stop.reason} and \code{runtime}.
#'
#' @seealso \code{\link{nmf.cox.inference}} for given-basis Wald inference,
#'   \code{\link{nmf.cox.cv}} for selecting \code{rank} and \code{X.L2.smooth},
#'   \code{\link{nmf.cox.cf}} for cross-fit inference on \eqn{\gamma}.
#'
#' @importFrom survival Surv coxph coxph.control
#' @examples
#' library(survival)
#' d <- survival::veteran
#' d$y <- survival::Surv(d$time, d$status)
#' fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                rank = 2, X.L2.smooth = 1000, verbose = FALSE)
#' coef(fit)
#' head(fit$beta.t)
#'
#' @export
nmf.cox <- function(formula, data, A, rank = 2, X.L2.smooth = 0,
                    epsilon = 1e-4, maxit = 50, verbose = TRUE, ...) {
  extra_args <- base::list(...)
  X.init        <- if (!is.null(extra_args$X.init))        extra_args$X.init        else "kmeans++"
  X.restriction <- if (!is.null(extra_args$X.restriction)) extra_args$X.restriction else "colSums"
  X.update      <- if (!is.null(extra_args$X.update))      extra_args$X.update      else "newton"
  X.maxit       <- if (!is.null(extra_args$X.maxit))       extra_args$X.maxit       else 100
  nstart        <- if (!is.null(extra_args$nstart))        extra_args$nstart        else 1
  seed          <- if (!is.null(extra_args$seed))          extra_args$seed          else 123
  prefix        <- if (!is.null(extra_args$prefix))        extra_args$prefix        else "Basis"
  smooth.time   <- if (!is.null(extra_args$smooth.time))   extra_args$smooth.time   else "gap"
  A.normalize   <- if (!is.null(extra_args$A.normalize))   extra_args$A.normalize   else "center"
  relax         <- if (!is.null(extra_args$relax))         extra_args$relax         else 1
  patience      <- if (!is.null(extra_args$patience))      extra_args$patience      else 5
  nonneg        <- if (!is.null(extra_args$nonneg))        extra_args$nonneg        else TRUE
  aic           <- if (!is.null(extra_args$aic))           extra_args$aic           else FALSE
  ties          <- if (!is.null(extra_args$ties))          extra_args$ties          else "breslow"
  A.normalize <- match.arg(A.normalize, c("center", "minmax", "none"))
  start.time <- Sys.time()
  X.restriction <- match.arg(X.restriction, c("colSums", "colSqSums", "totalSum", "none", "fixed"))
  if (!identical(ties, "breslow"))
    stop("nmf.cox: only ties='breslow' is supported (Efron is not implemented in the X-block).")
  Q <- rank
  .eps <- 1e-10

  mf <- stats::model.frame(formula, data)
  y  <- stats::model.response(mf)
  if (!inherits(y, "Surv")) stop("LHS must be a Surv() object.")
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2])
  N <- length(time)
  Z <- stats::model.matrix(attr(mf, "terms"), mf)
  Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  p <- ncol(Z); znames <- colnames(Z)

  if (inherits(A, "formula")) {
    Amf <- stats::model.frame(A, data)
    Amat <- stats::model.matrix(A, Amf)
    Amat <- Amat[, colnames(Amat) != "(Intercept)", drop = FALSE]
    A <- t(Amat)
  } else {
    A <- as.matrix(A)
    if (ncol(A) != N && nrow(A) == N) A <- t(A)
    if (ncol(A) != N) stop("A must have N columns (R x N).")
  }
  R <- nrow(A); anames <- rownames(A); if (is.null(anames)) anames <- paste0("A", 1:R)
  if (A.normalize == "center") {
    A.center <- rowMeans(A); A.scale <- apply(A, 1, sd); A.scale[A.scale < .eps] <- 1
  } else if (A.normalize == "minmax") {
    A.center <- apply(A, 1, min); A.rng <- apply(A, 1, function(v) diff(range(v)))
    A.scale <- ifelse(A.rng < .eps, 1, A.rng)
  } else {                                              # "none"
    A.center <- rep(0, R); A.scale <- rep(1, R)
  }
  A <- (A - A.center) / A.scale

  et <- sort(unique(time[status == 1])); P <- length(et)
  s_prev <- c(0, et[-P])
  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in seq_len(P)) {
    risk <- which(time >= et[j])
    m <- length(risk)
    rid   <- c(rid, risk)
    rj    <- c(rj, rep(j, m))
    rstart<- c(rstart, rep(s_prev[j], m))
    rstop <- c(rstop, rep(et[j], m))
    rev   <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j]))
  }
  rstart <- as.numeric(rstart); rstop <- as.numeric(rstop)
  Zexp <- Z[rid, , drop = FALSE]
  rows_by_j <- split(seq_along(rj), rj)

  xscale <- switch(X.restriction,
    colSums   = function(M) colSums(M) + .eps,
    colSqSums = function(M) sqrt(colSums(M * M)) + .eps,
    totalSum  = function(M) rep(sum(M) + .eps, ncol(M)),
    none      = function(M) rep(1, ncol(M)),
    fixed     = function(M) rep(1, ncol(M)))
  normX <- function(M) sweep(M, 2, pmax(xscale(M), .eps), "/")
  det_init <- function() {
    tt <- (et - min(et)) / (max(et) - min(et) + .eps)
    Xd <- matrix(0, P, Q); Xd[, 1] <- 1
    if (Q >= 2) for (q in 2:Q) Xd[, q] <- tt^(q - 1)
    Xd + 1e-3 * outer(tt, seq_len(Q), function(t, q) cos(q * pi * t))
  }
  if (is.matrix(X.init) || (is.numeric(X.init) && length(X.init) > 1)) {
    X <- as.matrix(X.init)
  } else if (is.character(X.init)) {
    a0 <- tryCatch(coef(survival::coxph(y ~ Z, ties = ties)), error = function(e) rep(0, p))
    sig <- matrix(0, P, R)
    for (j in seq_len(P)) {
      rows <- rows_by_j[[j]]; ids <- rid[rows]
      wj <- exp(as.vector(Zexp[rows, , drop = FALSE] %*% a0)); wj <- wj / sum(wj)
      abar <- as.vector(A[, ids, drop = FALSE] %*% wj)
      fid <- rid[rows[rev[rows] == 1]]
      sig[j, ] <- rowMeans(A[, fid, drop = FALSE]) - abar
    }
    Xk <- tryCatch(.init_X_method(X.init, sig, Q, seed = seed, nstart = nstart),
                   error = function(e) NULL)
    if (is.null(Xk)) Xk <- det_init()
    Xk[Xk < 0] <- 0
    bad <- colSums(Xk) < .eps; if (any(bad)) Xk[, bad] <- 1 / P
    X <- Xk
  } else X <- det_init()
  X <- normX(X)
  gamma <- tryCatch(coef(survival::coxph(y ~ Z, ties = ties)), error = function(e) rep(0, p))
  names(gamma) <- znames
  C_mat <- matrix(0, Q, R)

  sm.tt <- switch(smooth.time,
    uniform = as.numeric(seq_len(P)),
    gap     = (et - min(et)) / (max(et) - min(et) + .eps),
    log     = { lg <- log(pmax(et, .eps)); (lg - min(lg)) / (max(lg) - min(lg) + .eps) },
    stop("smooth.time must be 'gap', 'log', or 'uniform'"))
  sm.gap <- diff(sm.tt); if (length(sm.gap)) sm.gap <- pmax(sm.gap, max(sm.gap) * 1e-3)
  sm.w <- if (length(sm.gap)) { w <- 1 / sm.gap; w / mean(w) } else numeric(0)

  Sfun <- function(Th) Th %*% A
  negll_grad_X <- function(xvec, S, eta0, Msmooth) {
    Xc <- matrix(xvec, P, Q)
    negll <- 0; g <- matrix(0, P, Q)
    for (j in seq_len(P)) {
      rows <- rows_by_j[[j]]
      ids  <- rid[rows]
      Sj   <- S[, ids, drop = FALSE]
      eta  <- eta0[rows] + as.vector(crossprod(Sj, Xc[j, ]))  # length m
      mx <- max(eta); w <- exp(eta - mx); sw <- sum(w); wn <- w / sw
      dfail <- which(rev[rows] == 1); dj <- length(dfail)
      negll <- negll - (sum(eta[dfail]) - dj * (mx + log(sw)))
      sbar <- as.vector(Sj %*% wn)
      g[j, ] <- g[j, ] - (rowSums(Sj[, dfail, drop = FALSE]) - dj * sbar)
    }
    if (X.L2.smooth > 0) {
      D  <- diff(Xc)
      DM <- (D %*% Msmooth) * sm.w
      negll <- negll + 0.5 * X.L2.smooth * sum(DM * D)
      gp <- matrix(0, P, Q)
      gp[1:(P-1), ] <- gp[1:(P-1), ] - DM
      gp[2:P, ]     <- gp[2:P, ]     + DM
      g <- g + X.L2.smooth * gp
    }
    list(value = negll, grad = as.vector(g))
  }

  gamma.history <- gamma
  converged <- FALSE
  conv.msg  <- sprintf("reached maxit=%d without meeting tol/plateau", maxit)
  opt.conv  <- NA_integer_
  obj_old   <- NA_real_
  beta_old  <- NULL
  best.obj  <- -Inf
  best.gamma<- gamma; best.C_mat <- NULL; best.X <- X; best.it <- 0L
  stall     <- 0L
  objfunc.iter <- rep(NA_real_, maxit)
  for (it in seq_len(maxit)) {
    gamma_old <- gamma
    Fmat <- matrix(0, length(rj), Q * R)
    cn <- character(Q * R); col <- 1
    for (q in 1:Q) for (l in 1:R) {
      Fmat[, col] <- X[rj, q] * A[l, rid]
      cn[col] <- paste0("th_", q, "_", l); col <- col + 1
    }
    colnames(Fmat) <- cn
    dat <- data.frame(start = rstart, stop = rstop, ev = rev, id = rid)
    des <- cbind(Zexp, Fmat)
    ## ---- PENALIZED (gamma,C_mat)-block: maximise loglik - 0.5*lambda_X * sum_l theta_l' M_X theta_l ----
    ##   => both blocks optimise the SAME penalized partial likelihood (genuine block-coordinate ascent on eq:pll;
    ##      objfunc becomes monotone). Penalized Newton reusing coxph's score/information at the current iterate.
    MX <- if (X.L2.smooth > 0 && P >= 2) crossprod(diff(X) * sqrt(pmax(sm.w, 0))) else matrix(0, Q, Q)
    Ppad <- matrix(0, p + Q * R, p + Q * R)
    if (X.L2.smooth > 0) Ppad[(p + 1):(p + Q * R), (p + 1):(p + Q * R)] <- X.L2.smooth * kronecker(MX, diag(R))
    bcoef <- c(as.numeric(gamma), as.vector(t(C_mat)))          # vecR(C_mat): q outer, l inner (matches Fmat cols)
    okfit <- FALSE
    for (nb in 1:8) {
      f0 <- tryCatch(suppressWarnings(survival::coxph(survival::Surv(start, stop, ev) ~ des, data = dat, ties = ties,
                     init = bcoef, control = survival::coxph.control(iter.max = 0))), error = function(e) NULL)
      if (is.null(f0)) break
      U  <- tryCatch(colSums(residuals(f0, type = "score")), error = function(e) NULL)
      Hn <- tryCatch(solve(f0$var), error = function(e) NULL)
      if (is.null(U) || is.null(Hn) || any(!is.finite(U)) || any(!is.finite(Hn))) break
      g  <- U - as.vector(Ppad %*% bcoef)
      st <- tryCatch(solve(Hn + Ppad, g), error = function(e) NULL)
      if (is.null(st) || any(!is.finite(st))) break
      bcoef <- bcoef + st; okfit <- TRUE
      if (max(abs(st)) < 1e-7) break
    }
    if (!okfit) {                              # fallback: unpenalized coxph (heuristic behaviour)
      fit <- tryCatch(survival::coxph(survival::Surv(start, stop, ev) ~ des + cluster(id), data = dat, ties = ties), error = function(e) NULL)
      if (is.null(fit)) {
        if (it == 1)
          stop("nmf.cox: (gamma,C_mat)-block failed at the first iteration; cannot produce a solution.")
        conv.msg <- sprintf("(gamma,C_mat)-block failed at iter %d; returning last successful state", it)
        converged <- FALSE; break
      }
      cf <- coef(fit); cf[is.na(cf)] <- 0; bcoef <- cf
    }
    gamma_new <- bcoef[1:p]; names(gamma_new) <- znames
    C_mat <- matrix(bcoef[(p + 1):(p + Q * R)], Q, R, byrow = TRUE)  # rows q, cols l

    S <- Sfun(C_mat)
    eta0 <- as.vector(Zexp %*% gamma_new)
    Msm <- tcrossprod(C_mat)
    optval <- NA_real_; opt.conv <- 0L
    if (X.update == "MU" && nonneg) {
      fv  <- function(v) negll_grad_X(v, S, eta0, Msm)$value
      cur <- fv(as.vector(X)); stall.mu <- 0L
      for (mit in seq_len(X.maxit)) {
        G  <- matrix(negll_grad_X(as.vector(X), S, eta0, Msm)$grad, P, Q)
        gp <- pmax(G, 0); gm <- pmax(-G, 0)
        ratio <- pmin(pmax(gm / (gp + .eps), 1/5), 5)
        dX <- X * ratio - X
        t.step <- 1; fn <- fv(as.vector(X + dX))
        while (fn > cur && t.step > 1e-4) { t.step <- t.step / 2; fn <- fv(as.vector(X + t.step * dX)) }
        if (fn < cur - abs(cur) * 1e-8) { X <- X + t.step * dX; cur <- fn; stall.mu <- 0L }
        else { stall.mu <- stall.mu + 1L; if (stall.mu >= 3) break }
      }
      optval <- cur
    } else if (X.update == "poisson" && nonneg) {
      fv <- function(v) negll_grad_X(v, S, eta0, Msm)$value
      cur <- fv(as.vector(X)); st <- 0L
      for (mm in seq_len(X.maxit)) {
        G   <- matrix(negll_grad_X(as.vector(X), S, eta0, Msm)$grad, P, Q)
        dir <- matrix(0, P, Q)
        for (j in seq_len(P)) {
          rows <- rows_by_j[[j]]; ids <- rid[rows]
          Sj  <- S[, ids, drop = FALSE]
          eta <- eta0[rows] + as.vector(crossprod(Sj, X[j, ]))
          w   <- exp(eta - max(eta)); w <- w / sum(w)
          dj  <- sum(rev[rows] == 1)
          Hj  <- dj * ((Sj * rep(w, each = Q)) %*% t(Sj))
          dir[j, ] <- solve(Hj + diag(1e-6, Q), G[j, ])
        }
        t.step <- 1; Xn <- pmax(X - t.step * dir, 0); fn <- fv(as.vector(Xn))
        while (fn > cur && t.step > 1e-5) { t.step <- t.step / 2; Xn <- pmax(X - t.step * dir, 0); fn <- fv(as.vector(Xn)) }
        if (fn < cur - abs(cur) * 1e-10) { X <- Xn; cur <- fn; st <- 0L }
        else { st <- st + 1L; if (st >= 8) break }
      }
      optval <- cur
    } else if (X.update == "newton" && nonneg) {
      fv <- function(v) negll_grad_X(v, S, eta0, Msm)$value
      idxc <- function(a) (a - 1) * P + seq_len(P)
      Hsm <- matrix(0, P * Q, P * Q)
      if (X.L2.smooth > 0 && P >= 2) {
        Dm <- matrix(0, P - 1, P); for (j in 1:(P - 1)) { Dm[j, j] <- -1; Dm[j, j + 1] <- 1 }
        Lsm <- crossprod(Dm * sqrt(pmax(sm.w, 0)))
        Hsm <- X.L2.smooth * kronecker(Msm, Lsm)
      }
      for (nit in seq_len(X.maxit)) {
        cur <- fv(as.vector(X))
        if (!is.finite(cur)) break
        Gv <- negll_grad_X(as.vector(X), S, eta0, Msm)$grad
        if (any(!is.finite(Gv))) break
        Hd <- matrix(0, P * Q, P * Q)
        for (j in seq_len(P)) {
          rows <- rows_by_j[[j]]; ids <- rid[rows]; Sj <- S[, ids, drop = FALSE]
          eta <- pmax(pmin(eta0[rows] + as.vector(crossprod(Sj, X[j, ])), 30), -30)
          w <- exp(eta - max(eta)); w <- w / sum(w); dj <- sum(rev[rows] == 1)
          sbar <- as.vector(Sj %*% w)
          Vj <- dj * ((Sj * rep(w, each = Q)) %*% t(Sj) - outer(sbar, sbar))  # d_j Var_w(s)
          ii <- vapply(seq_len(Q), function(a) (a - 1) * P + j, numeric(1))
          Hd[ii, ii] <- Hd[ii, ii] + Vj
        }
        H <- Hd + Hsm
        Xv <- as.vector(X); active <- (Xv <= .eps) & (Gv > 0)      # active-set
        free <- which(!active); d <- numeric(P * Q)
        d[free] <- tryCatch(solve(H[free, free, drop = FALSE] + diag(1e-8, length(free)), Gv[free]),
                            error = function(e) Gv[free])
        if (any(!is.finite(d))) break
        slope <- sum(Gv[free] * d[free]); if (!is.finite(slope) || slope <= 0) slope <- 1e-12
        t.step <- 1; fn <- fv(pmax(Xv - t.step * d, 0))
        while (!isTRUE(fn <= cur - 1e-4 * t.step * slope) && t.step > 1e-6) {
          t.step <- t.step / 2; fn <- fv(pmax(Xv - t.step * d, 0))
        }
        if (!is.finite(fn) || fn >= cur) break
        X <- matrix(pmax(Xv - t.step * d, 0), P, Q)
        kkt <- max(abs(as.vector(X) - pmax(as.vector(X) - Gv, 0)))
        if (isTRUE(kkt < 1e-6 * max(1, max(abs(Gv))))) break
      }
      optval <- fv(as.vector(X)); if (!is.finite(optval)) optval <- 1e12
    } else {
      opt <- tryCatch(stats::optim(par = as.vector(X), method = "L-BFGS-B",
                   lower = if (nonneg) 0 else -Inf,
                   fn = function(v) negll_grad_X(v, S, eta0, Msm)$value,
                   gr = function(v) negll_grad_X(v, S, eta0, Msm)$grad,
                   control = list(maxit = 100)),
               error = function(e) NULL)
      if (is.null(opt)) {
        if (it == 1)
          stop("nmf.cox: X-block L-BFGS-B failed at the first iteration; cannot produce a solution.")
        conv.msg <- sprintf("X-block L-BFGS-B failed at iter %d; returning last successful state", it)
        converged <- FALSE; break
      }
      X <- matrix(opt$par, P, Q); opt.conv <- opt$convergence; optval <- opt$value
    }
    if (X.restriction != "fixed" && X.restriction != "none") {
      sc <- pmax(xscale(X), .eps)
      X <- sweep(X, 2, sc, "/"); C_mat <- sweep(C_mat, 1, sc, "*")
    }

    gamma <- (1 - relax) * gamma_old + relax * gamma_new
    gamma.history <- rbind(gamma.history, gamma)
    beta_new <- X %*% C_mat
    obj_new  <- -optval
    objfunc.iter[it] <- -obj_new
    dgamma <- sqrt(sum((gamma - gamma_old)^2)) / max(1, sqrt(sum(gamma_old^2)))
    dbeta  <- if (is.null(beta_old)) Inf else
              sqrt(sum((beta_new - beta_old)^2)) / max(1, sqrt(sum(beta_old^2)))
    dobj   <- if (is.na(obj_old)) Inf else abs(obj_new - obj_old) / max(1, abs(obj_old))
    beta_old <- beta_new; obj_old <- obj_new
    if (verbose) cat(sprintf("iter %2d: |d.gamma|=%.2e |d.beta|=%.2e |d.pll|=%.2e obj=%.4f\n",
                             it, dgamma, dbeta, dobj, obj_new))
    scale.obj <- if (is.finite(best.obj)) max(1, abs(best.obj)) else 1
    if (is.finite(obj_new) && obj_new > best.obj + scale.obj * epsilon) {
      best.obj <- obj_new; best.gamma <- gamma; best.C_mat <- C_mat; best.X <- X; best.it <- it; stall <- 0L
    } else stall <- stall + 1L
    if (is.finite(best.obj) && is.finite(obj_new) && obj_new < best.obj - max(1, abs(best.obj)) * 0.5) {
      conv.msg <- sprintf("halted: penalized loglik diverged at iter %d; returning best iterate (iter %d)", it, best.it)
      converged <- best.it > 0; break
    }
    if (max(dgamma, dbeta, dobj) <= epsilon) {
      converged <- TRUE
      conv.msg  <- "converged: max rel. change of (gamma, beta=X*C_mat, penalized loglik) <= tol"
      break
    }
    if (stall >= patience) {
      converged <- TRUE
      conv.msg  <- sprintf("converged: penalized loglik plateaued (no improvement > tol for %d iters; best iter %d)",
                           patience, best.it)
      break
    }
  }
  if (!is.null(best.C_mat)) { gamma <- best.gamma; C_mat <- best.C_mat; X <- best.X }
  if (Q > 1 && X.restriction != "fixed") {
    ord.b <- order(as.vector((seq_len(nrow(X)) / nrow(X)) %*% X))
    X <- X[, ord.b, drop = FALSE]; C_mat <- C_mat[ord.b, , drop = FALSE]
  }
  best.iter <- best.it
  if (converged && !is.na(opt.conv) && opt.conv != 0)
    conv.msg <- sprintf("%s; note: X-block optim$convergence=%d at some iter", conv.msg, opt.conv)
  if (!converged) warning("nmf.cox did NOT converge: ", conv.msg)

  loglik <- df.eff <- aic.value <- NA_real_; df.parts <- NULL
  if (isTRUE(aic)) {
    S.a <- Sfun(C_mat); eta0.a <- as.vector(Zexp %*% gamma); Msm.a <- tcrossprod(C_mat)
    ng  <- negll_grad_X(as.vector(X), S.a, eta0.a, Msm.a)$value
    pen <- if (X.L2.smooth > 0 && P >= 2) { Dd <- diff(X); 0.5 * X.L2.smooth * sum(((Dd %*% Msm.a) * sm.w) * Dd) } else 0
    loglik <- -(ng - pen)
    Hd <- matrix(0, P * Q, P * Q)
    for (j in seq_len(P)) {
      rows <- rows_by_j[[j]]; ids <- rid[rows]; Sj <- S.a[, ids, drop = FALSE]
      eta <- pmax(pmin(eta0.a[rows] + as.vector(crossprod(Sj, X[j, ])), 30), -30)
      w <- exp(eta - max(eta)); w <- w / sum(w); dj <- sum(rev[rows] == 1)
      sbar <- as.vector(Sj %*% w)
      Vj <- dj * ((Sj * rep(w, each = Q)) %*% t(Sj) - outer(sbar, sbar))
      ii <- (seq_len(Q) - 1) * P + j; Hd[ii, ii] <- Hd[ii, ii] + Vj
    }
    Hsm <- matrix(0, P * Q, P * Q)
    if (X.L2.smooth > 0 && P >= 2) {
      Dm <- matrix(0, P - 1, P); for (j in 1:(P - 1)) { Dm[j, j] <- -1; Dm[j, j + 1] <- 1 }
      Hsm <- X.L2.smooth * kronecker(Msm.a, crossprod(Dm * sqrt(pmax(sm.w, 0))))
    }
    df.X <- tryCatch(sum(diag(solve(Hd + Hsm + diag(1e-8, P * Q), Hd))), error = function(e) NA_real_)
    rot  <- min(Q * Q, Q * R); df.beta <- df.X + Q * R - rot; df.eff <- p + df.beta
    aic.value <- -2 * loglik + 2 * df.eff
    df.parts  <- c(gamma = p, df.X = df.X, theta = Q * R, rotation = rot, df.beta = df.beta, df.eff = df.eff)
    if (verbose) cat(sprintf("AIC: loglik=%.2f  df.eff=%.2f (df.X=%.2f)  AIC=%.2f\n", loglik, df.eff, df.X, aic.value))
  }

  o_final <- numeric(length(rj))
  S <- Sfun(C_mat)
  for (j in seq_len(P)) {
    rows <- rows_by_j[[j]]; ids <- rid[rows]
    o_final[rows] <- as.vector(crossprod(S[, ids, drop = FALSE], X[j, ]))
  }
  o_final <- pmax(pmin(o_final, 30), -30)
  dat2 <- data.frame(start = rstart, stop = rstop, ev = rev, id = rid, off = o_final)
  Zdf <- as.data.frame(Zexp); colnames(Zdf) <- znames
  dat2 <- cbind(dat2, Zdf)
  fml <- stats::as.formula(paste0("Surv(start,stop,ev) ~ ",
                           paste(znames, collapse = " + "),
                           " + offset(off) + cluster(id)"))
  cox.fit <- survival::coxph(fml, data = dat2, ties = ties)

  beta.t <- X %*% C_mat
  colnames(beta.t) <- anames
  rownames(beta.t) <- format(et)
  rownames(X) <- format(et); colnames(X) <- paste0(prefix, 1:Q)
  rownames(C_mat) <- paste0(prefix, 1:Q); colnames(C_mat) <- anames

  objfunc.iter <- objfunc.iter[seq_len(it)]
  objfunc <- -best.obj
  X.prob    <- sweep(X, 2, pmax(colSums(X), .eps), "/")
  X.cluster <- max.col(X.prob, ties.method = "first")
  dims    <- sprintf("Cox: N=%d, events=%d, P=%d, Q=%d, R=%d", N, sum(status == 1), P, Q, R)
  runtime <- as.numeric(difftime(Sys.time(), start.time, units = "sec"))

  structure(list(
    call = match.call(), dims = dims, runtime = runtime, method = X.update,
    rank = Q, objfunc = objfunc, objfunc.iter = objfunc.iter,
    X = X, X.prob = X.prob, X.cluster = X.cluster, X.restriction = X.restriction,
    gamma = coef(cox.fit)[znames], gamma.raw = gamma,
    C = C_mat, beta.t = beta.t,
    event.times = et,
    cox.fit = cox.fit, offset = o_final, gamma.history = gamma.history,
    X.L2.smooth = X.L2.smooth, iter = it, best.iter = best.iter, best.obj = best.obj,
    loglik = loglik, df.eff = df.eff, aic = aic.value, df.parts = df.parts,
    converged = converged, stop.reason = conv.msg, opt.convergence = opt.conv,
    ties = ties, Z = Z, A = A, A.center = A.center, A.scale = A.scale
  ), class = "nmf.cox")
}


# --------------------------------------------------------------------
# nmf.cox.inference : given-basis Wald inference
# --------------------------------------------------------------------

#' @title Statistical inference for an NMF-COX fit (given basis)
#' @description
#' \code{nmf.cox.inference} adds Wald inference to a fitted \code{\link{nmf.cox}}
#' object, mirroring the \code{\link{nmfkc}} / \code{\link{nmfkc.inference}}
#' split: a cluster-robust Cox refit with the interactions
#' \eqn{x_{pq} a_{rn}} as covariates yields a coefficients table for the
#' parameter matrix \eqn{C}, pointwise standard errors of
#' \eqn{\beta(t) = XC}, and per-covariate / global Wald tests of
#' \eqn{H_0: \beta_r(t) = 0} \strong{conditional on the estimated basis}
#' \eqn{\hat X}.
#'
#' For inference that accounts for basis selection, use the cross-fitted
#' \code{\link{nmf.cox.phtest}} (time-varying effects) and
#' \code{\link{nmf.cox.cf}} (proportional effect \eqn{\gamma}).
#'
#' @param object A fitted object of class \code{"nmf.cox"}.
#' @param formula,data,A The same arguments used in the \code{\link{nmf.cox}} call.
#' @param ... Not used.
#'
#' @return The input \code{object} with class
#'   \code{c("nmf.cox.inference", "nmf.cox")} and added components:
#'   \item{coefficients}{Data frame with columns \code{Basis}, \code{Covariate},
#'     \code{Estimate}, \code{SE}, \code{z_value}, \code{p_value} for the
#'     entries of \eqn{C} (cluster-robust, two-sided).}
#'   \item{se.beta.t}{Pointwise standard errors of \eqn{\beta(t)} (P x R).}
#'   \item{beta.vcov}{List of per-covariate covariance blocks of \eqn{C}.}
#'   \item{wald, wald.global}{Per-covariate and global Wald tests of
#'     \eqn{H_0: \beta_r(t) = 0} given the basis.}
#'   \item{inf.fit}{The underlying \code{\link[survival]{coxph}} refit.}
#'
#' @seealso \code{\link{nmf.cox}}, \code{\link{nmf.cox.phtest}},
#'   \code{\link{nmf.cox.cf}}
#' @examples
#' library(survival)
#' d <- survival::veteran
#' d$y <- survival::Surv(d$time, d$status)
#' fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                rank = 2, X.L2.smooth = 1000, verbose = FALSE)
#' fit <- nmf.cox.inference(fit, y ~ trt + age, data = d, A = ~ karno + celltype)
#' fit$coefficients
#' fit$wald
#' @export
nmf.cox.inference <- function(object, formula, data, A, ...) {
  if (!inherits(object, "nmf.cox")) stop("object must be of class 'nmf.cox'.")
  ties <- if (!is.null(object$ties)) object$ties else "breslow"
  X <- object$X; C_mat <- object$C
  Q <- ncol(X); et <- object$event.times; P <- length(et)

  mf <- stats::model.frame(formula, data); y <- stats::model.response(mf)
  if (!inherits(y, "Surv")) stop("LHS must be a Surv() object.")
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); N <- length(time)
  Z <- stats::model.matrix(attr(mf, "terms"), mf)
  Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  if (inherits(A, "formula")) {
    Am <- stats::model.matrix(A, stats::model.frame(A, data))
    Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]
    A <- t(Am)
  } else {
    A <- as.matrix(A)
    if (ncol(A) != N && nrow(A) == N) A <- t(A)
  }
  R <- nrow(A); anames <- colnames(C_mat)
  ## standardise A exactly as in the fit
  A <- (A - object$A.center) / object$A.scale

  ## risk-set expansion at the fit's event times
  s_prev <- c(0, et[-P])
  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in seq_len(P)) {
    risk <- which(time >= et[j]); m <- length(risk)
    rid   <- c(rid, risk); rj <- c(rj, rep(j, m))
    rstart<- c(rstart, rep(s_prev[j], m)); rstop <- c(rstop, rep(et[j], m))
    rev   <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j]))
  }
  rstart <- as.numeric(rstart); rstop <- as.numeric(rstop)
  Zexp <- Z[rid, , drop = FALSE]

  Finf <- matrix(0, length(rj), Q * R); ci <- 1
  for (q in 1:Q) for (r in 1:R) { Finf[, ci] <- X[rj, q] * A[r, rid]; ci <- ci + 1 }
  quad.wald <- function(theta, Vmat) {
    Vs <- (Vmat + t(Vmat)) / 2
    ee <- eigen(Vs, symmetric = TRUE)
    tol <- max(ee$values, 0) * 1e-8 * length(ee$values)
    keep <- ee$values > tol
    if (!any(keep)) return(c(chisq = NA_real_, df = 0))
    U <- ee$vectors[, keep, drop = FALSE]
    W <- sum((crossprod(U, theta)^2) / ee$values[keep])
    c(chisq = as.numeric(W), df = sum(keep))
  }
  av <- stats::setNames(rep(0, ncol(Zexp)), colnames(Zexp)); cc <- coef(object$cox.fit)
  av[names(cc)] <- ifelse(is.na(cc), 0, cc)
  init.inf <- c(as.numeric(av), as.vector(t(C_mat)))
  inf.fit <- tryCatch(survival::coxph(survival::Surv(rstart, rstop, rev) ~ Zexp + Finf + cluster(rid),
                            ties = ties, init = init.inf,
                            control = survival::coxph.control(iter.max = 200, eps = 1e-9)),
                      error = function(e) NULL)
  if (is.null(inf.fit)) stop("nmf.cox.inference: the inference coxph refit failed.")

  se.beta.t <- matrix(NA_real_, P, R); beta.vcov <- vector("list", R)
  V <- stats::vcov(inf.fit); cf <- coef(inf.fit); nz <- ncol(Zexp); idxTh <- (nz + 1):(nz + Q * R)
  wald <- data.frame(covariate = anames, df = Q, chisq = NA_real_, p.value = NA_real_,
                     stringsAsFactors = FALSE)
  for (r in 1:R) {
    posr <- idxTh[(0:(Q - 1)) * R + r]
    Vr <- V[posr, posr, drop = FALSE]; beta.vcov[[r]] <- Vr
    qw <- quad.wald(cf[posr], Vr)
    wald$chisq[r] <- qw["chisq"]; wald$df[r] <- qw["df"]
    wald$p.value[r] <- if (qw["df"] > 0) stats::pchisq(qw["chisq"], qw["df"], lower.tail = FALSE) else NA_real_
    for (j in 1:P) se.beta.t[j, r] <- sqrt(max(0, as.numeric(crossprod(X[j, ], Vr %*% X[j, ]))))
  }
  qwg <- quad.wald(cf[idxTh], V[idxTh, idxTh, drop = FALSE])
  wald.global <- list(chisq = as.numeric(qwg["chisq"]), df = as.integer(qwg["df"]),
                      p.value = if (qwg["df"] > 0) stats::pchisq(qwg["chisq"], qwg["df"], lower.tail = FALSE) else NA_real_)

  ## coefficients table for C (house style: Basis / Covariate columns)
  estTh <- cf[idxTh]; seTh <- sqrt(pmax(diag(V)[idxTh], 0))
  zTh <- estTh / seTh
  coefficients <- data.frame(
    Basis     = rep(rownames(C_mat), each = R),
    Covariate = rep(anames, Q),
    Estimate  = as.numeric(estTh),
    SE        = as.numeric(seTh),
    z_value   = as.numeric(zTh),
    p_value   = 2 * stats::pnorm(abs(zTh), lower.tail = FALSE),
    stringsAsFactors = FALSE, row.names = NULL)

  colnames(se.beta.t) <- anames; rownames(se.beta.t) <- format(et)
  object$coefficients <- coefficients
  object$se.beta.t <- se.beta.t; object$beta.vcov <- beta.vcov
  object$wald <- wald; object$wald.global <- wald.global
  object$inf.fit <- inf.fit
  if (!inherits(object, "nmf.cox.inference"))
    class(object) <- c("nmf.cox.inference", class(object))
  object
}


# --------------------------------------------------------------------
# nmf.cox.cv : rank / smoothing selection
# --------------------------------------------------------------------

#' @title Select the rank and smoothing parameter for NMF-COX
#' @description
#' \code{nmf.cox.cv} chooses the rank \eqn{Q} and the roughness penalty
#' \eqn{\lambda_X} by the cross-validated partial likelihood (CVPL) of Verweij and
#' van Houwelingen (1993), or by a partial-likelihood AIC / BIC with effective
#' degrees of freedom.
#'
#' @details
#' With \code{criterion = "cvpl"} the fold contribution is
#' \eqn{d_v = \ell(\hat\theta^{-v}) - \ell^{-v}(\hat\theta^{-v})} and CVPL is
#' \emph{maximised}. The returned \code{X.L2.smooth.best} is the CVPL argmax,
#' while \code{X.L2.smooth.best.1se} applies a one-standard-error rule (the most
#' parsimonious, i.e. smoothest, value whose CVPL is within one SE of the
#' maximum). Because the CVPL curve is typically flat near its optimum,
#' selecting within the 1-SE band is usually preferable to the raw maximiser.
#'
#' With \code{criterion = "aic"} or \code{"bic"} a single fit per grid point is used
#' and the criterion is \emph{minimised}; the effective degrees of freedom are a
#' heuristic.
#'
#' @param formula,data,A As in \code{\link{nmf.cox}}.
#' @param rank Integer vector of candidate ranks \eqn{Q} to evaluate. Default \code{2}.
#' @param X.L2.smooth Numeric vector of candidate roughness penalties \eqn{\lambda_X}.
#' @param criterion One of \code{"cvpl"} (default), \code{"aic"} or \code{"bic"}.
#' @param nfolds Number of folds for CVPL. Default \code{5}.
#' @param verbose Logical; print the criterion grid.
#' @param ... Additional arguments passed to \code{\link{nmf.cox}}
#'   (e.g. \code{X.update}, \code{X.restriction}, \code{X.init}, \code{nonneg},
#'   \code{nstart}). Also accepts: \code{seed} (fold assignment, default 123),
#'   \code{mc.cores} (cores for \code{\link[parallel]{mclapply}}, default 1),
#'   \code{maxit} (outer iterations per fit, default 30).
#'
#' @return An object of class \code{"nmf.cox.cv"}. For \code{"cvpl"} it contains
#'   the \code{cvpl} and \code{se} matrices (rank by smoothing),
#'   \code{rank.best}, \code{X.L2.smooth.best}, \code{rank.best.1se} and
#'   \code{X.L2.smooth.best.1se}; for \code{"aic"}/\code{"bic"} it contains
#'   \code{aic}, \code{bic}, \code{loglik}, \code{df.eff}, \code{nevents},
#'   \code{rank.best} and \code{X.L2.smooth.best}.
#'
#' @references
#' Verweij PJM, van Houwelingen HC (1993). Cross-validation in survival analysis.
#' \emph{Statistics in Medicine}, 12(24), 2305-2314.
#'
#' @seealso \code{\link{nmf.cox}}
#'
#' @examples
#' \donttest{
#' library(survival)
#' d <- survival::veteran
#' d$y <- survival::Surv(d$time, d$status)
#' cv <- nmf.cox.cv(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                  rank = 2, X.L2.smooth = c(300, 1000, 3000),
#'                  criterion = "cvpl", nfolds = 5, seed = 1, verbose = FALSE)
#' cv$X.L2.smooth.best
#' }
#'
#' @export
nmf.cox.cv <- function(formula, data, A, rank = 2,
                        X.L2.smooth = c(0, 30, 100, 300, 1000, 3000),
                        criterion = c("cvpl", "aic", "bic"),
                        nfolds = 5, verbose = TRUE, ...) {
  extra_args <- base::list(...)
  seed     <- if (!is.null(extra_args$seed))     extra_args$seed     else 123
  mc.cores <- if (!is.null(extra_args$mc.cores)) extra_args$mc.cores else 1
  maxit    <- if (!is.null(extra_args$maxit))    extra_args$maxit    else 30
  ## remaining ... are forwarded to nmf.cox (fitter-level options)
  fit_args <- extra_args[!names(extra_args) %in% c("seed", "mc.cores", "maxit")]
  criterion <- match.arg(criterion)
  rank.grid <- rank; X.L2.smooth.grid <- X.L2.smooth
  nQ <- length(rank.grid); nS <- length(X.L2.smooth.grid)

  if (criterion %in% c("aic", "bic")) {
    one <- function(qi, si) {
      f <- tryCatch(do.call(nmf.cox,
                    c(list(formula, data = data, A = A, rank = rank.grid[qi],
                           X.L2.smooth = X.L2.smooth.grid[si], verbose = FALSE,
                           aic = TRUE, maxit = maxit, seed = 1), fit_args)),
                    error = function(e) NULL)
      if (is.null(f) || is.null(f$aic) || !is.finite(f$aic)) return(c(NA_real_, NA_real_, NA_real_))
      c(f$aic, f$loglik, f$df.eff)
    }
    tasks <- expand.grid(qi = seq_len(nQ), si = seq_len(nS))
    m <- if (mc.cores > 1)
           do.call(rbind, parallel::mclapply(seq_len(nrow(tasks)),
                     function(r) one(tasks$qi[r], tasks$si[r]), mc.cores = mc.cores))
         else t(vapply(seq_len(nrow(tasks)), function(r) one(tasks$qi[r], tasks$si[r]), numeric(3)))
    aicM <- matrix(m[, 1], nQ, nS); llM <- matrix(m[, 2], nQ, nS); dfM <- matrix(m[, 3], nQ, nS)
    nev  <- sum(as.integer(stats::model.response(stats::model.frame(formula, data))[, 2]))
    bicM <- -2 * llM + log(nev) * dfM
    critM <- if (criterion == "bic") bicM else aicM
    dimnames(aicM) <- dimnames(bicM) <- dimnames(llM) <- dimnames(dfM) <- list(rank = rank.grid, smooth = X.L2.smooth.grid)
    bi <- which(critM == min(critM, na.rm = TRUE), arr.ind = TRUE)[1, ]
    best.rank <- rank.grid[bi[1]]; best.smooth <- X.L2.smooth.grid[bi[2]]
    if (verbose) {
      cat(sprintf("%s grid (rows=rank, cols=X.L2.smooth):\n", toupper(criterion))); print(round(critM, 1))
      cat(sprintf("best(min %s): rank=%g, X.L2.smooth=%g\n", toupper(criterion), best.rank, best.smooth))
    }
    return(structure(list(criterion = criterion, rank = rank.grid, X.L2.smooth = X.L2.smooth.grid,
                 aic = aicM, bic = bicM, loglik = llM, df.eff = dfM, nevents = nev,
                 rank.best = best.rank, X.L2.smooth.best = best.smooth), class = "nmf.cox.cv"))
  }

  ## ---- criterion = "cvpl": K-fold cross-validated partial likelihood ----
  mf <- stats::model.frame(formula, data); y <- stats::model.response(mf)
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); N <- length(time)
  Z <- stats::model.matrix(attr(mf, "terms"), mf); Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  if (inherits(A, "formula")) {
    Am <- stats::model.matrix(A, stats::model.frame(A, data)); Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]
    Araw <- t(Am)
  } else { Araw <- as.matrix(A); if (ncol(Araw) != N && nrow(Araw) == N) Araw <- t(Araw) }
  set.seed(seed)
  fold <- sample(rep_len(1:nfolds, N))

  one <- function(qi, si, k) {
    Q <- rank.grid[qi]; sm <- X.L2.smooth.grid[si]; tr <- which(fold != k)
    fit <- tryCatch(do.call(nmf.cox,
                  c(list(formula, data = data[tr, , drop = FALSE],
                         A = Araw[, tr, drop = FALSE], rank = Q, X.L2.smooth = sm,
                         verbose = FALSE, maxit = maxit, seed = 1), fit_args)),
                  error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    Astd <- (Araw - fit$A.center) / fit$A.scale
    etk <- fit$event.times; B <- fit$beta.t
    betaFun <- function(tj) vapply(seq_len(ncol(B)),
                  function(l) stats::approx(etk, B[, l], xout = tj, rule = 2)$y, numeric(1))
    al <- stats::setNames(rep(0, ncol(Z)), colnames(Z)); al[names(fit$gamma)] <- fit$gamma
    zA <- as.vector(Z %*% al); clip <- function(x) pmax(pmin(x, 30), -30)
    plik <- function(sub) {
      et <- sort(unique(time[sub][status[sub] == 1])); ll <- 0
      for (tj in et) {
        risk <- sub[time[sub] >= tj]; fail <- sub[status[sub] == 1 & time[sub] == tj]
        if (!length(fail) || !length(risk)) next
        bt <- betaFun(tj)
        eta <- clip(zA[risk] + as.vector(crossprod(Astd[, risk, drop = FALSE], bt)))
        m <- max(eta); lse <- m + log(sum(exp(eta - m)))
        ef <- clip(zA[fail] + as.vector(crossprod(Astd[, fail, drop = FALSE], bt)))
        ll <- ll + sum(ef - lse)
      }
      ll
    }
    plik(1:N) - plik(tr)
  }
  tasks <- expand.grid(qi = seq_len(nQ), si = seq_len(nS), k = 1:nfolds)
  vals <- if (mc.cores > 1)
            unlist(parallel::mclapply(seq_len(nrow(tasks)),
                     function(r) one(tasks$qi[r], tasks$si[r], tasks$k[r]), mc.cores = mc.cores))
          else vapply(seq_len(nrow(tasks)),
                     function(r) one(tasks$qi[r], tasks$si[r], tasks$k[r]), numeric(1))
  arr  <- array(vals, dim = c(nQ, nS, nfolds))
  cvpl <- apply(arr, c(1, 2), function(v) sum(v, na.rm = TRUE))
  se   <- apply(arr, c(1, 2), function(v) { v <- v[is.finite(v)]
             if (length(v) < 2) NA_real_ else stats::sd(v) * sqrt(length(v)) })
  dimnames(cvpl) <- dimnames(se) <- list(rank = rank.grid, smooth = X.L2.smooth.grid)
  bi <- which(cvpl == max(cvpl, na.rm = TRUE), arr.ind = TRUE)[1, ]
  best.rank <- rank.grid[bi[1]]; best.smooth <- X.L2.smooth.grid[bi[2]]
  thr <- cvpl[bi[1], bi[2]] - ifelse(is.na(se[bi[1], bi[2]]), 0, se[bi[1], bi[2]])
  elig <- which(cvpl >= thr, arr.ind = TRUE)
  ord  <- order(rank.grid[elig[, 1]], -X.L2.smooth.grid[elig[, 2]])
  best.rank.1se <- rank.grid[elig[ord[1], 1]]; best.smooth.1se <- X.L2.smooth.grid[elig[ord[1], 2]]
  if (verbose) {
    cat("CVPL grid (rows=rank, cols=X.L2.smooth):\n"); print(round(cvpl, 2))
    cat(sprintf("best(CVPL argmax): rank=%g, X.L2.smooth=%g\n", best.rank, best.smooth))
    cat(sprintf("best(1-SE, most parsimonious): rank=%g, X.L2.smooth=%g\n", best.rank.1se, best.smooth.1se))
  }
  structure(list(criterion = "cvpl", rank = rank.grid, X.L2.smooth = X.L2.smooth.grid,
                 cvpl = cvpl, se = se, rank.best = best.rank, X.L2.smooth.best = best.smooth,
                 rank.best.1se = best.rank.1se, X.L2.smooth.best.1se = best.smooth.1se,
                 nfolds = nfolds), class = "nmf.cox.cv")
}


# --------------------------------------------------------------------
# nmf.cox.cf : two-stage cross-fit estimator
# --------------------------------------------------------------------

#' @title Two-stage cross-fit estimator of the proportional-hazards effect
#' @description
#' \code{nmf.cox.cf} estimates the proportional-hazards estimand \eqn{\gamma} by a
#' two-stage cross-fit (double machine-learning) procedure: for each fold \eqn{k}
#' the time-varying nuisance \eqn{\beta^{(-k)}(t)=X^{(-k)}\Theta^{(-k)}} is learned
#' \emph{without} fold \eqn{k}, the out-of-fold offset
#' \eqn{\hat o^{(-k)}_n(t)=\bm{a}_n^\top\beta^{(-k)}(t)} is frozen, and \eqn{\gamma}
#' is estimated on the full sample by a fold-stratified Cox fit with that offset and
#' cluster-robust standard errors.
#'
#' @details
#' Because the offset is linear in \code{a} at each event time, the score for
#' \eqn{\gamma} is Neyman-orthogonal to the nuisance, which is what makes the
#' cross-fit construction remove the first-order effect of nuisance estimation
#' error. The target \eqn{\gamma} is estimated on the whole sample; only the
#' nuisance is cross-fitted.
#'
#' The accompanying oracle-equivalence result is a \emph{high-level conditional}
#' statement: it presumes the out-of-fold nuisance attains a uniform
#' \eqn{o_p(N^{-1/4})} rate, together with within-risk-set orthogonality and
#' fold-wise regularity. Whether the constrained (non-negative, low-rank, smoothed)
#' estimator attains that rate under a suitable scaling of \code{X.L2.smooth} is not
#' established here.
#'
#' @param formula,data,A As in \code{\link{nmf.cox}}.
#' @param rank Rank \eqn{Q} of the nuisance time basis. Default \code{2}.
#' @param X.L2.smooth Roughness penalty used for the nuisance fits. Default \code{0}.
#' @param nfolds Number of cross-fitting folds. Default \code{5}.
#' @param verbose Logical; report per-fold nuisance convergence.
#' @param ... Additional arguments passed to \code{\link{nmf.cox}}
#'   (e.g. \code{X.update}, \code{X.restriction}, \code{X.init}, \code{nonneg},
#'   \code{nstart}). Also accepts: \code{seed} (fold assignment, default 123)
#'   and \code{maxit} (outer iterations per nuisance fit, default 30).
#'
#' @return An object of class \code{"nmf.cox.cf"} with components \code{gamma}
#'   (cross-fit estimate), \code{se} (cluster-robust), \code{cox.fit},
#'   \code{fold} (fold assignment), \code{nuisance.converged} and \code{nfolds}.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C, Newey W, Robins J
#' (2018). Double/debiased machine learning for treatment and structural parameters.
#' \emph{The Econometrics Journal}, 21(1), C1-C68.
#'
#' @seealso \code{\link{nmf.cox}}, \code{\link{nmf.cox.phtest}}
#'
#' @examples
#' \donttest{
#' library(survival)
#' d <- survival::veteran
#' d$y <- survival::Surv(d$time, d$status)
#' cf <- nmf.cox.cf(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                  rank = 2, X.L2.smooth = 1000, nfolds = 5, seed = 1)
#' cbind(estimate = cf$gamma, se = cf$se)
#' }
#'
#' @export
nmf.cox.cf <- function(formula, data, A, rank = 2, X.L2.smooth = 0, nfolds = 5,
                        verbose = FALSE, ...) {
  extra_args <- base::list(...)
  seed  <- if (!is.null(extra_args$seed))  extra_args$seed  else 123
  maxit <- if (!is.null(extra_args$maxit)) extra_args$maxit else 30
  ties  <- if (!is.null(extra_args$ties))  extra_args$ties  else "breslow"
  fit_args <- extra_args[!names(extra_args) %in% c("seed", "maxit")]
  if (!identical(ties, "breslow")) stop("nmf.cox.cf: only ties='breslow' is supported.")
  mf <- stats::model.frame(formula, data); y <- stats::model.response(mf)
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); N <- length(time)
  Z <- stats::model.matrix(attr(mf, "terms"), mf); Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  znames <- colnames(Z)
  if (inherits(A, "formula")) {
    Am <- stats::model.matrix(A, stats::model.frame(A, data)); Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]; Araw <- t(Am)
  } else { Araw <- as.matrix(A); if (ncol(Araw) != N && nrow(Araw) == N) Araw <- t(Araw) }

  et <- sort(unique(time[status == 1])); P <- length(et); s_prev <- c(0, et[-P])
  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in seq_len(P)) { risk <- which(time >= et[j]); m <- length(risk)
    rid <- c(rid, risk); rj <- c(rj, rep(j, m)); rstart <- c(rstart, rep(s_prev[j], m))
    rstop <- c(rstop, rep(et[j], m)); rev <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j])) }
  rstart <- as.numeric(rstart); rstop <- as.numeric(rstop)

  set.seed(seed)
  fold <- sample(rep_len(1:nfolds, N))
  Ocross <- matrix(0, N, P); conv <- logical(nfolds)
  for (k in 1:nfolds) {
    tr <- which(fold != k); ho <- which(fold == k)
    fit <- tryCatch(do.call(nmf.cox,
                  c(list(formula, data = data[tr, , drop = FALSE],
                         A = Araw[, tr, drop = FALSE], rank = rank,
                         X.L2.smooth = X.L2.smooth, verbose = FALSE,
                         maxit = maxit, seed = 1), fit_args)),
                  error = function(e) NULL)
    if (is.null(fit))
      stop(sprintf("nmf.cox.cf: nuisance fit (trained without fold %d) failed; refusing zero-offset fold.", k))
    conv[k] <- isTRUE(fit$converged)
    Xf  <- sapply(seq_len(ncol(fit$X)),
                  function(q) stats::approx(fit$event.times, fit$X[, q], xout = et, rule = 2)$y)
    Aho <- (Araw[, ho, drop = FALSE] - fit$A.center) / fit$A.scale
    Sho <- fit$C %*% Aho
    Ocross[ho, ] <- t(Xf %*% Sho)
    if (verbose) cat(sprintf("  fold %d: nuisance converged=%s\n", k, conv[k]))
  }
  o_row <- pmax(pmin(Ocross[cbind(rid, rj)], 30), -30)
  dat <- data.frame(start = rstart, stop = rstop, ev = rev, id = rid, off = o_row, fold = fold[rid])
  Zdf <- as.data.frame(Z[rid, , drop = FALSE]); colnames(Zdf) <- znames
  dat <- cbind(dat, Zdf)
  fml <- stats::as.formula(paste0("Surv(start,stop,ev) ~ ", paste(znames, collapse = " + "),
                           " + offset(off) + strata(fold) + cluster(id)"))
  cox.fit <- survival::coxph(fml, data = dat, ties = ties)
  sm <- summary(cox.fit)$coefficients
  se.col <- if ("robust se" %in% colnames(sm)) "robust se" else "se(coef)"
  structure(list(gamma = coef(cox.fit)[znames], se = stats::setNames(sm[znames, se.col], znames),
                 cox.fit = cox.fit, offset = Ocross, fold = fold, nfolds = nfolds,
                 nuisance.converged = conv, rank = rank, X.L2.smooth = X.L2.smooth),
            class = "nmf.cox.cf")
}


# --------------------------------------------------------------------
# nmf.cox.phtest : cross-fitted Wald test
# --------------------------------------------------------------------

#' @title Cross-fitted Wald test for a time-varying effect
#' @description
#' \code{nmf.cox.phtest} tests \eqn{H_0:\beta_r(t)=0} for each non-proportional
#' covariate, using a basis learned out of fold so that the test is not
#' invalidated by having selected the basis on the same data.
#'
#' @details
#' For each fold the non-negative time basis is learned without that fold; the
#' resulting basis is then used to span \eqn{\beta_r(t)} in a Cox fit on the full
#' sample, and a cluster-robust Wald statistic is formed per covariate together
#' with a global test. Using an in-sample (non cross-fitted) basis would over-reject.
#'
#' @param formula,data,A As in \code{\link{nmf.cox}}.
#' @param rank Rank \eqn{Q} of the time basis. Default \code{2}.
#' @param X.L2.smooth Roughness penalty used for the basis fits. Default \code{0}.
#' @param nfolds Number of cross-fitting folds. Default \code{5}.
#' @param verbose Logical; print the Wald table.
#' @param ... Additional arguments passed to \code{\link{nmf.cox}}
#'   (e.g. \code{X.update}, \code{X.restriction}, \code{X.init}, \code{nonneg},
#'   \code{nstart}). Also accepts: \code{seed} (fold assignment, default 123)
#'   and \code{maxit} (outer iterations per basis fit, default 30).
#'
#' @return An object of class \code{"nmf.cox.phtest"} containing the per-covariate
#'   Wald table (\code{wald}), the global test (\code{wald.global}) and the fold
#'   assignment.
#'
#' @seealso \code{\link{nmf.cox}}, \code{\link{nmf.cox.cf}}
#'
#' @examples
#' \donttest{
#' library(survival)
#' d <- survival::veteran
#' d$y <- survival::Surv(d$time, d$status)
#' ph <- nmf.cox.phtest(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                      rank = 2, X.L2.smooth = 1000, nfolds = 5, seed = 1)
#' ph$wald
#' }
#'
#' @export
nmf.cox.phtest <- function(formula, data, A, rank = 2, X.L2.smooth = 0,
                            nfolds = 5, verbose = FALSE, ...) {
  extra_args <- base::list(...)
  seed  <- if (!is.null(extra_args$seed))  extra_args$seed  else 123
  maxit <- if (!is.null(extra_args$maxit)) extra_args$maxit else 30
  ties  <- if (!is.null(extra_args$ties))  extra_args$ties  else "breslow"
  fit_args <- extra_args[!names(extra_args) %in% c("seed", "maxit")]
  if (!identical(ties, "breslow")) stop("nmf.cox.phtest: only ties='breslow' is supported.")
  Q <- rank
  mf <- stats::model.frame(formula, data); y <- stats::model.response(mf)
  time <- as.numeric(y[, 1]); status <- as.integer(y[, 2]); N <- length(time)
  Z <- stats::model.matrix(attr(mf, "terms"), mf); Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  znames <- colnames(Z)
  if (inherits(A, "formula")) {
    Am <- stats::model.matrix(A, stats::model.frame(A, data)); Am <- Am[, colnames(Am) != "(Intercept)", drop = FALSE]
    Araw <- t(Am)
  } else { Araw <- as.matrix(A); if (ncol(Araw) != N && nrow(Araw) == N) Araw <- t(Araw) }
  R <- nrow(Araw); anames <- rownames(Araw); if (is.null(anames)) anames <- paste0("A", 1:R)
  Ac <- rowMeans(Araw); As <- apply(Araw, 1, sd); As[As < 1e-10] <- 1
  Astd <- (Araw - Ac) / As
  et <- sort(unique(time[status == 1])); P <- length(et); s_prev <- c(0, et[-P])
  set.seed(seed)
  fold <- sample(rep_len(1:nfolds, N))

  Xcf <- vector("list", nfolds); conv <- logical(nfolds)
  for (k in 1:nfolds) {
    tr <- which(fold != k)
    fit <- tryCatch(do.call(nmf.cox,
                    c(list(formula, data = data[tr, , drop = FALSE],
                           A = Araw[, tr, drop = FALSE], rank = Q,
                           X.L2.smooth = X.L2.smooth, verbose = FALSE,
                           maxit = maxit, seed = 1), fit_args)),
                    error = function(e) NULL)
    if (is.null(fit)) stop(sprintf("nmf.cox.phtest: basis fit (trained without fold %d) failed.", k))
    conv[k] <- isTRUE(fit$converged)
    Xcf[[k]] <- sapply(1:ncol(fit$X), function(q) stats::approx(fit$event.times, fit$X[, q], xout = et, rule = 2)$y)
  }

  rid <- rj <- rstart <- rstop <- rev <- integer(0)
  for (j in 1:P) { risk <- which(time >= et[j]); m <- length(risk)
    rid <- c(rid, risk); rj <- c(rj, rep(j, m)); rstart <- c(rstart, rep(s_prev[j], m))
    rstop <- c(rstop, rep(et[j], m)); rev <- c(rev, as.integer(status[risk] == 1 & time[risk] == et[j])) }
  Zexp <- Z[rid, , drop = FALSE]

  frow <- fold[rid]; basisByRow <- matrix(0, length(rj), Q)
  for (k in 1:nfolds) { rk <- which(frow == k); if (length(rk)) basisByRow[rk, ] <- Xcf[[k]][rj[rk], , drop = FALSE] }
  Finf <- matrix(0, length(rj), Q * R); ci <- 1
  for (q in 1:Q) for (l in 1:R) { Finf[, ci] <- basisByRow[, q] * Astd[l, rid]; ci <- ci + 1 }

  inf.fit <- tryCatch(survival::coxph(survival::Surv(rstart, rstop, rev) ~ Zexp + Finf + cluster(rid), ties = ties,
                            control = survival::coxph.control(iter.max = 200, eps = 1e-9)), error = function(e) NULL)
  if (is.null(inf.fit)) return(NULL)
  V <- stats::vcov(inf.fit); cf <- coef(inf.fit); nz <- ncol(Zexp); idxTh <- (nz + 1):(nz + Q * R)
  quad.wald <- function(theta, Vmat) {
    Vs <- (Vmat + t(Vmat)) / 2; ee <- eigen(Vs, symmetric = TRUE)
    tol <- max(ee$values, 0) * 1e-8 * length(ee$values); keep <- ee$values > tol
    if (!any(keep)) return(c(chisq = NA_real_, df = 0))
    U <- ee$vectors[, keep, drop = FALSE]
    c(chisq = as.numeric(sum((crossprod(U, theta)^2) / ee$values[keep])), df = sum(keep))
  }
  wald <- data.frame(covariate = anames, df = NA_integer_, chisq = NA_real_, p.value = NA_real_, stringsAsFactors = FALSE)
  beta.vcov <- vector("list", R)
  for (l in 1:R) {
    posl <- idxTh[(0:(Q - 1)) * R + l]
    Vl <- V[posl, posl, drop = FALSE]; beta.vcov[[l]] <- Vl
    qw <- quad.wald(cf[posl], Vl)
    wald$chisq[l] <- qw["chisq"]; wald$df[l] <- qw["df"]
    wald$p.value[l] <- if (qw["df"] > 0) stats::pchisq(qw["chisq"], qw["df"], lower.tail = FALSE) else NA_real_
  }
  qwg <- quad.wald(cf[idxTh], V[idxTh, idxTh, drop = FALSE])
  wald.global <- list(chisq = as.numeric(qwg["chisq"]), df = as.integer(qwg["df"]),
                      p.value = if (qwg["df"] > 0) stats::pchisq(qwg["chisq"], qwg["df"], lower.tail = FALSE) else NA_real_)
  Xbar <- Reduce(`+`, Xcf) / nfolds
  C.cf <- matrix(cf[idxTh], nrow = Q, ncol = R, byrow = TRUE)
  beta.t <- Xbar %*% C.cf; se.beta.t <- matrix(NA_real_, P, R)
  for (l in 1:R) for (j in 1:P)
    se.beta.t[j, l] <- sqrt(max(0, as.numeric(crossprod(Xbar[j, ], beta.vcov[[l]] %*% Xbar[j, ]))))
  colnames(beta.t) <- colnames(se.beta.t) <- anames
  if (verbose) { cat("Cross-fitted Wald test  H0: beta_r(t)=0\n"); print(wald, row.names = FALSE) }
  structure(list(wald = wald, wald.global = wald.global, nfolds = nfolds, nuisance.converged = conv,
                 rank = Q, X.L2.smooth = X.L2.smooth, X = Xbar, C = C.cf,
                 beta.t = beta.t, se.beta.t = se.beta.t, beta.vcov = beta.vcov, event.times = et),
            class = "nmf.cox.phtest")
}

# --------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------

#' @title Extract the proportional-hazards coefficients from an NMF-COX fit
#' @param object An object of class \code{"nmf.cox"} or \code{"nmf.cox.cf"}.
#' @param ... Ignored.
#' @return A named numeric vector of the estimated \eqn{\gamma}.
#' @examples
#' library(survival)
#' d <- survival::veteran; d$y <- survival::Surv(d$time, d$status)
#' fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                rank = 2, X.L2.smooth = 1000, verbose = FALSE)
#' coef(fit)
#' @export
coef.nmf.cox <- function(object, ...) object$gamma

#' @rdname coef.nmf.cox
#' @export
coef.nmf.cox.cf <- function(object, ...) object$gamma

#' @title Time-varying coefficients implied by an NMF-COX fit
#' @description Returns \eqn{\beta(t)=XC}, the time-varying coefficient of the
#'   non-proportional covariates, evaluated at the event times
#'   (\code{object$event.times}); rows are event times, columns are covariates.
#' @param object An object of class \code{"nmf.cox"}.
#' @param ... Ignored.
#' @return A numeric matrix of per-SD log hazard ratios.
#' @examples
#' library(survival)
#' d <- survival::veteran; d$y <- survival::Surv(d$time, d$status)
#' fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                rank = 2, X.L2.smooth = 1000, verbose = FALSE)
#' head(fitted(fit))
#' @export
fitted.nmf.cox <- function(object, ...) object$beta.t

#' @title Print an NMF-COX fit
#' @param x An object of class \code{"nmf.cox"}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @examples
#' library(survival)
#' d <- survival::veteran; d$y <- survival::Surv(d$time, d$status)
#' fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                rank = 2, X.L2.smooth = 1000, verbose = FALSE)
#' print(fit)
#' @export
print.nmf.cox <- function(x, ...) {
  cat("NMF-COX: Cox model with a low-rank non-negative time-varying offset\n")
  cat(x$dims, "\n")
  cat(sprintf("rank=%d, X.L2.smooth=%g, X.restriction=%s, solver=%s\n",
              x$rank, x$X.L2.smooth, x$X.restriction, x$method))
  cat(sprintf("iterations=%d, converged=%s\n",
              x$iter, ifelse(isTRUE(x$converged), "TRUE", "FALSE")))
  if (!isTRUE(x$converged) && !is.null(x$stop.reason))
    cat("  ** ", x$stop.reason, " **\n", sep = "")
  cat("\nProportional-hazards coefficients (gamma), net of the time-varying effect of a:\n")
  if (!is.null(x$cox.fit)) print(round(coef(summary(x$cox.fit)), 4)) else print(round(x$gamma, 4))
  if (!is.null(x$wald)) {
    cat("\nWald test for a time-varying effect  H0: beta_r(t)=0  (given basis X):\n")
    w <- x$wald
    w$chisq <- round(w$chisq, 3)
    w$p.value <- format.pval(w$p.value, digits = 3, eps = 1e-4)
    print(w, row.names = FALSE)
    if (!is.null(x$wald.global))
      cat(sprintf("Global (all covariates):  chisq=%.3f, df=%d, p=%s\n",
                  x$wald.global$chisq, x$wald.global$df,
                  format.pval(x$wald.global$p.value, digits = 3, eps = 1e-4)))
  } else {
    cat("\n(run nmf.cox.inference() for SEs of beta(t), the C coefficients table and Wald tests)\n")
  }
  invisible(x)
}

#' @title Summary method for nmf.cox objects
#' @description
#' Summarises an NMF-COX fit: dimensions, convergence, the proportional-hazards
#' coefficients \eqn{\gamma}, and (after \code{\link{nmf.cox.inference}}) the
#' coefficients table of \eqn{C} and the given-basis Wald tests.
#' @param object An object of class \code{"nmf.cox"} (or \code{"nmf.cox.inference"}).
#' @param ... Not used.
#' @return An object of class \code{"summary.nmf.cox"}.
#' @examples
#' library(survival)
#' d <- survival::veteran; d$y <- survival::Surv(d$time, d$status)
#' fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                rank = 2, X.L2.smooth = 1000, verbose = FALSE)
#' summary(fit)
#' @export
summary.nmf.cox <- function(object, ...) {
  ans <- list(call = object$call, dims = object$dims, rank = object$rank,
              X.L2.smooth = object$X.L2.smooth, X.restriction = object$X.restriction,
              method = object$method, iter = object$iter,
              converged = object$converged, stop.reason = object$stop.reason,
              runtime = object$runtime, objfunc = object$objfunc,
              cox.fit = object$cox.fit, gamma = object$gamma,
              coefficients = object$coefficients,
              wald = object$wald, wald.global = object$wald.global)
  class(ans) <- "summary.nmf.cox"
  ans
}

#' @title Print method for summary.nmf.cox objects
#' @param x An object of class \code{"summary.nmf.cox"}.
#' @param digits Minimum number of significant digits.
#' @param ... Not used.
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @export
print.summary.nmf.cox <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Dimensions:", x$dims, "\n")
  cat(sprintf("Settings:   rank=%d, X.L2.smooth=%g, X.restriction=%s, solver=%s\n",
              x$rank, x$X.L2.smooth, x$X.restriction, x$method))
  cat(sprintf("Convergence: %s (iter=%d), objfunc=%.4f, runtime=%.1f sec\n",
              ifelse(isTRUE(x$converged), "converged", "NOT converged"),
              x$iter, x$objfunc, x$runtime))
  cat("\nProportional-hazards coefficients (gamma):\n")
  if (!is.null(x$cox.fit)) print(round(coef(summary(x$cox.fit)), 4)) else print(round(x$gamma, 4))
  if (!is.null(x$coefficients)) {
    cf <- x$coefficients
    stars <- ifelse(!is.finite(cf$p_value), " ",
             ifelse(cf$p_value < 0.001, "***",
             ifelse(cf$p_value < 0.01, "**",
             ifelse(cf$p_value < 0.05, "*",
             ifelse(cf$p_value < 0.1, ".", " ")))))
    cat("\nCoefficients C (conditional on the estimated basis X):\n")
    out <- data.frame(cf[, c("Basis", "Covariate")],
                      Estimate = round(cf$Estimate, 4), SE = round(cf$SE, 4),
                      z_value = round(cf$z_value, 2),
                      p_value = format.pval(cf$p_value, digits = 3, eps = 1e-4),
                      sig = stars)
    print(out, row.names = FALSE)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  if (!is.null(x$wald)) {
    cat("\nWald test  H0: beta_r(t)=0  (given basis X):\n")
    w <- x$wald
    w$chisq <- round(w$chisq, 3)
    w$p.value <- format.pval(w$p.value, digits = 3, eps = 1e-4)
    print(w, row.names = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' @title Plot method for nmf.cox objects (convergence)
#' @description Plots the penalised partial log-likelihood objective over
#'   iterations (house-style convergence diagnostic).
#' @param x An object of class \code{"nmf.cox"}.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @return Invisible \code{NULL}.
#' @examples
#' library(survival)
#' d <- survival::veteran; d$y <- survival::Surv(d$time, d$status)
#' fit <- nmf.cox(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                rank = 2, X.L2.smooth = 1000, verbose = FALSE)
#' plot(fit)
#' @export
plot.nmf.cox <- function(x, ...) {
  extra_args <- list(...)
  args <- list(x = seq_along(x$objfunc.iter), y = x$objfunc.iter)
  if (is.null(extra_args$type)) args$type <- "l"
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  if (is.null(extra_args$main)) args$main <- "nmf.cox convergence"
  do.call("plot", c(args, extra_args))
  invisible(NULL)
}

#' @title Print a cross-fit NMF-COX estimate
#' @param x An object of class \code{"nmf.cox.cf"}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @examples
#' \donttest{
#' library(survival)
#' d <- survival::veteran; d$y <- survival::Surv(d$time, d$status)
#' cf <- nmf.cox.cf(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                  rank = 2, X.L2.smooth = 1000, nfolds = 5, seed = 1)
#' print(cf)
#' }
#' @export
print.nmf.cox.cf <- function(x, ...) {
  cat("NMF-COX two-stage cross-fit estimate of gamma\n")
  cat(sprintf("rank=%d, X.L2.smooth=%g, folds=%d, nuisance converged in %d/%d folds\n",
              x$rank, x$X.L2.smooth, x$nfolds,
              sum(isTRUE(x$nuisance.converged) | x$nuisance.converged), x$nfolds))
  est <- x$gamma; se <- x$se
  tab <- cbind(estimate = est, se = se,
               lower = est - 1.96 * se, upper = est + 1.96 * se)
  cat("\nCluster-robust 95% confidence intervals:\n")
  print(round(tab, 4))
  invisible(x)
}

#' @title Print a cross-fitted time-varying-effect test
#' @param x An object of class \code{"nmf.cox.phtest"}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @examples
#' \donttest{
#' library(survival)
#' d <- survival::veteran; d$y <- survival::Surv(d$time, d$status)
#' ph <- nmf.cox.phtest(y ~ trt + age, data = d, A = ~ karno + celltype,
#'                      rank = 2, X.L2.smooth = 1000, nfolds = 5, seed = 1)
#' print(ph)
#' }
#' @export
print.nmf.cox.phtest <- function(x, ...) {
  cat("NMF-COX cross-fitted Wald test  H0: beta_r(t)=0\n")
  cat(sprintf("rank=%d, X.L2.smooth=%g, folds=%d\n", x$rank, x$X.L2.smooth, x$nfolds))
  w <- x$wald
  w$chisq <- round(w$chisq, 3)
  w$p.value <- format.pval(w$p.value, digits = 3, eps = 1e-4)
  print(w, row.names = FALSE)
  if (!is.null(x$wald.global))
    cat(sprintf("Global (all covariates):  chisq=%.3f, df=%d, p=%s\n",
                x$wald.global$chisq, x$wald.global$df,
                format.pval(x$wald.global$p.value, digits = 3, eps = 1e-4)))
  invisible(x)
}

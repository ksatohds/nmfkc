# =====================================================================
#  R/nmf.gmm.R -- NMF-GMM: a Gaussian-mixture latent-class extension of
#                 NMF with covariates (Satoh 2026).
#
#  Model:  b_n | (z_n = k) ~ N_Q(C a_n + mu_k, Sigma_k),  y_n = X b_n + eps_n,
#          X >= 0, colSums(X) = 1, eps ~ N(0, sigma2 I_P), sum_k xi_k mu_k = 0.
#  C (Q x R) is the covariate-coefficient matrix (Theta in the paper); mu
#  (Q x K) are the class means; the mixture on the scores turns NMF-RE's
#  arg-max hard clustering into model-based soft clustering via the posterior
#  responsibilities. K = 1 (tied) reduces exactly to NMF-RE (nmfre()).
#
#  House-style API (optimization / inference split):
#    nmf.gmm           EM fit (optimization only)
#    nmf.gmm.inference given-X Wald inference on C (sandwich / outer-product
#                      information + wild bootstrap), tied covariance only
#    nmf.gmm.select    choose K by BIC / ICL (optionally ARI vs known labels)
#
#  The internal .nmfgmm.* engine below is the verbatim research EM core
#  (estep / mstep / inner ANCOVA solve / multistart driver); the public
#  wrappers only rename arguments/returns to the package conventions
#  (Theta -> C, x0 -> X.init, ...) and add BIC/ICL, labels and S3 support.
#  Reference: Satoh K. (2026), NMF-GMM.
# =====================================================================

.nmfgmm.WOODBURY_P <- 100L   # switch to the Woodbury E-step at P >= this

.nmfgmm.pos <- function(M) pmax(M, 0)
.nmfgmm.neg <- function(M) pmax(-M, 0)

## k-means++ seeding (rows of Xr are points)
.nmfgmm.kmpp <- function(Xr, K, seed) {
  set.seed(seed); n <- nrow(Xr); ce <- matrix(NA, K, ncol(Xr))
  ce[1, ] <- Xr[sample(n, 1), ]
  d2 <- colSums((t(Xr) - ce[1, ])^2)
  for (k in 2:K) {
    i <- sample(n, 1, prob = d2 / sum(d2)); ce[k, ] <- Xr[i, ]
    d2 <- pmin(d2, colSums((t(Xr) - ce[k, ])^2))
  }
  ce
}

## ---------------------------------------------------------------------
## E-step
## ---------------------------------------------------------------------
.nmfgmm.estep <- function(par, Y, A, cov = "tied") {
  X <- par$X; Theta <- par$Theta; mu <- par$mu; tau2 <- par$tau2
  sigma2 <- par$sigma2; xi <- par$xi
  P <- nrow(Y); N <- ncol(Y); Q <- ncol(X); K <- length(xi)
  XtX <- crossprod(X); M <- Theta %*% A; XtY <- crossprod(X, Y)

  if (cov != "free") {                             # "tied" or "scalar": shared diagonal
    tau2 <- as.numeric(tau2)                       # length-Q vector (all equal if scalar)
    Dinv <- diag(1 / tau2, Q)
    Omega <- solve(Dinv + XtX / sigma2)            # Q x Q, both paths
    use_wood <- (P >= .nmfgmm.WOODBURY_P)
    if (use_wood) {
      W <- Omega / sigma2
      logdetV <- P * log(sigma2) +
        as.numeric(determinant(diag(Q) + XtX %*% diag(tau2, Q) / sigma2,
                               logarithm = TRUE)$modulus)
      YtY <- colSums(Y^2)
    } else {
      V <- X %*% (tau2 * t(X)) + sigma2 * diag(P)
      L <- chol(V); logdetV <- 2 * sum(log(diag(L))); Vinv <- chol2inv(L)
    }
    logcomp <- matrix(0, N, K); bhat <- vector("list", K)
    for (k in 1:K) {
      prior <- M + matrix(mu[, k], Q, N)
      if (use_wood) {
        Xprior <- X %*% prior; Xtr <- XtY - XtX %*% prior
        r_sq <- YtY - 2 * colSums(Y * Xprior) + colSums(Xprior^2)
        quad <- (r_sq - colSums(Xtr * (W %*% Xtr))) / sigma2
      } else {
        Resid <- Y - X %*% prior; quad <- colSums(Resid * (Vinv %*% Resid))
      }
      logcomp[, k] <- log(xi[k]) - 0.5 * (P * log(2 * pi) + logdetV + quad)
      bhat[[k]] <- Omega %*% (Dinv %*% prior + XtY / sigma2)
    }
    Omega_out <- Omega                             # single, class-free
  } else {                                         # cov == "free": per-class
    logcomp <- matrix(0, N, K); bhat <- vector("list", K)
    Omega_out <- vector("list", K)
    for (k in 1:K) {
      Dinv <- diag(1 / tau2[, k], Q)
      Ok <- solve(Dinv + XtX / sigma2)
      Vk <- X %*% (tau2[, k] * t(X)) + sigma2 * diag(P)
      Lk <- chol(Vk); logdetV <- 2 * sum(log(diag(Lk))); Vinv <- chol2inv(Lk)
      prior <- M + matrix(mu[, k], Q, N)
      Resid <- Y - X %*% prior; quad <- colSums(Resid * (Vinv %*% Resid))
      logcomp[, k] <- log(xi[k]) - 0.5 * (P * log(2 * pi) + logdetV + quad)
      bhat[[k]] <- Ok %*% (Dinv %*% prior + XtY / sigma2)
      Omega_out[[k]] <- Ok
    }
  }

  ## Vectorized softmax over rows of logcomp (bit-identical to the per-row loop:
  ## m/s broadcast down rows, gamma keeps the exp(lc - s) form, and Reduce("+")
  ## reproduces the loop's sequential double-precision loglik accumulation).
  m <- apply(logcomp, 1, max)
  s <- m + log(rowSums(exp(logcomp - m)))
  gamma <- exp(logcomp - s)
  ll <- Reduce("+", s)
  list(gamma = gamma, bhat = bhat, Omega = Omega_out, loglik = ll,
       Nj = colSums(gamma))
}

## ---------------------------------------------------------------------
## Inner Theta/mu alternation: the ANCOVA normal equations, solved to
## convergence (rel. change < inner_tol) with a cap (inner_max).  Returns
## the pair (Theta, mu) satisfying  E A' = 0  and  sum_n gamma_nk e_nk = 0
## to machine precision (paper rem:ancova).  Used by both .nmfgmm.mstep and the
## final polish in .nmfgmm.fit.
## ---------------------------------------------------------------------
.nmfgmm.solve_theta_mu <- function(Theta, mu, bhat, gamma, Nj, A, intercept,
                           inner_tol = 1e-12, inner_max = 50) {
  Q <- nrow(mu); N <- ncol(A); K <- ncol(mu)
  xi <- Nj / sum(Nj); AAt_inv <- solve(tcrossprod(A))
  prev <- c(Theta, mu)
  for (it in 1:inner_max) {
    bbar <- matrix(0, Q, N)
    for (k in 1:K) bbar <- bbar + (bhat[[k]] - mu[, k]) * rep(gamma[, k], each = Q)
    Theta <- (bbar %*% t(A)) %*% AAt_inv; M <- Theta %*% A
    for (k in 1:K) mu[, k] <- rowSums((bhat[[k]] - M) * rep(gamma[, k], each = Q)) / Nj[k]
    mubar <- as.numeric(mu %*% xi); mu <- mu - mubar
    Theta[, intercept] <- Theta[, intercept] + mubar
    cur <- c(Theta, mu)
    if (sqrt(sum((cur - prev)^2)) / (sqrt(sum(cur^2)) + 1e-12) < inner_tol) break
    prev <- cur
  }
  list(Theta = Theta, mu = mu)
}

## ---------------------------------------------------------------------
## M-step
## ---------------------------------------------------------------------
.nmfgmm.mstep <- function(par, es, Y, A, intercept = 1, cov = "tied",
                  inner_tol = 1e-12, inner_max = 50) {
  X <- par$X; tau2 <- par$tau2; sigma2 <- par$sigma2
  P <- nrow(Y); N <- ncol(Y); Q <- ncol(X); K <- ncol(par$mu)
  gamma <- es$gamma; bhat <- es$bhat; Omega <- es$Omega; Nj <- es$Nj
  xi <- Nj / N

  ## --- inner Theta/mu alternation (ANCOVA normal equations) ---
  tm <- .nmfgmm.solve_theta_mu(par$Theta, par$mu, bhat, gamma, Nj, A, intercept,
                       inner_tol, inner_max)
  Theta <- tm$Theta; mu <- tm$mu

  ## --- variances ---
  M <- Theta %*% A
  if (cov != "free") {                               # "tied" or "scalar": one shared Omega
    num <- rep(0, Q)
    for (k in 1:K) {
      resid <- (bhat[[k]] - M) - mu[, k]
      num <- num + rowSums(resid^2 * rep(gamma[, k], each = Q))
    }
    tau2 <- pmax((num + N * diag(Omega)) / N, 1e-6)
  } else {
    for (k in 1:K) {
      resid <- (bhat[[k]] - M) - mu[, k]
      tau2[, k] <- pmax((rowSums(resid^2 * rep(gamma[, k], each = Q)) +
                         Nj[k] * diag(Omega[[k]])) / Nj[k], 1e-6)
    }
  }

  ## --- basis: responsibility-weighted semi-NMF update ---
  Cstat <- matrix(0, P, Q); Gstat <- matrix(0, Q, Q)
  for (k in 1:K) {
    Cstat <- Cstat + (Y * rep(gamma[, k], each = P)) %*% t(bhat[[k]])
    Gstat <- Gstat + (bhat[[k]] * rep(gamma[, k], each = Q)) %*% t(bhat[[k]])
    if (cov == "free") Gstat <- Gstat + Nj[k] * Omega[[k]]
  }
  if (cov != "free") Gstat <- Gstat + N * Omega
  Xn <- X * sqrt((.nmfgmm.pos(Cstat) + X %*% .nmfgmm.neg(Gstat)) /
                 (.nmfgmm.neg(Cstat) + X %*% .nmfgmm.pos(Gstat) + 1e-12))

  ## --- noise variance ---
  ss <- 0; tr <- 0
  for (k in 1:K) {
    rY <- Y - Xn %*% bhat[[k]]; ss <- ss + sum(colSums(rY^2) * gamma[, k])
    if (cov == "free") tr <- tr + Nj[k] * sum((Xn %*% Omega[[k]]) * Xn)
  }
  if (cov != "free") tr <- N * sum((Xn %*% Omega) * Xn)
  sigma2 <- max((ss + tr) / (P * N), 1e-6)

  ## --- column-normalize X and rescale ---
  D <- colSums(Xn); Xn <- Xn / rep(D, each = P)
  Theta <- Theta * D; mu <- mu * D
  tau2 <- tau2 * D^2   # length-Q D^2 recycles down cols of the Q x K free matrix; elementwise for tied/scalar
  if (cov == "scalar") tau2 <- rep(mean(tau2), Q)     # isotropic: pool to one variance

  list(X = Xn, Theta = Theta, mu = mu, tau2 = tau2, sigma2 = sigma2, xi = xi)
}

## ---------------------------------------------------------------------
## Initialization  (NMF-RE fit for K=1, then k-means++ on residual scores)
## ---------------------------------------------------------------------
.nmfgmm.init_par <- function(Y, A, X0, K, sigma2_0, intercept = 1, seed = 1, cov = "tied") {
  Q <- ncol(X0); N <- ncol(Y)
  bls <- solve(crossprod(X0), crossprod(X0, Y))
  Theta <- (bls %*% t(A)) %*% solve(tcrossprod(A))
  r <- bls - Theta %*% A
  ms <- function(M) pmax(apply(M, 1, function(v) mean(v^2)), 1e-3)
  if (K == 1) {
    tau2 <- if (cov == "free") matrix(ms(r), Q, 1)
            else if (cov == "scalar") rep(mean(ms(r)), Q) else ms(r)
    return(list(X = X0, Theta = Theta, mu = matrix(0, Q, 1),
                tau2 = tau2, sigma2 = sigma2_0, xi = 1))
  }
  rt <- t(r)
  km <- stats::kmeans(rt, centers = .nmfgmm.kmpp(rt, K, seed), iter.max = 50)
  xi <- as.numeric(table(factor(km$cluster, levels = 1:K))) / N
  mu <- matrix(t(km$centers), Q, K)
  if (cov != "free") {
    resid <- r - mu[, km$cluster]; tau2 <- pmax(rowMeans(resid^2), 1e-3)
    if (cov == "scalar") tau2 <- rep(mean(tau2), Q)
  } else {
    tau2 <- sapply(1:K, function(k) { idx <- km$cluster == k
      if (sum(idx) > 1) ms(r[, idx, drop = FALSE] - mu[, k]) else ms(r) })
    tau2 <- matrix(tau2, Q, K)
  }
  mubar <- as.numeric(mu %*% xi); mu <- mu - mubar
  Theta[, intercept] <- Theta[, intercept] + mubar
  list(X = X0, Theta = Theta, mu = mu, tau2 = tau2, sigma2 = sigma2_0, xi = xi)
}

## ---------------------------------------------------------------------
## Multistart EM driver
## ---------------------------------------------------------------------
.nmfgmm.fit <- function(Y, A, X0, K, cov = "tied", intercept = 1,
                       maxit = 500, tol = 1e-7,
                       nstart = if (K == 1) 1 else 8,
                       inner_tol = 1e-12, inner_max = 50, cores = 1L) {
  bls <- solve(crossprod(X0), crossprod(X0, Y))
  s2_0 <- max(mean((Y - X0 %*% bls)^2), 1e-3)
  ## One EM restart. Each start re-seeds from `s` (via .nmfgmm.init_par) and the
  ## EM loop is RNG-free, so run_start(s) is a deterministic function of s alone
  ## -- independent of the other starts and of execution order.
  run_start <- function(s) {
    par <- .nmfgmm.init_par(Y, A, X0, K, s2_0, intercept, seed = s, cov = cov)
    ll_old <- -Inf; hist <- numeric(0)
    for (it in 1:maxit) {
      es <- .nmfgmm.estep(par, Y, A, cov = cov); hist <- c(hist, es$loglik)
      if (abs(es$loglik - ll_old) / (abs(es$loglik) + 1) < tol) break
      ll_old <- es$loglik
      par <- .nmfgmm.mstep(par, es, Y, A, intercept, cov = cov,
                   inner_tol = inner_tol, inner_max = inner_max)
    }
    es <- .nmfgmm.estep(par, Y, A, cov = cov)
    ## polish: re-solve Theta/mu at the final responsibilities so the returned
    ## pair satisfies the ANCOVA normal equations (rem:ancova) to machine
    ## precision, removing the one-step outer-loop lag.  Holds X, variances,
    ## and the loglik fixed to outer-tol; .nmfgmm.ARI/BIC unaffected.
    tm <- .nmfgmm.solve_theta_mu(par$Theta, par$mu, es$bhat, es$gamma, es$Nj, A,
                         intercept, inner_tol, inner_max)
    par$Theta <- tm$Theta; par$mu <- tm$mu
    list(par = par, es = es, loglik = es$loglik, iter = it, hist = hist)
  }
  ## .nmfkc.parlapply preserves input order, so which.max() picks the same
  ## first-best restart as the sequential `> best$loglik` scan -- the returned
  ## fit is identical for any `cores`.
  starts <- .nmfkc.parlapply(seq_len(nstart), run_start,
                             cores = if (nstart > 1L) cores else 1L)
  starts[[which.max(vapply(starts, function(r) r$loglik, numeric(1)))]]
}

## ---------------------------------------------------------------------
## Utilities
## ---------------------------------------------------------------------
## free (per-class) parameter count uses Q*K variances; tied uses Q.
.nmfgmm.pcount <- function(Q, P, R, K, cov = "tied") {
  vars <- switch(cov, scalar = 1L, tied = Q, free = Q * K)
  Q * (P - 1) + Q * R + Q * (K - 1) + vars + (K - 1) + 1
}

.nmfgmm.ARI <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); cc <- function(x) sum(choose(x, 2))
  ai <- rowSums(tab); bj <- colSums(tab)
  idx <- cc(as.vector(tab)); ei <- cc(ai) * cc(bj) / choose(n, 2)
  mi <- (cc(ai) + cc(bj)) / 2
  if (mi - ei == 0) return(0); (idx - ei) / (mi - ei)
}

.nmfgmm.hardclass <- function(f) max.col(f$es$gamma, ties.method = "first")

# =====================================================================
#  Public API
# =====================================================================

#' @title Fit NMF-GMM: a Gaussian-mixture latent-class extension of NMF with covariates
#' @description
#' \code{nmf.gmm} fits the model
#' \deqn{\bm b_n\mid(z_n=k)\sim N_Q(C\bm a_n+\bm\mu_k,\Sigma_k),\qquad
#'       \bm y_n = X\bm b_n+\bm\varepsilon_n,\quad X\ge 0,}
#' by a generalized EM algorithm: a \eqn{K}-component Gaussian mixture is placed
#' on the latent NMF scores \eqn{\bm b_n}, whose class-conditional mean combines
#' a covariate regression \eqn{C\bm a_n} (\eqn{C=\Theta} in the paper, the
#' \eqn{Q\times R} coefficient matrix) with a class shift \eqn{\bm\mu_k}. The
#' shared basis \eqn{X} stays non-negative and column-normalized. Clustering is
#' model based, through the posterior responsibilities. \eqn{K=1} (tied
#' covariance) reduces exactly to NMF with random effects (\code{\link{nmfre}}).
#'
#' Like the rest of the package, \code{nmf.gmm} performs \strong{optimization
#' only}; use \code{\link{nmf.gmm.inference}} for standard errors and tests of
#' \eqn{C}, and \code{\link{nmf.gmm.select}} to choose \eqn{K}.
#'
#' @param Y Data matrix \eqn{Y} (P x N).
#' @param A Covariate matrix \eqn{A} (R x N). Default \code{NULL} uses an
#'   intercept-only \eqn{1\times N} matrix (a plain Gaussian mixture on the
#'   scores, with no covariate adjustment). A row of ones should be present for
#'   the identifiability constraint \eqn{\sum_k \xi_k\bm\mu_k = 0} (see
#'   \code{intercept}).
#' @param rank Integer rank \eqn{Q} of the basis.
#' @param K Integer number of mixture components. Default 1 (= NMF-RE).
#' @param ... Additional arguments:
#'   \itemize{
#'     \item \code{cov}: score-covariance structure. \code{"tied"} (default,
#'       shared diagonal \eqn{\Sigma_k=\mathrm{diag}(\tau^2_1,\dots,\tau^2_Q)},
#'       \eqn{Q} variances); \code{"free"} (per-class diagonal, \eqn{Q\times K});
#'       or \code{"scalar"} (isotropic \eqn{\Sigma_k=\tau^2 I_Q}, a single
#'       variance --- the most parsimonious variant; at \eqn{K=1} it is the
#'       \code{\link{nmfre}} model, and \code{tau2} is returned as one number).
#'     \item \code{X.init}: initial basis. A \eqn{P\times Q} non-negative matrix,
#'       or an initialization method name passed to the shared initializer
#'       (\code{"nndsvd"} default, \code{"kmeans++"}, \code{"kmeans"}, ...).
#'     \item \code{intercept}: index of the intercept row of \code{A} (default 1);
#'       the class-mean average is absorbed into that column of \eqn{C}.
#'     \item \code{nstart}: EM restarts (default 1 for \code{K=1}, else 8).
#'     \item \code{cores}: run the \code{nstart} restarts in parallel (default
#'       \code{getOption("mc.cores", 1L)}; PSOCK cluster on Windows, forking
#'       elsewhere). Each restart is an independent self-seeded EM and results
#'       are combined in order, so the returned fit is identical for any
#'       \code{cores}. (\code{\link{nmf.gmm.select}} parallelizes over \code{K}
#'       instead and runs each inner fit with \code{cores = 1}.)
#'     \item \code{maxit}, \code{tol}: outer EM cap / tolerance (500, 1e-7).
#'     \item \code{seed}: RNG seed (default 1). \code{prefix}: basis-name prefix
#'       (default \code{"Basis"}).
#'   }
#'
#' @return An object of class \code{"nmf.gmm"}: a list with \code{X}
#'   (basis), \code{C} (\eqn{=\Theta}, Q x R), \code{mu} (class means, Q x K),
#'   \code{tau2}, \code{sigma2}, \code{xi} (mixing proportions), \code{gamma}
#'   (responsibilities, N x K), \code{cluster} (hard labels), \code{loglik},
#'   \code{BIC}, \code{ICL}, \code{n.params}, \code{entropy}, \code{Yhat},
#'   \code{objfunc(.iter)}, \code{iter}, \code{converged}, \code{K},
#'   \code{rank}, \code{cov}, \code{dims} and \code{runtime}.
#'
#' @seealso \code{\link{nmf.gmm.inference}}, \code{\link{nmf.gmm.select}},
#'   \code{\link{nmfre}} (the \eqn{K=1} special case).
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(abs(rnorm(20 * 60)) + 1, 20, 60)
#' A <- rbind(1, rnorm(60))                 # intercept + one covariate
#' fit <- nmf.gmm(Y, A, rank = 2, K = 2)
#' fit$BIC
#' table(fit$cluster)
#' }
#' @export
nmf.gmm <- function(Y, A = NULL, rank, K = 1, ...) {
  extra <- base::list(...)
  getopt <- function(nm, default) if (!is.null(extra[[nm]])) extra[[nm]] else default
  cov       <- getopt("cov", "tied")
  X.init    <- extra$X.init
  intercept <- getopt("intercept", 1L)
  maxit     <- getopt("maxit", 500L)
  tol       <- getopt("tol", 1e-7)
  seed      <- getopt("seed", 1L)
  prefix    <- getopt("prefix", "Basis")
  inner_tol <- getopt("inner_tol", 1e-12)
  inner_max <- getopt("inner_max", 50L)
  cov <- match.arg(cov, c("tied", "free", "scalar"))

  t0 <- proc.time()
  Y <- as.matrix(Y); P <- nrow(Y); N <- ncol(Y); Q <- as.integer(rank)
  if (is.null(A)) A <- matrix(1, 1, N)
  A <- as.matrix(A); if (ncol(A) != N && nrow(A) == N) A <- t(A)
  if (ncol(A) != N) stop("A must have N columns (R x N).")
  R <- nrow(A)
  nstart <- getopt("nstart", if (K == 1) 1L else 8L)
  cores  <- getopt("cores", getOption("mc.cores", 1L))

  ## --- initial non-negative, column-normalized basis X0 ---
  if (is.matrix(X.init) || (is.numeric(X.init) && length(X.init) > 1)) {
    X0 <- as.matrix(X.init)
  } else {
    method <- if (is.character(X.init)) X.init else "nndsvd"
    X0 <- .init_X_method(method, Y, Q, seed = seed, nstart = max(nstart, 1L))
  }
  X0 <- X0 / rep(pmax(colSums(X0), 1e-12), each = nrow(X0))

  fit <- .nmfgmm.fit(Y, A, X0, K, cov = cov, intercept = intercept,
                     maxit = maxit, tol = tol, nstart = nstart,
                     inner_tol = inner_tol, inner_max = inner_max, cores = cores)
  par <- fit$par; es <- fit$es

  ## --- labels (house style) ---
  X <- par$X; C <- par$Theta; mu <- par$mu
  blab <- paste0(prefix, seq_len(Q)); clab <- paste0("Class", seq_len(K))
  alab <- rownames(A); if (is.null(alab)) alab <- paste0("Cov", seq_len(R))
  dimnames(X) <- list(rownames(Y), blab)
  dimnames(C) <- list(blab, alab)
  dimnames(mu) <- list(blab, clab)
  gamma <- es$gamma; colnames(gamma) <- clab
  cluster <- max.col(gamma, ties.method = "first")
  tau2 <- par$tau2
  if (cov == "scalar") tau2 <- as.numeric(tau2)[1]   # single isotropic variance
  else if (cov == "tied") { tau2 <- as.numeric(tau2); names(tau2) <- blab }
  else dimnames(tau2) <- list(blab, clab)

  ## --- BIC / ICL ---
  n.params <- .nmfgmm.pcount(Q, P, R, K, cov = cov)
  bic <- -2 * fit$loglik + n.params * log(N)
  ent <- -sum(ifelse(gamma > 1e-300, gamma * log(gamma), 0))
  icl <- bic + 2 * ent

  ## --- least-squares scores (Q x N): bls = (X'X)^{-1} X'Y ---
  ## The geometric projection of Y onto the estimated basis X, used by
  ## plot(type = "scores"). This matches the paper's score-space figure: the
  ## covariate mean C %*% A is removed at plot time to show the size-adjusted
  ## residual scores. (Deliberately NOT the responsibility-averaged posterior
  ## mean sum_k gamma_nk E[b_n | z=k], which shrinks points toward their class
  ## centroid and would make the clusters look artificially separated.)
  scores <- solve(crossprod(X), crossprod(X, Y))
  dimnames(scores) <- list(blab, colnames(Y))

  ## --- responsibility-averaged reconstruction Yhat = X (C A + mu gamma') ---
  Yhat <- X %*% (C %*% A + mu %*% t(gamma))
  dimnames(Yhat) <- dimnames(Y)

  dims <- sprintf("Y(%d,%d) ~ X(%d,%d) [K=%d, %s]", P, N, P, Q, K, cov)
  runtime <- as.numeric((proc.time() - t0)["elapsed"])

  structure(list(
    call = match.call(), dims = dims, runtime = runtime,
    rank = Q, K = K, cov = cov, intercept = intercept,
    X = X, C = C, mu = mu, tau2 = tau2, sigma2 = par$sigma2, xi = par$xi,
    gamma = gamma, cluster = cluster, scores = scores, Yhat = Yhat,
    loglik = fit$loglik, BIC = bic, ICL = icl, entropy = ent, n.params = n.params,
    objfunc = -fit$loglik, objfunc.iter = -fit$hist, iter = fit$iter,
    converged = fit$iter < maxit, A = A
  ), class = "nmf.gmm")
}


#' @title Statistical inference for an NMF-GMM fit (given basis)
#' @description
#' \code{nmf.gmm.inference} adds Wald inference on the covariate-coefficient
#' matrix \eqn{C} (\eqn{=\Theta}) of a fitted \code{\link{nmf.gmm}} object,
#' conditional on the estimated basis \eqn{\hat X} and mixture. It reports the
#' analytic outer-product (mixture-information) standard error --- the
#' coverage-calibrated choice at \eqn{K>1} --- together with a
#' \strong{wild-bootstrap} standard error and percentile confidence interval,
#' mirroring \code{\link{nmfre.inference}} / \code{\link{nmfkc.inference}}.
#' Supported for the tied-covariance variant.
#'
#' @param object A fitted \code{"nmf.gmm"} object (with \code{cov = "tied"}).
#' @param Y The data matrix used in the fit.
#' @param A The covariate matrix used in the fit. Default \code{object$A}.
#' @param ... Additional arguments: \code{boot} (bootstrap replicates, default
#'   2000), \code{seed} (default 123), \code{ci.level} (default 0.95).
#'
#' @return \code{object} with class \code{c("nmf.gmm.inference", "nmf.gmm")}
#'   and added components: \code{coefficients} (data frame with \code{Basis},
#'   \code{Covariate}, \code{Estimate}, \code{SE}, \code{BSE}, \code{z_value},
#'   \code{p_value}, \code{CI_low}, \code{CI_high}), \code{C.se.outer},
#'   \code{C.se.sandwich}, \code{C.se.boot} (Q x R matrices) and \code{Xi}
#'   (\eqn{= X C}).
#' @seealso \code{\link{nmf.gmm}}, \code{\link{nmf.gmm.select}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(abs(rnorm(20 * 60)) + 1, 20, 60); A <- rbind(1, rnorm(60))
#' fit <- nmf.gmm(Y, A, rank = 2, K = 2)
#' fit <- nmf.gmm.inference(fit, Y, A, boot = 300)
#' head(fit$coefficients)
#' }
#' @export
nmf.gmm.inference <- function(object, Y, A = object$A, ...) {
  if (!inherits(object, "nmf.gmm")) stop("object must be of class 'nmf.gmm'.")
  if (!identical(object$cov, "tied"))
    stop("nmf.gmm.inference currently supports cov = 'tied' only.")
  extra <- base::list(...)
  B     <- if (!is.null(extra$boot))     extra$boot     else 2000L
  seed  <- if (!is.null(extra$seed))     extra$seed     else 123L
  level <- if (!is.null(extra$ci.level)) extra$ci.level else 0.95

  X <- object$X; C <- object$C; mu <- object$mu
  tau2 <- as.numeric(object$tau2); s2 <- object$sigma2; gamma <- object$gamma
  Y <- as.matrix(Y); A <- as.matrix(A)
  P <- nrow(Y); N <- ncol(Y); Q <- ncol(X); R <- nrow(A)

  V     <- X %*% (tau2 * t(X)) + s2 * diag(P)
  Vinv  <- chol2inv(chol(V))
  XtVX  <- crossprod(X, Vinv %*% X)
  Info  <- kronecker(tcrossprod(A), XtVX)               # vec(C) column-stacked
  Hinv  <- solve(Info)

  mubar <- mu %*% t(gamma)                              # Q x N
  Rm    <- Y - X %*% (C %*% A + mubar)
  G     <- crossprod(X, Vinv %*% Rm)                    # Q x N
  score <- matrix(0, Q * R, N)
  for (n in 1:N) score[, n] <- as.vector(tcrossprod(G[, n], A[, n]))

  J      <- tcrossprod(score) * (N / (N - 1))
  V_sand <- Hinv %*% J %*% Hinv
  se_sand <- matrix(sqrt(pmax(diag(V_sand), 0)), Q, R)
  SS     <- tcrossprod(score)
  V_J    <- tryCatch(solve(SS),
                     error = function(e) solve(SS + diag(1e-8 * mean(diag(SS)), nrow(SS))))
  se_J   <- matrix(sqrt(pmax(diag(V_J), 0)), Q, R)

  ## wild bootstrap on vec(C)
  set.seed(seed)
  th <- as.vector(C); boot <- matrix(NA_real_, Q * R, B)
  for (b in 1:B) {
    w <- stats::rexp(N, rate = 1) - 1
    boot[, b] <- th + as.vector(Hinv %*% (score %*% w))
  }
  a <- 1 - level
  lo  <- matrix(apply(boot, 1, stats::quantile, probs = a/2,     names = FALSE), Q, R)
  hi  <- matrix(apply(boot, 1, stats::quantile, probs = 1 - a/2, names = FALSE), Q, R)
  bsd <- matrix(apply(boot, 1, stats::sd), Q, R)

  z <- C / se_J
  coefficients <- data.frame(
    Basis     = rep(rownames(C), times = R),
    Covariate = rep(colnames(C), each = Q),
    Estimate  = as.vector(C),
    SE        = as.vector(se_J),
    BSE       = as.vector(bsd),
    z_value   = as.vector(z),
    p_value   = as.vector(2 * stats::pnorm(-abs(z))),
    CI_low    = as.vector(lo),
    CI_high   = as.vector(hi),
    stringsAsFactors = FALSE, row.names = NULL)

  Xi <- X %*% C; dimnames(Xi) <- list(rownames(X), colnames(C))
  object$coefficients <- coefficients
  object$C.se.outer <- se_J; object$C.se.sandwich <- se_sand; object$C.se.boot <- bsd
  dimnames(object$C.se.outer) <- dimnames(object$C.se.sandwich) <-
    dimnames(object$C.se.boot) <- dimnames(C)
  object$Xi <- Xi; object$C.p.side <- "two.sided"
  if (!inherits(object, "nmf.gmm.inference"))
    class(object) <- c("nmf.gmm.inference", class(object))
  object
}


#' @title Choose the number of mixture components K for NMF-GMM
#' @description
#' \code{nmf.gmm.select} fits \code{\link{nmf.gmm}} over a vector of \eqn{K}
#' values and reports the log-likelihood, BIC and ICL for each, selecting
#' \code{K.best} by BIC (\code{K.best.icl} by ICL). If a vector of known labels
#' is supplied via \code{truth}, the adjusted Rand index of the hard clustering
#' is added.
#'
#' @param Y,A,rank As in \code{\link{nmf.gmm}}.
#' @param K Integer vector of candidate component counts. Default \code{1:5}.
#' @param ... Additional arguments passed to \code{\link{nmf.gmm}} (e.g.
#'   \code{cov}, \code{X.init}, \code{nstart}, \code{maxit}, \code{seed}). Also
#'   accepts \code{truth} (a length-N vector of known class labels for the ARI
#'   column), \code{verbose} (logical, print the table; default \code{TRUE}), and
#'   \code{cores} (evaluate the \code{K} candidates in parallel; default
#'   \code{getOption("mc.cores", 1L)}). Parallelism uses a PSOCK cluster on
#'   Windows and forking elsewhere; because each \code{K} is an independent
#'   self-seeded fit and results are returned in order, the table and the
#'   selected \code{K} are identical for any \code{cores}.
#'
#' @return An object of class \code{"nmf.gmm.select"}: a list with the
#'   \code{table} (data frame over \code{K}), \code{K.best} (min BIC),
#'   \code{K.best.icl} (min ICL) and \code{fits} (the fitted objects).
#' @seealso \code{\link{nmf.gmm}}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(abs(rnorm(20 * 80)) + 1, 20, 80); A <- rbind(1, rnorm(80))
#' sel <- nmf.gmm.select(Y, A, rank = 2, K = 1:4, verbose = FALSE)
#' sel$K.best
#' }
#' @export
nmf.gmm.select <- function(Y, A = NULL, rank, K = 1:5, ...) {
  extra <- base::list(...)
  truth   <- extra$truth
  verbose <- if (!is.null(extra$verbose)) extra$verbose else TRUE
  cores   <- if (!is.null(extra$cores)) extra$cores else getOption("mc.cores", 1L)
  fit_args <- extra[!names(extra) %in% c("truth", "verbose", "cores")]

  ## Each K is an independent, self-seeded nmf.gmm() fit; .nmfkc.parlapply
  ## preserves input order so the table and selection are identical to the
  ## sequential loop for any `cores`. The inner fit runs its nstart restarts
  ## sequentially (cores = 1) so the K-level parallelism here does not nest.
  fit_one <- function(k)
    do.call(nmf.gmm, c(list(Y = Y, A = A, rank = rank, K = k, cores = 1L), fit_args))
  fits <- .nmfkc.parlapply(K, fit_one, cores = cores)
  rows <- lapply(seq_along(K), function(i) {
    f <- fits[[i]]
    ari <- if (!is.null(truth) && K[i] > 1) .nmfgmm.ARI(f$cluster, truth) else NA_real_
    data.frame(K = K[i], logLik = f$loglik, n.params = f$n.params,
               BIC = f$BIC, ICL = f$ICL, ARI = ari)
  })
  tab <- do.call(rbind, rows)
  K.best     <- tab$K[which.min(tab$BIC)]
  K.best.icl <- tab$K[which.min(tab$ICL)]
  if (isTRUE(verbose)) {
    show <- tab; show$logLik <- round(show$logLik, 1)
    show$BIC <- round(show$BIC, 1); show$ICL <- round(show$ICL, 1)
    show$ARI <- round(show$ARI, 3)
    if (all(is.na(show$ARI))) show$ARI <- NULL
    print(show, row.names = FALSE)
    cat(sprintf("BIC selects K=%d; ICL selects K=%d\n", K.best, K.best.icl))
  }
  structure(list(table = tab, K.best = K.best, K.best.icl = K.best.icl,
                 fits = fits), class = "nmf.gmm.select")
}


# =====================================================================
#  S3 methods
# =====================================================================

#' @title Extract the covariate-coefficient matrix from an NMF-GMM fit
#' @param object An object of class \code{"nmf.gmm"}.
#' @param ... Ignored.
#' @return The \eqn{Q\times R} coefficient matrix \eqn{C} (\eqn{=\Theta}).
#' @export
coef.nmf.gmm <- function(object, ...) object$C

#' @title Fitted (responsibility-averaged) reconstruction of an NMF-GMM fit
#' @param object An object of class \code{"nmf.gmm"}.
#' @param ... Ignored.
#' @return The \eqn{P\times N} fitted matrix \eqn{\hat Y = X(CA + \mu\gamma^\top)}.
#' @export
fitted.nmf.gmm <- function(object, ...) object$Yhat

#' @title Class assignments / responsibilities from an NMF-GMM fit
#' @param object An object of class \code{"nmf.gmm"}.
#' @param type \code{"class"} (default) returns the hard cluster labels;
#'   \code{"responsibility"} returns the N x K posterior-probability matrix.
#' @param ... Ignored.
#' @return A vector of labels or the responsibility matrix.
#' @export
predict.nmf.gmm <- function(object, type = c("class", "responsibility"), ...) {
  type <- match.arg(type)
  if (type == "class") object$cluster else object$gamma
}

#' @title Print an NMF-GMM fit
#' @param x An object of class \code{"nmf.gmm"}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.nmf.gmm <- function(x, ...) {
  cat("NMF-GMM: Gaussian-mixture latent-class NMF with covariates\n")
  cat(x$dims, "\n")
  cat(sprintf("K=%d components, rank Q=%d, covariance=%s\n", x$K, x$rank, x$cov))
  cat(sprintf("logLik=%.2f, BIC=%.2f, ICL=%.2f (params=%d), converged=%s\n",
              x$loglik, x$BIC, x$ICL, x$n.params, isTRUE(x$converged)))
  cat("Mixing proportions (xi):", paste(round(x$xi, 3), collapse = " "), "\n")
  cat("Cluster sizes:\n"); print(table(cluster = x$cluster))
  if (!is.null(x$coefficients))
    cat("\n(coefficients table for C available in $coefficients)\n")
  invisible(x)
}

#' @title Summary of an NMF-GMM fit
#' @param object An object of class \code{"nmf.gmm"}.
#' @param ... Ignored.
#' @return An object of class \code{"summary.nmf.gmm"}.
#' @export
summary.nmf.gmm <- function(object, ...) {
  ans <- list(dims = object$dims, K = object$K, rank = object$rank,
              cov = object$cov, loglik = object$loglik, BIC = object$BIC,
              ICL = object$ICL, n.params = object$n.params, xi = object$xi,
              sizes = table(object$cluster), converged = object$converged,
              runtime = object$runtime, coefficients = object$coefficients,
              C.p.side = object$C.p.side)
  class(ans) <- "summary.nmf.gmm"
  ans
}

#' @title Print method for summary.nmf.gmm objects
#' @param x A \code{"summary.nmf.gmm"} object.
#' @param digits Minimum significant digits.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.summary.nmf.gmm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("NMF-GMM fit\n"); cat(x$dims, "\n")
  cat(sprintf("K=%d, rank=%d, cov=%s | logLik=%.2f, BIC=%.2f, ICL=%.2f, params=%d\n",
              x$K, x$rank, x$cov, x$loglik, x$BIC, x$ICL, x$n.params))
  cat("Mixing proportions:", paste(round(x$xi, 3), collapse = " "), "\n")
  cat("Cluster sizes:\n"); print(x$sizes)
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
                      BSE = round(cf$BSE, 4), z = round(cf$z_value, 2),
                      p_value = format.pval(cf$p_value, digits = 3, eps = 1e-4),
                      sig = stars)
    print(out, row.names = FALSE)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    cat("(SE = outer-product mixture-information; BSE = wild bootstrap.)\n")
  }
  invisible(x)
}

#' @title Plot method for nmf.gmm objects
#' @description
#' \code{type = "convergence"} (default) plots the EM objective
#' (\eqn{-\log L}) over iterations. \code{type = "scores"} draws a cluster
#' scatterplot: the covariate-adjusted least-squares scores
#' (\eqn{(X'X)^{-1}X'Y - C A}) projected onto their first two principal
#' components, coloured by the hard cluster assignment (or by \code{group}).
#' @param x An object of class \code{"nmf.gmm"}.
#' @param type \code{"convergence"} or \code{"scores"}.
#' @param group Optional length-N factor/vector to colour the score plot by
#'   (e.g. known labels). Default \code{NULL} colours by \code{x$cluster}.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @return Invisible \code{NULL}.
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(abs(rnorm(20 * 80)) + 1, 20, 80); A <- rbind(1, rnorm(80))
#' fit <- nmf.gmm(Y, A, rank = 2, K = 3)
#' plot(fit, type = "scores")
#' }
#' @export
plot.nmf.gmm <- function(x, type = c("convergence", "scores"), group = NULL, ...) {
  type <- match.arg(type)
  extra <- list(...)
  if (type == "convergence") {
    args <- list(x = seq_along(x$objfunc.iter), y = x$objfunc.iter)
    if (is.null(extra$type)) args$type <- "l"
    if (is.null(extra$xlab)) args$xlab <- "EM iteration"
    if (is.null(extra$ylab)) args$ylab <- "-loglik"
    if (is.null(extra$main)) args$main <- "nmf.gmm convergence"
    do.call("plot", c(args, extra))
    return(invisible(NULL))
  }
  ## type == "scores": covariate-adjusted scores, PCA to 2D, coloured by cluster
  if (is.null(x$scores)) stop("no $scores stored; refit with nmf.gmm().")
  adj <- x$scores - x$C %*% x$A                      # remove the covariate mean
  g   <- if (!is.null(group)) as.factor(group) else as.factor(x$cluster)
  pal <- grDevices::hcl.colors(max(nlevels(g), 2), "Dark 3")
  col <- pal[as.integer(g)]
  if (nrow(adj) >= 2) {
    pc <- stats::prcomp(t(adj))$x[, 1:2, drop = FALSE]
    xl <- "adjusted score PC1"; yl <- "adjusted score PC2"
  } else {                                            # rank 1: strip plot
    pc <- cbind(as.numeric(adj), stats::runif(ncol(adj), -1, 1))
    xl <- "adjusted score"; yl <- ""
  }
  args <- list(x = pc[, 1], y = pc[, 2], col = col)
  if (is.null(extra$pch))  args$pch  <- 19
  if (is.null(extra$xlab)) args$xlab <- xl
  if (is.null(extra$ylab)) args$ylab <- yl
  if (is.null(extra$main)) args$main <- sprintf("nmf.gmm scores (K=%d)", x$K)
  do.call("plot", c(args, extra))
  graphics::legend("topright", legend = levels(g), col = pal[seq_len(nlevels(g))],
                   pch = 19, bty = "n", cex = 0.8)
  invisible(NULL)
}

#' @title Print method for nmf.gmm.select objects
#' @param x An object of class \code{"nmf.gmm.select"}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.nmf.gmm.select <- function(x, ...) {
  tab <- x$table; tab$logLik <- round(tab$logLik, 1)
  tab$BIC <- round(tab$BIC, 1); tab$ICL <- round(tab$ICL, 1); tab$ARI <- round(tab$ARI, 3)
  if (all(is.na(tab$ARI))) tab$ARI <- NULL
  print(tab, row.names = FALSE)
  cat(sprintf("BIC selects K=%d; ICL selects K=%d\n", x$K.best, x$K.best.icl))
  invisible(x)
}

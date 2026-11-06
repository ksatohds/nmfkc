## nmfkc.ident(): identifiability of a covariate-NMF fit Y ~ X C A
## (X >= 0 with columns summing to 1, C >= 0, A known).
##
## 1. test: the anchor test H(pi, rho) on the OLS estimate of Xi = X C (HC1, intersection-union over
##    candidates, Bonferroni over the 2Q conditions of a candidate); pi.max and rho.max.
## 2. identification region (R = Q): the maximum over the model, Xi.hat (row-wise non-negative least
##    squares in the metric AA'; the group means when A is one-hot), and the range of every entry of C and
##    X over all exact factorizations X = Xt N^-1, C = N S of Xi.hat.
## 3. repair: nmfkc re-run from a point of the region (the one nearest to object's C, or the center); the
##    start is already optimal, so the re-fit stops at once at the exact optimum.
## 4. intervals (ci = TRUE): profile-likelihood intervals of C (Gaussian errors with one variance common to
##    all entries of Y), a confidence interval of the identification region.
##
## Internal helpers carry the prefix .nid. (".ident." is taken by nmfkc.ar.R).

#' @title Identifiability of an nmfkc fit: test, identification region, repair and intervals
#' @description
#' \code{nmfkc.ident} examines whether the factorization \eqn{Y \approx X C A}
#' of a fitted \code{\link{nmfkc}} model is unique, and what can be said about
#' \eqn{C} (\eqn{\Theta}) when it is not.  For the model
#' \eqn{Y = X C A + \varepsilon} with \eqn{X \ge 0} (columns summing to 1),
#' \eqn{C \ge 0} and \eqn{A} known, any column-stochastic \eqn{T} with
#' \eqn{XT \ge 0} and \eqn{T^{-1}C \ge 0} gives the same fit, so the bases and
#' the coefficients are determined only up to such a remixing.  The function
#' \enumerate{
#'   \item tests the anchor condition that makes the factorization unique
#'     (every basis has an anchor row of \eqn{X} with membership probability at
#'     least \code{pi} and size at least \code{rho}, and an anchor column of
#'     \eqn{C}), from the least-squares estimate of \eqn{\Xi = XC} only, and
#'     reports the largest \code{pi} and \code{rho} not rejected;
#'   \item computes the identification region: the range of every entry of
#'     \eqn{C} and \eqn{X} over all exact factorizations of the maximum
#'     \eqn{\hat\Xi};
#'   \item re-optimizes the fit from a point of that region (\code{fit}), which
#'     repairs a fit that stopped short of the optimum;
#'   \item optionally (\code{ci = TRUE}) computes profile-likelihood intervals
#'     for \eqn{C}, which cover the identification region with the stated
#'     confidence.
#' }
#' Steps 2--4 need as many covariates as bases (\code{nrow(A) == rank}), e.g.
#' group indicators or \code{rank} continuous covariates; for
#' \code{nrow(A) > rank} (or \code{region = FALSE}) only the test is computed.
#'
#' @param object An object of class \code{"nmfkc"} fitted with covariates
#'   (\code{X.restriction = "colSums"}, the default).
#' @param Y Observation matrix (P x N) used in the fit.
#' @param A Covariate matrix (R x N) used in the fit, of full row rank with
#'   \code{R < N}.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{pi}}{Membership-probability threshold of the anchor rows at
#'       which the test decision is reported. Default \code{1} (exact anchors).}
#'     \item{\code{rho}}{Size threshold of the anchor rows (value within a
#'       column of \eqn{X}). Default \code{0.5 / P}, half the value of a uniform
#'       basis.}
#'     \item{\code{alpha}}{Level of the test. Default \code{0.05}.}
#'     \item{\code{start}}{Starting point of the repair: \code{"nearest"}
#'       (default; the point of the region closest to \code{object}'s \eqn{C})
#'       or \code{"center"} (the point closest to the midpoints of the ranges
#'       of \eqn{C}).}
#'     \item{\code{region}}{Logical. Compute the identification region and
#'       the repaired fit (steps 2 and 3). Default \code{TRUE}; \code{FALSE}
#'       returns the test only, in a second or so.}
#'     \item{\code{ci}}{Logical. Compute profile-likelihood intervals for
#'       \eqn{C}. Default \code{FALSE} (it takes minutes for moderate
#'       \code{P}; see \code{cores}).}
#'     \item{\code{ci.level}}{Confidence level of the intervals. Default 0.95.}
#'     \item{\code{ci.step}}{Grid step of the profile, as a fraction of the
#'       column sum of \eqn{\hat\Xi}. Default 0.02. The interval ends are refined
#'       by bisection between grid points.}
#'     \item{\code{region.n}}{Number of random factorizations drawn to explore
#'       the region. Default 50000. Drawing them is cheap (about 2 seconds for
#'       50000 at P = 67); fewer than a few thousand leave the region visibly
#'       too narrow.}
#'     \item{\code{region.nstart}}{Number of the most extreme random
#'       factorizations from which each end of the region is refined by
#'       Nelder--Mead. Default 3. This is where the time goes (about 20
#'       seconds at 3, a minute at 10, for P = 67 and rank 3).}
#'     \item{\code{epsilon}, \code{maxit}}{Tolerance and iteration cap of the
#'       repair re-fit. Defaults \code{1e-8} and \code{1e5}.}
#'     \item{\code{seed}}{Seed. Default 123. The caller's random stream is
#'       restored on exit.}
#'     \item{\code{cores}}{Number of cores for the profile intervals. Default
#'       \code{getOption("mc.cores", 1L)}; results do not depend on it.}
#'   }
#'
#' @details
#' \strong{Test.} \eqn{\hat\Xi = Y A^\top (AA^\top)^{-1}} with HC1 standard
#' errors. A candidate is an ordered set of \eqn{Q} columns \eqn{L} of
#' \eqn{\Xi} (the anchor columns of \eqn{C}) and \eqn{Q} distinct rows (the
#' anchor rows of \eqn{X}); in the representative where the columns \eqn{L}
#' are pure, \eqn{X = \Xi_L S_L^{-1}}, so the membership probability and the
#' size of a row are smooth functions of \eqn{\hat\Xi}. A candidate survives if
#' none of its \eqn{2Q} one-sided \eqn{t} statistics falls below
#' \eqn{-t_{1-\alpha/(2Q)}(N-R)}; the null hypothesis is rejected if no
#' candidate survives (intersection--union test). \code{pi.max} is the largest
#' \code{pi} not rejected at \code{rho}, an upper confidence bound of the best
#' attainable membership probability of the anchor rows; \code{rho.max} is the
#' largest \code{rho} not rejected at \code{pi = 1}. Non-rejection is not
#' evidence of identification, as in a goodness-of-fit test. For
#' \code{nrow(A) > rank} the condition that \eqn{C} is non-negative outside the
#' anchor columns is not tested.
#'
#' \strong{Identification region (R = Q).} With as many covariates as bases,
#' every non-negative \eqn{\Xi} factorizes exactly, so the maximum of the
#' likelihood over the model is the row-wise non-negative least squares
#' \eqn{\hat\Xi} in the metric \eqn{AA^\top} (for one-hot \eqn{A}, the group
#' means). All exact factorizations are \eqn{X = \tilde X N^{-1}},
#' \eqn{C = N S} with \eqn{\tilde X} the column-normalized \eqn{\hat\Xi},
#' \eqn{S} its column sums and \eqn{N \ge 0} column-stochastic with
#' \eqn{\tilde X N^{-1} \ge 0}; the bases are labelled by matching the columns
#' of \eqn{X} to those of \eqn{\tilde X} (so that relabellings are excluded).
#' The range of each entry of \eqn{C} is found from \code{region.n} random
#' factorizations and refined by Nelder--Mead from the \code{region.nstart} most extreme of them; the ranges of \eqn{X} come from
#' the random factorizations only. Both are inner approximations.
#' \code{identified} is \code{TRUE} when every range of \eqn{C} is shorter than
#' \code{1e-6} times the largest column sum.
#'
#' \strong{Intervals.} Gaussian errors with one variance common to all entries
#' of \eqn{Y} (the working model whose maximum is nmfkc's least squares),
#' profiled over \eqn{X}, the other entries of \eqn{C} and the variance. The
#' profile is flat over the identification region, so the interval contains
#' it; with the \eqn{\chi^2_1} cut each end of the region is missed with
#' probability about \eqn{(1-\code{ci.level})/2}, so the interval is a
#' confidence interval of the identification region (and a conservative one
#' for the true value of the coefficient). The asymptotics need the ends of
#' the region to move smoothly with \eqn{\Xi}; an exact zero in
#' \eqn{\hat\Xi} that pins a column breaks this. The profile is computed by
#' continuation from the ends of the region with local re-fits, so the
#' intervals may be slightly too short.
#'
#' @return An object of class \code{"nmfkc.ident"}:
#' \item{test}{List: \code{pi}, \code{rho}, \code{alpha}, \code{reject}
#'   (decision at \code{pi}, \code{rho}), \code{pi.max}, \code{rho.max},
#'   \code{pi.point} (point estimate of the attainable membership
#'   probability), \code{crit} (critical value), \code{candidate} (rows and
#'   columns of a surviving candidate, or \code{NULL}).}
#' \item{Xi.hat}{The maximum over the model (P x R); \code{NULL} if
#'   \code{R > rank}.}
#' \item{C.region.lower, C.region.upper}{Identification region of \eqn{C}
#'   (Q x R), bases in \code{object}'s order.}
#' \item{X.region.lower, X.region.upper}{Identification region of \eqn{X}
#'   (P x Q).}
#' \item{identified}{Logical: the region of \eqn{C} is a single point.}
#' \item{X, C}{The repaired basis and coefficients (bases in \code{object}'s
#'   order).}
#' \item{fit}{The repaired \code{"nmfkc"} object (bases in \code{object}'s
#'   order), usable with \code{\link{nmfkc.inference}} and the other methods.}
#' \item{repair}{List: \code{start}, \code{iter}, \code{change} (largest change
#'   of \eqn{C} during the re-fit), \code{gap} (largest
#'   \eqn{|XC - \hat\Xi|}, relative to the largest entry of \eqn{\hat\Xi}).}
#' \item{C.ci.lower, C.ci.upper}{Profile-likelihood interval of \eqn{C}
#'   (\code{NULL} unless \code{ci = TRUE}).}
#' \item{ci.level, rank, call}{As named.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.inference}},
#'   \code{\link{plot.nmfkc.ident}}
#' @references
#' Imbens, G. W. and Manski, C. F. (2004). Confidence intervals for partially
#' identified parameters. \emph{Econometrica}, 72(6), 1845--1857.
#' \doi{10.1111/j.1468-0262.2004.00555.x}
#' @export
#' @examples
#' ## two groups, two bases; variables 1 and 5 are anchors (multiplicative noise keeps their zeros)
#' set.seed(1)
#' X <- cbind(c(0.4, 0.3, 0.2, 0.1, 0), c(0, 0.1, 0.2, 0.3, 0.4))
#' A <- rbind(g1 = rep(c(1, 0), each = 20), g2 = rep(c(0, 1), each = 20))
#' Y <- (X %*% diag(c(10, 8)) %*% A) * exp(matrix(rnorm(200, sd = 0.1), 5, 40))
#' fit <- nmfkc(Y, A, rank = 2, epsilon = 1e-8, verbose = FALSE)
#' id <- nmfkc.ident(fit, Y, A, region.n = 5000)
#' id
#' \donttest{
#' id2 <- nmfkc.ident(fit, Y, A, region.n = 5000, ci = TRUE)  # profile intervals
#' plot(id2)
#' }
nmfkc.ident <- function(object, Y, A, ...) {
  if (!inherits(object, "nmfkc")) stop("object must be of class 'nmfkc'")
  ex <- base::list(...)
  opt <- function(nm, def) if (!base::is.null(ex[[nm]])) ex[[nm]] else def
  Y <- base::as.matrix(Y); A <- base::as.matrix(A)
  P <- base::nrow(Y); N <- base::ncol(Y); R <- base::nrow(A); Q <- object$rank
  pi0      <- opt("pi", 1)
  rho0     <- opt("rho", 0.5 / P)
  alpha    <- opt("alpha", 0.05)
  start    <- base::match.arg(opt("start", "nearest"), c("nearest", "center"))
  region   <- !base::isFALSE(opt("region", TRUE))
  ci       <- base::isTRUE(opt("ci", FALSE))
  ci.level <- opt("ci.level", 0.95)
  ci.step  <- opt("ci.step", 0.02)
  region.n <- opt("region.n", 50000)
  region.nstart <- opt("region.nstart", 3)
  epsilon  <- opt("epsilon", 1e-8)
  maxit    <- opt("maxit", 1e5)
  seed     <- opt("seed", 123)
  cores    <- opt("cores", base::getOption("mc.cores", 1L))
  .rng <- .nmfkc.rng.save(seed)
  base::on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
  base::set.seed(seed)

  if (base::ncol(A) != N) stop("ncol(A) must equal ncol(Y)")
  if (base::anyNA(Y) || base::anyNA(A)) stop("Y and A must not contain missing values")
  if (!base::identical(object$X.restriction, "colSums"))
    stop("nmfkc.ident() needs a fit with X.restriction = \"colSums\"")
  if (R >= N) stop("nrow(A) must be smaller than ncol(Y)")
  if (base::qr(A)$rank < R) stop("A must have full row rank")
  if (Q > R) stop("the rank must not exceed nrow(A)")
  if (!base::all(base::dim(object$C) == base::c(Q, R)))
    stop("object does not match A (dim(object$C) must be rank x nrow(A))")

  ## ---- 1. test ----
  f <- .nid.ols(Y, A)
  cval <- stats::qt(1 - alpha / (2 * Q), N - R)
  dec <- .nid.decide(f, Q, pi0, rho0, cval)
  test <- base::list(pi = pi0, rho = rho0, alpha = alpha, reject = dec$reject,
                     pi.max = .nid.max.pi(f, Q, rho0, cval),
                     rho.max = .nid.max.rho(f, Q, cval),
                     pi.point = .nid.pi.point(f, Q, rho0)$value, crit = cval,
                     candidate = if (dec$reject) NULL else base::list(rows = dec$J, columns = dec$L))

  out <- base::list(test = test, Xi.hat = NULL, C.region.lower = NULL, C.region.upper = NULL,
                    X.region.lower = NULL, X.region.upper = NULL, identified = NA,
                    X = object$X, C = object$C, fit = object, repair = NULL,
                    C.ci.lower = NULL, C.ci.upper = NULL, ci.level = ci.level, rank = Q,
                    call = base::match.call())
  if (R > Q || !region) {
    base::class(out) <- "nmfkc.ident"
    return(out)
  }

  ## ---- 2. identification region ----
  s <- .nid.suff(Y, A)
  Xi <- .nid.xi.hat(s)
  dimnames(Xi) <- base::list(base::rownames(Y), base::rownames(A))
  s0 <- base::colSums(Xi)
  if (base::any(s0 <= 0)) stop("a column of the maximum Xi.hat is zero; the region is not defined")
  Xt <- base::sweep(Xi, 2, s0, "/")
  reg <- .nid.region(Xt, s0, region.n, region.nstart)
  ord <- .nid.match(object$X, Xt)              # ord[q] = internal basis of object's basis q
  Cn <- base::dimnames(object$C); Xn <- base::dimnames(object$X)
  put_C <- function(M) { M <- M[ord, , drop = FALSE]; dimnames(M) <- Cn; M }
  put_X <- function(M) { M <- M[, ord, drop = FALSE]; dimnames(M) <- Xn; M }
  out$Xi.hat <- Xi
  out$C.region.lower <- put_C(reg$C.lo); out$C.region.upper <- put_C(reg$C.hi)
  out$X.region.lower <- put_X(reg$X.lo); out$X.region.upper <- put_X(reg$X.hi)
  out$identified <- base::max(reg$C.hi - reg$C.lo) <= 1e-6 * base::max(s0)

  ## ---- 3. repair ----
  target <- if (start == "nearest") {
    Ct <- base::matrix(0, Q, R); Ct[ord, ] <- object$C; Ct
  } else (reg$C.lo + reg$C.hi) / 2
  Nst <- .nid.nearest(target, reg, Xt, s0)
  X0 <- base::pmax(Xt %*% base::solve(Nst), 0); C0 <- base::pmax(Nst %*% base::diag(s0, Q), 0)
  X0 <- base::sweep(X0, 2, base::colSums(X0), "/")
  fit <- nmfkc(Y, A, rank = Q, epsilon = epsilon, maxit = maxit, verbose = FALSE,
               X.init = X0[, ord, drop = FALSE], C.init = C0[ord, , drop = FALSE])
  idx <- .nid.match(X0[, ord, drop = FALSE], fit$X)   # idx[q] = fit's basis for object's basis q
  fit <- .nid.reorder(fit, idx, Y, A, object)
  out$X <- fit$X; out$C <- fit$C; out$fit <- fit
  out$repair <- base::list(start = start, iter = fit$iter,
                           change = base::max(base::abs(fit$C - C0[ord, , drop = FALSE])),
                           gap = base::max(base::abs(fit$X %*% fit$C - Xi)) / base::max(Xi))

  ## ---- 4. intervals ----
  if (ci) {
    pr <- .nid.profile(s, Xt, s0, reg, Q, ci.level, ci.step, seed, cores)
    out$C.ci.lower <- put_C(pr$lo); out$C.ci.upper <- put_C(pr$hi)
    out$ci.dip <- pr$dip
  }
  base::class(out) <- "nmfkc.ident"
  out
}

## ===================== test =====================

.nid.perms <- function(n) {
  if (n == 1) return(base::matrix(1L, 1, 1))
  p <- .nid.perms(n - 1)
  base::do.call(base::rbind, base::lapply(base::seq_len(n), function(i) base::cbind(i, base::ifelse(p >= i, p + 1L, p))))
}

## OLS of Xi and the pieces of the HC1 delta method
.nid.ols <- function(Y, A) {
  N <- base::ncol(Y); R <- base::nrow(A)
  W  <- base::solve(A %*% base::t(A), A)
  Xi <- Y %*% base::t(W); E <- Y - Xi %*% A
  base::list(Xi = Xi, s = base::colSums(Xi), W = W, E = E, e1 = base::colSums(E), N = N, R = R,
             hc = N / (N - R))
}

## one-sided t of the membership contrast (1 - pi) u_jq - pi sum_{q' != q} u_jq', rows j, columns L
.nid.purity.t <- function(f, L, pi) {
  Q <- base::length(L); P <- base::nrow(f$Xi)
  sL <- f$s[L]; XL <- f$Xi[, L, drop = FALSE]; U <- base::sweep(XL, 2, sL, "/")
  WL <- f$W[L, , drop = FALSE]
  tt <- base::matrix(NA_real_, P, Q)
  for (q in base::seq_len(Q)) {
    wq <- base::ifelse(base::seq_len(Q) == q, 1 - pi, -pi)
    g <- base::as.vector(U %*% wq)
    a <- wq / sL
    B <- -base::sweep(XL, 2, wq / sL^2, "*")
    psi <- base::sweep(f$E, 2, base::as.vector(base::crossprod(a, WL)), "*") +
      base::sweep(B %*% WL, 2, f$e1, "*")
    se <- base::sqrt(base::rowSums(psi^2) * f$hc)
    tt[, q] <- base::ifelse(se > 0, g / se, base::ifelse(g >= 0, Inf, -Inf))
  }
  tt
}

## one-sided t of the size contrast Xi[j, l_q] - rho s_{l_q}
.nid.size.t <- function(f, L, rho) {
  Q <- base::length(L); P <- base::nrow(f$Xi); tt <- base::matrix(NA_real_, P, Q)
  for (q in base::seq_len(Q)) {
    l <- L[q]; g <- f$Xi[, l] - rho * f$s[l]
    psi <- base::sweep(f$E, 2, f$W[l, ], "*") - rho * base::outer(base::rep(1, P), f$e1 * f$W[l, ])
    se <- base::sqrt(base::rowSums(psi^2) * f$hc)
    tt[, q] <- base::ifelse(se > 0, g / se, base::ifelse(g >= 0, Inf, -Inf))
  }
  tt
}

## distinct rows j_q in S_q (q = 1..Q), or NULL
.nid.matching <- function(S) {
  if (base::any(base::lengths(S) == 0)) return(NULL)
  rec <- function(q, used) {
    if (q > base::length(S)) return(integer(0))
    for (j in S[[q]]) if (!(j %in% used)) {
      r <- rec(q + 1, base::c(used, j)); if (!base::is.null(r)) return(base::c(j, r))
    }
    NULL
  }
  rec(1, integer(0))
}

.nid.all.L <- function(R, Q) {
  cs <- utils::combn(R, Q); pm <- .nid.perms(Q)
  base::do.call(base::rbind, base::lapply(base::seq_len(base::ncol(cs)), function(k)
    base::t(base::apply(pm, 1, function(p) cs[p, k]))))
}

.nid.decide <- function(f, Q, pi, rho, cval) {
  Ls <- .nid.all.L(f$R, Q)
  for (k in base::seq_len(base::nrow(Ls))) {
    L <- Ls[k, ]; pt <- .nid.purity.t(f, L, pi); st <- .nid.size.t(f, L, rho)
    S <- base::lapply(base::seq_len(Q), function(q) base::which(pt[, q] >= -cval & st[, q] >= -cval))
    J <- .nid.matching(S)
    if (!base::is.null(J)) return(base::list(reject = FALSE, L = L, J = J))
  }
  base::list(reject = TRUE, L = NULL, J = NULL)
}

## largest pi not rejected at rho (bisection); -Inf if rejected at pi = 1/Q
.nid.max.pi <- function(f, Q, rho, cval) {
  lo <- 1 / Q
  if (.nid.decide(f, Q, lo, rho, cval)$reject) return(-Inf)
  if (!.nid.decide(f, Q, 1, rho, cval)$reject) return(1)
  hi <- 1
  for (it in 1:30) {
    mid <- (lo + hi) / 2
    if (.nid.decide(f, Q, mid, rho, cval)$reject) hi <- mid else lo <- mid
  }
  lo
}

## largest rho not rejected at pi = 1 (bisection on [0, 1]); -Inf if rejected at rho = 0
.nid.max.rho <- function(f, Q, cval) {
  if (.nid.decide(f, Q, 1, 0, cval)$reject) return(-Inf)
  lo <- 0; hi <- 1
  if (!.nid.decide(f, Q, 1, hi, cval)$reject) return(1)
  for (it in 1:30) {
    mid <- (lo + hi) / 2
    if (.nid.decide(f, Q, 1, mid, cval)$reject) hi <- mid else lo <- mid
  }
  lo
}

## point estimate of the attainable membership probability (negative OLS entries counted as 0)
.nid.pi.point <- function(f, Q, rho) {
  Ls <- .nid.all.L(f$R, Q); best <- -Inf; arg <- NULL
  for (k in base::seq_len(base::nrow(Ls))) {
    L <- Ls[k, ]; U <- base::sweep(f$Xi[, L, drop = FALSE], 2, f$s[L], "/")
    Up <- base::pmax(U, 0); sh <- Up / base::pmax(base::rowSums(Up), 1e-300); ok <- U >= rho
    m <- base::vapply(base::seq_len(Q), function(q) if (base::any(ok[, q])) base::max(sh[ok[, q], q]) else -Inf, 0)
    jj <- base::vapply(base::seq_len(Q), function(q)
      if (base::any(ok[, q])) base::which(ok[, q])[base::which.max(sh[ok[, q], q])] else NA_integer_, 0L)
    if (base::min(m) > best && !base::anyDuplicated(jj)) { best <- base::min(m); arg <- base::list(L = L, J = jj) }
  }
  base::c(base::list(value = best), arg)
}

## ===================== maximum over the model (R = Q) =====================

## sufficient statistics of the row-wise regressions
.nid.suff <- function(Y, A) {
  G <- A %*% base::t(A); B <- Y %*% base::t(A) %*% base::solve(G)
  W <- base::rowSums((Y - B %*% A)^2)
  base::list(B = B, G = G, Wsum = base::sum(W), N = base::ncol(Y), P = base::nrow(Y), R = base::nrow(A))
}

## Lawson-Hanson non-negative least squares: argmin_{x >= 0} ||M x - b||
.nid.nnls <- function(M, b, tol = 1e-10) {
  n <- base::ncol(M); x <- base::numeric(n); Pset <- base::logical(n)
  w <- base::as.vector(base::crossprod(M, b - M %*% x))
  it <- 0
  while (base::any(!Pset & w > tol) && it < 3 * n + 30) {
    it <- it + 1
    j <- base::which(!Pset)[base::which.max(w[!Pset])]; Pset[j] <- TRUE
    repeat {
      z <- base::numeric(n)
      z[Pset] <- base::qr.solve(M[, Pset, drop = FALSE], b)
      if (base::all(z[Pset] > 0)) { x <- z; break }
      neg <- Pset & z <= 0
      a <- base::min(x[neg] / (x[neg] - z[neg]))
      x <- x + a * (z - x)
      Pset <- Pset & x > tol
      x[!Pset] <- 0
    }
    w <- base::as.vector(base::crossprod(M, b - M %*% x))
  }
  x
}

## row-wise NNLS in the metric G: the maximum over the model (R = Q) for both the common-variance
## likelihood and nmfkc's least squares (both separate over rows)
.nid.xi.hat <- function(s) {
  L <- base::chol(s$G)
  base::t(base::apply(s$B, 1, function(b) .nid.nnls(L, base::as.vector(L %*% b))))
}

## -2 log L (common variance, profiled out) and its gradient in Xi
.nid.f <- function(Xi, s) {
  D <- Xi - s$B; RSS <- s$Wsum + base::sum((D %*% s$G) * D)
  s$N * s$P * base::log(RSS / (s$N * s$P))
}
.nid.g <- function(Xi, s) {
  D <- Xi - s$B; RSS <- s$Wsum + base::sum((D %*% s$G) * D)
  2 * s$N * s$P / RSS * (D %*% s$G)
}

## ===================== identification region (R = Q) =====================

.nid.cos <- function(M) base::sweep(M, 2, base::sqrt(base::pmax(base::colSums(M^2), 1e-300)), "/")

## best matching of the columns of Xa to those of Xb (cosine): m[q] = column of Xb for column q of Xa
.nid.match <- function(Xa, Xb) {
  Q <- base::ncol(Xb); C <- base::crossprod(.nid.cos(base::pmax(Xa, 0)), .nid.cos(Xb))
  if (Q <= 7) {
    pm <- .nid.perms(Q)
    return(pm[base::which.max(base::apply(pm, 1, function(p) base::sum(C[base::cbind(base::seq_len(Q), p)]))), ])
  }
  m <- base::integer(Q); free_a <- base::seq_len(Q); free_b <- base::seq_len(Q)   # greedy for large Q
  for (it in base::seq_len(Q)) {
    sub <- C[free_a, free_b, drop = FALSE]; w <- base::which(sub == base::max(sub), arr.ind = TRUE)[1, ]
    m[free_a[w[1]]] <- free_b[w[2]]; free_a <- free_a[-w[1]]; free_b <- free_b[-w[2]]
  }
  m
}
## labels kept: the matching of the columns of X to those of Xt is the identity
.nid.label.ok <- function(X, Xt) base::all(.nid.match(X, Xt) == base::seq_len(base::ncol(Xt)))

.nid.feas <- function(N, Xt, tol = 1e-10) {
  if (base::any(N < -tol)) return(FALSE)
  Ni <- base::tryCatch(base::solve(N), error = function(e) NULL)
  if (base::is.null(Ni)) return(FALSE)
  X <- Xt %*% Ni
  !base::any(X < -tol) && .nid.label.ok(X, Xt)
}
.nid.viol <- function(N, Xt) {
  v <- base::sum(base::pmax(-N, 0))
  Ni <- base::tryCatch(base::solve(N), error = function(e) NULL)
  if (base::is.null(Ni)) return(1e3)
  X <- Xt %*% Ni; v <- v + 100 * base::sum(base::pmax(-X, 0))
  if (!.nid.label.ok(X, Xt)) v <- v + 1
  v
}
.nid.N <- function(par, Q) { M <- base::matrix(par, Q - 1, Q); base::rbind(M, 1 - base::colSums(M)) }

## random feasible factorizations: moves that keep the column sums at 1, in both signs (an exact zero in
## Xi.hat makes every one-signed mixing of I infeasible) and sparse positive (N = I is a corner of N >= 0)
.nid.sample <- function(Xt, n) {
  Q <- base::ncol(Xt); keep <- base::list(base::diag(Q))
  rE <- function() { E <- base::matrix(stats::rnorm(Q * Q), Q); base::sweep(E, 2, base::colMeans(E)) }
  rP <- function() {
    E <- base::matrix(base::abs(stats::rnorm(Q * Q)) * (stats::runif(Q * Q) < 0.5), Q)
    base::diag(E) <- 0; E - base::diag(base::colSums(E), Q)
  }
  for (i in base::seq_len(n)) {
    u <- stats::runif(1)
    N <- if (u < 0.4) keep[[base::sample.int(base::length(keep), 1)]] + 10^stats::runif(1, -3, -0.3) * rE()
         else if (u < 0.8) keep[[base::sample.int(base::length(keep), 1)]] + 10^stats::runif(1, -3, -0.3) * rP()
         else base::diag(Q) + 10^stats::runif(1, -3, 0) * rP()
    if (.nid.feas(N, Xt)) keep[[base::length(keep) + 1]] <- N
  }
  base::array(base::unlist(keep), base::c(Q, Q, base::length(keep)))
}

.nid.region <- function(Xt, s0, n, nstart = 3) {
  Q <- base::ncol(Xt); P <- base::nrow(Xt)
  arr <- .nid.sample(Xt, n)
  C.lo <- C.hi <- base::matrix(NA_real_, Q, Q); arg <- base::list()
  for (q in base::seq_len(Q)) for (k in base::seq_len(Q)) for (dir in base::c("lo", "hi")) {
    sg <- if (dir == "lo") 1 else -1
    v <- sg * arr[q, k, ]; starts <- base::order(v)[base::seq_len(base::min(nstart, base::length(v)))]
    ob <- function(par) { N <- .nid.N(par, Q); sg * N[q, k] + 1e4 * .nid.viol(N, Xt) }
    Nb <- arr[, , starts[1]]
    for (j in starts) {                       # refine from the most extreme samples; keep the best feasible
      par <- base::as.vector(arr[-Q, , j]); val <- ob(par)
      for (rs in 1:4) {                       # Nelder-Mead stalls; restart it from where it stopped
        op <- stats::optim(par, ob, method = "Nelder-Mead", control = base::list(maxit = 4000, reltol = 1e-12))
        if (op$value >= val - 1e-12) break
        par <- op$par; val <- op$value
      }
      N1 <- .nid.N(par, Q)
      if (.nid.feas(N1, Xt, tol = 1e-9) && sg * N1[q, k] < sg * Nb[q, k]) Nb <- N1
    }
    if (dir == "lo") C.lo[q, k] <- Nb[q, k] * s0[k] else C.hi[q, k] <- Nb[q, k] * s0[k]
    arg[[base::paste(q, k, dir)]] <- Nb
  }
  ## X ranges from the samples and the extreme points (inner approximation)
  pts <- base::c(base::lapply(base::seq_len(base::dim(arr)[3]), function(i) arr[, , i]), arg)
  Xs <- base::vapply(pts, function(N) Xt %*% base::solve(N), base::matrix(0, P, Q))
  base::list(C.lo = C.lo, C.hi = C.hi, X.lo = base::pmax(base::apply(Xs, c(1, 2), base::min), 0),
             X.hi = base::apply(Xs, c(1, 2), base::max), argN = arg, samples = arr, n = base::dim(arr)[3])
}

## the feasible N whose C = N S is closest (least squares) to a target C (internal basis order)
.nid.nearest <- function(target, reg, Xt, s0) {
  Q <- base::ncol(Xt)
  gap <- function(N) base::sum((N %*% base::diag(s0, Q) - target)^2)
  cand <- base::c(base::lapply(base::seq_len(base::dim(reg$samples)[3]), function(i) reg$samples[, , i]), reg$argN)
  N1 <- cand[[base::which.min(base::vapply(cand, gap, 0))]]
  ob <- function(par) { N <- .nid.N(par, Q); gap(N) + 1e6 * .nid.viol(N, Xt) }
  op <- stats::optim(base::as.vector(N1[-Q, ]), ob, method = "Nelder-Mead",
                     control = base::list(maxit = 20000, reltol = 1e-14))
  N2 <- .nid.N(op$par, Q)
  if (.nid.feas(N2, Xt, tol = 1e-10) && gap(N2) < gap(N1)) N2 else N1
}

## put the bases of an nmfkc fit in a given order (idx[q] = fit's basis for basis q) and refresh the
## quantities derived from X, B and C
.nid.reorder <- function(fit, idx, Y, A, object) {
  if (base::identical(base::as.integer(idx), base::seq_along(idx))) return(fit)
  fit$X <- fit$X[, idx, drop = FALSE]; fit$C <- fit$C[idx, , drop = FALSE]; fit$B <- fit$B[idx, , drop = FALSE]
  base::dimnames(fit$X) <- base::dimnames(object$X); base::dimnames(fit$C) <- base::dimnames(object$C)
  base::rownames(fit$B) <- base::rownames(object$C)
  .eps <- 1e-10
  fit$B.prob <- base::t(base::t(fit$B) / (base::colSums(fit$B) + .eps))
  fit$B.cluster <- base::apply(fit$B.prob, 2, base::which.max)
  fit$B.cluster[base::colSums(fit$B.prob) == 0] <- NA
  fit$X.prob <- fit$X / (base::rowSums(fit$X) + .eps)
  fit$X.cluster <- base::apply(fit$X.prob, 1, base::which.max)
  if (!base::is.null(fit$anchors)) fit$anchors <- .nmfkc.anchor(fit$X, fit$C, A)
  fit
}

## ===================== profile-likelihood intervals (R = Q) =====================

.nid.softmax <- function(al) {
  e <- base::exp(base::sweep(al, 2, base::apply(al, 2, base::max))); base::sweep(e, 2, base::colSums(e), "/")
}

## minimize -2 log L over X (column softmax) and C >= 0 with C[fix] = t
.nid.fit.fixed <- function(s, Q, X0, C0, fix, t, pscale) {
  P <- s$P; R <- s$R
  free <- base::rep(TRUE, Q * R); free[(fix[2] - 1) * Q + fix[1]] <- FALSE
  unpack <- function(par) {
    al <- base::matrix(par[base::seq_len(P * Q)], P, Q); th <- base::numeric(Q * R)
    th[free] <- par[-base::seq_len(P * Q)]; th[!free] <- t
    base::list(X = .nid.softmax(al), C = base::matrix(th, Q, R))
  }
  fn <- function(par) { u <- unpack(par); .nid.f(u$X %*% u$C, s) }
  gr <- function(par) {
    u <- unpack(par); G <- .nid.g(u$X %*% u$C, s)
    GX <- G %*% base::t(u$C)
    Gal <- u$X * base::sweep(GX, 2, base::colSums(u$X * GX))
    base::c(base::as.vector(Gal), base::as.vector(base::crossprod(u$X, G))[free])
  }
  X0 <- base::pmax(X0, 1e-6); X0 <- base::sweep(X0, 2, base::colSums(X0), "/"); C0[fix[1], fix[2]] <- t
  par0 <- base::c(base::as.vector(base::log(X0)), base::as.vector(C0)[free]); nf <- base::sum(free)
  op <- stats::optim(par0, fn, gr, method = "L-BFGS-B",
                     lower = base::c(base::rep(-40, P * Q), base::rep(0, nf)),
                     upper = base::c(base::rep(40, P * Q), base::rep(Inf, nf)),
                     control = base::list(maxit = 20000, factr = 1e3, pgtol = 0,
                                          parscale = base::c(base::rep(1, P * Q), base::rep(pscale, nf))))
  u <- unpack(op$par)
  base::list(f = .nid.f(u$X %*% u$C, s), X = u$X, C = u$C)
}

.nid.perturb <- function(X, C, eps) {
  Q <- base::ncol(X); Tm <- base::diag(Q) + eps * (base::matrix(stats::runif(Q * Q), Q) - 1 / Q)
  Tm <- base::sweep(Tm, 2, base::colSums(Tm), "/")
  X1 <- base::pmax(X %*% Tm, 1e-6); X1 <- base::sweep(X1, 2, base::colSums(X1), "/")
  base::list(X = X1, C = base::pmax(base::solve(Tm, C), 0))
}

.nid.profile <- function(s, Xt, s0, reg, Q, level, step, seed, cores) {
  fmax <- .nid.f(Xt %*% base::diag(s0, Q), s); cut <- stats::qchisq(level, 1)
  pscale <- base::mean(s0) / Q
  jobs <- base::expand.grid(dir = base::c("lo", "hi"), q = base::seq_len(Q), k = base::seq_len(Q),
                            stringsAsFactors = FALSE)
  run_one <- function(i) {
    base::set.seed(seed + i)
    jb <- jobs[i, ]; Nb <- reg$argN[[base::paste(jb$q, jb$k, jb$dir)]]
    cur <- base::list(X = Xt %*% base::solve(Nb), C = Nb %*% base::diag(s0, Q))
    t_b <- if (jb$dir == "lo") reg$C.lo[jb$q, jb$k] else reg$C.hi[jb$q, jb$k]
    sg <- if (jb$dir == "hi") 1 else -1; h <- step * s0[jb$k]
    tt_all <- t_b; ff <- fmax; states <- base::list(cur); at0 <- FALSE
    try_at <- function(t, base_state) {
      cand <- base::c(base::list(base_state), base::lapply(1:2, function(r) .nid.perturb(base_state$X, base_state$C, 0.05)))
      rs <- base::lapply(cand, function(st) .nid.fit.fixed(s, Q, st$X, st$C, base::c(jb$q, jb$k), t, pscale))
      fv <- base::vapply(rs, function(r) if (.nid.label.ok(r$X, Xt)) r$f else Inf, 0)
      if (!base::any(base::is.finite(fv))) return(NULL)
      rs[[base::which.min(fv)]]
    }
    for (j in base::seq_len(80)) {
      t <- t_b + sg * j * h
      if (t > 2 * s0[jb$k]) break
      if (t < 0) { if (t_b + sg * (j - 1) * h > 1e-8) { t <- 0; at0 <- TRUE } else break }
      r <- try_at(t, cur); if (base::is.null(r)) break
      cur <- r; tt_all <- base::c(tt_all, t); ff <- base::c(ff, r$f); states[[base::length(states) + 1]] <- r
      if (r$f - fmax > cut + 2 || at0) break
    }
    g <- ff - fmax; ib <- base::which(g > cut)[1]
    end <- if (base::is.na(ib)) (if (sg > 0) Inf else base::max(0, base::min(tt_all))) else {
      ta <- tt_all[ib - 1]; tb <- tt_all[ib]; sa <- states[[ib - 1]]
      for (hh in 1:8) {
        tm <- (ta + tb) / 2; r <- try_at(tm, sa)
        if (!base::is.null(r) && r$f - fmax <= cut) { ta <- tm; sa <- r } else tb <- tm
      }
      (ta + tb) / 2
    }
    dip <- if (base::length(ff) > 1) base::max(0, -base::diff(ff)) else 0
    base::list(end = end, dip = dip)
  }
  res <- .nmfkc.parlapply(base::seq_len(base::nrow(jobs)), run_one, cores = cores)
  lo <- hi <- base::matrix(NA_real_, Q, Q)
  for (i in base::seq_len(base::nrow(jobs))) {
    jb <- jobs[i, ]
    if (jb$dir == "lo") lo[jb$q, jb$k] <- base::max(0, res[[i]]$end) else hi[jb$q, jb$k] <- res[[i]]$end
  }
  base::list(lo = lo, hi = hi, dip = base::max(base::vapply(res, `[[`, 0, "dip")))
}

## ===================== methods =====================

#' @export
print.nmfkc.ident <- function(x, ...) {
  t <- x$test
  base::cat("Identifiability of an nmfkc fit (rank ", x$rank, ")\n", sep = "")
  base::cat(base::sprintf("  anchor test H(pi = %g, rho = %.3g): %s (alpha = %g)\n", t$pi, t$rho,
                          if (t$reject) "rejected" else "not rejected", t$alpha))
  base::cat(base::sprintf("  pi.max = %s  (point estimate %s);  rho.max = %s\n",
                          base::formatC(t$pi.max, digits = 3, format = "f"),
                          base::formatC(t$pi.point, digits = 3, format = "f"),
                          base::formatC(t$rho.max, digits = 3, format = "f")))
  if (base::is.null(x$Xi.hat)) {
    base::cat("  identification region: not computed (needs nrow(A) == rank)\n")
    return(base::invisible(x))
  }
  base::cat(base::sprintf("  identification region of C: %s\n",
                          if (base::isTRUE(x$identified)) "a single point (identified)" else "has width (not identified)"))
  base::cat(base::sprintf("  repaired fit: start = %s, %d iterations, largest change of C %.2g\n",
                          x$repair$start, x$repair$iter, x$repair$change))
  base::print(summary(x)$coefficients, row.names = FALSE, digits = 4)
  base::invisible(x)
}

#' @export
summary.nmfkc.ident <- function(object, ...) {
  x <- object
  co <- NULL
  if (!base::is.null(x$Xi.hat)) {
    C <- x$C
    co <- base::data.frame(Basis = base::rep(base::rownames(C), base::ncol(C)),
                           Covariate = base::rep(base::colnames(C), each = base::nrow(C)),
                           Estimate = base::as.vector(C),
                           region_low = base::as.vector(x$C.region.lower),
                           region_high = base::as.vector(x$C.region.upper),
                           stringsAsFactors = FALSE)
    if (!base::is.null(x$C.ci.lower)) { co$CI_low <- base::as.vector(x$C.ci.lower); co$CI_high <- base::as.vector(x$C.ci.upper) }
  }
  base::structure(base::list(test = x$test, identified = x$identified, coefficients = co,
                             ci.level = x$ci.level, rank = x$rank), class = "summary.nmfkc.ident")
}

#' @export
print.summary.nmfkc.ident <- function(x, ...) {
  base::cat(base::sprintf("pi.max = %.3f, rho.max = %.3f, identified = %s\n", x$test$pi.max, x$test$rho.max,
                          base::format(x$identified)))
  if (!base::is.null(x$coefficients)) base::print(x$coefficients, row.names = FALSE, digits = 4)
  base::invisible(x)
}

#' @title Plot the identification region and intervals of an nmfkc.ident object
#' @description One panel per entry of \eqn{C} (rows = bases, columns =
#'   covariates): the identification region (black), the profile-likelihood
#'   interval if computed (gray) and the repaired estimate (red cross).
#' @param x An object of class \code{"nmfkc.ident"}.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @export
plot.nmfkc.ident <- function(x, ...) {
  if (base::is.null(x$Xi.hat)) stop("no identification region (needs nrow(A) == rank)")
  Q <- base::nrow(x$C); R <- base::ncol(x$C)
  hi <- if (!base::is.null(x$C.ci.upper)) x$C.ci.upper else x$C.region.upper
  xmax <- base::max(hi[base::is.finite(hi)], x$C) * 1.05
  op <- graphics::par(mfrow = base::c(Q, R), mar = base::c(2, 2, 1.6, 0.5))
  base::on.exit(graphics::par(op))
  for (q in base::seq_len(Q)) for (k in base::seq_len(R)) {
    graphics::plot(NA, xlim = base::c(0, xmax), ylim = base::c(0.5, 1.5), yaxt = "n", xlab = "", ylab = "",
                   main = base::paste(base::rownames(x$C)[q], "x", base::colnames(x$C)[k]), cex.main = 0.8)
    if (!base::is.null(x$C.ci.lower))
      graphics::segments(x$C.ci.lower[q, k], 1, x$C.ci.upper[q, k], 1, lwd = 2.5, col = "gray55")
    graphics::segments(x$C.region.lower[q, k], 1, x$C.region.upper[q, k], 1, lwd = 6, lend = 1)
    graphics::points(x$C.region.lower[q, k], 1, pch = 19, cex = 0.6)
    graphics::points(x$C[q, k], 1, pch = 4, col = "red3", cex = 1.6, lwd = 2)
  }
  base::invisible(x)
}

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
.nmfre.seminmf.X <- function(X, Y, B, S = NULL, pen.num = 0, pen.den = 0, eps = 1e-10) {
  YBt <- tcrossprod(Y, B)
  G   <- B %*% t(B)   # keep the gemm form: one-arg tcrossprod (dsyrk) is not bitwise-identical
  if (!is.null(S)) G <- G + S
  pos <- function(M) pmax(M, 0)
  neg <- function(M) pmax(-M, 0)
  ## pen.num / pen.den: non-negative MAP-penalty contributions on X
  ## (ortho / smooth), added to the numerator / denominator before the sqrt.
  numer <- pos(YBt) + X %*% neg(G) + pen.num
  denom <- neg(YBt) + X %*% pos(G) + pen.den
  X * sqrt(.nmfre.safe.div(numer, denom, eps = eps))
}

#' Exact non-negative least squares in normal-equation form
#'
#' Minimizes \eqn{\tfrac12 x^\top H x - g^\top x} subject to \eqn{x\ge0}
#' (\eqn{H} symmetric positive definite) by the active-set method of Lawson
#' and Hanson (1974, ch. 23), which ends at the exact solution after finitely
#' many steps.  Used for the \eqn{C}-step of \code{nmfre(C.signed = FALSE)},
#' where \eqn{H=(AA^\top)\otimes(X^\top X)} is only \eqn{QK\times QK}.
#' @keywords internal
#' @noRd
.nmfre.nnls <- function(H, g, tol = 1e-10, maxit = 30L * length(g)) {
  n <- length(g); x <- numeric(n); P <- rep(FALSE, n)
  sc <- max(1, max(abs(g)))
  w <- as.vector(g - H %*% x); it <- 0L
  while (it < maxit && any(!P & w > tol * sc)) {
    it <- it + 1L
    j <- which(!P)[which.max(w[!P])]; P[j] <- TRUE
    repeat {
      z <- numeric(n); z[P] <- solve(H[P, P, drop = FALSE], g[P])
      if (all(z[P] > 0)) { x <- z; break }
      ## step back to the boundary and drop the variables that reach zero
      Qs <- P & z <= 0
      alpha <- min(x[Qs] / (x[Qs] - z[Qs]))
      x <- x + alpha * (z - x)
      P <- P & x > tol * max(1, max(abs(x)))
      x[!P] <- 0
      it <- it + 1L
      if (it >= maxit) break
    }
    w <- as.vector(g - H %*% x)
  }
  x
}

#' Row-wise non-negative least squares for many rows at once
#'
#' For each row \eqn{c_p} of \code{Cm} solves
#' \eqn{\min_{x\ge0,\ x_j=0\ (j\notin \mathrm{allow}_p)} \tfrac12 x^\top G x - c_p^\top x}
#' (\eqn{G} positive definite).  For small \eqn{Q} every free set is tried for all
#' rows at once and each row keeps the one meeting the KKT conditions (unique,
#' the problem being strictly convex); rows left over (ties at the tolerance) and
#' large \eqn{Q} use \code{.nmfre.nnls()} one row at a time.
#' @keywords internal
#' @noRd
.nmfre.rows.nnls <- function(G, Cm, allow = NULL, tol = 1e-12) {
  P <- nrow(Cm); Q <- ncol(Cm)
  if (is.null(allow)) allow <- matrix(TRUE, P, Q)
  X <- matrix(0, P, Q)
  one.by.one <- function(rr) for (p in rr) {
    a <- allow[p, ]
    if (any(a)) X[p, a] <<- .nmfre.nnls(G[a, a, drop = FALSE], Cm[p, a])
  }
  if (Q > 8L) { one.by.one(seq_len(P)); return(X) }
  sc <- tol * max(1, max(abs(Cm)))
  done <- rep(FALSE, P)
  pats <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), Q)))
  pats <- pats[order(rowSums(pats)), , drop = FALSE]
  for (k in seq_len(nrow(pats))) {
    f <- pats[k, ]
    rr <- which(!done & rowSums(allow[, f, drop = FALSE]) == sum(f))   # the pattern is allowed for the row
    if (!length(rr)) next
    Xc <- matrix(0, length(rr), Q)
    if (any(f)) Xc[, f] <- Cm[rr, f, drop = FALSE] %*% solve(G[f, f, drop = FALSE])
    gr <- Xc %*% G - Cm[rr, , drop = FALSE]                             # gradient
    ok <- rowSums(Xc[, f, drop = FALSE] < -sc) == 0 &
          rowSums((gr < -sc) & !matrix(f, length(rr), Q, byrow = TRUE) & allow[rr, , drop = FALSE]) == 0
    if (any(ok)) { X[rr[ok], ] <- pmax(Xc[ok, , drop = FALSE], 0); done[rr[ok]] <- TRUE }
    if (all(done)) break
  }
  if (any(!done)) one.by.one(which(!done))
  X
}

#' Exact X-step on the column simplices
#'
#' Minimizes \eqn{\mathrm{tr}(XGX^\top) - 2\,\mathrm{tr}(X^\top C)} over
#' \eqn{X\ge0} with \eqn{\mathbf 1^\top X=\mathbf 1^\top} (and \eqn{X=0} where
#' \code{mask} is \code{TRUE}), a convex quadratic program.  Given the
#' multipliers \eqn{\mu} of the column sums the rows decouple into
#' \eqn{Q}-variable non-negative least squares (\code{.nmfre.nnls()}), and
#' \eqn{\mu} solves \eqn{\mathbf 1^\top X(\mu)=\mathbf 1^\top} by a damped
#' semismooth Newton method on the concave dual.  With
#' \eqn{G=BB^\top+S} and \eqn{C=YB^\top} this is the \eqn{X}-step of
#' \code{nmfre(X.update = "exact")}: the block is minimized exactly under the
#' column-sum constraint, so no renormalization is needed.
#' @keywords internal
#' @noRd
.nmfre.simplex.X <- function(G, CY, mu = NULL, mask = NULL, tol = 1e-10, maxit = 100L) {
  P <- nrow(CY); Q <- ncol(CY)
  if (is.null(mu) || length(mu) != Q || any(!is.finite(mu))) mu <- numeric(Q)
  allow <- if (is.null(mask)) matrix(TRUE, P, Q) else !mask
  rows <- function(mu) .nmfre.rows.nnls(G, sweep(CY, 2, mu), allow)
  dual <- function(X, mu) sum((X %*% G) * X) - 2 * sum(X * sweep(CY, 2, mu)) - 2 * sum(mu)
  X <- rows(mu); Fv <- colSums(X) - 1; ph <- dual(X, mu); it <- 0L
  while (max(abs(Fv)) > tol && it < maxit) {
    it <- it + 1L
    ## generalized Jacobian of colSums(X(mu)) (negated): sum over rows of G_FF^{-1},
    ## accumulated once per free-set pattern
    Mm <- matrix(0, Q, Q)
    code <- as.vector((X > 0) %*% 2^(seq_len(Q) - 1))
    for (cd in unique(code[code > 0])) {
      f <- bitwAnd(cd, 2^(seq_len(Q) - 1)) > 0
      Mm[f, f] <- Mm[f, f] + sum(code == cd) * solve(G[f, f, drop = FALSE])
    }
    dg <- diag(Mm); fill <- if (any(dg > 0)) max(dg) else 1 / max(diag(G))
    diag(Mm)[dg <= 0] <- fill                      # a column with no free entry: a scaled gradient step
    d <- solve(Mm + diag(1e-10 * max(diag(Mm)), Q), Fv)
    a <- 1; ok <- FALSE; F2 <- sum(Fv^2)
    repeat {
      mu1 <- mu + a * d; X1 <- rows(mu1); ph1 <- dual(X1, mu1); Fv1 <- colSums(X1) - 1
      ## Armijo ascent on the concave dual, or -- once the dual's change is below
      ## rounding -- a decrease of the column-sum violation itself
      if (ph1 >= ph + 2e-4 * a * sum(Fv * d) || sum(Fv1^2) < (1 - 1e-4 * a) * F2) { ok <- TRUE; break }
      if (a < 1e-8) break
      a <- a / 2
    }
    if (!ok) break                                   # no progress possible: keep the best point
    mu <- mu1; X <- X1; ph <- ph1; Fv <- Fv1
  }
  cs <- colSums(X); X <- sweep(X, 2, ifelse(cs > 0, cs, 1), "/")   # remove the last ~tol of infeasibility
  list(X = X, mu = mu, iter = it, viol = max(abs(Fv)))
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
.nmfre.marginal.nll <- function(Y, X, C, A, sigma2, tau2, eps = 1e-12, CA = NULL) {
  P <- nrow(Y); N <- ncol(Y); Q <- ncol(X)
  s2 <- max(sigma2, eps); t2 <- max(tau2, eps)
  lambda <- s2 / t2
  if (is.null(CA)) CA <- C %*% A                    # caller may pass C %*% A
  Rm   <- Y - X %*% CA                              # P x N fixed-effects residual
  XtRm <- crossprod(X, Rm)                          # Q x N
  XtX  <- crossprod(X)                              # Q x Q (computed once)
  Mll  <- XtX + diag(pmax(lambda, 1e-12), Q)
  quad <- tryCatch((sum(Rm^2) - sum(XtRm * solve(Mll, XtRm))) / s2,
                   error = function(e) NA_real_)
  logdetQ <- tryCatch(
    as.numeric(determinant(diag(Q) + (t2 / s2) * XtX, logarithm = TRUE)$modulus),
    error = function(e) NA_real_)
  if (!is.finite(quad) || !is.finite(logdetQ)) return(NA_real_)
  0.5 * (N * (P * log(2 * pi) + P * log(s2) + logdetQ) + quad)
}

#' Degeneracy diagnostic of an NMF-RE fit (two bases merging, C diverging)
#'
#' With \eqn{Q \ge 2} and a signed \eqn{C}, a fit can drift along a path on which
#' two columns of \eqn{X} approach each other while the corresponding rows of
#' \eqn{C} grow with opposite signs: the difference \eqn{(x_q - x_{q'}) c_q} then
#' represents a direction outside the non-negative cone, and the fit does not
#' settle.  Two bases suffice, so this is not tied to the rotational
#' indeterminacy of \eqn{Q \ge 3} (seen at \eqn{Q = 2} on Orthodont, a rank the
#' data do not support, and on test data).  In simulations such fits usually
#' had a clearly LOWER marginal likelihood than a non-degenerate fit, i.e. the
#' algorithm had been trapped -- in the HS-based simulation of the NMF-RE paper
#' (\eqn{Q = 3}) mostly by the old variance start \eqn{\sigma^2 = \tau^2 = 1},
#' from which every data-driven basis start fell in; refit from another start
#' and keep the larger likelihood.  Three continuous measures (columns of
#' \eqn{X} sum to one, so \eqn{\|x_q\|_1 = 1}):
#' \code{max.cosine}, the largest cosine similarity between two columns of
#' \eqn{X} (and that \code{pair}); \code{cancellation}, the largest over
#' covariates \eqn{k} of \eqn{\sum_q |c_{qk}| / \|X c_{\cdot k}\|_1}, i.e. how far
#' the composite effect is a difference of large opposite terms; and
#' \code{kappa.X}, the condition number of \eqn{X}.  \code{flag} is
#' \code{max.cosine >= cos.tol || cancellation >= canc.tol}.  The defaults were
#' set on two simulation designs (non-degenerate fits: cosine 0.55--0.83,
#' cancellation 1.9--3.8; degenerate: cosine 0.89--0.995, cancellation
#' 5.5--680); read the measures, not only the flag.
#' @keywords internal
#' @noRd
.nmfre.degeneracy <- function(X, C, cos.tol = 0.9, canc.tol = 5) {
  Q <- ncol(X)
  if (Q < 2) return(list(max.cosine = NA_real_, pair = c(NA_integer_, NA_integer_),
                         cancellation = NA_real_, kappa.X = NA_real_, flag = FALSE,
                         cos.tol = cos.tol, canc.tol = canc.tol))
  nrm <- sqrt(colSums(X^2)); nrm[nrm == 0] <- 1
  cs <- crossprod(sweep(X, 2, nrm, "/")); diag(cs) <- -Inf
  ij <- which(cs == max(cs), arr.ind = TRUE)[1, ]
  Xi <- X %*% C
  den <- colSums(abs(Xi)); den[den == 0] <- NA
  canc <- suppressWarnings(max(colSums(abs(C)) / den, na.rm = TRUE))
  kx <- tryCatch(kappa(X, exact = TRUE), error = function(e) NA_real_)
  canc <- if (is.finite(canc)) canc else NA_real_
  list(max.cosine = max(cs), pair = sort(as.integer(ij)), cancellation = canc,
       kappa.X = kx, flag = isTRUE(max(cs) >= cos.tol) || isTRUE(canc >= canc.tol),
       cos.tol = cos.tol, canc.tol = canc.tol)
}

#' Variance start of nmfre() from the residual covariance (Internal)
#'
#' Given \eqn{A}, the residual \eqn{R = Y(I - H_A)} has covariance
#' \eqn{\tau^2 XX' + \sigma^2 I}: its eigenvalues beyond the \eqn{Q}-th are
#' \eqn{\sigma^2} and the top \eqn{Q} exceed it by \eqn{\tau^2 d_q} (\eqn{d_q} the
#' eigenvalues of \eqn{X'X}), as in probabilistic PCA (Tipping and Bishop, 1999,
#' JRSS B 61, 611--622, eqs. 7--8).  So \eqn{\sigma^2} starts at the mean of the
#' eigenvalues beyond \eqn{Q} (within the rank of \eqn{R}) and \eqn{\tau^2} at
#' U's share \eqn{\sum_{q \le Q} (l_q - \sigma^2)_+} over \eqn{\mathrm{tr}(X'X)}
#' of the start basis.  In the HS simulation of the NMF-RE paper (2800
#' samples, true \eqn{\lambda} about 0.06) the start \eqn{\sigma^2 = \tau^2 = 1}
#' led every data-driven basis start (k-means, SPA, NNDSVD, random) into a
#' degenerate fit in 509 samples; from this start none did, and the other
#' samples reached the same fits.
#' @return \code{c(sigma2, tau2)}, or \code{NULL} when \eqn{Y} has missing
#'   values or the residual has no eigenvalues beyond the \eqn{Q}-th.
#' @keywords internal
#' @noRd
.nmfre.variance.start <- function(Y, A, X) {
  Q <- ncol(X); N <- ncol(Y); K <- nrow(A)
  if (anyNA(Y)) return(NULL)
  R <- Y - (Y %*% t(A)) %*% solve(tcrossprod(A) + diag(1e-10, K), A)
  r <- min(nrow(Y), N - K)                         # rank of the residual
  if (r <= Q) return(NULL)
  l <- eigen(tcrossprod(R) / (N - K), symmetric = TRUE, only.values = TRUE)$values[seq_len(r)]
  s2 <- mean(l[(Q + 1):r])
  if (!is.finite(s2) || s2 <= 0) return(NULL)
  share <- max(sum(pmax(l[seq_len(Q)] - s2, 0)), 1e-3 * sum(l))
  c(sigma2 = s2, tau2 = share / sum(X^2))          # tr(X'X) = sum(X^2)
}

#' Resolve fewer than Q anchor rows for nmfre() (Internal)
#'
#' \code{X.anchor} with fewer rows than bases anchors the first bases only:
#' row \code{J[q]} carries basis \code{q} alone, the other bases are free.  For
#' \eqn{Q = 3} one exact anchor row already pins the rotation that leaves the
#' likelihood unchanged, so an anchor known from the design can be imposed
#' without imposing anchors the data do not have.  The checks are those of
#' \code{.nmfkc_anchor_rows()} for the rows given (no SPA here: "spa" picks
#' all \code{Q} rows).
#' @param Y The response as the initialisation sees it (shifted to be non-negative).
#' @return Integer row indices, one per anchored basis.
#' @keywords internal
#' @noRd
.nmfre.partial.anchor.rows <- function(X.anchor, Y, Q) {
  P <- nrow(Y)
  if (is.character(X.anchor)) {
    if (is.null(rownames(Y)))
      stop("nmfre(): X.anchor gives row names, but Y has none.  Give row indices instead.",
           call. = FALSE)
    J <- match(X.anchor, rownames(Y))
    if (anyNA(J))
      stop("nmfre(): X.anchor names row(s) not in Y: ",
           paste(X.anchor[is.na(J)], collapse = ", "), ".", call. = FALSE)
  } else if (is.numeric(X.anchor)) {
    if (anyNA(X.anchor) || any(X.anchor != round(X.anchor)) ||
        any(X.anchor < 1) || any(X.anchor > P))
      stop(sprintf("nmfre(): X.anchor must be row indices between 1 and nrow(Y) = %d.", P),
           call. = FALSE)
    J <- as.integer(X.anchor)
  } else {
    stop("nmfre(): X.anchor must be NULL, \"spa\", or anchor rows of Y (indices or names).",
         call. = FALSE)
  }
  if (anyDuplicated(J))
    stop("nmfre(): X.anchor repeats a row; each basis needs its own anchor row.", call. = FALSE)
  empty <- J[rowSums(Y[J, , drop = FALSE] > 0) == 0]
  if (length(empty))
    stop(sprintf(paste("nmfre(): anchor row(s) %s of Y have no positive entry.  An anchor row",
                       "carries its basis alone, so that basis would be zero."),
                 paste(empty, collapse = ", ")), call. = FALSE)
  if (qr(Y[J, , drop = FALSE])$rank < length(J))
    stop("nmfre(): the anchor rows of Y are linearly dependent.  Choose other rows.", call. = FALSE)
  J
}

#' nmfre() with fewer anchor rows than bases: assignment chosen by likelihood (Internal)
#'
#' Which basis an anchor row belongs to is part of the estimate, so it is
#' chosen by the marginal likelihood: (1) a free fit -- from \code{X.init}
#' (by default \code{"mix2"}: the better of the k-means and SPA starts);
#' (2) for every ordered choice of distinct bases for the anchor rows, a refit
#' started from the free fit with each anchor row made pure on its basis (so the
#' start is already on the right side of every other local solution); (3) the
#' refit with the smallest marginal NLL.  In an HS1939-based simulation
#' (\code{rank = 3}, one true anchor row) a single fit from the k-means start put
#' the anchor on the wrong basis and degenerated in 9 of 24 samples; this
#' procedure gave the true basis in 24 of 24 at a cost of 0--2.4 in NLL.
#' @return The chosen fit, with \code{X.anchor.search} (one row per assignment:
#'   \code{bases}, the column of the free fit each anchor row was made pure on;
#'   \code{nll}; \code{degenerate}; \code{converged}; \code{chosen}) and
#'   \code{nll.free}, the marginal NLL of the free fit.
#' @keywords internal
#' @noRd
.nmfre.partial.anchor.fit <- function(Y, A, Q, C.signed, epsilon, maxit, extra_args, call,
                                      X.update = "exact") {
  ea <- extra_args; ea$Q <- NULL; ea$X.anchor <- NULL; ea$.anchor.direct <- NULL
  ## the search compares fits as they stop; the result is oriented once, at the end
  orient <- if (!is.null(ea[["orient"]])) match.arg(ea[["orient"]], c("none", "geomin")) else "none"
  orient.eps <- if (!is.null(ea$orient.eps)) ea$orient.eps else 0.01
  ident <- .nmfre.ident.args(ea)
  ea$orient <- "none"; ea$ident.search <- FALSE
  Ym <- as.matrix(Y)
  Am <- if (is.null(A)) matrix(1, 1, ncol(Ym)) else as.matrix(A)
  signed <- !(identical(C.signed, FALSE) || identical(C.signed, "nonneg"))
  shift <- if (signed && min(Ym, na.rm = TRUE) < 0) -min(Ym, na.rm = TRUE) else 0
  J <- .nmfre.partial.anchor.rows(extra_args$X.anchor, Ym + shift, Q)
  fit1 <- function(args) do.call(nmfre, c(list(Y = Y, A = A, rank = Q, C.signed = C.signed,
                                              X.update = X.update, epsilon = epsilon, maxit = maxit), args))
  nll <- function(f) .nmfre.marginal.nll(Ym, f$X, f$C, Am, f$sigma2, f$tau2)
  ## (1) the free fit (by default X.init = "mix2": the k-means and SPA starts)
  f0 <- tryCatch(fit1(ea), error = function(e) NULL)
  if (is.null(f0)) stop("nmfre(): the free fit before imposing X.anchor failed.", call. = FALSE)
  ## (2) every ordered choice of distinct bases for the anchor rows
  asg <- as.matrix(expand.grid(rep(list(seq_len(Q)), length(J))))
  asg <- asg[apply(asg, 1, function(a) !anyDuplicated(a)), , drop = FALSE]
  fits <- lapply(seq_len(nrow(asg)), function(i) {
    X0 <- f0$X
    for (k in seq_along(J)) {
      s <- sum(f0$X[J[k], ]); X0[J[k], ] <- 0; X0[J[k], asg[i, k]] <- s
    }
    a <- utils::modifyList(ea, list(X.anchor = J, X.init = X0, C.init = f0$C, U.init = f0$U,
                                    sigma2 = f0$sigma2, tau2 = f0$tau2, .anchor.direct = TRUE))
    tryCatch(fit1(a), error = function(e) NULL)
  })
  ok <- !vapply(fits, is.null, logical(1))
  if (!any(ok)) stop("nmfre(): every fit with X.anchor failed.", call. = FALSE)
  ## (3) the smallest marginal NLL
  nl <- rep(Inf, length(fits)); nl[ok] <- vapply(fits[ok], nll, numeric(1))
  b <- which.min(nl)
  best <- fits[[b]]
  best$X.anchor.search <- data.frame(
    bases = apply(asg, 1, paste, collapse = ","), nll = nl,
    degenerate = vapply(fits, function(f) if (is.null(f)) NA else isTRUE(f$degeneracy$flag), logical(1)),
    converged = vapply(fits, function(f) if (is.null(f)) NA else isTRUE(f$converged), logical(1)),
    chosen = seq_along(fits) == b)
  best$nll.free <- nll(f0)
  best$call <- call
  .nmfre.orient.fit(best, nonneg = !signed, orient = orient, eps = orient.eps, ident = ident)
}

## Evaluate a fit, holding back its messages and warnings (and an error) so the
## caller can pass on those of the fit it keeps only.
.nmfre.quietly <- function(expr) {
  msg <- character(0); wrn <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { wrn <<- c(wrn, conditionMessage(w)); invokeRestart("muffleWarning") })
  err <- inherits(val, "error")
  list(fit = if (err) NULL else val, error = if (err) conditionMessage(val) else NULL,
       messages = msg, warnings = wrn)
}

#' nmfre() from the k-means and SPA starts (X.init = "mix2", the default) (Internal)
#'
#' (1) Fit from the k-means start (\code{seed}) and from the SPA start and keep
#' the smaller marginal NLL (the k-means fit on a tie).  (2) While the kept fit
#' is degenerate (\code{.nmfre.degeneracy}), refit from k-means starts with
#' seeds \code{seed + 1}, \code{seed + 2}, \ldots (at most \code{mix.restarts},
#' default 20) and keep the smallest NLL.  No early stop: in the R2 simulation
#' (2800 samples) stopping after two restarts in a row that reached the kept fit
#' again left 24 samples at a degenerate fit up to 328 worse in 2 x NLL, since a
#' later seed still found the non-degenerate one.  A degenerate fit is a trapped one,
#' with a clearly lower likelihood than a fit from another start as a rule, so
#' the NLL alone sets it aside.  This is the procedure the R2 simulations of the
#' NMF-RE paper used.  Only the kept fit's messages and warnings are passed on.
#' @return The kept fit, with \code{X.init.search} (one row per fit:
#'   \code{start}, \code{seed}, \code{nll}, \code{degenerate}, \code{converged},
#'   \code{chosen}).
#' @keywords internal
#' @noRd
.nmfre.mix2.fit <- function(Y, A, Q, C.signed, X.update, epsilon, maxit, extra_args, call) {
  ea <- extra_args; ea$Q <- NULL; ea$X.init <- NULL; ea$mix.restarts <- NULL
  restarts <- if (!is.null(extra_args$mix.restarts)) as.integer(extra_args$mix.restarts) else 20L
  seed <- if (!is.null(ea$seed)) ea$seed else 1
  ## the fits are compared as they stop; the kept one is oriented once, at the end
  orient <- if (!is.null(ea[["orient"]])) match.arg(ea[["orient"]], c("none", "geomin")) else "none"
  orient.eps <- if (!is.null(ea$orient.eps)) ea$orient.eps else 0.01
  ident <- .nmfre.ident.args(ea)
  ea$orient <- "none"; ea$ident.search <- FALSE
  Ym <- as.matrix(Y)
  Am <- if (is.null(A)) matrix(1, 1, ncol(Ym)) else as.matrix(A)
  signed <- !(identical(C.signed, FALSE) || identical(C.signed, "nonneg"))
  fit1 <- function(init, s) {
    args <- c(list(Y = Y, A = A, rank = Q, C.signed = C.signed, epsilon = epsilon, maxit = maxit,
                   X.init = init), utils::modifyList(ea, list(seed = s)))
    if (!is.null(X.update)) args$X.update <- X.update   # only when the caller gave it
    .nmfre.quietly(do.call(nmfre, args))
  }
  nll <- function(r) {
    if (is.null(r$fit)) return(Inf)
    v <- .nmfre.marginal.nll(Ym, r$fit$X, r$fit$C, Am, r$fit$sigma2, r$fit$tau2)
    if (is.finite(v)) v else Inf
  }
  ## (1) the k-means and SPA starts
  runs <- list(fit1("kmeans", seed), fit1("spa", seed))
  start <- c("kmeans", "spa"); seeds <- c(seed, seed)
  nl <- vapply(runs, nll, numeric(1))
  b <- which.min(nl)
  if (!is.finite(nl[b]))
    stop("nmfre(): the fit failed from both the k-means and the SPA start: ",
         if (!is.null(runs[[1]]$error)) runs[[1]]$error else "no finite likelihood", call. = FALSE)
  ## (2) while the kept fit is degenerate, more k-means starts (no early stop:
  ## a seed after several that returned to the degenerate fit can still escape)
  n.extra <- 0L
  while (isTRUE(runs[[b]]$fit$degeneracy$flag) && n.extra < restarts) {
    n.extra <- n.extra + 1L
    runs[[length(runs) + 1L]] <- fit1("kmeans", seed + n.extra)
    start <- c(start, "kmeans"); seeds <- c(seeds, seed + n.extra)
    nl <- c(nl, nll(runs[[length(runs)]]))
    if (nl[length(nl)] < nl[b]) b <- length(nl)
  }
  best <- runs[[b]]$fit
  for (m in runs[[b]]$messages) message(m, appendLF = FALSE)
  for (w in runs[[b]]$warnings) warning(w, call. = FALSE)
  if (n.extra > 0L && isTRUE(best$degeneracy$flag))
    message(sprintf(paste("nmfre(): the fit is still degenerate after %d more k-means starts;",
                          "its two merging bases are not interpretable, or the rank may be larger",
                          "than the data support (see $degeneracy)."), n.extra))
  best$X.init.search <- data.frame(
    start = start, seed = seeds, nll = nl,
    degenerate = vapply(runs, function(r) if (is.null(r$fit)) NA else isTRUE(r$fit$degeneracy$flag), logical(1)),
    converged = vapply(runs, function(r) if (is.null(r$fit)) NA else isTRUE(r$fit$converged), logical(1)),
    chosen = seq_along(runs) == b)
  best$call <- call
  .nmfre.orient.fit(best, nonneg = !signed, orient = orient, eps = orient.eps, ident = ident)
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
#' ridge/BLUP update for \eqn{U}, an update for \eqn{X} (by default the exact
#' minimization under the column sums; see \code{X.update}), and
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
#'   with a two-sided test (interior null). \code{FALSE} constrains
#'   \eqn{C \ge 0}, updated by an exact non-negative least squares (see
#'   \code{C.nonneg.update} under \code{...}), with a one-sided test (boundary
#'   null). The basis \eqn{X} is always non-negative (see \code{X.update}). A
#'   character value (\code{"signed"} / \code{"nonneg"}) is also accepted for
#'   backward compatibility.
#' @param X.update The \eqn{X}-step.  \code{"exact"} (default): the
#'   \eqn{X}-block is minimized exactly under \eqn{X\ge0} and
#'   \eqn{\mathbf 1^\top X=\mathbf 1^\top} (a convex quadratic program, solved
#'   through the multipliers of the column sums), with no renormalization; with
#'   the exact \eqn{U}- and \eqn{C}-steps every block update is an exact
#'   minimization, the objective decreases monotonically and the fit stops at a
#'   KKT point of the likelihood.  \code{"semi"}: the complete-EM semi-NMF
#'   multiplicative step (Ding, Li and Jordan, 2010) followed by the column
#'   renormalization \eqn{X\to XD^{-1}}, \eqn{(C,U)\to(DC,DU)}, which keeps the
#'   fitted mean but not the likelihood, so the fit can stop short of a
#'   stationary point; the default of nmfkc 1.0.3 and earlier, kept to reproduce
#'   those fits.  On the package's examples \code{"exact"} took as long as
#'   \code{"semi"} or less and reached the same best point from every start.  A
#'   \code{"semi"} fit checks how far it is from a stationary point
#'   (\code{kkt.check} under \code{...}).  \code{"exact"} is not available with
#'   \code{X.L2.ortho} or \code{X.L2.smooth}.
#' @param epsilon Convergence tolerance for relative change in objective (default 1e-5).
#' @param maxit Maximum number of iterations.  Default \code{5000}
#'   (matches \code{\link{nmfkc}} and the other MU functions in the
#'   package).  When the cap is hit without meeting the relative-
#'   tolerance criterion, a \code{"maximum iterations (...) reached..."}
#'   warning is emitted so users notice unconverged fits.
#' @param ... Additional arguments for initialization, variance control, dfU control,
#'   optimization, and inference settings.
#'   \itemize{
#'     \item \code{X.init}: Start of the basis.  \code{"mix2"} (default; also
#'       \code{NULL}): fit from the k-means start and from the SPA start and
#'       keep the smaller marginal NLL; while the kept fit is degenerate (see
#'       \code{degeneracy} in the return value), refit from k-means starts with
#'       seeds \code{seed + 1}, \code{seed + 2}, \ldots (at most
#'       \code{mix.restarts}, default 20) and keep the smallest NLL.  The fits
#'       are returned as \code{X.init.search}; only the kept fit's messages and
#'       warnings are shown.  This takes about twice the time of one fit, and
#'       up to 22 times when every start ends degenerate (e.g. a rank the data
#'       do not support).
#'       A method name (\code{"kmeans"}, \code{"spa"}, \code{"nndsvd"},
#'       \code{"kmeans++"}, \ldots) fits once from that start of
#'       \code{\link{nmfkc}} (\code{"kmeans"} is the single start of nmfkc
#'       1.0.3 and earlier); a basis matrix (P x Q) is used as given.  With
#'       \code{X.anchor} the anchors build the start instead.
#'     \item \code{mix.restarts}: Most k-means restarts after a degenerate fit
#'       with \code{X.init = "mix2"} (default 20; 0 for none).
#'     \item \code{X.anchor}: anchor rows of \eqn{X}, imposed as in
#'       \code{\link{nmfkc}}: \code{NULL} (default), \code{"spa"}, or one row
#'       of \code{Y} per basis.  The zeros are kept for the whole fit, including
#'       through the floor this function applies to \eqn{X}.  There is no
#'       \code{rank <= nrow(A)} limit here: the random effects give the scores
#'       full rank.  The rows used are returned as \code{X.anchor}.
#'       Fewer rows than bases (indices or row names) anchor the first bases
#'       only: row \code{X.anchor[q]} carries basis \code{q} alone and the
#'       other bases are free.  Which basis an anchor row belongs to is then
#'       chosen by likelihood: a free fit (from \code{X.init}, by default
#'       \code{"mix2"}), a refit from
#'       it for every assignment of the anchor rows to bases (\code{rank} refits
#'       for one anchor row), and the refit with the smallest marginal NLL.  The
#'       search is returned as \code{X.anchor.search} and the free fit's
#'       marginal NLL as \code{nll.free}; their difference is the cost of the
#'       anchors.  For \code{rank = 3} one exact anchor row already fixes the
#'       rotation that leaves the likelihood unchanged, so impose only anchors
#'       the design supports.  With \code{C.signed = TRUE} and a \code{Y} with
#'       negative entries, \code{"spa"} chooses the rows on \code{Y - min(Y)},
#'       the same shifted \code{Y} the initialisation sees (as
#'       \code{X.init = "spa"} does); the shift is a device for the start, so
#'       the rows chosen are then only approximately those of the
#'       non-negative case.
#'     \item \code{C.init}: Initial coefficient matrix (Q x K), or \code{NULL}.
#'       When \code{NULL}, \code{\link{nmfkc}} is called internally to generate initial values.
#'     \item \code{U.init}: Initial random effects matrix (Q x N), or \code{NULL} (all zeros).
#'     \item \code{prefix}: Prefix for basis names (default \code{"Basis"}).
#'     \item \code{sigma2}: Initial residual variance.  By default (when
#'       \eqn{\sigma^2} is updated) it starts from the residual covariance of
#'       \eqn{Y} on \eqn{A}, \eqn{\tau^2 XX' + \sigma^2 I}, as in probabilistic
#'       PCA: \eqn{\sigma^2} at the mean of its eigenvalues beyond the
#'       \eqn{Q}-th, \eqn{\tau^2} at U's share (the top \eqn{Q} eigenvalues in
#'       excess of \eqn{\sigma^2}) over \eqn{\mathrm{tr}(X'X)} of the start
#'       basis.  Returned as \code{variance.start}.  nmfkc 1.0.3 and earlier
#'       started at \eqn{\sigma^2 = \tau^2 = 1} (\eqn{\lambda = 1} whatever the
#'       scale of \eqn{Y}), from which data-driven starts could fall into a
#'       degenerate fit; give \code{sigma2 = 1, tau2 = 1} to reproduce them.
#'       With \code{sigma2.update = FALSE} the default stays 1.
#'     \item \code{sigma2.update}: Logical. Update \eqn{\sigma^2} during iterations (default \code{TRUE}).
#'     \item \code{tau2}: Initial random effect variance (default: from the
#'       residual covariance, see \code{sigma2}; 1 with \code{tau2.update = FALSE}).
#'     \item \code{tau2.update}: Logical. Update \eqn{\tau^2} by moment matching (default \code{TRUE}).
#'     \item \code{x.postvar}: Logical. Include the posterior-variance term in
#'       the \eqn{X}-step (default \code{TRUE}; advanced). Applies to every
#'       \eqn{X}-step except the positive-part multiplicative one used with
#'       \code{C.nonneg.update = "mu"} and \code{X.update = "semi"}.
#'     \item \code{kkt.check}: Logical (default \code{TRUE}).  With
#'       \code{X.update = "semi"}, run 50 iterations of \code{X.update = "exact"}
#'       from the fit and record the decrease of \eqn{2\times}NLL as
#'       \code{kkt$improvement}; a message suggests \code{X.update = "exact"}
#'       when it exceeds \code{kkt.tol} (default 0.1).  Skipped with anchor rows
#'       and with the \eqn{X} penalties.
#'     \item \code{C.nonneg.update}: Updates used with \code{C.signed = FALSE}.
#'       \code{"exact"} (default): the \eqn{C}-step is an exact non-negative
#'       least squares (Lawson and Hanson's active-set method) and the
#'       \eqn{X}-step is that of the sign-free path (\code{X.update}), so every block
#'       update is exact or monotone and the fit stops at a KKT point of the
#'       likelihood.  \code{"mu"}: the positive-part multiplicative updates of
#'       nmfkc 1.0.3 and earlier, kept to reproduce those fits; their fixed
#'       points need not be stationary.
#'     \item \code{X.L2.smooth}: Non-negative row-smoothness penalty on the
#'       basis \eqn{X} (default 0), adding
#'       \eqn{(\lambda/2)\,\mathrm{tr}(X^\top L X)} with \eqn{L} the path-graph
#'       Laplacian over rows (squared adjacent-row differences).  Well suited to
#'       ordered rows (e.g.\ longitudinal / spatial bases).  Injected into the
#'       \eqn{X}-step for both \code{C.signed = TRUE/FALSE}; leaves the
#'       random-effect variance updates unchanged.
#'     \item \code{X.L2.ortho}: Non-negative column-orthogonality penalty on
#'       \eqn{X} (default 0), penalizing
#'       \eqn{(\lambda/2)\lVert\mathrm{offdiag}(X^\top X)\rVert^2}.
#'     \item \code{C.L2}: Non-negative ridge penalty on \eqn{\Theta = C}
#'       (default 0), adding \eqn{\lambda\lVert C\rVert^2}.  A Gaussian prior on
#'       the fixed effects; for \code{C.signed = TRUE} the \eqn{C}-step stays a
#'       closed-form Sylvester ridge-LS solve, for \code{C.signed = FALSE} it is
#'       added to the non-negative least squares (or, with
#'       \code{C.nonneg.update = "mu"}, to the multiplicative-update denominator).
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
#'     \item \code{orient}: Representation reported within the identification
#'       region (see the return value): \code{"none"} (default; the point where
#'       the algorithm stopped, as up to 1.0.2) or \code{"geomin"}.
#'     \item \code{orient.eps}: \eqn{\varepsilon} of the geomin criterion (default 0.01).
#'     \item \code{ident.search}: For \eqn{Q\ge4}, search the identification
#'       region numerically when, after the bases pinned by anchors are
#'       removed, between 4 and 7 remain (default \code{TRUE}; see the return
#'       value).
#'     \item \code{ident.dirs}, \code{ident.draws}: Directions followed from the
#'       fit and random draws over the whole family in that search (defaults
#'       \code{300} and \code{20000}).
#'     \item \code{ident.zero.tol}: Entries of \eqn{X} (and of \eqn{C} for
#'       \code{C.signed = FALSE}) at most this value count as zeros when
#'       looking for anchors (default \code{0}: exact zeros only, which the
#'       exact updates return for active constraints).  A positive value makes
#'       the pinning a numerical diagnostic, not an exact result.
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
#'   \strong{Degeneracy diagnostic}
#'   \describe{
#'     \item{\code{degeneracy}}{A list for detecting a degenerate fit, in which
#'       two columns of \eqn{X} approach each other while the corresponding rows
#'       of \eqn{C} grow with opposite signs (with \eqn{Q \ge 2} and a signed
#'       \eqn{C}: two bases suffice, so it can happen at \eqn{Q = 2} as well,
#'       e.g. when the rank exceeds what the data support; the fit does not
#'       settle and those rows of \eqn{C} are not interpretable).  In simulations such fits usually had a clearly lower
#'       marginal likelihood than a fit from another start, so the default
#'       \code{X.init = "mix2"} refits from other starts while \code{flag} is
#'       \code{TRUE} and keeps the fit with the larger likelihood:
#'       \code{max.cosine} (largest cosine similarity between two columns of
#'       \eqn{X}) and its \code{pair}; \code{cancellation}
#'       (\eqn{\max_k \sum_q |c_{qk}| / \|X c_{\cdot k}\|_1}, large when a composite
#'       effect is a difference of large opposite terms); \code{kappa.X} (condition
#'       number of \eqn{X}); and \code{flag} (\code{max.cosine >= 0.9} or
#'       \code{cancellation >= 5}).  The thresholds were set on two simulation
#'       designs; read the measures as well.  \code{summary()} prints a line when
#'       \code{flag} is \code{TRUE}.}
#'     \item{\code{X.init.search}}{With \code{X.init = "mix2"}: one row per fit
#'       (\code{start}, \code{seed}, marginal \code{nll}, \code{degenerate},
#'       \code{converged}, \code{chosen}).  \code{summary()} prints the start
#'       kept and how far the others stopped from it.}
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
#'   \strong{Identification region and orientation}
#'
#'   Under the working model the fits giving the same distribution are
#'   \eqn{(XO, O^\top\Theta)} with \eqn{O} orthogonal, \eqn{O\mathbf 1=\mathbf 1}
#'   and \eqn{XO\ge0} (and \eqn{O^\top\Theta\ge0} for \code{C.signed = FALSE}):
#'   the identification region.  \eqn{X\Theta} is constant on it.  For
#'   \eqn{Q=3} the \eqn{O} form a circle, and the arc through the fit is the
#'   local identification region.  The ends of the arc are found numerically
#'   (a grid of 7200 angles refined by root finding, constraints met to
#'   \eqn{10^{-10}}; an arc shorter than \eqn{10^{-7}} is reported as a point);
#'   the range of each entry of \eqn{\Theta} and of \eqn{X} over the arc is
#'   then computed in closed form.  A width of 0 means the fit cannot be moved
#'   within the region, to this resolution (a row of \eqn{X} on one basis only
#'   pins it).  For \eqn{Q\le2} the region is a point up to relabelling.  For
#'   \eqn{Q\ge4} the bases pinned by anchors are removed one at a time: if a row
#'   of \eqn{X} is non-zero, outside the bases already removed, on one basis
#'   only (and, for \code{C.signed = FALSE}, if a column of \eqn{C} is), every
#'   admissible \eqn{O} keeps that basis, because a non-negative row of \eqn{O}
#'   with sum 1 and length 1 is a unit vector.  Zeros are exact zeros unless
#'   \code{ident.zero.tol} is set.  If at most two bases remain the region is a
#'   point up to relabelling; if three remain their local region (an arc, a
#'   point, or the whole circle) is computed as for \eqn{Q=3}; if 4 to 7 remain
#'   it is searched numerically (\code{ident.search}): the ranges then come from
#'   points reached along paths from the fit that are feasible at each step,
#'   and are approximate, not certified bounds (finding no other
#'   representation is evidence, not a proof).  Feasible transformations found
#'   elsewhere in the family show that the fit is not unique globally; they are
#'   counted (\code{identification$search$global.nonperm}) but not used in the
#'   ranges, since they need not belong to the local region.  The geomin
#'   orientation and the profile intervals' use of the region are for
#'   \eqn{Q=3} only.  By default (\code{orient = "none"}) the
#'   representation returned is the point where the algorithm stopped.  With
#'   \code{orient = "geomin"} it is the member of the local region that minimizes
#'   the geomin criterion on the row shares of \eqn{X} (rows of zeros left out;
#'   as rotated loadings in exploratory factor analysis); nothing moves when the
#'   region is a point.  The search is local: it finds the minimizer over the
#'   whole identification region only when the local region contains it.
#'   The likelihood and \eqn{X\Theta} do not change.
#'   \describe{
#'     \item{\code{C.ident.lower}, \code{C.ident.upper}, \code{C.ident.width}}{Range
#'       of each entry of \eqn{C} over the local identification region and its
#'       width (Q x K matrices; absent when the region is not computed).}
#'     \item{\code{X.ident.lower}, \code{X.ident.upper}, \code{X.ident.width}}{The
#'       same for \eqn{X} (P x Q matrices).}
#'     \item{\code{identification}}{List: \code{width} (largest entry of
#'       \code{C.ident.width}; \code{NA} when not computed), \code{arc} and
#'       \code{full.circle} of the region, \code{method} (\code{"exact"},
#'       \code{"anchors"}, \code{"anchors+exact"}, \code{"search"},
#'       \code{"anchors+search"} or \code{"none"}; \code{"exact"} means the arc
#'       computed as above, in closed form over an arc whose ends are found
#'       numerically, not a certified computation), \code{pinned} (bases pinned
#'       by anchors, \eqn{Q\ge4}), \code{Q.rest} (bases left after them),
#'       \code{search} (counts of the numerical search), and \code{orient}, the orientation
#'       applied (\code{rule}, \code{angle}, \code{O}, \code{moved}: largest
#'       change of \eqn{C}, the geomin \code{criterion} before and after,
#'       \code{eps}, \code{note}).  \code{\link{nmfre.inference}} adds its
#'       identifiability diagnostic to this list.}
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
nmfre <- function(Y, A = NULL, rank = 2, C.signed = TRUE, X.update = c("exact", "semi"),
                  epsilon = 1e-5, maxit = 5000, ...) {
  cl <- base::match.call()
  X.update.given <- !base::missing(X.update)
  X.update <- base::match.arg(X.update)
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
  ## The X-step is set by X.update (exact by default) for a sign-free C and for
  ## C >= 0 with the exact updates (the default); with X.update = "semi" it is
  ## the complete-EM semi-NMF step (Ding-Li-Jordan 2010).  Only the legacy
  ## C.nonneg.update = "mu" keeps the positive-part multiplicative step.
  ## Complete-EM X-step: add the posterior-variance term tr(X S X'),
  ## S = N*sigma2*(X'X+lambda I)^{-1}; FALSE recovers the conditional X-step.
  x.postvar <- if (!is.null(extra_args$x.postvar)) isTRUE(extra_args$x.postvar) else TRUE
  ## With C >= 0, "exact" (default) solves the C-step as an exact non-negative
  ## least squares and uses the X-step of the sign-free path (X.update), so each block update is exact or monotone and the ECM stops at a KKT
  ## point of the likelihood.  "mu" is the positive-part multiplicative update
  ## of nmfkc <= 1.0.3 (kept to reproduce those fits): it acts on [Y - XU]+ and
  ## [CA + U]+, so its fixed points need not be stationary.
  C.nonneg.update <- if (!is.null(extra_args$C.nonneg.update))
    base::match.arg(extra_args$C.nonneg.update, c("exact", "mu")) else "exact"
  nn.exact <- C.mode == "nonneg" && C.nonneg.update == "exact"
  ## X-step (argument X.update).  "semi": the (complete-EM) semi-NMF
  ## multiplicative step followed by the column renormalization X -> X D^-1,
  ## (C, U) -> (D C, D U), which keeps the fitted mean but not the likelihood,
  ## so the ECM may stop short of a stationary point.  "exact" (default): the X-block is
  ## minimized exactly under 1'X = 1' (.nmfre.simplex.X()), no renormalization;
  ## with the exact U- and C-steps every block of the inner loop is then
  ## minimized exactly.
  X.exact <- X.update == "exact"

  ## Fewer anchor rows than bases (e.g. X.anchor = "x7" for rank = 3): which
  ## basis an anchor row belongs to is part of the estimate.  Fit freely first,
  ## refit from that fit with each assignment of the anchor rows to bases, and
  ## keep the smallest marginal NLL (.nmfre.partial.anchor.fit()).
  X.anchor.arg <- extra_args$X.anchor
  if (!is.null(X.anchor.arg) && !identical(X.anchor.arg, "spa") &&
      length(X.anchor.arg) < Q && !isTRUE(extra_args$.anchor.direct))
    return(.nmfre.partial.anchor.fit(Y, A, Q = Q, C.signed = C.signed, epsilon = epsilon,
                                     maxit = maxit, extra_args = extra_args, call = cl,
                                     X.update = X.update))
  ## X.init = "mix2" (the default): the better of the k-means and SPA starts,
  ## with more k-means starts while the kept fit is degenerate
  ## (.nmfre.mix2.fit()).  Anchor rows build their own start, so it does not
  ## apply to them.
  if (is.null(X.init) || identical(X.init, "mix2")) {
    if (is.null(X.anchor.arg))
      return(.nmfre.mix2.fit(Y, A, Q = Q, C.signed = C.signed,
                             X.update = if (X.update.given) X.update else NULL,
                             epsilon = epsilon, maxit = maxit, extra_args = extra_args, call = cl))
    X.init <- NULL
  }

  ## MAP penalties (default off).  These act on X (basis) and C (= Theta) as
  ## Gaussian priors and are orthogonal to the random-effect machinery
  ## (U, lambda, sigma2, tau2), so the variance M-steps are unchanged.
  ##   X.L2.ortho : (l/2)||offdiag(X'X)||^2   -- column orthogonality of X
  ##   X.L2.smooth: (l/2) tr(X' L X)          -- path-graph row smoothness of X
  ##   C.L2       : l * ||C||^2               -- ridge on Theta (keeps the signed
  ##                C-step in closed form: a Sylvester ridge-LS solve)
  X.L2.ortho  <- if (!is.null(extra_args$X.L2.ortho))  extra_args$X.L2.ortho  else 0
  X.L2.smooth <- if (!is.null(extra_args$X.L2.smooth)) extra_args$X.L2.smooth else 0
  if (X.exact && (X.L2.ortho > 0 || X.L2.smooth > 0)) {
    if (X.update.given)
      stop("X.update = \"exact\" does not support X.L2.ortho or X.L2.smooth; use X.update = \"semi\".")
    ## the default cannot take the X penalties: use the semi-NMF step, as before
    message("nmfre(): X.L2.ortho / X.L2.smooth need X.update = \"semi\"; using it.")
    X.update <- "semi"; X.exact <- FALSE
  }
  C.L2        <- if (!is.null(extra_args$C.L2))        extra_args$C.L2        else 0

  # variance handling
  sigma2        <- if (!is.null(extra_args[["sigma2"]])) extra_args[["sigma2"]] else 1
  sigma2.update <- if (!is.null(extra_args$sigma2.update)) extra_args$sigma2.update else TRUE
  tau2          <- if (!is.null(extra_args[["tau2"]])) extra_args[["tau2"]] else 1
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
  ## Keep the self-seeding of this fit out of the caller's random stream.
  .rng <- .nmfkc.rng.save(seed)
  on.exit(.nmfkc.rng.restore(.rng), add = TRUE)

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
  ## nmfkc() needs Y >= 0, but the sign-free fit does not: its semi-NMF X-step
  ## splits the signs.  With C.signed = TRUE a Y with negative entries is
  ## shifted to start the basis only, and C is then started by the exact
  ## sign-free least squares on the original Y (below).  Before 1.0.3 such
  ## a Y stopped in the initialisation ("The matrix Y should be non-negative").
  Y.shift <- if (C.mode == "signed" && min(Y, na.rm = TRUE) < 0) -min(Y, na.rm = TRUE) else 0

  ## X.anchor: anchor rows of X, kept as exact zeros for the whole fit.  The
  ## start built from them goes to the nmfkc() initialisation as a matrix,
  ## whose multiplicative updates keep the zeros, and the floor pmax(X, .eps)
  ## below puts them back wherever it would lift them.  No rank check as in
  ## nmfkc(): with the random effects U the scores have full rank anyway.
  ## The rows are chosen, and the start is built, on the same shifted Y as the
  ## nmfkc() initialisation (as X.init = "spa" does); before 1.0.3
  ## X.anchor = "spa" stopped on a Y with negative entries.  Fewer than Q
  ## rows anchor the first bases only (see .nmfre.partial.anchor.rows()).
  X.anchor <- extra_args$X.anchor
  anchor.rows <- NULL; anchor.mask <- NULL
  anchor.partial <- !is.null(X.anchor) && !identical(X.anchor, "spa") &&
    length(X.anchor) < Q
  if (anchor.partial) {
    anchor.rows <- .nmfre.partial.anchor.rows(X.anchor, Y + Y.shift, Q)
    anchor.mask <- matrix(FALSE, P, Q)
    for (q in seq_along(anchor.rows)) anchor.mask[anchor.rows[q], -q] <- TRUE
  } else if (!is.null(X.anchor)) {
    an <- .nmfkc_anchor_setup(X.anchor, Y + Y.shift, Q, D_A = NULL, fun = "nmfre",
                              seed = seed)
    anchor.rows <- an$rows; anchor.mask <- an$mask
    if (!is.null(X.init) && !identical(X.init, "spa"))
      message("nmfre(): X.anchor builds its own starting basis, so X.init was not used.")
    X.init <- NULL
  }
  anchor.full <- !is.null(anchor.rows) && !anchor.partial

  C.from.init <- is.null(C.init)
  if (is.null(X.init) || is.null(C.init)) {
    init_args <- list(Y + Y.shift, A, Q = Q, epsilon = epsilon, seed = seed,
                      nstart = nstart, print.trace = print.trace,
                      print.dims = FALSE)
    if (anchor.full) {
      init_args$X.init <- an$X                    # the anchored start
    } else if (is.character(X.init)) {
      init_args$X.init <- X.init                  # forward the named method
    } else if (is.matrix(X.init) || is.data.frame(X.init)) {
      init_args$X.init <- as.matrix(X.init)       # fix the supplied basis
      init_args$X.restriction <- "fixed"
    }
    res0 <- do.call(nmfkc, init_args)
    if (is.null(X.init) || is.character(X.init)) X.init <- res0$X
    if (is.null(C.init)) C.init <- res0$C
    if (anchor.full) {
      ## nmfkc() reorders its bases: put basis q back on anchor row q
      k <- apply(X.init[anchor.rows, , drop = FALSE] > 0, 1, which)
      X.init <- X.init[, k, drop = FALSE]
      if (C.from.init) C.init <- C.init[k, , drop = FALSE]
    }
  }
  if (anchor.partial) {
    ## the free start gives every row some share of each basis: anchor row q
    ## goes to the basis with its largest share among those not yet taken, and
    ## the anchored bases are moved to the front (basis q on anchor row q)
    X0 <- as.matrix(X.init)
    sh <- X0 / pmax(rowSums(X0), .eps)
    k <- integer(0)
    for (j in anchor.rows) {
      cand <- setdiff(seq_len(Q), k)
      k <- c(k, cand[which.max(sh[j, cand])])
    }
    ord <- c(k, setdiff(seq_len(Q), k))
    X.init <- X0[, ord, drop = FALSE]
    C.init <- as.matrix(C.init)[ord, , drop = FALSE]
    if (!is.null(U.init)) U.init <- as.matrix(U.init)[ord, , drop = FALSE]
  }
  if (Y.shift > 0 && C.from.init) {
    ## C from the shifted Y would carry the shift; restart it by the sign-free
    ## least squares C = (X'X)^{-1} X' Y A' (A A')^{-1} on the original Y.
    X0 <- as.matrix(X.init)
    C.init <- solve(crossprod(X0) + diag(1e-10, Q), crossprod(X0, Y)) %*%
      t(A) %*% solve(tcrossprod(A) + diag(1e-10, K))
    dimnames(C.init) <- list(colnames(X0), rownames(A))
  }
  ## pmax(X, .eps), except at the anchors' zeros
  .floorX <- function(X) {
    X <- pmax(X, .eps)
    if (!is.null(anchor.mask)) X[anchor.mask] <- 0
    X
  }

  X <- X.init
  stopifnot(nrow(X) == P, ncol(X) == Q)
  if (!is.null(anchor.mask)) X[anchor.mask] <- 0   # a free start (partial anchors) has them non-zero

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
  ## Clip C to +eps only in nonneg mode: under C.signed = TRUE a warm-start
  ## C.init may legitimately carry negative entries (e.g. full-refit
  ## bootstrap), and clipping here would destroy their sign at every refit.
  ## This mirrors the in-loop rule `if (C.mode == "nonneg") pmax(...)`.
  C0 <- if (C.mode == "nonneg") pmax(C_mat, .eps) else C_mat
  normed <- .nmfre.normalize.X(.floorX(X), C0, U)
  X <- normed$X
  C_mat <- normed$C
  U <- normed$U

  ## Variance start, when a variance is updated and not given: from the
  ## eigenvalues of the residual covariance of Y on A (.nmfre.variance.start()).
  ## sigma2 = tau2 = 1 (nmfkc 1.0.3 and earlier) put lambda at 1 whatever the
  ## scale of Y, and from there every data-driven basis start could fall into a
  ## degenerate fit.
  var.start <- list(rule = "given", sigma2 = sigma2, tau2 = tau2)
  s2.auto <- sigma2.update && is.null(extra_args[["sigma2"]])
  t2.auto <- tau2.update && is.null(extra_args[["tau2"]])
  if (s2.auto || t2.auto) {
    vs <- .nmfre.variance.start(Y, A, X)
    if (!is.null(vs)) {
      if (s2.auto) sigma2 <- vs[["sigma2"]]
      if (t2.auto) tau2 <- vs[["tau2"]]
      var.start <- list(rule = "residual eigenvalues", sigma2 = sigma2, tau2 = tau2)
    } else var.start$rule <- "1 (the residual has no eigenvalues beyond rank)"
  }

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

  ## X-penalty MU contributions (num / den) and the penalty value added to the
  ## inner objective.  ortho -> denominator only; smooth -> +W X (num) / +D X (den).
  .xpen <- function(X) {
    num <- 0; den <- 0
    if (X.L2.ortho > 0) { G <- crossprod(X); diag(G) <- 0; den <- den + X.L2.ortho * (X %*% G) }
    if (X.L2.smooth > 0 && nrow(X) >= 2) {
      Pr <- nrow(X); WX <- X * 0
      WX[-Pr, ] <- WX[-Pr, ] + X[-1, , drop = FALSE]
      WX[-1, ]  <- WX[-1, ]  + X[-Pr, , drop = FALSE]
      num <- num + X.L2.smooth * WX
      den <- den + X.L2.smooth * (c(1, rep(2, Pr - 2), 1) * X)
    }
    list(num = num, den = den)
  }
  .xpen.val <- function(X) {
    v <- 0
    if (X.L2.ortho > 0) { G <- crossprod(X); diag(G) <- 0; v <- v + (X.L2.ortho / 2) * sum(G^2) }
    if (X.L2.smooth > 0 && nrow(X) >= 2)
      v <- v + (X.L2.smooth / 2) * sum((X[-1, , drop = FALSE] -
                                        X[-nrow(X), , drop = FALSE])^2)
    v
  }

  ## ---- loop-invariant precomputations (A and Y never change in the ECM) ----
  tA <- t(A)                          # N x K, reused in every C-step
  if (C.mode == "signed") {
    if (C.L2 > 0) {
      ea  <- eigen(tcrossprod(A), symmetric = TRUE)   # eigenbasis of AA' (invariant)
      tea <- t(ea$vectors)
    } else {
      AAt_t <- tcrossprod(A) + diag(1e-10, K)
    }
  } else {
    Y_tilde <- pmax(Y, 0)             # nonneg-path stabilized Y (invariant)
    AAt_mm  <- A %*% t(A)             # keep the gemm form (one-arg tcrossprod/dsyrk is not bitwise-identical)
  }

  outer_done <- 0L
  mu_X <- NULL                        # column-sum multipliers of the exact X-step (warm start)
  for (outer in 1:outer.maxit) {
    outer_done <- outer

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
      tcholM <- t(cholM)                          # hoist transpose out of the loop
      RHS <- matrix(0, Q, N)                      # per-column crossprod keeps the exact FP ops
      for (n in 1:N) RHS[, n] <- crossprod(X, R0[, n, drop = FALSE])
      U[] <- backsolve(cholM, forwardsolve(tcholM, RHS))  # batched column-independent triangular solves
      ## NOTE: do NOT row-center U here. The (U, Theta) indeterminacy is
      ## U -> U + Delta A, Theta -> Theta - Delta, so the identification
      ## condition is U A' = 0, not U 1 = 0. The alternating fixed point
      ## satisfies U A' = 0 automatically; imposing row-centering on top
      ## breaks it whenever 1' is not in rowspace(A).

      # (2) X-step: X >= 0 with sign-free score matrix B = CA + U
      B_sem <- CA + U
      xp <- .xpen(X)
      if (X.exact) {
        ## exact minimization of ||Y - X B||^2 + tr(X S X') over X >= 0, 1'X = 1'
        ## (anchor zeros kept); C and U are not rescaled
        Gx <- B_sem %*% t(B_sem)
        if (isTRUE(x.postvar)) Gx <- Gx + N * clip_val(sigma2, sigma2.min, sigma2.max) * chol2inv(cholM)
        sx <- .nmfre.simplex.X(Gx, Y %*% t(B_sem), mu = mu_X, mask = anchor.mask)
        X <- sx$X; mu_X <- sx$mu
      } else if (C.mode == "signed" || nn.exact) {
        ## complete-EM X-step: add posterior-variance S = N*sigma2*(X'X+lambda I)^{-1}
        S_pv <- if (isTRUE(x.postvar))
          N * clip_val(sigma2, sigma2.min, sigma2.max) * chol2inv(cholM) else NULL
        X <- .nmfre.seminmf.X(X, Y, B_sem, S = S_pv,
                              pen.num = xp$num, pen.den = xp$den, eps = .eps)
      } else {
        B_pos <- pmax(B_sem, 0)       # Y_tilde hoisted above the outer loop
        numX <- tcrossprod(Y_tilde, B_pos) + xp$num
        denX <- X %*% (B_pos %*% t(B_pos)) + xp$den
        X <- X * .nmfre.safe.div(numX, denX, eps = .eps)
      }
      if (!X.exact) {
        X <- .floorX(X)
        normed <- .nmfre.normalize.X(X, C_mat, U)
        X     <- .floorX(normed$X)
        C_mat <- if (C.mode == "nonneg") pmax(normed$C, if (nn.exact) 0 else .eps) else normed$C
        U     <- normed$U
      }

      # (3) C-step
      if (C.mode == "signed") {
        ## sign-free least squares.  With C.L2 = 0 this is the exact separable
        ## solution C = (X'X)^{-1} X'(Y - X U) A'(A A')^{-1}.  With C.L2 > 0 the
        ## normal equations become the Sylvester ridge system
        ##   (X'X) C (A A') + C.L2 * C = X'(Y - X U) A',
        ## solved in closed form via the eigenbases of X'X and A A'.
        Y_star <- Y - X %*% U
        rhs    <- crossprod(X, Y_star) %*% tA
        if (C.L2 > 0) {
          ## ea / tea (eigen of AA') hoisted above the outer loop (A invariant)
          ex <- eigen(crossprod(X),  symmetric = TRUE)
          Mt <- crossprod(ex$vectors, rhs) %*% ea$vectors
          Z  <- Mt / (outer(ex$values, ea$values) + C.L2)
          C_mat <- ex$vectors %*% Z %*% tea
        } else {
          ## AAt_t hoisted above the outer loop (A invariant)
          XtX_t <- crossprod(X) + diag(1e-10, Q)
          C_mat <- solve(XtX_t, rhs)
          C_mat <- t(solve(AAt_t, t(C_mat)))
        }
      } else if (nn.exact) {
        ## exact non-negative least squares on vec(C):
        ##   min_{C >= 0} ||Y - X U - X C A||^2 + C.L2 ||C||^2,
        ## H = (A A') (x) (X'X) + C.L2 I,  g = vec(X'(Y - X U) A')
        Hn <- kronecker(AAt_mm, crossprod(X))
        if (C.L2 > 0) Hn <- Hn + C.L2 * diag(Q * K)
        gn <- as.vector(crossprod(X, Y - X %*% U) %*% tA)
        C_mat[] <- .nmfre.nnls(Hn, gn)
      } else {
        ## non-negative MU (Ding 2006), positive-part stabilization
        ## (Y_tilde, tA, AAt_mm hoisted above the outer loop; Y and A invariant)
        Y_star <- Y_tilde - X %*% U; Y_star_pos <- pmax(Y_star, 0)
        numC <- crossprod(X, Y_star_pos) %*% tA
        denC <- (crossprod(X, X) %*% C_mat) %*% AAt_mm
        if (C.L2 > 0) denC <- denC + C.L2 * C_mat
        C_mat <- C_mat * .nmfre.safe.div(numC, denC, eps = .eps)
        C_mat <- pmax(C_mat, .eps)
      }

      # (4) objective at the FIXED lambda (what the inner loop minimizes)
      CA  <- C_mat %*% A
      R   <- Y - X %*% (CA + U)
      obj <- sum(R^2) + lambda * sum(U^2) +
             .xpen.val(X) + if (C.L2 > 0) C.L2 * sum(C_mat^2) else 0
      if (total_iter <= maxit) {
        obj_trace[total_iter] <- obj
        rss_trace[total_iter] <- sum(R^2)
        ## marginal NLL (U integrated out): the ECM-monotone diagnostic
        nll_trace[total_iter] <- .nmfre.marginal.nll(Y, X, C_mat, A, sigma2, tau2, CA = CA)
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

  ## The outer ECM loop running to its own cap is a real stopping reason, not an
  ## unknown one.  stop_reason starts at "maxit", and if the outer loop finishes
  ## without any break firing it is still "maxit" while total_iter < maxit --
  ## which is what used to be relabelled "unknown_break".  Name it: this is the
  ## only path that reached that branch (verified -- the reported iteration
  ## count tracks outer.maxit exactly: 50 -> 793, 200 -> 1093, 1000 -> 2693).
  if (identical(stop_reason, "maxit") && outer_done >= outer.maxit &&
      iter_done < maxit) {
    stop_reason <- "outer_maxit"
    converged <- FALSE
  } else if (iter_done < maxit && identical(stop_reason, "maxit")) {
    ## Kept as a genuine catch-all; no known path leads here.
    stop_reason <- "unknown_break"
    converged <- FALSE
  }

  ## Warn when the MU loop exhausted maxit without meeting the
  ## relative-tolerance criterion (matches nmfkc() / nmf.ffb()
  ## convention).  The non-finite-objective break path issues its
  ## own diagnostic via stop_reason; here we only fire on maxit.
  if (identical(stop_reason, "maxit") && iter_done == maxit)
    warning(paste0("maximum iterations (", maxit, ") reached..."))

  # ---- reorder basis ----
  ## (not when the anchor rows were given: basis q is the one anchored at the
  ## q-th row, and reordering would break that)
  anchor.keep.order <- !is.null(anchor.rows) && !identical(X.anchor, "spa")
  if (ncol(X) > 1 && !anchor.keep.order) {
    w_ord <- matrix((1:P) / P, nrow = 1)
    score <- as.numeric(w_ord %*% X)
    index <- order(score)
    X <- X[, index, drop = FALSE]
    C_mat <- C_mat[index, , drop = FALSE]
    U <- U[index, , drop = FALSE]
    if (!is.null(anchor.rows)) anchor.rows <- anchor.rows[index]
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

    ## Header for print.nmf(), matching the other fitters.
    call = cl,
    dims = if (base::is.null(A))
             base::sprintf("Y(%d,%d)~X(%d,%d)B(%d,%d)",
                           base::nrow(Y), base::ncol(Y), base::nrow(X), Q,
                           Q, base::ncol(Y))
           else base::sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)",
                              base::nrow(Y), base::ncol(Y), base::nrow(X), Q,
                              Q, base::nrow(A), base::nrow(A), base::ncol(A)),

    # convergence diagnostics
    converged = converged,
    stop.reason = stop_reason,
    iter = iter_done,
    ## The outer ECM loop has its own cap, and hitting it is the commonest
    ## reason a run stops short of `maxit`.  Without these two the stop reason
    ## could not be told apart from the inner loop's, which is how the
    ## "unknown_break" label came about.
    outer.iter = outer_done,
    outer.maxit = outer.maxit,
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
    C.signed = (C.mode == "signed"),
    C.nonneg.update = if (C.mode == "nonneg") C.nonneg.update else NULL,
    X.update = X.update,

    ## two bases merging while C diverges (see .nmfre.degeneracy)
    degeneracy = .nmfre.degeneracy(X, C_mat),
    ## where the variances started (see .nmfre.variance.start)
    variance.start = var.start
  )

  ## Added only when anchors were used, so every other fit keeps its fields.
  if (!is.null(anchor.rows))
    out$X.anchor <- stats::setNames(anchor.rows,
      (if (!is.null(colnames(out$X))) colnames(out$X) else
         paste0("Basis", seq_len(Q)))[seq_along(anchor.rows)])
  ## Has the "semi" fit stopped short of a stationary point?  Run a few exact
  ## iterations from it: they reach most of what "exact" would gain (85-98% in
  ## 50 iterations on the examples) at a small cost, and say so if it matters.
  kkt.check <- if (!is.null(extra_args$kkt.check)) isTRUE(extra_args$kkt.check) else TRUE
  kkt.tol   <- if (!is.null(extra_args$kkt.tol)) extra_args$kkt.tol else 0.1
  if (X.update == "semi" && kkt.check && !isTRUE(extra_args$.kkt.inner) && is.null(anchor.rows) &&
      X.L2.ortho == 0 && X.L2.smooth == 0) {
    ea <- extra_args
    ea[c("X.init", "C.init", "U.init", "sigma2", "tau2", "kkt.check", "kkt.tol", ".kkt.inner", "orient",
         "orient.eps", "print.trace", "verbose", "seed", "nstart", "Q", "ident.search")] <- NULL
    fe <- tryCatch(suppressMessages(suppressWarnings(do.call(nmfre, c(list(
      Y = Y, A = A, rank = Q, C.signed = C.signed, X.update = "exact", epsilon = epsilon, maxit = 50L,
      X.init = out$X, C.init = out$C, U.init = out$U, sigma2 = out$sigma2, tau2 = out$tau2,
      verbose = FALSE, print.trace = FALSE, .kkt.inner = TRUE, ident.search = FALSE), ea)))), error = function(e) NULL)
    if (!is.null(fe)) {
      imp <- 2 * (.nmfre.marginal.nll(Y, out$X, out$C, A, out$sigma2, out$tau2) -
                  .nmfre.marginal.nll(Y, fe$X, fe$C, A, fe$sigma2, fe$tau2))
      out$kkt <- list(improvement = imp, iterations = 50L, tol = kkt.tol)
      if (is.finite(imp) && imp > kkt.tol)
        message(sprintf(paste0("nmfre(): the fit stopped short of a stationary point of the likelihood: ",
                               "50 iterations of X.update = \"exact\" from it lower 2*NLL by %.3g.\n",
                               "  For results to report, refit with X.update = \"exact\" (it can start from ",
                               "this fit: X.init, C.init, U.init, sigma2, tau2)."), imp))
    }
  }
  ## identification region; the representation within it (as stopped by default)
  orient <- if (!is.null(extra_args[["orient"]])) base::match.arg(extra_args[["orient"]], c("none", "geomin")) else "none"
  out <- .nmfre.orient.fit(out, nonneg = (C.mode == "nonneg"), orient = orient,
                           eps = if (!is.null(extra_args$orient.eps)) extra_args$orient.eps else 0.01,
                           ident = .nmfre.ident.args(extra_args))
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
#' @param ci.show Logical. If \code{TRUE}, show confidence interval columns
#'   (default \code{FALSE}). Named object-first to match the package style
#'   (\code{C.signed}, \code{X.init}, ...). Legacy \code{show_ci} is accepted
#'   via \code{...}.
#' @param max.coef Largest number of coefficient rows to print (default 20).
#'   The table has one row per (basis, covariate) pair, so it grows as
#'   \eqn{Q\times K}; beyond the cap the significant rows are preferred and the
#'   omission is reported.
#' @param by Grouping order of the coefficient rows: \code{"covariate"}
#'   (default; list all bases for each covariate) or \code{"basis"} (list all
#'   covariates for each basis).
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
summary.nmfre <- function(object, ci.show = FALSE, max.coef = 20,
                          by = c("covariate", "basis"), ...) {
  ## Every other summary.* in the package builds a summary.* object that a
  ## paired print method renders.  This one printed as a side effect and
  ## returned the fit itself, so `s <- summary(fit)` printed unexpectedly and
  ## `print(s)` showed the fit header instead of the summary.  Worse, once
  ## nmfre.inference() started tagging its result "nmf.inference",
  ## print.nmf.inference's `print(summary(x))` re-entered this function
  ## forever.  Tagging the returned object breaks the cycle and gives the
  ## class the same shape as its siblings.
  .summary.nmfre.body(object, ci.show = ci.show, max.coef = max.coef,
                      by = match.arg(by), print = FALSE)
}

#' @rdname summary.nmfre
#' @param x An object of class \code{"summary.nmfre"}.
#' @export
print.summary.nmfre <- function(x, ...) {
  .summary.nmfre.body(x$object, ci.show = x$ci.show, max.coef = x$max.coef,
                      by = x$by, print = TRUE)
  invisible(x)
}

## One block for summary(): the identification region nmfre() computed.
.nmfre.print.ident <- function(x, Q) {
  idn <- x$identification
  if (!is.list(idn) || (is.null(idn$width) && is.null(idn$orient))) return(invisible(NULL))
  w <- idn$width; o <- idn$orient
  pad <- "                      "
  srch <- !is.null(idn$method) && grepl("search", idn$method)
  np <- length(idn$pinned)
  pin <- if (np == 1L) sprintf("basis %d pinned by an anchor", idn$pinned)
         else if (np > 1L) sprintf("bases %s pinned by anchors", paste(idn$pinned, collapse = ", ")) else ""
  if (Q <= 2L) {
    cat("Identification region: a point up to relabelling of the bases (Q <= 2)\n")
  } else if (is.null(w) || !is.finite(w)) {
    cat(sprintf("Identification region: not computed (%s%d bases left; at most 7 are searched)\n",
                if (np) paste0(pin, "; ") else "", if (is.null(idn$Q.rest)) Q else idn$Q.rest))
  } else if (w <= 1e-8 && identical(idn$method, "anchors")) {
    cat(sprintf("Identification region: a point up to relabelling (%s)\n", pin))
  } else if (w <= 1e-8 && srch) {
    cat(sprintf("Identification region: no other representation found near the fit by a numerical search%s\n",
                if (np) paste0(" (", pin, ")") else ""))
  } else if (w <= 1e-8) {
    cat("Identification region: a point (width 0); Theta is locally unique (no other representation nearby)\n")
  } else {
    if (srch)
      cat(sprintf("Identification region: width of Theta about %.3g (numerical search, not a certified bound%s)\n", w,
                  if (np) paste0("; ", pin) else ""))
    else
      cat(sprintf("Identification region: width of Theta up to %.3g (arc %.1f deg%s)%s\n", w,
                  diff(idn$arc) * 180 / pi, if (isTRUE(idn$full.circle)) ", full circle" else "",
                  if (np) paste0(" (", pin, ")") else ""))
    if (!is.null(o))
      cat(pad, " reported: ",
          if (identical(o$rule, "geomin")) sprintf("the geomin member (Theta moved by %.3g)", o$moved)
          else "the fit as it stopped (orient = \"none\")", "\n", sep = "")
    if (!is.null(x$C.ident.lower)) {
      both <- which(x$C.ident.lower < 0 & x$C.ident.upper > 0)
      if (length(both)) {
        rn <- if (!is.null(rownames(x$C))) rownames(x$C) else paste0("Basis", seq_len(nrow(x$C)))
        cn <- if (!is.null(colnames(x$C))) colnames(x$C) else paste0("Cov", seq_len(ncol(x$C)))
        lab <- paste0(rep(cn, each = length(rn)), ":", rep(rn, times = length(cn)))[both]
        cat(pad, if (srch) " sign not determined (both signs found by the search): "
                 else " sign not determined by the data: ", paste(lab, collapse = ", "), "\n", sep = "")
      }
    }
    cat(pad, " ranges: C.ident.lower / C.ident.upper (X: X.ident.*)\n", sep = "")
  }
  invisible(NULL)
}

## Shared body: prints when `print = TRUE`, otherwise just packages the
## arguments so print.summary.nmfre can render later.
.summary.nmfre.body <- function(object, ci.show = FALSE, max.coef = 20,
                                by = c("covariate", "basis"), print = TRUE,
                                ...) {
  by <- match.arg(by)
  if (!isTRUE(print))
    return(structure(list(object = object, ci.show = ci.show,
                          max.coef = max.coef, by = by),
                     class = "summary.nmfre"))

  ## legacy alias (original snake_case name)
  .extra <- list(...)
  if (!is.null(.extra$show_ci)) ci.show <- .extra$show_ci

  x <- object
  P <- nrow(x$X); Q <- ncol(x$X)
  K <- ncol(x$C); N <- ncol(x$U)

  # ---- header ----
  cat(sprintf("NMF-RE: Y(%d,%d) = X(%d,%d) [C(%d,%d) A + U(%d,%d)]\n",
              P, N, P, Q, Q, K, Q, N))

  # ---- convergence ----
  ## Shared formatter, so print() and every summary() word this identically.
  .print.convergence(x)

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

  # ---- degeneracy (two bases merging, C diverging) ----
  dg <- x$degeneracy
  if (!is.null(dg) && isTRUE(dg$flag))
    cat(sprintf(paste0("Degeneracy:           bases %d and %d nearly coincide (cosine %.3f);",
                       " C cancellation %.1f.\n",
                       "                      Their rows of C are not interpretable; refit from another\n",
                       "                      start (e.g. another seed) and keep the larger likelihood.\n"),
                dg$pair[1], dg$pair[2], dg$max.cosine, dg$cancellation))

  # ---- starts (X.init = "mix2") ----
  st <- x$X.init.search
  if (!is.null(st) && any(st$chosen)) {
    ch <- which(st$chosen)
    gap <- 2 * (st$nll[-ch] - st$nll[ch]); gap <- gap[is.finite(gap)]
    cat(sprintf("Starts (mix2):        kept %s of %d fits (k-means, SPA%s); %s\n",
                if (st$start[ch] == "spa") "SPA" else sprintf("k-means seed %d", st$seed[ch]), nrow(st),
                if (nrow(st) > 2) sprintf(", %d k-means after a degenerate fit", nrow(st) - 2) else "",
                if (!length(gap)) "the other fits failed"
                else if (max(gap) < 1e-3) sprintf("%s the same likelihood",
                                                  if (length(gap) == 1) "the other reaches" else "the others reach")
                else if (length(gap) == 1) sprintf("the other is higher in 2*NLL by %.3g", gap)
                else sprintf("the others are higher in 2*NLL by %.3g to %.3g", min(gap), max(gap))))
  }

  # ---- effective rank (BLUP scores B = Theta A + U, Q x N) ----
  eff <- if (!is.null(x$B.blup)) .effective.rank(x$B.blup) else NA_real_
  if (is.finite(eff)) cat(sprintf("Effective Rank:       %.2f / %d  (%.1f%%)\n",
                                  eff, Q, 100 * eff / Q))

  # ---- identification region (from nmfre(); Q = 3) ----
  .nmfre.print.ident(x, Q)
  # ---- stationarity check of a "semi" fit ----
  if (!is.null(x$kkt) && is.finite(x$kkt$improvement))
    cat(sprintf("Stationarity check:   %d exact X-steps lower 2*NLL by %.3g%s\n", x$kkt$iterations,
                x$kkt$improvement, if (x$kkt$improvement > x$kkt$tol) " -- refit with X.update = \"exact\""
                else " (negligible)"))

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
    cat("\n")
    if (identical(x$identified.approx, FALSE) && !is.null(x$identification$note)) {
      cat("  Identification: ", x$identification$note, "\n", sep = "")
    } else .ident.print.line(x$identified.approx, x$identification, Q, what = "C")
    cat("Coefficients:\n")

    cf <- x$coefficients
    cf <- cf[.coef.order.by(cf, by), , drop = FALSE]   # grouping order (by)
    ## Same Q x K cap the nmfkc / nmfae printers apply; this copy predates it.
    .n.all <- nrow(cf)
    .sh <- .coef.show.idx(cf, max.coef)
    cf <- cf[.sh$idx, , drop = FALSE]
    if (.sh$truncated || nrow(cf) < .n.all)
      cat(sprintf("(showing %d of %d coefficients; use x$coefficients for all)\n",
                  nrow(cf), .n.all))

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
    if (isTRUE(ci.show) && all(c("CI_low", "CI_high") %in% names(cf))) {
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
# predict.nmfre() - fixed-effect prediction
# ============================================================

#' @title Predict method for nmfre objects
#' @description
#' Predicts the response from a fitted NMF-RE model. With \code{newA} supplied,
#' returns the \strong{fixed-effect} prediction \eqn{X\,\Theta\,A_{new}} (new
#' units carry no estimated random effect \eqn{U}, so only the population mean
#' is predicted). Without \code{newA}, returns the in-sample BLUP fit
#' \eqn{X(\Theta A + U)} (\code{object$XB.blup}), matching \code{fitted()}.
#' @param object An object of class \code{"nmfre"} from \code{\link{nmfre}}.
#' @param newA Optional covariate matrix (K x M) for new units. If \code{NULL}
#'   (default), the in-sample BLUP fit is returned.
#' @param ... Not used.
#' @return A matrix of predicted (P x M) values.
#' @seealso \code{\link{nmfre}}, \code{\link{fitted.nmf}}
#' @export
predict.nmfre <- function(object, newA = NULL, ...) {
  if (is.null(newA)) {
    object$XB.blup
  } else {
    newA <- base::as.matrix(newA)
    object$X %*% (object$C %*% newA)
  }
}


# ============================================================
# Identification region of C under the isotropic random effects
# ============================================================

## Under the working model, (X, C) and (XO, O'C) give the same distribution
## exactly when O is orthogonal with O1 = 1 and XO >= 0 (and O'C >= 0 for the
## non-negative variant).  These (X, C) form the identification region.  For
## Q = 3 the O lie on one circle O(phi) about 1; the feasible arc through
## phi = 0 is the local identification region.  Each entry of O(phi)'C is
## a + b cos(phi) + c sin(phi), so its range over the arc is exact; only the
## ends of the arc are found numerically.  Q = 1, 2: the region is a point up
## to relabelling.  Q >= 4: components pinned by anchors are removed one at a
## time (.nmfre.ident.pin()); the Q' that remain are a point (Q' <= 2), an arc
## computed as for Q = 3 (Q' = 3), or searched numerically (4 <= Q' <= 7).
.nmfre.O3.basis <- function() base::qr.Q(base::qr(base::cbind(1, base::diag(3))))[, 2:3, drop = FALSE]

.nmfre.O3 <- function(phi, V) {
  R <- base::matrix(c(base::cos(phi), base::sin(phi), -base::sin(phi), base::cos(phi)), 2, 2)
  base::matrix(1 / 3, 3, 3) + V %*% R %*% base::t(V)
}

.nmfre.ident.region <- function(X, C, nonneg = FALSE, tol = 1e-10, ngrid = 7200L,
                                ident = base::list(search = TRUE)) {
  Q <- base::ncol(X)
  if (Q <= 2L)
    return(base::list(Q = Q, computed = TRUE, method = "exact", full = FALSE, arc = c(0, 0), lower = C,
                      upper = C, width = C * 0, phi.lower = C * 0, phi.upper = C * 0,
                      X.lower = X, X.upper = X, X.width = X * 0))
  if (Q > 3L) return(.nmfre.ident.reduce(X, C, nonneg = nonneg, tol = tol, ngrid = ngrid, ident = ident))
  V <- .nmfre.O3.basis()
  ctol <- tol * base::max(1, base::max(base::abs(C)))
  g <- function(phi) {
    O <- .nmfre.O3(phi, V)
    m <- base::min(X %*% O)
    if (nonneg) m <- base::min(m, base::min(base::crossprod(O, C)) * tol / ctol)
    m
  }
  ph <- base::seq(-pi, pi, length.out = ngrid + 1L)
  ok <- base::vapply(ph, g, 0) >= -tol
  i0 <- base::which.min(base::abs(ph)); ok[i0] <- TRUE
  full <- base::all(ok)
  if (full) {
    lo <- -pi; hi <- pi
  } else {
    n <- base::length(ph); h <- ph[2] - ph[1]
    a <- i0; while (a > 1L && ok[a - 1L]) a <- a - 1L
    b <- i0; while (b < n && ok[b + 1L]) b <- b + 1L
    ## the grid's two ends are the same point of the circle: wrap around
    if (a == 1L) { a2 <- n; while (ok[a2 - 1L]) a2 <- a2 - 1L; lo.in <- ph[a2] - 2 * pi } else lo.in <- ph[a]
    if (b == n)  { b2 <- 1L; while (ok[b2 + 1L]) b2 <- b2 + 1L; hi.in <- ph[b2] + 2 * pi } else hi.in <- ph[b]
    gt <- function(p) g(p) + tol
    lo <- if (gt(lo.in - h) < 0) stats::uniroot(gt, c(lo.in - h, lo.in), tol = 1e-12)$root else lo.in
    hi <- if (gt(hi.in + h) < 0) stats::uniroot(gt, c(hi.in, hi.in + h), tol = 1e-12)$root else hi.in
    lo <- base::min(lo, 0); hi <- base::max(hi, 0)
    ## an arc this short is the tolerance, not the region: a point (e.g. an anchor row)
    if (hi - lo < 1e-7) lo <- hi <- 0
  }
  ## O(phi)'C = A0 + cos(phi) Bm + sin(phi) Cm;  X O(phi) likewise
  Jt <- base::matrix(c(0, -1, 1, 0), 2, 2)               # R(phi)' = cos I + sin Jt
  rc <- .nmfre.arc.range(base::matrix(base::colSums(C) / 3, 3, base::ncol(C), byrow = TRUE),
                         V %*% base::crossprod(V, C), V %*% Jt %*% base::crossprod(V, C), lo, hi)
  rx <- .nmfre.arc.range(base::matrix(base::rowSums(X) / 3, base::nrow(X), 3),
                         X %*% V %*% base::t(V), X %*% V %*% base::t(Jt) %*% base::t(V), lo, hi)
  base::list(Q = 3L, computed = TRUE, method = "exact", full = full, arc = c(lo, hi), V = V,
             lower = rc$lower, upper = rc$upper, width = rc$upper - rc$lower,
             phi.lower = rc$phi.lower, phi.upper = rc$phi.upper,
             X.lower = rx$lower, X.upper = rx$upper, X.width = rx$upper - rx$lower)
}

## Exact range over phi in [lo, hi] of each entry of A0 + cos(phi) Bm + sin(phi) Cm:
## the extremes are at the ends or where the derivative vanishes.
.nmfre.arc.range <- function(A0, Bm, Cm, lo, hi) {
  lower <- upper <- A0; phl <- phu <- A0 * 0
  for (j in base::seq_along(A0)) {
    fj <- function(p) A0[j] + base::cos(p) * Bm[j] + base::sin(p) * Cm[j]
    p0 <- base::atan2(Cm[j], Bm[j])
    cand <- c(lo, hi, base::as.vector(base::outer(c(p0, p0 + pi), 2 * pi * (-2:2), "+")))
    cand <- cand[cand >= lo - 1e-12 & cand <= hi + 1e-12]
    v <- base::vapply(cand, fj, 0)
    lower[j] <- base::min(v); upper[j] <- base::max(v)
    phl[j] <- cand[base::which.min(v)]; phu[j] <- cand[base::which.max(v)]
  }
  base::list(lower = lower, upper = upper, phi.lower = phl, phi.upper = phu)
}

## Options of the identification computation from the arguments of nmfre()
## (exact matching: `ident.search` must not pick up another argument).
.nmfre.ident.args <- function(ea) base::list(
  search   = if (!is.null(ea[["ident.search"]])) base::isTRUE(ea[["ident.search"]]) else TRUE,
  n.dir    = if (!is.null(ea[["ident.dirs"]])) base::as.integer(ea[["ident.dirs"]]) else 300L,
  n.global = if (!is.null(ea[["ident.draws"]])) base::as.integer(ea[["ident.draws"]]) else 20000L,
  zero.tol = if (!is.null(ea[["ident.zero.tol"]])) ea[["ident.zero.tol"]] else 0,
  seed     = if (!is.null(ea[["seed"]])) ea[["seed"]] else 1L)

## Components pinned by anchors, found one at a time.  If a row of X is
## non-zero, outside the components already pinned, on component q alone, then
## row q of every admissible O is non-negative (XO >= 0 on that row).  A
## non-negative row with sum 1 (O1 = 1) and length 1 (O orthogonal) is a unit
## vector, so O keeps component q (up to relabelling) and is the identity on
## it; what is left is the same problem for the other components.  Under the
## non-negative variant a column of C non-zero on one remaining component pins
## it in the same way (O'C >= 0).  The result is about exact zeros: by default
## (zero.tol = 0) only entries equal to 0 count, which the exact X- and C-steps
## return for active constraints; a positive zero.tol (an absolute threshold)
## makes the pinning a numerical diagnostic, not a proof.
.nmfre.ident.pin <- function(X, C, nonneg = FALSE, zero.tol = 0) {
  Q <- base::ncol(X); S <- base::integer(0)
  nzx <- X > zero.tol
  nzc <- base::abs(C) > zero.tol
  repeat {
    rest <- base::setdiff(base::seq_len(Q), S)
    if (base::length(rest) <= 1L) break
    r <- nzx[, rest, drop = FALSE]
    new <- base::integer(0)
    for (i in base::which(base::rowSums(r) == 1L)) new <- base::c(new, rest[base::which(r[i, ])])
    if (nonneg) {
      cc <- nzc[rest, , drop = FALSE]
      for (k in base::which(base::colSums(cc) == 1L)) new <- base::c(new, rest[base::which(cc[, k])])
    }
    new <- base::setdiff(base::unique(new), S)
    if (!base::length(new)) break
    S <- base::sort(base::union(S, new))
  }
  S
}

## Q >= 4: remove the pinned components, then treat the Q' that remain.
.nmfre.ident.reduce <- function(X, C, nonneg = FALSE, tol = 1e-10, ngrid = 7200L,
                                ident = base::list(search = TRUE)) {
  Q <- base::ncol(X)
  zt <- if (is.null(ident$zero.tol)) 0 else ident$zero.tol
  S <- .nmfre.ident.pin(X, C, nonneg = nonneg, zero.tol = zt)
  rest <- base::setdiff(base::seq_len(Q), S); Qr <- base::length(rest)
  base.out <- base::list(Q = Q, computed = TRUE, full = FALSE, arc = NULL, pinned = S, Q.rest = Qr, zero.tol = zt,
                         lower = C, upper = C, width = C * 0 + 0, X.lower = X, X.upper = X, X.width = X * 0 + 0)
  if (Qr <= 2L) return(base::c(base.out, base::list(method = "anchors")))
  if (Qr == 3L) {
    r <- .nmfre.ident.region(X[, rest, drop = FALSE], C[rest, , drop = FALSE], nonneg = nonneg,
                             tol = tol, ngrid = ngrid)
    method <- if (base::length(S)) "anchors+exact" else "exact"
  } else if (base::isTRUE(ident$search) && Qr <= 7L) {
    r <- .nmfre.ident.search(X[, rest, drop = FALSE], C[rest, , drop = FALSE], nonneg = nonneg, tol = tol,
                             n.dir = if (is.null(ident$n.dir)) 300L else ident$n.dir,
                             n.global = if (is.null(ident$n.global)) 20000L else ident$n.global,
                             seed = if (is.null(ident$seed)) 1L else ident$seed)
    method <- if (base::length(S)) "anchors+search" else "search"
  } else return(base::list(Q = Q, computed = FALSE, pinned = S, Q.rest = Qr, zero.tol = zt, method = "none"))
  out <- base.out
  out$lower[rest, ] <- r$lower; out$upper[rest, ] <- r$upper; out$width <- out$upper - out$lower
  out$X.lower[, rest] <- r$X.lower; out$X.upper[, rest] <- r$X.upper; out$X.width <- out$X.upper - out$X.lower
  out$method <- method
  if (!is.null(r$arc)) out$arc <- r$arc
  if (!is.null(r$search)) out$search <- r$search
  out
}

## Numerical search of the region for 4 <= Q <= 7; approximate, without a
## certificate.  The O with O1 = 1 are O = 11'/Q + V R V' with R orthogonal of
## size Q - 1 (a family of dimension (Q-1)(Q-2)/2).  Local part: from O = I,
## random directions (Cayley curves R(t) = cayley(t A), A skew-symmetric) are
## walked outward in small steps and stopped at the first infeasible step,
## whose position is then refined by bisection, so every point used is reached
## along a path feasible at each step.  The ranges over these points
## approximate the local region from inside; a violation between two steps can
## be missed, so they are not a certified bound.  Global part: random draws over
## the whole family.  A feasible O other than a relabelling shows that the
## representation is not unique globally, but it need not lie in the local
## region (another connected component), so these draws are counted and kept
## out of the ranges.
.nmfre.ident.search <- function(X, C, nonneg = FALSE, tol = 1e-10, n.dir = 300L, n.global = 20000L,
                                seed = 1L, step = 1.05) {
  Q <- base::ncol(X); m <- Q - 1L
  .rng <- .nmfkc.rng.save(seed)
  on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
  V <- base::qr.Q(base::qr(base::cbind(1, base::diag(Q))))[, 2:Q, drop = FALSE]
  J <- base::matrix(1 / Q, Q, Q)
  Omat <- function(R) J + V %*% R %*% base::t(V)
  tolX <- tol * base::max(1, base::max(X)); tolC <- tol * base::max(1, base::max(base::abs(C)))
  feas <- function(O) base::min(X %*% O) >= -tolX && (!nonneg || base::min(base::crossprod(O, C)) >= -tolC)
  cayley <- function(A) base::solve(base::diag(m) + A, base::diag(m) - A)
  rskew <- function() {
    A <- base::matrix(0, m, m); A[base::upper.tri(A)] <- stats::rnorm(m * (m - 1L) / 2)
    A <- A - base::t(A); A / base::sqrt(base::sum(A^2))
  }
  found <- base::list(base::diag(Q)); n.loc <- 0L
  for (d in base::seq_len(n.dir)) {
    A <- rskew()
    t <- 1e-6
    if (!feas(Omat(cayley(t * A)))) next
    lo <- t; hi <- NA_real_
    while (t < 1e3) {                                     # walk out in small steps to the first violation
      t <- t * step
      if (feas(Omat(cayley(t * A)))) lo <- t else { hi <- t; break }
    }
    if (!base::is.na(hi)) for (it in 1:30) { mid <- (lo + hi) / 2; if (feas(Omat(cayley(mid * A)))) lo <- mid else hi <- mid }
    if (lo < 1e-5) next                                   # a segment this short is the tolerance
    n.loc <- n.loc + 1L
    found[[base::length(found) + 1L]] <- Omat(cayley(lo * A))
  }
  n.glob <- 0L
  for (b in base::seq_len(n.global)) {
    qz <- base::qr(base::matrix(stats::rnorm(m * m), m, m))
    O <- Omat(base::qr.Q(qz) %*% base::diag(base::sign(base::diag(base::qr.R(qz))), m))
    if (!feas(O)) next
    jm <- base::apply(O, 2, base::which.max)
    if (base::anyDuplicated(jm) || base::max(base::abs(O - base::diag(Q)[, jm])) > 1e-6) n.glob <- n.glob + 1L
  }
  Cs <- base::matrix(base::vapply(found, function(O) base::as.vector(base::crossprod(O, C)), base::numeric(base::length(C))),
                     ncol = base::length(found))
  Xs <- base::matrix(base::vapply(found, function(O) base::as.vector(X %*% O), base::numeric(base::length(X))),
                     ncol = base::length(found))
  lower <- upper <- C; lower[] <- base::apply(Cs, 1, base::min); upper[] <- base::apply(Cs, 1, base::max)
  Xl <- Xu <- X; Xl[] <- base::pmax(base::apply(Xs, 1, base::min), 0); Xu[] <- base::apply(Xs, 1, base::max)
  base::list(Q = Q, computed = TRUE, method = "search", full = FALSE, arc = NULL,
             lower = lower, upper = upper, width = upper - lower, X.lower = Xl, X.upper = Xu, X.width = Xu - Xl,
             search = base::list(dim = m * (m - 1L) / 2, n.dir = n.dir, n.global = n.global,
                                 local.feasible = n.loc, global.nonperm = n.glob))
}

## geomin criterion on the row shares of X (Yates 1987; Browne 2001)
.nmfre.geomin.crit <- function(X, eps = 0.01) {
  rs <- base::rowSums(X); keep <- rs > 0
  S <- X[keep, , drop = FALSE] / rs[keep]
  base::sum(base::apply(S^2 + eps, 1, base::prod)^(1 / base::ncol(X)))
}

## The member of the local identification region that minimizes geomin.  On a
## full circle the cyclic permutations O(+-2pi/3) are in the region too, so the
## search is kept to the arc between them that contains phi = 0.
.nmfre.orient.geomin <- function(X, region, eps = 0.01, ngrid = 2000L) {
  V <- region$V
  crit <- function(p) .nmfre.geomin.crit(X %*% .nmfre.O3(p, V), eps)
  arc <- if (region$full) c(-pi / 3, pi / 3) else region$arc
  ph <- base::seq(arc[1], arc[2], length.out = ngrid + 1L)
  v <- base::vapply(ph, crit, 0)
  j <- base::which.min(v); h <- ph[2] - ph[1]
  o <- stats::optimize(crit, c(base::max(arc[1], ph[j] - h), base::min(arc[2], ph[j] + h)), tol = 1e-10)
  p <- if (o$objective <= v[j]) o$minimum else ph[j]
  base::list(phi = p, O = .nmfre.O3(p, V), criterion = base::min(o$objective, v[j]), criterion.start = crit(0))
}

## Rotate a fit within its identification region: (X, C, U) -> (XO, O'C, O'U),
## with the fields derived from them recomputed.  XB, XB.blup and the fit
## statistics do not change.
.nmfre.rotate <- function(object, O) {
  colnorm_prob <- function(M, eps = 1e-12) base::sweep(M, 2, base::pmax(base::colSums(M), eps), "/")
  dX <- base::dimnames(object$X); dC <- base::dimnames(object$C); dU <- base::dimnames(object$U)
  object$X <- base::pmax(object$X %*% O, 0); base::dimnames(object$X) <- dX
  object$C <- base::crossprod(O, object$C); base::dimnames(object$C) <- dC
  object$U <- base::crossprod(O, object$U); base::dimnames(object$U) <- dU
  if (!is.null(object$B)) {
    dB <- base::dimnames(object$B); object$B <- base::crossprod(O, object$B); base::dimnames(object$B) <- dB
    object$B.prob <- colnorm_prob(base::pmax(object$B, 0))
  }
  if (!is.null(object$B.blup)) {
    dB <- base::dimnames(object$B.blup); object$B.blup <- base::crossprod(O, object$B.blup)
    base::dimnames(object$B.blup) <- dB
    object$B.blup.pos <- base::pmax(object$B.blup, 0)
    object$B.blup.prob <- colnorm_prob(object$B.blup.pos)
  }
  object$X.prob <- object$X / (base::rowSums(object$X) + 1e-10)
  object$X.cluster <- base::apply(object$X.prob, 1, base::which.max)
  object$degeneracy <- .nmfre.degeneracy(object$X, object$C)
  object
}

## Identification region and orientation of a fit, as returned by nmfre():
## orient the representation within its local region (geomin, or none) and
## add the ranges of C and X over the region and the identification list.
.nmfre.orient.fit <- function(out, nonneg = FALSE, orient = "none", eps = 0.01,
                              ident = base::list(search = TRUE)) {
  Q <- base::ncol(out$X)
  region <- .nmfre.ident.region(out$X, out$C, nonneg = nonneg, ident = ident)
  orientation <- base::list(rule = orient, angle = 0, O = base::diag(Q), moved = 0, eps = eps)
  if (orient == "geomin" && isTRUE(region$computed) && region$Q == 3L && base::diff(region$arc) > 1e-9) {
    og <- .nmfre.orient.geomin(out$X, region, eps = eps)
    C.before <- out$C
    out <- .nmfre.rotate(out, og$O)
    orientation <- base::list(rule = "geomin", angle = og$phi, O = og$O,
                              moved = base::max(base::abs(out$C - C.before)),
                              criterion = og$criterion, criterion.start = og$criterion.start, eps = eps)
    region <- .nmfre.ident.region(out$X, out$C, nonneg = nonneg, ident = ident)
  }
  w.max <- if (isTRUE(region$computed)) base::max(region$width) else NA_real_
  orientation$note <- if (!isTRUE(region$computed))
      base::sprintf("identification region not computed (%d components left after the anchors%s)",
                    region$Q.rest, if (isTRUE(ident$search)) "; searched for at most 7" else "; ident.search = FALSE")
    else if (region$Q <= 2L) "Q <= 2: the region is a point up to relabelling"
    else if (region$Q >= 4L && w.max <= 1e-8 && identical(region$method, "anchors"))
      "Q >= 4: every basis pinned by an anchor, one at a time: the region is a point up to relabelling"
    else if (region$Q >= 4L && grepl("search", region$method))
      base::paste0("Q >= 4: region searched numerically", if (base::length(region$pinned)) " after the anchors" else "", "; ",
                   if (w.max <= 1e-8) "no other representation found near the fit (evidence, not a proof)"
                   else "the ranges are approximate (along feasible paths from the fit, not certified bounds)",
                   if (isTRUE(region$search$global.nonperm > 0))
                     base::sprintf("; %d feasible transformations other than relabellings found elsewhere (not unique globally)",
                                   region$search$global.nonperm) else "")
    else if (region$Q >= 4L && w.max > 1e-8) "Q >= 4: orientation is computed for Q = 3 only"
    else if (region$Q >= 4L) "Q >= 4: the region is a point after the anchors"
    else if (base::diff(region$arc) <= 1e-9 && orientation$angle == 0) "the region is a point: nothing to orient"
    else NULL
  if (isTRUE(region$computed)) {
    out$C.ident.lower <- region$lower; out$C.ident.upper <- region$upper; out$C.ident.width <- region$width
    base::dimnames(out$C.ident.lower) <- base::dimnames(out$C.ident.upper) <-
      base::dimnames(out$C.ident.width) <- base::dimnames(out$C)
    out$X.ident.lower <- region$X.lower; out$X.ident.upper <- region$X.upper; out$X.ident.width <- region$X.width
    base::dimnames(out$X.ident.lower) <- base::dimnames(out$X.ident.upper) <-
      base::dimnames(out$X.ident.width) <- base::dimnames(out$X)
  }
  out$identification <- base::list(
    width = w.max,
    arc = if (isTRUE(region$computed)) region$arc else NULL,
    full.circle = isTRUE(region$full), method = region$method,
    pinned = region$pinned, Q.rest = region$Q.rest, search = region$search, orient = orientation)
  out
}

## Profile marginal likelihood intervals for entries of C (sign-free C).  For
## each entry the profile is followed outward from the member of the local
## identification region attaining the end (from the estimate when the region
## is not computed), in steps of `step` standard errors, each step warm-started
## from the previous one.  The optimization is local (a continuation), which is
## meant to keep the profile in the component of the estimate, away from other
## components and relabelings that also attain the maximum; it guarantees
## neither this nor a global maximum.  `below.fit` is how far a step went below
## the fit's 2 x NLL (a higher likelihood than the fit's); it is not returned.
.nmfre.profile.ci <- function(object, Y, A, se, region, which, level = 0.95, step = 0.25,
                              maxstep = 24L) {
  X <- object$X; C <- object$C; P <- base::nrow(Y); N <- base::ncol(Y)
  Q <- base::ncol(X); K <- base::nrow(A); crit <- stats::qchisq(level, 1)
  nll <- function(Xm, Cm, s2, t2) {
    S <- t2 * base::tcrossprod(Xm) + s2 * base::diag(P); R <- Y - Xm %*% Cm %*% A
    if (!base::all(base::is.finite(S)) || !base::all(base::is.finite(R))) return(1e10)
    ch <- base::tryCatch(base::chol(S), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    0.5 * (N * (P * base::log(2 * pi) + 2 * base::sum(base::log(base::diag(ch)))) +
             base::sum(base::backsolve(ch, R, transpose = TRUE)^2))
  }
  unpack <- function(par, idx, cval) {
    z <- base::matrix(par[1:(P * Q)], P, Q); Xm <- z^2; Xm <- base::sweep(Xm, 2, base::colSums(Xm), "/")
    th <- base::numeric(Q * K); th[-idx] <- par[P * Q + base::seq_len(Q * K - 1L)]; th[idx] <- cval
    base::list(X = Xm, C = base::matrix(th, Q, K), s2 = base::exp(par[P * Q + Q * K]),
               t2 = base::exp(par[P * Q + Q * K + 1L]))
  }
  pack <- function(Xm, Cm, idx) c(base::as.vector(base::sqrt(base::pmax(Xm, 1e-6))), base::as.vector(Cm)[-idx],
                                  base::log(object$sigma2), base::log(object$tau2))
  nll0 <- nll(X, C, object$sigma2, object$tau2)
  have.region <- isTRUE(region$computed) && region$Q == 3L
  out <- base::lapply(which, function(idx) {
    sd <- base::max(se[idx], 1e-8); est <- C[idx]
    side <- function(dir) {
      if (have.region) {
        p <- if (dir > 0) region$phi.upper[idx] else region$phi.lower[idx]
        O <- .nmfre.O3(p, region$V); X0 <- X %*% O; C0 <- base::crossprod(O, C)
        c0 <- if (dir > 0) region$upper[idx] else region$lower[idx]
      } else { X0 <- X; C0 <- C; c0 <- est }
      par <- pack(X0, C0, idx); dprev <- 0; cprev <- c0; low <- 0
      for (s in base::seq_len(maxstep)) {
        cval <- c0 + dir * step * s * sd
        o <- stats::optim(par, function(pp) { u <- unpack(pp, idx, cval); nll(u$X, u$C, u$s2, u$t2) },
                          method = "BFGS", control = base::list(maxit = 5000, reltol = 1e-12))
        par <- o$par; d <- 2 * (o$value - nll0); low <- base::min(low, d)
        if (d > crit)
          return(base::list(end = cprev + (crit - dprev) / (d - dprev) * (cval - cprev), reached = TRUE, low = low))
        dprev <- base::max(d, 0); cprev <- cval
      }
      base::list(end = c0 + dir * step * maxstep * sd, reached = FALSE, low = low)
    }
    up <- side(1); dn <- side(-1)
    base::data.frame(idx = idx, Lower = dn$end, Upper = up$end, reached = up$reached && dn$reached,
                     below.fit = -base::min(up$low, dn$low, 0))
  })
  base::do.call(base::rbind, out)
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
#'     \item{\code{profile}}{Logical. If \code{TRUE}, also computes profile-likelihood
#'       intervals for the entries of \eqn{C} (sign-free \eqn{C} only). Default
#'       \code{FALSE}: each entry needs a sequence of constrained refits.}
#'     \item{\code{profile.which}}{Entries of \eqn{vec(C)} to profile. Default all.}
#'     \item{\code{profile.level}}{Level of the profile intervals. Default 0.95.}
#'     \item{\code{profile.step}, \code{profile.maxstep}}{Step, in standard errors, and
#'       number of steps of the outward search. Defaults 0.25 and 24.}
#'   }
#'
#' @details
#' Inference is on the representation returned by \code{\link{nmfre}}
#' (by default the point where its algorithm stopped; with
#' \code{orient = "geomin"} the geomin member of its identification region; see
#' there).  Its \code{C.ident.*} and \code{X.ident.*} ranges are kept.
#'
#' \strong{Profile intervals.} With \code{profile = TRUE}, for each entry the
#' profile marginal likelihood is followed outward from the end of its
#' identification region until twice the drop exceeds the
#' \eqn{\chi^2_1} quantile; the interval is the level set.  The profile is
#' computed by continuation: the entry moves outward in steps of
#' \code{profile.step} standard errors, and each constrained maximization
#' (BFGS) starts from the previous solution.  The search is local, which is
#' meant to keep it away from other components of the region and from
#' relabellings (they attain the same maximum), but it neither certifies a
#' global maximum nor excludes other branches.  \code{C.profile.reached}
#' records whether the level was crossed on both sides within
#' \code{profile.maxstep} steps.
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
#' \item{identified.approx}{Whether \eqn{(X, C)} is unique up to a permutation
#'   of the bases, by the data.  With \code{C.signed = FALSE} it is the
#'   diagnostic of \code{\link{nmfkc.inference}} on \eqn{(X, C)}.  With a signed
#'   \eqn{C} (the default) non-negativity cannot pin the rotation, so it is
#'   \code{NA}: only the working model's isotropic random effects
#'   \eqn{U\sim N(0,\tau^2 I)} do, by allowing orthogonal \eqn{T} alone, which
#'   with the column sums of \eqn{X} fixed is a permutation when \eqn{Q=2} or
#'   every basis has an anchor row.  That is identification by assumption.
#'   An algebraic check on the point estimate: no re-fits.}
#' \item{identification}{The parts of the diagnostic, plus \code{by.assumption}
#'   (the isotropic argument above holds, with \eqn{\tau^2>0}: \eqn{Q\le2}, an
#'   anchor row on every basis, or, for \eqn{Q=3}, an identification region
#'   computed by \code{\link{nmfre}} that is a point) and \code{note}: for a
#'   signed \eqn{C}, and for \eqn{C\ge0} when non-negativity alone does not pin
#'   the factorization but the isotropic random effects do.  The identification-region entries that
#'   \code{\link{nmfre}} put in this list (\code{width}, \code{arc},
#'   \code{full.circle}, \code{orient}) are kept.}
#' \item{C.profile.lower, C.profile.upper, C.profile.reached}{With
#'   \code{profile = TRUE}: profile-likelihood interval of each entry of
#'   \eqn{C} (Q x K matrices, \code{NA} for entries not profiled) and whether
#'   the level was crossed on both sides.}
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
  profile <- base::isTRUE(extra_args$profile)
  wild.B      <- if (!is.null(extra_args$wild.B))      extra_args$wild.B      else 500
  wild.seed   <- if (!is.null(extra_args$wild.seed))   extra_args$wild.seed   else 123
  ## Keep our own seeding out of the caller's random stream.
  .rng <- .nmfkc.rng.save(wild.seed)
  on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
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

  # Sandwich covariance.  The per-observation scores are also exactly what the
  # wild bootstrap needs, so compute score_mat once and share it.
  V_sand <- NULL
  Xt <- base::t(X)
  s2_floor <- base::max(sigma2_used, 1e-12)
  score_mat <- base::matrix(0, Q * K, N)
  for (n in 1:N) {
    a_n <- A[, n, drop = FALSE]
    r_n <- R_C[, n, drop = FALSE]
    g_n <- Xt %*% r_n
    S_n <- -base::tcrossprod(g_n, a_n) / s2_floor
    score_mat[, n] <- base::as.vector(S_n)
  }
  J <- base::matrix(0, Q * K, Q * K)
  for (n in 1:N) J <- J + base::tcrossprod(score_mat[, n])
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
    ## score_mat computed once above (shared with the sandwich-J accumulation)
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
  object$C.se.sandwich <- C.se   # analytic SE (house-style name, matches nmfkc.inference)
  object$C.se.boot    <- C.se.boot
  object$C.ci.lower   <- C.ci.lower
  object$C.ci.upper   <- C.ci.upper
  object$coefficients <- coefficients
  object$C.p.side     <- C.p.side

  ## ---- optional profile-likelihood intervals (sign-free C) ----
  if (base::isTRUE(profile)) {
    if (C.mode != "signed") {
      base::warning("profile = TRUE is implemented for a sign-free C (C.signed = TRUE) only; skipped.")
    } else {
      p.which <- if (!is.null(extra_args$profile.which)) extra_args$profile.which else base::seq_along(C_mat)
      region <- .nmfre.ident.region(X, C_mat, ident = base::list(search = FALSE))   # used for Q = 3 only
      pr <- .nmfre.profile.ci(object, Y, A, se = se_vec, region = region, which = p.which,
              level = if (!is.null(extra_args$profile.level)) extra_args$profile.level else 0.95,
              step = if (!is.null(extra_args$profile.step)) extra_args$profile.step else 0.25,
              maxstep = if (!is.null(extra_args$profile.maxstep)) extra_args$profile.maxstep else 24L)
      plo <- pup <- C_mat * NA_real_; prc <- base::matrix(NA, Q, K, dimnames = base::dimnames(C_mat))
      plo[pr$idx] <- pr$Lower; pup[pr$idx] <- pr$Upper; prc[pr$idx] <- pr$reached
      object$C.profile.lower <- plo; object$C.profile.upper <- pup; object$C.profile.reached <- prc
    }
  }
  ## Is (X, C) unique up to a permutation?  Y ~ X (C A + U): X >= 0 with
  ## columns summing to one; U and B = C A + U are sign-free.  With C >= 0
  ## (C.signed = FALSE) the diagnostic of nmfkc.inference() applies to (X, C).
  ## With a signed C (the default) non-negativity cannot pin T^-1 at all;
  ## only the working model's isotropic random effects U ~ N(0, tau2 I) do --
  ## they allow orthogonal T alone, which with 1'T = 1' is a permutation when
  ## Q = 2 or every basis has an anchor row.  That is identification by
  ## assumption, not by the data, so identified.approx is NA there and
  ## identification$by.assumption records it.  No re-fits.
  sep.tol <- if (!is.null(extra_args$sep.tol)) extra_args$sep.tol else 0.999
  pin.tol <- if (!is.null(extra_args$pin.tol)) extra_args$pin.tol else 0.01
  v <- .ident.verdict(X, C_mat, A, scale.fixed = TRUE, sep.tol = sep.tol, pin.tol = pin.tol)
  ## For Q = 3 nmfre() computed the identification region under those random
  ## effects (the arc; its ends numerically); an anchor row on every basis is only a sufficient
  ## condition for it to be a point, so the region decides when it is there.
  w.region <- if (base::is.list(object$identification)) object$identification$width else NULL
  m.region <- if (base::is.list(object$identification)) object$identification$method else NULL
  ## computed, not searched: Q = 3, or Q >= 4 with the anchors and an arc for what remains
  ## (a numerical search is evidence, not a proof, and does not count here)
  region.known <- !is.null(w.region) && base::is.finite(w.region) &&
    (Q == 3L || (Q >= 4L && !is.null(m.region) && m.region %in% c("anchors", "anchors+exact", "exact")))
  v$identification$by.assumption <- isTRUE(tau2 > 0) &&
    (Q <= 2L || all(v$identification$anchor) || (region.known && w.region <= 1e-8))
  if (C.mode == "signed") {
    v$identified.approx <- NA
    v$identification$note <- paste0(
      "not by the data (C is signed); ",
      if (v$identification$by.assumption) {
        if (Q <= 2L || all(v$identification$anchor))
          "unique up to a permutation under the isotropic random effects"
        else "locally unique under the isotropic random effects (identification region a point)"
      } else if (region.known && isTRUE(tau2 > 0)) {
        base::sprintf("under the isotropic random effects, an identification region of width %.3g (C.ident.lower, C.ident.upper)",
                      w.region)
      } else "not established under the isotropic random effects")
  } else if (identical(v$identified.approx, FALSE) && v$identification$by.assumption) {
    ## C >= 0: non-negativity alone (any invertible T) does not pin the
    ## factorization, but the random effects (orthogonal T) do -- say both,
    ## so the line does not contradict the identification region
    v$identification$note <- paste0(
      "not by non-negativity alone; ",
      if (Q <= 2L || all(v$identification$anchor))
        "unique up to a permutation under the isotropic random effects"
      else "locally unique under the isotropic random effects (identification region a point)")
  }
  object$identified.approx <- v$identified.approx
  ## keep the identification-region entries nmfre() put in this list
  id <- if (base::is.list(object$identification)) object$identification else base::list()
  for (nm in base::names(v$identification)) id[[nm]] <- v$identification[[nm]]
  object$identification    <- id
  ## Tag it as an inference result, as nmfkc.inference / nmf.rrr.inference do.
  ## Without this the class stayed c("nmfre", "nmf") and print() fell to
  ## print.nmf, showing the ten-line fit header instead of the coefficient
  ## table that print.nmf.inference draws.
  if (!inherits(object, "nmf.inference"))
    class(object) <- c("nmf.inference", class(object))
  return(object)
}


# ============================================================
# nmfre.ecv() - NMF-RE-native element-wise cross-validation
# ============================================================

#' @title NMF-RE element-wise cross-validation for rank selection
#' @description
#' Selects the basis rank \eqn{Q} for \code{\link{nmfre}} by Wold-style
#' element-wise (entry-holdout) cross-validation that \strong{exercises the
#' random effects}. For each fold the held-out entries are hidden and filled by
#' iterative imputation: fit the NMF-RE model on the current matrix, replace the
#' held-out entries with the BLUP prediction \eqn{X(\Theta A + U)}, and repeat.
#' The score is the held-out prediction RMSE (\code{sigma}); the selected
#' rank minimizes it.
#'
#' Because the held-out entries of a column are predicted using the random
#' effect \eqn{u_n} fitted from that column's \emph{retained} entries, this
#' evaluates the full NMF-RE model (including \eqn{U}) --- unlike
#' \code{\link{nmfkc.ecv}}, which masks entries by zero weight and predicts from
#' the fixed-effect fit \eqn{X\Theta A} only. The two scores are therefore
#' \emph{not} directly comparable.
#'
#' @param Y Observation matrix (P x N), non-negative.
#' @param A Covariate matrix (K x N). Default is a row of ones (intercept only).
#' @param rank Integer vector of ranks \eqn{Q} to evaluate (default \code{1:3}).
#' @param C.signed Logical. Sign convention for \eqn{\Theta} passed to each
#'   \code{\link{nmfre}} fit (\code{TRUE} = sign-free, default; \code{FALSE} =
#'   non-negative). The basis update rule follows this choice automatically.
#' @param ... Additional arguments:
#'   \itemize{
#'     \item \code{nfolds}: Number of folds (default 5; legacy \code{nfold} also accepted).
#'     \item \code{rounds}: Iterative-imputation rounds per fold (default 4).
#'     \item \code{seed}: RNG seed for the fold assignment (default 1).
#'     \item \code{cores}: Integer; number of parallel workers for the rank
#'       sweep (default \code{getOption("mc.cores", 1L)}). Because the fold
#'       matrix is fixed before the sweep and each \code{\link{nmfre}} fit
#'       self-seeds, results are identical to the sequential run for any
#'       \code{cores}.
#'     \item \code{print.trace}: Logical; print per-rank scores (default \code{FALSE}).
#'     \item Convergence controls forwarded to \code{\link{nmfre}}:
#'       \code{epsilon} (default \code{1e-5}), \code{epsilon.outer}
#'       (default \code{1e-3}), \code{inner.maxit} (default \code{1500}),
#'       \code{outer.maxit} (default \code{80}), \code{maxit}
#'       (default \code{40000}). CV does not need the tight tolerances of a
#'       final fit, so these are loosened by default.
#'   }
#' @return A list of class \code{"nmfre.ecv"} with components:
#'   \describe{
#'     \item{\code{rank}}{The ranks evaluated.}
#'     \item{\code{sigma}}{Named numeric vector of held-out RMSE per rank
#'       (same field name as \code{\link{nmfkc.ecv}}).}
#'     \item{\code{best}}{The rank minimizing \code{sigma}.}
#'     \item{\code{nfolds}, \code{rounds}, \code{C.signed}}{Settings used.}
#'   }
#' @seealso \code{\link{nmfre}}, \code{\link{nmfre.inference}}, \code{\link{nmfkc.ecv}}
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("nlme", quietly = TRUE)) {
#'   Y <- matrix(nlme::Orthodont$distance, 4, 27)
#'   male <- ifelse(nlme::Orthodont$Sex[seq(1, 108, 4)] == "Male", 1, 0)
#'   A <- rbind(intercept = 1, male = male)
#'   cv <- nmfre.ecv(Y, A, rank = 1:3)
#'   cv$best
#'   plot(cv)
#' }
#' }
nmfre.ecv <- function(Y, A = NULL, rank = 1:3, C.signed = TRUE, ...) {
  extra <- base::list(...)
  if (!is.null(extra$Q)) rank <- extra$Q
  ## nfolds is the house-style name (matches nmfkc.ecv etc.); legacy nfold via ...
  nfolds <- if (!is.null(extra$nfolds)) extra$nfolds
            else if (!is.null(extra$nfold)) extra$nfold else 5L
  rounds <- if (!is.null(extra$rounds)) extra$rounds else 4L
  seed   <- if (!is.null(extra$seed))   extra$seed   else 1L
  ## Keep our own seeding out of the caller's random stream.
  .rng <- .nmfkc.rng.save(seed)
  on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
  cores  <- if (!is.null(extra$cores))  extra$cores  else getOption("mc.cores", 1L)
  print.trace <- isTRUE(extra$print.trace)
  ## loosened CV tolerances (overridable via ...)
  epsilon       <- if (!is.null(extra$epsilon))       extra$epsilon       else 1e-5
  epsilon.outer <- if (!is.null(extra$epsilon.outer)) extra$epsilon.outer else 1e-3
  inner.maxit   <- if (!is.null(extra$inner.maxit))   extra$inner.maxit   else 1500L
  outer.maxit   <- if (!is.null(extra$outer.maxit))   extra$outer.maxit   else 80L
  maxit         <- if (!is.null(extra$maxit))         extra$maxit         else 40000L

  Y <- base::as.matrix(Y)
  P <- base::nrow(Y); N <- base::ncol(Y)
  if (is.null(A)) A <- base::matrix(1, nrow = 1, ncol = N)
  A <- base::as.matrix(A)

  base::set.seed(seed)
  fold <- base::matrix(base::sample(base::rep(1:nfolds, length.out = P * N)), P, N)  # fixed across ranks
  rowm <- base::rowMeans(Y, na.rm = TRUE)
  ridx <- base::row(Y)

  ## one warm-started, output-silenced nmfre fit
  fit_quiet <- function(fill, Q, Xi, Ci, Ui) {
    tryCatch(
      suppressMessages({
        utils::capture.output(
          ## first round from the k-means start alone, as before "mix2"
          f <- nmfre(fill, A, rank = Q, C.signed = C.signed,
                     X.init = if (is.null(Xi)) "kmeans" else Xi, C.init = Ci, U.init = Ui,
                     sigma2 = 1, tau2 = 1,
                     epsilon = epsilon, epsilon.outer = epsilon.outer,
                     inner.maxit = inner.maxit, outer.maxit = outer.maxit,
                     maxit = maxit, ident.search = FALSE))
        f
      }),
      error = function(e) {
        if (print.trace) base::message("   fit err Q=", Q, ": ", base::conditionMessage(e))
        NULL
      })
  }

  ## Per-rank task. The fold matrix is fixed above and the only RNG consumer
  ## here is nmfre() (self-seeds), so ranks are independent and deterministic;
  ## .nmfkc.parlapply returns results in input (rank) order, so cores > 1 is
  ## bit-identical to the sequential sweep.
  rank_fun <- function(Q) {
    sse <- 0; nt <- 0
    for (k in 1:nfolds) {
      mask <- (fold == k)
      fill <- Y; fill[mask] <- rowm[ridx[mask]]        # warm fill by row mean
      Xi <- NULL; Ci <- NULL; Ui <- NULL; pred <- NULL
      for (it in 1:rounds) {
        f <- fit_quiet(fill, Q, Xi, Ci, Ui)
        if (is.null(f)) break
        pred <- f$XB.blup
        fill[mask] <- pred[mask]                        # re-impute held-out
        Xi <- f$X; Ci <- f$C; Ui <- f$U                 # warm-start next round
      }
      if (!is.null(pred)) {
        sse <- sse + base::sum((Y[mask] - pred[mask])^2)
        nt  <- nt  + base::sum(mask)
      }
    }
    s <- if (nt > 0) base::sqrt(sse / nt) else NA_real_
    if (print.trace) base::cat(sprintf("  Q=%d  sigma=%.4f\n", Q, s))
    s
  }
  sig_list <- .nmfkc.parlapply(base::as.list(rank), rank_fun,
                               cores = cores, envir = environment())
  sig <- stats::setNames(base::unlist(sig_list, use.names = FALSE),
                         base::as.character(rank))

  best <- rank[base::which.min(sig)]
  out <- base::list(rank = rank, sigma = sig, best = best,
                    nfolds = nfolds, rounds = rounds, C.signed = C.signed)
  base::class(out) <- "nmfre.ecv"
  out
}

#' @title Plot method for nmfre.ecv objects
#' @description Draws the element-wise CV score (held-out RMSE) against the rank,
#'   marking the minimizing rank as \dQuote{Best (Min)}.
#' @param x An object of class \code{"nmfre.ecv"}.
#' @param main Plot title.
#' @param xlab,ylab Axis labels.
#' @param ... Passed to \code{\link[graphics]{plot}}.
#' @return The input object, invisibly.
#' @seealso \code{\link{nmfre.ecv}}
#' @export
plot.nmfre.ecv <- function(x, main = "NMF-RE element-wise CV",
                           xlab = "Rank (Q)",
                           ylab = "sigma.ecv (held-out RMSE)", ...) {
  rk <- x$rank; s <- base::as.numeric(x$sigma)
  graphics::plot(rk, s, type = "b", pch = 19, xlab = xlab, ylab = ylab,
                 main = main, ...)
  if (base::any(is.finite(s))) {
    bi <- base::which.min(s)
    graphics::points(rk[bi], s[bi], pch = 19, col = "red", cex = 1.7)
    ## place the label away from the plot edge; xpd=NA avoids clipping
    lab.pos <- if (bi == 1L) 4 else if (bi == length(rk)) 2 else 3
    graphics::text(rk[bi], s[bi], labels = "Best (Min)", pos = lab.pos,
                   col = "red", xpd = NA)
  }
  invisible(x)
}

#' @title Print method for nmfre.ecv objects
#' @param x An object of class \code{"nmfre.ecv"}.
#' @param ... Not used.
#' @return The input object, invisibly.
#' @export
print.nmfre.ecv <- function(x, ...) {
  base::cat(sprintf("NMF-RE element-wise CV (%d-fold, %d imputation rounds, C.signed=%s)\n",
                    x$nfolds, x$rounds, base::as.character(x$C.signed)))
  tab <- base::data.frame(rank = x$rank, sigma = base::as.numeric(x$sigma))
  base::print(tab, row.names = FALSE)
  base::cat(sprintf("Best rank (min sigma): %d\n", x$best))
  invisible(x)
}

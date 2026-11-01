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
#' With \eqn{Q \ge 3} and a signed \eqn{C}, a fit can drift along a path on which
#' two columns of \eqn{X} approach each other while the corresponding rows of
#' \eqn{C} grow with opposite signs: the difference \eqn{(x_q - x_{q'}) c_q} then
#' represents a direction outside the non-negative cone, and the fit does not
#' settle.  In simulations such fits usually had a clearly LOWER marginal
#' likelihood than a non-degenerate fit from another start (e.g.
#' \code{X.init = "spa"}), i.e. the algorithm had been trapped; refit from another
#' start and keep the larger likelihood.  Three continuous measures (columns of
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
#' chosen by the marginal likelihood: (1) a free fit -- from \code{X.init} if
#' given, else the better of the default start and \code{X.init = "spa"};
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
.nmfre.partial.anchor.fit <- function(Y, A, Q, C.signed, epsilon, maxit, extra_args, call) {
  ea <- extra_args; ea$Q <- NULL; ea$X.anchor <- NULL; ea$.anchor.direct <- NULL
  ## the search compares fits as they stop; the result is oriented once, at the end
  orient <- if (!is.null(ea$orient)) match.arg(ea$orient, c("geomin", "none")) else "geomin"
  orient.eps <- if (!is.null(ea$orient.eps)) ea$orient.eps else 0.01
  ea$orient <- "none"
  Ym <- as.matrix(Y)
  Am <- if (is.null(A)) matrix(1, 1, ncol(Ym)) else as.matrix(A)
  signed <- !(identical(C.signed, FALSE) || identical(C.signed, "nonneg"))
  shift <- if (signed && min(Ym, na.rm = TRUE) < 0) -min(Ym, na.rm = TRUE) else 0
  J <- .nmfre.partial.anchor.rows(extra_args$X.anchor, Ym + shift, Q)
  fit1 <- function(args) do.call(nmfre, c(list(Y = Y, A = A, rank = Q, C.signed = C.signed,
                                              epsilon = epsilon, maxit = maxit), args))
  nll <- function(f) .nmfre.marginal.nll(Ym, f$X, f$C, Am, f$sigma2, f$tau2)
  ## (1) the free fit
  starts <- if (!is.null(ea$X.init)) list(ea) else list(ea, utils::modifyList(ea, list(X.init = "spa")))
  frees <- Filter(Negate(is.null), lapply(starts, function(a) tryCatch(fit1(a), error = function(e) NULL)))
  if (!length(frees)) stop("nmfre(): the free fit before imposing X.anchor failed.", call. = FALSE)
  f0 <- frees[[which.min(vapply(frees, nll, numeric(1)))]]
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
  .nmfre.orient.fit(best, nonneg = !signed, orient = orient, eps = orient.eps)
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
#'       A method name (\code{"kmeans"}, \code{"nndsvd"}, \code{"spa"}, \ldots)
#'       is forwarded to that \code{nmfkc()} initialization.
#'     \item \code{X.anchor}: anchor rows of \eqn{X}, imposed as in
#'       \code{\link{nmfkc}}: \code{NULL} (default), \code{"spa"}, or one row
#'       of \code{Y} per basis.  The zeros are kept for the whole fit, including
#'       through the floor this function applies to \eqn{X}.  There is no
#'       \code{rank <= nrow(A)} limit here: the random effects give the scores
#'       full rank.  The rows used are returned as \code{X.anchor}.
#'       Fewer rows than bases (indices or row names) anchor the first bases
#'       only: row \code{X.anchor[q]} carries basis \code{q} alone and the
#'       other bases are free.  Which basis an anchor row belongs to is then
#'       chosen by likelihood: a free fit (from \code{X.init} if given, else the
#'       better of the default start and \code{X.init = "spa"}), a refit from
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
#'     \item \code{sigma2}: Initial residual variance (default 1).
#'     \item \code{sigma2.update}: Logical. Update \eqn{\sigma^2} during iterations (default \code{TRUE}).
#'     \item \code{tau2}: Initial random effect variance (default 1).
#'     \item \code{tau2.update}: Logical. Update \eqn{\tau^2} by moment matching (default \code{TRUE}).
#'     \item \code{x.postvar}: Logical. Include the posterior-variance term in
#'       the semi-NMF \eqn{X}-step (default \code{TRUE}; advanced). Applies only
#'       to the sign-free / semi-NMF path (\code{C.signed = TRUE}).
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
#'       added to the multiplicative-update denominator.
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
#'       region (see the return value): \code{"geomin"} (default) or \code{"none"}
#'       (the point where the algorithm stopped).
#'     \item \code{orient.eps}: \eqn{\varepsilon} of the geomin criterion (default 0.01).
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
#'       of \eqn{C} grow with opposite signs (with \eqn{Q \ge 3} and a signed
#'       \eqn{C}; the fit does not settle and those rows of \eqn{C} are not
#'       interpretable).  In simulations such fits usually had a clearly lower
#'       marginal likelihood than a fit from another start, so when \code{flag} is
#'       \code{TRUE}, refit with e.g. \code{X.init = "spa"} and keep the fit with the
#'       larger likelihood:
#'       \code{max.cosine} (largest cosine similarity between two columns of
#'       \eqn{X}) and its \code{pair}; \code{cancellation}
#'       (\eqn{\max_k \sum_q |c_{qk}| / \|X c_{\cdot k}\|_1}, large when a composite
#'       effect is a difference of large opposite terms); \code{kappa.X} (condition
#'       number of \eqn{X}); and \code{flag} (\code{max.cosine >= 0.9} or
#'       \code{cancellation >= 5}).  The thresholds were set on two simulation
#'       designs; read the measures as well.  \code{summary()} prints a line when
#'       \code{flag} is \code{TRUE}.}
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
#'   local identification region; the range of each entry of \eqn{\Theta} and
#'   of \eqn{X} over the arc is computed exactly.  A width of 0 means the fit
#'   cannot be moved within the region (a row of \eqn{X} on one basis only pins
#'   it).  For \eqn{Q\le2} the region is a point up to relabelling; for
#'   \eqn{Q\ge4} it is not computed.  With \code{orient = "geomin"} the
#'   representation returned is the member of the local region that minimizes
#'   the geomin criterion on the row shares of \eqn{X} (as rotated loadings in
#'   exploratory factor analysis); nothing moves when the region is a point.
#'   The likelihood and \eqn{X\Theta} do not change.
#'   \describe{
#'     \item{\code{C.ident.lower}, \code{C.ident.upper}, \code{C.ident.width}}{Range
#'       of each entry of \eqn{C} over the local identification region and its
#'       width (Q x K matrices; absent for \eqn{Q\ge4}).}
#'     \item{\code{X.ident.lower}, \code{X.ident.upper}, \code{X.ident.width}}{The
#'       same for \eqn{X} (P x Q matrices).}
#'     \item{\code{identification}}{List: \code{width} (largest entry of
#'       \code{C.ident.width}; \code{NA} for \eqn{Q\ge4}), \code{arc} and
#'       \code{full.circle} of the region, and \code{orient}, the orientation
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
nmfre <- function(Y, A = NULL, rank = 2, C.signed = TRUE,
                  epsilon = 1e-5, maxit = 5000, ...) {
  cl <- base::match.call()
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

  ## Fewer anchor rows than bases (e.g. X.anchor = "x7" for rank = 3): which
  ## basis an anchor row belongs to is part of the estimate.  Fit freely first,
  ## refit from that fit with each assignment of the anchor rows to bases, and
  ## keep the smallest marginal NLL (.nmfre.partial.anchor.fit()).
  X.anchor.arg <- extra_args$X.anchor
  if (!is.null(X.anchor.arg) && !identical(X.anchor.arg, "spa") &&
      length(X.anchor.arg) < Q && !isTRUE(extra_args$.anchor.direct))
    return(.nmfre.partial.anchor.fit(Y, A, Q = Q, C.signed = C.signed, epsilon = epsilon,
                                     maxit = maxit, extra_args = extra_args, call = cl))

  ## MAP penalties (default off).  These act on X (basis) and C (= Theta) as
  ## Gaussian priors and are orthogonal to the random-effect machinery
  ## (U, lambda, sigma2, tau2), so the variance M-steps are unchanged.
  ##   X.L2.ortho : (l/2)||offdiag(X'X)||^2   -- column orthogonality of X
  ##   X.L2.smooth: (l/2) tr(X' L X)          -- path-graph row smoothness of X
  ##   C.L2       : l * ||C||^2               -- ridge on Theta (keeps the signed
  ##                C-step in closed form: a Sylvester ridge-LS solve)
  X.L2.ortho  <- if (!is.null(extra_args$X.L2.ortho))  extra_args$X.L2.ortho  else 0
  X.L2.smooth <- if (!is.null(extra_args$X.L2.smooth)) extra_args$X.L2.smooth else 0
  C.L2        <- if (!is.null(extra_args$C.L2))        extra_args$C.L2        else 0

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
  ## sign-free least squares on the original Y (below).  Before 1.0.2.9000 such
  ## a Y stopped in the initialisation ("The matrix Y should be non-negative").
  Y.shift <- if (C.mode == "signed" && min(Y, na.rm = TRUE) < 0) -min(Y, na.rm = TRUE) else 0

  ## X.anchor: anchor rows of X, kept as exact zeros for the whole fit.  The
  ## start built from them goes to the nmfkc() initialisation as a matrix,
  ## whose multiplicative updates keep the zeros, and the floor pmax(X, .eps)
  ## below puts them back wherever it would lift them.  No rank check as in
  ## nmfkc(): with the random effects U the scores have full rank anyway.
  ## The rows are chosen, and the start is built, on the same shifted Y as the
  ## nmfkc() initialisation (as X.init = "spa" does); before 1.0.2.9000
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
      if (C.mode == "signed") {
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
      X <- .floorX(X)
      normed <- .nmfre.normalize.X(X, C_mat, U)
      X     <- .floorX(normed$X)
      C_mat <- if (C.mode == "nonneg") pmax(normed$C, .eps) else normed$C
      U     <- normed$U

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

    ## two bases merging while C diverges (see .nmfre.degeneracy)
    degeneracy = .nmfre.degeneracy(X, C_mat)
  )

  ## Added only when anchors were used, so every other fit keeps its fields.
  if (!is.null(anchor.rows))
    out$X.anchor <- stats::setNames(anchor.rows,
      (if (!is.null(colnames(out$X))) colnames(out$X) else
         paste0("Basis", seq_len(Q)))[seq_along(anchor.rows)])
  ## identification region; the representation within it (geomin by default)
  orient <- if (!is.null(extra_args$orient)) base::match.arg(extra_args$orient, c("geomin", "none")) else "geomin"
  out <- .nmfre.orient.fit(out, nonneg = (C.mode == "nonneg"), orient = orient,
                           eps = if (!is.null(extra_args$orient.eps)) extra_args$orient.eps else 0.01)
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
  if (Q <= 2L) {
    cat("Identification region: a point up to relabelling of the bases (Q <= 2)\n")
  } else if (is.null(w) || !is.finite(w)) {
    cat("Identification region: not computed (Q >= 4)\n")
  } else if (w <= 1e-8) {
    cat("Identification region: a point (width 0); Theta is unique up to relabelling\n")
  } else {
    cat(sprintf("Identification region: width of Theta up to %.3g (arc %.1f deg%s)\n", w,
                diff(idn$arc) * 180 / pi, if (isTRUE(idn$full.circle)) ", full circle" else ""))
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
        cat(pad, " sign not determined by the data: ", paste(lab, collapse = ", "), "\n", sep = "")
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
                       "                      start (e.g. X.init = \"spa\") and keep the larger likelihood.\n"),
                dg$pair[1], dg$pair[2], dg$max.cosine, dg$cancellation))

  # ---- effective rank (BLUP scores B = Theta A + U, Q x N) ----
  eff <- if (!is.null(x$B.blup)) .effective.rank(x$B.blup) else NA_real_
  if (is.finite(eff)) cat(sprintf("Effective Rank:       %.2f / %d  (%.1f%%)\n",
                                  eff, Q, 100 * eff / Q))

  # ---- identification region (from nmfre(); Q = 3) ----
  .nmfre.print.ident(x, Q)

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
    .ident.print.line(x$identified.approx, x$identification, Q, what = "C")
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
## to relabelling.  Q >= 4: not computed (the region is multi-dimensional).
.nmfre.O3.basis <- function() base::qr.Q(base::qr(base::cbind(1, base::diag(3))))[, 2:3, drop = FALSE]

.nmfre.O3 <- function(phi, V) {
  R <- base::matrix(c(base::cos(phi), base::sin(phi), -base::sin(phi), base::cos(phi)), 2, 2)
  base::matrix(1 / 3, 3, 3) + V %*% R %*% base::t(V)
}

.nmfre.ident.region <- function(X, C, nonneg = FALSE, tol = 1e-10, ngrid = 7200L) {
  Q <- base::ncol(X)
  if (Q <= 2L)
    return(base::list(Q = Q, computed = TRUE, full = FALSE, arc = c(0, 0), lower = C, upper = C,
                      width = C * 0, phi.lower = C * 0, phi.upper = C * 0,
                      X.lower = X, X.upper = X, X.width = X * 0))
  if (Q > 3L) return(base::list(Q = Q, computed = FALSE))
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
  base::list(Q = 3L, computed = TRUE, full = full, arc = c(lo, hi), V = V,
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
.nmfre.orient.fit <- function(out, nonneg = FALSE, orient = "geomin", eps = 0.01) {
  Q <- base::ncol(out$X)
  region <- .nmfre.ident.region(out$X, out$C, nonneg = nonneg)
  orientation <- base::list(rule = orient, angle = 0, O = base::diag(Q), moved = 0, eps = eps)
  if (orient == "geomin" && isTRUE(region$computed) && region$Q == 3L && base::diff(region$arc) > 1e-9) {
    og <- .nmfre.orient.geomin(out$X, region, eps = eps)
    C.before <- out$C
    out <- .nmfre.rotate(out, og$O)
    orientation <- base::list(rule = "geomin", angle = og$phi, O = og$O,
                              moved = base::max(base::abs(out$C - C.before)),
                              criterion = og$criterion, criterion.start = og$criterion.start, eps = eps)
    region <- .nmfre.ident.region(out$X, out$C, nonneg = nonneg)
  }
  orientation$note <- if (!isTRUE(region$computed)) "identification region not computed for Q >= 4"
    else if (region$Q <= 2L) "Q <= 2: the region is a point up to relabelling"
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
    width = if (isTRUE(region$computed)) base::max(region$width) else NA_real_,
    arc = if (isTRUE(region$computed)) region$arc else NULL,
    full.circle = isTRUE(region$full), orient = orientation)
  out
}

## Profile marginal likelihood intervals for entries of C (sign-free C).  For
## each entry the profile is followed outward from the member of the local
## identification region attaining the end (from the estimate when the region
## is not computed), in steps of `step` standard errors, each step warm-started
## from the previous one.  The optimization is local, so the profile stays in
## the component of the estimate (other components and relabelings, which
## also attain the maximum, are not visited).
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
#' (by default oriented by geomin within its identification region; see
#' there).  Its \code{C.ident.*} and \code{X.ident.*} ranges are kept.
#'
#' \strong{Profile intervals.} With \code{profile = TRUE}, for each entry the
#' profile marginal likelihood is followed outward from the end of its
#' identification region until twice the drop exceeds the
#' \eqn{\chi^2_1} quantile; the interval is the level set.  The search is
#' local: other components of the region and relabellings, which attain the
#' same maximum, are not visited.
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
#'   computed by \code{\link{nmfre}} that is a point) and, for a signed
#'   \eqn{C}, \code{note}.  The identification-region entries that
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
      region <- .nmfre.ident.region(X, C_mat)
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
  ## effects exactly; an anchor row on every basis is only a sufficient
  ## condition for it to be a point, so the region decides when it is there.
  w.region <- if (base::is.list(object$identification)) object$identification$width else NULL
  region.known <- Q == 3L && !is.null(w.region) && base::is.finite(w.region)
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
          f <- nmfre(fill, A, rank = Q, C.signed = C.signed,
                     X.init = Xi, C.init = Ci, U.init = Ui,
                     sigma2 = 1, tau2 = 1,
                     epsilon = epsilon, epsilon.outer = epsilon.outer,
                     inner.maxit = inner.maxit, outer.maxit = outer.maxit,
                     maxit = maxit))
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

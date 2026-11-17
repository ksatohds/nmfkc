## nmfkc.ident(): identifiability of a covariate-NMF fit Y ~ X C A
## (X >= 0 with columns summing to 1, C >= 0, A known).
##
## 1. test: the anchor test H(rho) on the OLS estimate of Xi = X C (HC1, intersection-union over candidates):
##    the Q(Q-1) off-diagonal entries and the Q sizes of a candidate are tested one by one (Bonferroni over
##    Q^2), with no ratios; rho.max.  pi.max (reported only) uses the membership contrasts with the column sums
##    of Xi held fixed (Bonferroni over 2Q); pi.max = 1 exactly when H(rho) is not rejected.
## 2. identification region (R = Q): the maximum over the model, Xi.hat (row-wise non-negative least
##    squares in the metric AA'; the group means when A is one-hot), and the range of every entry of C and
##    X over all exact factorizations X = Xt N^-1, C = N S of Xi.hat.
## 3. repair: nmfkc re-run from a point of the region (the one nearest to object's C, or the center); the
##    start is already optimal, so the re-fit stops at once at the exact optimum.
## 4. intervals (ci = TRUE): profile-likelihood intervals of C (Gaussian errors with one variance common to
##    all entries of Y); they contain the identification region and are meant to cover each point of it
##    (the whole region can be covered less often; coverage not yet studied by simulation).
##
## 0. algebraic check (any A, also none; needs only the fitted X and C): sufficient condition, the
##    first-order cone of remixings decided by NNLS (Farkas), and the extent along pairwise, cone and
##    random directions (inner bounds of the region).  Computed last, with its own seed.
##
## Internal helpers carry the prefix .nid. (".ident." is taken by nmfkc.ar.R).

#' @title Identifiability of an nmfkc fit: test, identification region, repair and intervals
#' @description
#' \code{nmfkc.ident} examines whether the factorization \eqn{Y \approx X C A}
#' of a fitted \code{\link{nmfkc}} model is unique, and what can be said about
#' \eqn{C} (\eqn{\Theta}) when it is not.  For the model
#' \eqn{Y = X C A + \varepsilon} with \eqn{X \ge 0} (columns summing to 1),
#' \eqn{C \ge 0} and \eqn{A} known, any \eqn{T} with columns summing to 1,
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
#'     for \eqn{C} that contain the identification region (approximate; see
#'     Intervals below).
#' }
#' Steps 2--4 need as many covariates as bases (\code{nrow(A) == rank}), e.g.
#' group indicators or \code{rank} continuous covariates; for
#' \code{nrow(A) > rank} (or \code{region = FALSE}) only the test is computed.
#'
#' In every case, with or without covariates, the function also runs an
#' \strong{algebraic check} on the fitted \eqn{X} and \eqn{C} themselves
#' (component \code{algebraic}): whether any remixing moves them, and how far.
#' It needs neither \code{Y} nor \code{A}, so \code{nmfkc.ident(fit)} works
#' for a fit without covariates (\eqn{Y \approx XB}, where the test is not
#' available) and for one with covariates alike.
#'
#' @param object An object of class \code{"nmfkc"}
#'   (\code{X.restriction = "colSums"}, the default), with or without
#'   covariates.
#' @param Y Observation matrix (P x N) used in the fit. Needed, with \code{A},
#'   for the test, the region, the repair and the intervals; not needed for the
#'   algebraic check.
#' @param A Covariate matrix (R x N) used in the fit, of full row rank with
#'   \code{R < N}. \code{NULL} (default): only the algebraic check is computed
#'   (for a fit without covariates, or when only the fitted matrices are at
#'   hand).
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
#'     \item{\code{zero.tol}}{Algebraic check: an entry of \eqn{X} (\eqn{C})
#'       at most \code{zero.tol} times its column (row) maximum counts as 0.
#'       Default \code{1e-3}. The threshold is relative; the sufficient
#'       condition does not use it but keeps the rule of \code{summary.nmfkc()}
#'       (an entry at most 1e-3 in absolute terms, or at most 1\% of the
#'       entry it is compared with, counts as 0).}
#'     \item{\code{move.tol}}{Algebraic check: a remixing counts as moving
#'       the fit when it changes an entry of \eqn{X} or an entry of \eqn{C}
#'       relative to its column sum by more than this. Default \code{1e-3}.}
#'     \item{\code{ident.n}}{Algebraic check: number of random directions
#'       added to the pairwise and cone directions. Default 2000.}
#'   }
#'
#' @section Algebraic check:
#' Every \eqn{T} with columns summing to 1, \eqn{XT \ge 0} and
#' \eqn{T^{-1}C \ge 0} gives the same fit; with \eqn{A} of full row rank the
#' non-negativity is that of \eqn{C} itself (of \eqn{B = C} without
#' covariates), so the check needs only \eqn{X} and \eqn{C}. Write
#' \eqn{T = I + tD} with \eqn{D} of zero column sums (\eqn{Q(Q-1)} free
#' entries). Entries at most \code{zero.tol} times their column (row) maximum
#' are set to 0 first, since the multiplicative updates leave small positive
#' values where the fit is at the boundary.
#' \enumerate{
#'   \item \emph{Sufficient condition}: two-sided anchors, or the peeling of
#'     \code{summary.nmfkc()}, with its own rule for zeros (not
#'     \code{zero.tol}); when it holds the fit is identified.
#'   \item \emph{Local cone}: the directions \eqn{D} that keep both
#'     inequalities to first order, \eqn{(XD)_{ij} \ge 0} where
#'     \eqn{X_{ij} = 0} and \eqn{(DC)_{kj} \le 0} where \eqn{C_{kj} = 0}. Whether
#'     it contains a direction with a positive (negative) value of each
#'     off-diagonal entry is decided exactly by a non-negative least-squares
#'     problem (Farkas' lemma), whose residual is such a direction. If none
#'     does, no remixing moves the fit: it is locally identified.
#'   \item \emph{Extent}: along the pairwise remixings
#'     \eqn{x_q \to (1-a)x_q + a x_r}, the cone directions and \code{ident.n}
#'     random directions, the largest step that keeps \eqn{XT \ge 0} (exact,
#'     linear in \eqn{t}) and \eqn{T^{-1}C \ge 0} (grid, then bisection) is
#'     found. The ranges of \eqn{X} and \eqn{C} over these factorizations are
#'     inner bounds of the identification region, and \code{pairs} gives the
#'     exact interval of \eqn{a} for each pair of bases.
#' }
#' Random directions alone miss narrow sets of remixings (a single possible
#' direction is hit with probability 0), which is why the pairwise and cone
#' directions are always included. The check concerns the fitted matrices: a
#' zero in \eqn{\hat X} or \eqn{\hat C} that pins the fit is an estimate and
#' may not survive sampling variation (use the test with \code{Y} and \code{A}
#' for that), and it examines only remixings connected to the fit, not
#' factorizations far away.
#'
#' @details
#' \strong{Test.} \eqn{\hat\Xi = Y A^\top (AA^\top)^{-1}} with HC1 standard
#' errors. A candidate is an ordered set of \eqn{Q} columns \eqn{L} of
#' \eqn{\Xi} (the anchor columns of \eqn{C}) and \eqn{Q} distinct rows \eqn{J}
#' (the anchor rows of \eqn{X}). The null hypothesis \eqn{H(\rho)} (at
#' \code{pi = 1}) says that some candidate has zero off-diagonal entries
#' \eqn{\Xi_{j_q l_{q'}}} (\eqn{q \ne q'}) and sizes
#' \eqn{\Xi_{j_q l_q} \ge \rho \sum_p \Xi_{p l_q}}. Each of these \eqn{Q^2}
#' conditions is tested on the entries of \eqn{\hat\Xi} directly by a one-sided
#' \eqn{t} test at \eqn{\alpha/Q^2} (\eqn{t_{N-R}} critical value; no ratios of
#' estimates enter the standard errors); a candidate survives if none of them
#' rejects, and \eqn{H(\rho)} is rejected if no candidate survives
#' (intersection--union test). \code{rho.max} is the largest \code{rho} not
#' rejected, an upper confidence bound of the size of the anchors.
#'
#' The membership probability of a row (its share in the representative where
#' the columns \eqn{L} are pure, \eqn{X = \Xi_L S_L^{-1}}) is reported, not
#' used for the decision at \code{pi = 1}: \code{pi.max} is 1 when
#' \eqn{H(\rho)} is not rejected, and otherwise the largest \code{pi} for which
#' some candidate passes the \eqn{2Q} one-sided tests of the membership
#' contrasts \eqn{(1-\pi) u_{j l_q} - \pi \sum_{q' \ne q} u_{j l_{q'}} \ge 0}
#' and the sizes, at \eqn{\alpha/(2Q)}, with \eqn{u = \Xi / s} and the column
#' sums \eqn{s} held fixed in the standard errors (letting them vary, with the
#' estimate plugged into the gradient, made the test liberal when \eqn{P} is
#' small). A decision at \code{pi < 1} uses the same contrasts. Non-rejection
#' is not evidence of identification, as in a goodness-of-fit test. For
#' \code{nrow(A) > rank} the condition that \eqn{C} is non-negative outside the
#' anchor columns is not tested.
#'
#' \strong{Identification region (R = Q).} With as many covariates as bases,
#' every non-negative \eqn{\Xi} factorizes exactly, so the maximum of the
#' likelihood over the model is the row-wise non-negative least squares
#' \eqn{\hat\Xi} in the metric \eqn{AA^\top} (for one-hot \eqn{A}, the group
#' means). All exact factorizations are \eqn{X = \tilde X T},
#' \eqn{C = T^{-1} S} with \eqn{\tilde X} the column-normalized \eqn{\hat\Xi},
#' \eqn{S} its column sums and \eqn{T} the remixing above (columns summing to
#' 1, \eqn{\tilde X T \ge 0}, \eqn{T^{-1} \ge 0}); \eqn{T = I} is the
#' representative in which every covariate has its own basis. The bases are
#' labelled by matching the columns of \eqn{X} to those of \eqn{\tilde X} (so
#' that relabellings are excluded). The region is explored from \eqn{T = I}
#' by \code{region.n} random proposals, each a small move of a factorization
#' already found that is kept only if it satisfies these conditions; each end
#' of the range of an entry of \eqn{C} is then refined by Nelder--Mead (with a
#' penalty for negative entries) from the \code{region.nstart} most extreme
#' factorizations. The ranges of \eqn{X} come from the factorizations found.
#' Both are inner approximations: when the region is thin (for example, when
#' only some columns of \eqn{T} can move), the ends can be found noticeably
#' inside, and a larger \code{region.nstart} helps.
#' \code{identified} is \code{TRUE} when every range of \eqn{C} is shorter than
#' \code{1e-6} times the largest column sum.
#'
#' \strong{Intervals.} Gaussian errors with one variance common to all entries
#' of \eqn{Y} (the working model whose maximum is nmfkc's least squares),
#' profiled over \eqn{X}, the other entries of \eqn{C} and the variance. The
#' profile is flat over the identification region, so the interval contains
#' it. The \eqn{\chi^2_1} cut is meant to cover each point of the region with
#' probability about \code{ci.level}, as for a point of an identified set
#' (Imbens and Manski, 2004; Chernozhukov, Hong and Tamer, 2007); the true
#' coefficient is one such point. The probability of covering the whole
#' region at once can be lower: the ends of the region are where an entry of
#' \eqn{X} or \eqn{C} reaches 0, so they are minima and maxima rather than
#' smooth functions of \eqn{\Xi} and need not be asymptotically normal, in
#' particular when the width of the region is comparable to its standard
#' error or when an exact zero in \eqn{\hat\Xi} pins a column. The coverage
#' has not yet been studied by simulation in this setting, so read the
#' intervals as approximate. The profile is computed by continuation from the
#' ends of the region with local re-fits, so the intervals may be slightly
#' too short.
#'
#' @return An object of class \code{"nmfkc.ident"}:
#' \item{algebraic}{List: \code{verdict} (text), \code{sufficient} (two-sided
#'   anchors / peeling), \code{A.full.rank}, \code{cone.trivial} (no remixing
#'   moves the fit to first order), \code{pairs} (data frame \code{basis},
#'   \code{from}, \code{a.min}, \code{a.max}: column \code{basis} of \eqn{X}
#'   can become \eqn{(1-a)x_{basis} + a x_{from}}), \code{X.lower},
#'   \code{X.upper}, \code{C.lower}, \code{C.upper} (inner bounds of the
#'   identification region, containing the fit), \code{width.X},
#'   \code{width.C} (largest width; for \eqn{C} relative to its column sums,
#'   which do not change), \code{moved} (number of directions that moved the
#'   fit, by kind), \code{zero.tol}, \code{move.tol}, \code{n.dir}.}
#' \item{test}{\code{NULL} without \code{A}; otherwise a list: \code{pi}, \code{rho}, \code{alpha}, \code{reject}
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
#' \item{C.ci.lower, C.ci.upper}{Profile-likelihood interval of \eqn{C}, which
#'   contains the identification region (\code{NULL} unless \code{ci = TRUE}).}
#' \item{ci.level, rank, call}{As named; \code{dims} holds \code{P}, \code{N}, \code{R}.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.inference}},
#'   \code{\link{summary.nmfkc.ident}}, \code{\link{plot.nmfkc.ident}}
#' @references
#' Imbens, G. W. and Manski, C. F. (2004). Confidence intervals for partially
#' identified parameters. \emph{Econometrica}, 72(6), 1845--1857.
#' \doi{10.1111/j.1468-0262.2004.00555.x}
#'
#' Chernozhukov, V., Hong, H. and Tamer, E. (2007). Estimation and confidence
#' regions for parameter sets in econometric models. \emph{Econometrica},
#' 75(5), 1243--1284. \doi{10.1111/j.1468-0262.2007.00794.x}
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
#'
#' ## without covariates: the algebraic check only
#' fit0 <- nmfkc(t(as.matrix(iris[, 1:4])), rank = 2, epsilon = 1e-8, verbose = FALSE)
#' nmfkc.ident(fit0)
#' \donttest{
#' id2 <- nmfkc.ident(fit, Y, A, region.n = 5000, ci = TRUE)  # profile intervals
#' summary(id2)
#' plot(id2)
#' plot(id2, type = "X")
#' plot(id2, col = c(estimate = "blue", shade = NA))  # colours; NA leaves a part out
#' plot(id2, type = "X", cex = c(estimate = 0.5), lwd = c(range = 2))  # sizes, widths
#' }
nmfkc.ident <- function(object, Y = NULL, A = NULL, ...) {
  if (!inherits(object, "nmfkc")) stop("object must be of class 'nmfkc'")
  ex <- base::list(...)
  opt <- function(nm, def) if (!base::is.null(ex[[nm]])) ex[[nm]] else def
  if (!base::identical(object$X.restriction, "colSums"))
    stop("nmfkc.ident() needs a fit with X.restriction = \"colSums\"")
  alg.args <- base::list(zero.tol = opt("zero.tol", 1e-3), move.tol = opt("move.tol", 1e-3),
                         n.dir = opt("ident.n", 2000), seed = opt("seed", 123))
  ## Without A (a fit without covariates, or no A given), only the algebraic
  ## check on the fitted X and C is possible: the test needs A's residual df.
  if (base::is.null(A)) {
    if (!base::is.null(object$criterion$anchors$A.full.rank) && base::isFALSE(object$criterion$anchors$A.full.rank))
      alg.args$A.full.rank <- FALSE
    out <- base::list(test = NULL, Xi.hat = NULL, C.region.lower = NULL, C.region.upper = NULL,
                      X.region.lower = NULL, X.region.upper = NULL, identified = NA,
                      X = object$X, C = object$C, fit = object, repair = NULL,
                      C.ci.lower = NULL, C.ci.upper = NULL, ci.level = NA, rank = object$rank,
                      dims = base::c(P = base::nrow(object$X), N = NA, R = NA),
                      call = base::match.call())
    out$algebraic <- base::do.call(.nid.algebraic, base::c(base::list(object$X, object$C), alg.args))
    base::class(out) <- "nmfkc.ident"
    return(out)
  }
  if (base::is.null(Y)) stop("Y is needed when A is given (the test uses both)")
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
  if (R >= N) stop("nrow(A) must be smaller than ncol(Y)")
  if (base::qr(A)$rank < R) stop("A must have full row rank")
  if (Q > R) stop("the rank must not exceed nrow(A)")
  if (!base::all(base::dim(object$C) == base::c(Q, R)))
    stop("object does not match A (dim(object$C) must be rank x nrow(A))")

  ## ---- 1. test ----
  f <- .nid.ols(Y, A)
  crit.rho <- stats::qt(1 - alpha / Q^2, N - R)       # H(rho): Q(Q-1) off-diagonal entries and Q sizes
  crit.pi  <- stats::qt(1 - alpha / (2 * Q), N - R)   # membership contrasts (pi.max, pi < 1): 2Q per candidate
  dec <- if (pi0 >= 1) .nid.decide.rho(f, Q, rho0, crit.rho) else .nid.decide(f, Q, pi0, rho0, crit.pi)
  test <- base::list(pi = pi0, rho = rho0, alpha = alpha, reject = dec$reject,
                     pi.max = .nid.max.pi(f, Q, rho0, crit.rho, crit.pi),
                     rho.max = .nid.max.rho(f, Q, crit.rho),
                     pi.point = .nid.pi.point(f, Q, rho0)$value,
                     crit = if (pi0 >= 1) crit.rho else crit.pi,
                     candidate = if (dec$reject) NULL else
                       base::list(rows = dec$J, columns = dec$L, row.names = base::rownames(Y)[dec$J],
                                  column.names = base::rownames(A)[dec$L]))

  out <- base::list(test = test, Xi.hat = NULL, C.region.lower = NULL, C.region.upper = NULL,
                    X.region.lower = NULL, X.region.upper = NULL, identified = NA,
                    X = object$X, C = object$C, fit = object, repair = NULL,
                    C.ci.lower = NULL, C.ci.upper = NULL, ci.level = ci.level, rank = Q,
                    dims = base::c(P = P, N = N, R = R),
                    call = base::match.call())
  ## the algebraic check comes last: it draws random directions after a set.seed()
  ## of its own, so the test, region, repair and intervals are as without it
  alg <- function(out) {
    out$algebraic <- base::do.call(.nid.algebraic, base::c(base::list(object$X, object$C, A = A), alg.args))
    base::class(out) <- "nmfkc.ident"
    out
  }
  if (R > Q || !region) return(alg(out))

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
  alg(out)
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

## t of every entry of Xi.hat (HC1): the off-diagonal entries of H(rho) are tested one by one
.nid.entry.t <- function(f) {
  se <- base::sqrt(base::sapply(base::seq_len(f$R), function(l)
    base::rowSums(base::sweep(f$E, 2, f$W[l, ], "*")^2)) * f$hc)
  se <- base::matrix(se, base::nrow(f$Xi))
  base::ifelse(se > 0, f$Xi / se, base::ifelse(f$Xi > 0, Inf, base::ifelse(f$Xi < 0, -Inf, 0)))
}

## one-sided t of the membership contrast (1 - pi) u_jq - pi sum_{q' != q} u_jq', rows j, columns L.
## The column sums s are held fixed in the standard error: with s random, the gradient in s carries the
## estimate Xi.hat, which shrinks the standard error exactly when Xi.hat is large (liberal for small P).
.nid.purity.t <- function(f, L, pi) {
  Q <- base::length(L); P <- base::nrow(f$Xi)
  sL <- f$s[L]; XL <- f$Xi[, L, drop = FALSE]; U <- base::sweep(XL, 2, sL, "/")
  WL <- f$W[L, , drop = FALSE]
  tt <- base::matrix(NA_real_, P, Q)
  for (q in base::seq_len(Q)) {
    wq <- base::ifelse(base::seq_len(Q) == q, 1 - pi, -pi)
    g <- base::as.vector(U %*% wq)
    psi <- base::sweep(f$E, 2, base::as.vector(base::crossprod(wq / sL, WL)), "*")
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

## Candidates for the decisions: the column sets in increasing order only. Reordering L permutes the
## row sets S_q and leaves the existence of a matching unchanged, and the identity order was the first
## one tried among the Q! orders, so the decision and the reported candidate are the same.
.nid.comb.L <- function(R, Q) base::t(utils::combn(R, Q))

.nid.decide <- function(f, Q, pi, rho, cval) {
  Ls <- .nid.comb.L(f$R, Q); stR <- .nid.size.t(f, base::seq_len(f$R), rho)
  for (k in base::seq_len(base::nrow(Ls))) {
    L <- Ls[k, ]; pt <- .nid.purity.t(f, L, pi); st <- stR[, L, drop = FALSE]
    S <- base::lapply(base::seq_len(Q), function(q) base::which(pt[, q] >= -cval & st[, q] >= -cval))
    J <- .nid.matching(S)
    if (!base::is.null(J)) return(base::list(reject = FALSE, L = L, J = J))
  }
  base::list(reject = TRUE, L = NULL, J = NULL)
}

## H(rho): a candidate survives when none of its Q(Q-1) off-diagonal entries is significantly positive and
## none of its Q sizes Xi[j, l_q] - rho s_{l_q} is significantly negative (one-sided t at cval each)
.nid.decide.rho <- function(f, Q, rho, cval, tE = .nid.entry.t(f)) {
  Ls <- .nid.comb.L(f$R, Q); P <- base::nrow(f$Xi); stR <- .nid.size.t(f, base::seq_len(f$R), rho)
  for (k in base::seq_len(base::nrow(Ls))) {
    L <- Ls[k, ]; st <- stR[, L, drop = FALSE]
    S <- base::lapply(base::seq_len(Q), function(q) {
      off <- if (Q > 1) base::apply(tE[, L[-q], drop = FALSE], 1, base::max) else base::rep(-Inf, P)
      base::which(off <= cval & st[, q] >= -cval)
    })
    J <- .nid.matching(S)
    if (!base::is.null(J)) return(base::list(reject = FALSE, L = L, J = J))
  }
  base::list(reject = TRUE, L = NULL, J = NULL)
}

## pi.max (reported): 1 when H(rho) is not rejected; otherwise the largest pi not rejected by the membership
## contrasts (bisection on [1/Q, 1)); -Inf if rejected at pi = 1/Q
.nid.max.pi <- function(f, Q, rho, cval.rho, cval.pi) {
  if (!.nid.decide.rho(f, Q, rho, cval.rho)$reject) return(1)
  lo <- 1 / Q
  if (.nid.decide(f, Q, lo, rho, cval.pi)$reject) return(-Inf)
  hi <- 1
  for (it in 1:30) {
    mid <- (lo + hi) / 2
    if (.nid.decide(f, Q, mid, rho, cval.pi)$reject) hi <- mid else lo <- mid
  }
  lo
}

## largest rho not rejected by H(rho) (bisection on [0, 1]); -Inf if rejected at rho = 0
.nid.max.rho <- function(f, Q, cval) {
  tE <- .nid.entry.t(f)
  if (.nid.decide.rho(f, Q, 0, cval, tE)$reject) return(-Inf)
  lo <- 0; hi <- 1
  if (!.nid.decide.rho(f, Q, hi, cval, tE)$reject) return(1)
  for (it in 1:30) {
    mid <- (lo + hi) / 2
    if (.nid.decide.rho(f, Q, mid, cval, tE)$reject) hi <- mid else lo <- mid
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

## put the bases of the refit in the object's order (idx[q] = fit's basis for basis q) and give
## them the object's names -- also when the order is already the same: the refit is named by
## nmfkc()'s default prefix ("Basis"), the object may carry other names (prefix, or set by hand)
.nid.reorder <- function(fit, idx, Y, A, object) {
  if (!base::identical(base::as.integer(idx), base::seq_along(idx))) {
    fit$X <- fit$X[, idx, drop = FALSE]; fit$C <- fit$C[idx, , drop = FALSE]; fit$B <- fit$B[idx, , drop = FALSE]
    .eps <- 1e-10
    fit$B.prob <- base::t(base::t(fit$B) / (base::colSums(fit$B) + .eps))
    fit$B.cluster <- base::apply(fit$B.prob, 2, base::which.max)
    fit$B.cluster[base::colSums(fit$B.prob) == 0] <- NA
    fit$X.prob <- fit$X / (base::rowSums(fit$X) + .eps)
    fit$X.cluster <- base::apply(fit$X.prob, 1, base::which.max)
  }
  base::dimnames(fit$X) <- base::dimnames(object$X); base::dimnames(fit$C) <- base::dimnames(object$C)
  base::rownames(fit$B) <- base::rownames(object$C)
  if (base::is.matrix(fit$B.prob)) base::rownames(fit$B.prob) <- base::rownames(object$C)
  if (base::is.matrix(fit$X.prob)) base::colnames(fit$X.prob) <- base::colnames(object$X)
  ## anchor purity (named by basis; nmfkc() computes it with the defaults, as here)
  if (!base::is.null(fit$criterion$anchors)) fit$criterion$anchors <- .nmfkc.anchor(fit$X, fit$C, A)
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

## ===================== algebraic check (any A, also none) =====================
## Every T with 1'T = 1', X T >= 0 and T^-1 C >= 0 gives the same fit; with A of
## full row rank the constraint is on C itself (on B = C when there is no A), so
## the check needs only the fitted X and C.  T = I + t D with D of zero column
## sums, parametrised by its Q(Q-1) off-diagonal entries d (D_qq = -sum_r D_rq).

## off-diagonal positions (r, q), column-major
.nid.alg.ij <- function(Q) base::which(base::diag(Q) == 0, arr.ind = TRUE)

.nid.alg.D <- function(d, Q) {
  D <- base::matrix(0, Q, Q); D[.nid.alg.ij(Q)] <- d; base::diag(D) <- -base::colSums(D); D
}

## entries at most zero.tol times their column maximum (X) or row maximum (C)
## count as zero: multiplicative updates leave 1e-10..1e-5 at the boundary
.nid.alg.snap <- function(X, C, zero.tol) {
  X[X <= zero.tol * base::rep(base::apply(X, 2, base::max), each = base::nrow(X))] <- 0
  C[C <= zero.tol * base::apply(C, 1, base::max)] <- 0
  base::list(X = X, C = C)
}

## first-order constraints on d, as rows of M with M d >= 0:
## (XD)_ij = sum_{q != j} X_iq D_qj >= 0 where X_ij = 0, and
## (DC)_kj = sum_{q != k} D_kq C_qj <= 0 where C_kj = 0
.nid.alg.M <- function(X, C) {
  ij <- .nid.alg.ij(base::ncol(X)); m <- base::nrow(ij)
  zx <- base::which(X == 0, arr.ind = TRUE); zc <- base::which(C == 0, arr.ind = TRUE)
  Mx <- base::matrix(0, base::nrow(zx), m); Mc <- base::matrix(0, base::nrow(zc), m)
  for (k in base::seq_len(m)) {
    s <- zx[, 2] == ij[k, 2]
    if (base::any(s)) Mx[s, k] <- X[base::cbind(zx[s, 1], ij[k, 1])]
    s <- zc[, 1] == ij[k, 1]
    if (base::any(s)) Mc[s, k] <- -C[base::cbind(base::rep(ij[k, 2], base::sum(s)), zc[s, 2])]
  }
  M <- base::rbind(Mx, Mc)
  nr <- base::sqrt(base::rowSums(M^2)); M <- M[nr > 0, , drop = FALSE] / nr[nr > 0]
  base::unique(base::round(M, 12))
}

## The cone {d : M d >= 0}.  By Farkas' lemma there is d in it with s d_k > 0
## iff M'y = -s e_k has no solution y >= 0; the residual of the non-negative
## least squares min ||M'y + s e_k|| is then such a d (the KKT conditions give
## M r >= 0 and s r_k = ||r||^2 > 0), so each of the 2 Q(Q-1) problems both
## decides and, when the answer is yes, returns a direction.
.nid.alg.cone <- function(M, m, tol = 1e-8) {
  feas <- base::matrix(FALSE, m, 2, dimnames = base::list(NULL, c("plus", "minus")))
  rays <- base::list(); err <- FALSE
  Mt <- base::t(M)
  for (k in base::seq_len(m)) for (si in 1:2) {
    s <- c(1, -1)[si]; b <- base::numeric(m); b[k] <- -s
    r <- if (base::nrow(M) == 0) -b else base::tryCatch({
      y <- .nid.nnls(Mt, b); base::as.vector(Mt %*% y) - b
    }, error = function(e) NULL)
    if (base::is.null(r)) { err <- TRUE; next }
    nr <- base::sqrt(base::sum(r^2))
    if (nr > tol && s * r[k] > 0 && (base::nrow(M) == 0 || base::all(M %*% r >= -1e-8 * nr))) {
      feas[k, si] <- TRUE; rays[[base::length(rays) + 1L]] <- r / nr
    }
  }
  base::list(trivial = !base::any(feas) && !err, feasible = feas, rays = rays, error = err)
}

## largest t >= 0 (first exit, capped) with X (I + tD) >= 0 and (I + tD)^-1 C >= 0:
## the X side is linear in t (exact bound); the C side is checked on a grid on
## [0, t_X] and the first exit bisected
.nid.alg.line <- function(X, C, D, t.cap = 10, grid = 20, bis = 30, ftol = 1e-12) {
  XD <- X %*% D
  neg <- XD < -ftol
  tX <- base::min(if (base::any(neg)) base::min(X[neg] / -XD[neg]) else t.cap, t.cap)
  if (tX <= 0) return(0)
  DC <- D %*% C
  if (base::any(C == 0 & DC > ftol * base::max(1, base::abs(DC)))) return(0)   # blocked at first order
  Q <- base::nrow(D); I <- base::diag(Q); cs <- base::max(C)
  ok <- function(t) {
    Tm <- I + t * D
    base::det(Tm) > 1e-12 && base::all(base::solve(Tm, C) >= -ftol * cs)
  }
  ts <- tX * base::seq_len(grid) / grid
  for (g in base::seq_along(ts)) if (!ok(ts[g])) {
    lo <- if (g == 1) 0 else ts[g - 1]; hi <- ts[g]
    for (i in base::seq_len(bis)) { mid <- (lo + hi) / 2; if (ok(mid)) lo <- mid else hi <- mid }
    return(lo)
  }
  tX
}

.nid.algebraic <- function(X, C, A = NULL, A.full.rank = NULL, zero.tol = 1e-3, move.tol = 1e-3,
                           n.dir = 2000, seed = 123, t.cap = 10) {
  X <- base::as.matrix(X); C <- base::as.matrix(C); Q <- base::ncol(X)
  if (base::is.null(A.full.rank)) A.full.rank <- base::is.null(A) || base::qr(A)$rank == base::nrow(A)
  suff <- base::isTRUE(.ident.verdict(X, C, A = A)$identified.approx) && base::isTRUE(A.full.rank)
  res <- base::list(verdict = NA_character_, sufficient = suff, A.full.rank = A.full.rank,
                    cone.trivial = NA, pairs = NULL, X.lower = X, X.upper = X,
                    C.lower = C, C.upper = C, width.X = 0, width.C = 0, moved = NULL,
                    zero.tol = zero.tol, move.tol = move.tol, n.dir = n.dir)
  if (!base::isTRUE(A.full.rank)) {
    res$verdict <- "not identified: A is not of full row rank, so C is not unique even given X"
    return(res)
  }
  if (Q == 1L) { res$verdict <- "identified (rank 1)"; res$cone.trivial <- TRUE; return(res) }
  st <- .nmfkc.rng.save(seed); base::on.exit(.nmfkc.rng.restore(st), add = TRUE)
  base::set.seed(seed)
  sn <- .nid.alg.snap(X, C, zero.tol); X0 <- sn$X; C0 <- sn$C
  ij <- .nid.alg.ij(Q); m <- base::nrow(ij)
  cone <- .nid.alg.cone(.nid.alg.M(X0, C0), m)
  ## directions: the pairwise remixings (unnormalised, so t is the share a), the
  ## cone directions from the Farkas certificates, and n.dir random ones
  dirs <- base::list(); type <- base::character(0)
  for (k in base::seq_len(m)) for (s in c(1, -1)) {
    E <- base::matrix(0, Q, Q); E[ij[k, 1], ij[k, 2]] <- s; E[ij[k, 2], ij[k, 2]] <- -s
    dirs[[base::length(dirs) + 1L]] <- E; type <- c(type, "pair")
  }
  for (r in cone$rays) { dirs[[base::length(dirs) + 1L]] <- .nid.alg.D(r, Q); type <- c(type, "cone") }
  for (i in base::seq_len(n.dir)) {
    D <- base::matrix(stats::rnorm(Q * Q), Q, Q); D <- base::sweep(D, 2, base::colMeans(D))
    dirs[[base::length(dirs) + 1L]] <- D / base::sqrt(base::sum(D^2)); type <- c(type, "random")
  }
  tt <- base::vapply(dirs, function(D) .nid.alg.line(X0, C0, D, t.cap), 0)
  cs <- base::pmax(base::colSums(C0), 1e-300)
  Xlo <- Xhi <- X0; Clo <- Chi <- C0; dX <- dC <- base::numeric(base::length(dirs))
  for (i in base::which(tt > 0)) {
    Tm <- base::diag(Q) + tt[i] * dirs[[i]]
    Xn <- X0 %*% Tm; Cn <- base::solve(Tm, C0)
    dX[i] <- base::max(base::abs(Xn - X0)); dC[i] <- base::max(base::sweep(base::abs(Cn - C0), 2, cs, "/"))
    if (base::max(dX[i], dC[i]) > move.tol) {
      Xlo <- base::pmin(Xlo, Xn); Xhi <- base::pmax(Xhi, Xn); Clo <- base::pmin(Clo, Cn); Chi <- base::pmax(Chi, Cn)
    }
  }
  moved <- base::pmax(dX, dC) > move.tol
  ## pairwise table: column q of X becomes (1 - a) x_q + a x_r, a in [a.min, a.max]
  pt <- tt[type == "pair"]
  bn <- if (!base::is.null(base::colnames(X))) base::colnames(X) else base::paste0("Basis", base::seq_len(Q))
  res$pairs <- base::data.frame(basis = bn[ij[, 2]], from = bn[ij[, 1]],
                                a.min = -pt[c(FALSE, TRUE)], a.max = pt[c(TRUE, FALSE)],
                                stringsAsFactors = FALSE)
  res$cone.trivial <- cone$trivial
  res$X.lower <- base::pmin(Xlo, X); res$X.upper <- base::pmax(Xhi, X)
  res$C.lower <- base::pmin(Clo, C); res$C.upper <- base::pmax(Chi, C)
  base::dimnames(res$X.lower) <- base::dimnames(res$X.upper) <- base::dimnames(X)
  base::dimnames(res$C.lower) <- base::dimnames(res$C.upper) <- base::dimnames(C)
  res$width.X <- base::max(Xhi - Xlo)
  res$width.C <- base::max(base::sweep(Chi - Clo, 2, cs, "/"))
  res$moved <- base::c(pair = base::sum(moved[type == "pair"]), cone = base::sum(moved[type == "cone"]),
                       random = base::sum(moved[type == "random"]))
  res$verdict <- if (suff) "identified (two-sided anchors / peeling, a sufficient condition)"
    else if (base::any(moved)) "not identified: another factorization gives the same fit"
    else if (base::isTRUE(cone$trivial)) "locally identified: no remixing moves the fit (fits far away are not examined)"
    else "undetermined: some remixings are possible but move X and C by less than move.tol"
  res
}

## ===================== methods =====================

#' @export
print.nmfkc.ident <- function(x, ...) {
  base::print(summary(x), ...)
  base::invisible(x)
}

#' @title Summary of an nmfkc.ident object
#' @description The anchor test, the identification region, the repaired fit
#'   and, when computed, the profile-likelihood intervals, with a coefficient
#'   table for \eqn{C}.
#' @param object An object of class \code{"nmfkc.ident"}.
#' @param ... Not used.
#' @return An object of class \code{"summary.nmfkc.ident"}: \code{test},
#'   \code{dims}, \code{rank}, \code{identified}, \code{repair},
#'   \code{algebraic} (as in \code{\link{nmfkc.ident}}), \code{X.dim},
#'   \code{C.dim},
#'   \code{region.share} (largest share of an interval taken by the
#'   identification region), \code{ci.level}, and \code{coefficients}, a data
#'   frame with columns \code{Basis}, \code{Covariate}, \code{Estimate}
#'   (repaired fit), \code{region_low}, \code{region_high}, and with intervals
#'   \code{CI_low}, \code{CI_high} and \code{nonzero} (the interval excludes 0:
#'   positive even allowing for the non-uniqueness).
#' @seealso \code{\link{nmfkc.ident}}, \code{\link{plot.nmfkc.ident}}
#' @export
summary.nmfkc.ident <- function(object, ...) {
  x <- object
  co <- NULL; share <- NA_real_
  if (!base::is.null(x$Xi.hat)) {
    C <- x$C
    co <- base::data.frame(Basis = base::rep(base::rownames(C), base::ncol(C)),
                           Covariate = base::rep(base::colnames(C), each = base::nrow(C)),
                           Estimate = base::as.vector(C),
                           region_low = base::as.vector(x$C.region.lower),
                           region_high = base::as.vector(x$C.region.upper),
                           stringsAsFactors = FALSE)
    if (!base::is.null(x$C.ci.lower)) {
      co$CI_low <- base::as.vector(x$C.ci.lower); co$CI_high <- base::as.vector(x$C.ci.upper)
      co$nonzero <- co$CI_low > 0
      w <- co$CI_high - co$CI_low
      share <- base::max(base::ifelse(w > 0, (co$region_high - co$region_low) / w, 0))
    }
  }
  base::structure(base::list(test = x$test, dims = x$dims, rank = x$rank, identified = x$identified,
                             repair = x$repair, region.share = share, ci.level = x$ci.level,
                             coefficients = co, call = x$call, algebraic = x$algebraic,
                             X.dim = base::dim(x$X), C.dim = base::dim(x$C)),
                  class = "summary.nmfkc.ident")
}

## the algebraic-check block of print.summary.nmfkc.ident()
.nid.print.algebraic <- function(a, digits = 4) {
  if (base::is.null(a)) return(base::invisible(NULL))
  base::cat(base::sprintf("\nAlgebraic check on the fitted X and C  [remixings XT, T^-1 C with 1'T = 1'; entries <= %g x their column (X) / row (C) maximum count as 0]\n",
                          a$zero.tol))
  base::cat("  ", a$verdict, "\n", sep = "")
  yn <- function(v) if (base::isTRUE(v)) "yes" else if (base::isFALSE(v)) "no" else "not computed"
  base::cat(base::sprintf("  sufficient condition (two-sided anchors / peeling): %s;  first-order cone of remixings is {0}: %s\n",
                          if (a$sufficient) "met" else "not met", yn(a$cone.trivial)))
  if (!base::is.null(a$moved))
    base::cat(base::sprintf("  remixings that move the fit (> %g): %d pairwise, %d cone, %d of %d random\n",
                            a$move.tol, a$moved[["pair"]], a$moved[["cone"]], a$moved[["random"]], a$n.dir))
  if (a$width.X > 0 || a$width.C > 0)
    base::cat(base::sprintf("  widths found: X up to %.4f;  C up to %.1f%% of a column sum (inner bounds: the region is at least this wide)\n",
                            a$width.X, 100 * a$width.C))
  pr <- a$pairs
  if (!base::is.null(pr) && base::nrow(pr) <= 30 && base::any(pr$a.max - pr$a.min > 0)) {
    base::cat("  pairwise remixings, column 'basis' -> (1 - a) x_basis + a x_from  (a = 0: no change):\n")
    pr$a.min <- base::signif(pr$a.min, digits); pr$a.max <- base::signif(pr$a.max, digits)
    l <- utils::capture.output(base::print(pr, row.names = FALSE))
    base::cat(base::paste0("    ", l), sep = "\n")
  }
  base::invisible(NULL)
}

#' @export
print.summary.nmfkc.ident <- function(x, digits = 4, ...) {
  t <- x$test; d <- x$dims
  fmt <- function(v, k = 3) if (base::is.finite(v)) base::formatC(v, digits = k, format = "f") else "none (rejected at 0)"
  base::cat("Identifiability of an nmfkc fit\n")
  if (base::is.null(t)) {
    ## no A: the algebraic check only
    if (!base::is.null(x$X.dim))
      base::cat(base::sprintf("  X: %d x %d,  C: %d x %d,  rank %d  (no A given: the test, region and intervals need Y and A)\n",
                              x$X.dim[1], x$X.dim[2], x$C.dim[1], x$C.dim[2], x$rank))
    .nid.print.algebraic(x$algebraic, digits)
    return(base::invisible(x))
  }
  if (!base::is.null(d))
    base::cat(base::sprintf("  Y: %d x %d,  A: %d covariates,  rank %d\n", d[["P"]], d[["N"]], d[["R"]], x$rank))
  .nid.print.algebraic(x$algebraic, digits)
  base::cat(if (t$pi >= 1)
    "\nAnchor test H(rho)  [OLS of Xi = XC, HC1, entries tested one by one, intersection-union over candidates]\n" else
    "\nAnchor test H(pi, rho)  [OLS of Xi = XC, HC1, membership contrasts, intersection-union over candidates]\n")
  base::cat(base::sprintf("  at pi = %g, rho = %.3g:  %s   (alpha = %g, critical t = %.3f)\n", t$pi, t$rho,
                          if (t$reject) "rejected" else "not rejected", t$alpha, t$crit))
  base::cat(base::sprintf("  pi.max  = %s   (largest membership probability not rejected; point estimate %s)\n",
                          fmt(t$pi.max), fmt(t$pi.point)))
  rm <- if (base::is.finite(t$rho.max) && !base::is.null(d))
    base::sprintf("%s  (P * rho.max = %.2f; 1 = a uniform basis)", fmt(t$rho.max), d[["P"]] * t$rho.max) else fmt(t$rho.max)
  base::cat(base::sprintf("  rho.max = %s\n", rm))
  if (!base::is.null(t$candidate)) {
    cr <- if (!base::is.null(t$candidate$row.names)) t$candidate$row.names else t$candidate$rows
    cc <- if (!base::is.null(t$candidate$column.names)) t$candidate$column.names else t$candidate$columns
    base::cat(base::sprintf("  surviving candidate: rows %s  x  covariates %s\n",
                            base::paste(cr, collapse = ", "), base::paste(cc, collapse = ", ")))
  }
  base::cat("  (non-rejection is not evidence of uniqueness, as in a goodness-of-fit test)\n")
  if (base::is.null(x$coefficients)) {
    base::cat("\nIdentification region: not computed (needs nrow(A) == rank and region = TRUE)\n")
    return(base::invisible(x))
  }
  base::cat("\nIdentification region of C\n")
  base::cat(base::sprintf("  %s\n", if (base::isTRUE(x$identified)) "a single point: the factorization is unique"
                          else "has width: C is not unique; read it from the region and its intervals"))
  if (base::is.finite(x$region.share))
    base::cat(base::sprintf("  largest share of an interval taken by the region: %.0f%%\n", 100 * x$region.share))
  r <- x$repair
  base::cat(base::sprintf("  repaired fit (start = %s): %d iterations, |XC - Xi.hat| / max(Xi.hat) = %.1e\n",
                          r$start, r$iter, r$gap))
  co <- x$coefficients
  base::cat("\nCoefficients (Estimate = repaired fit; region = identification region")
  if (!base::is.null(co$CI_low)) base::cat(base::sprintf(";\n  CI = %g%% profile-likelihood interval containing the region (approximate), one error variance for all of Y",
                                                         100 * x$ci.level))
  base::cat(")\n")
  show <- co
  tol <- 1e-6 * base::max(base::abs(base::c(co$region_low, co$region_high)), 1)
  near <- show$Estimate >= show$region_low - tol & show$Estimate <= show$region_high + tol
  show$Estimate[near] <- base::pmin(base::pmax(show$Estimate[near], show$region_low[near]), show$region_high[near])
  if (!base::is.null(show$nonzero)) { show$nonzero <- base::ifelse(show$nonzero, "*", ""); base::names(show)[base::names(show) == "nonzero"] <- "" }
  num <- base::vapply(show, base::is.numeric, TRUE)
  show[num] <- base::lapply(show[num], function(v) base::signif(v, digits))
  base::print(show, row.names = FALSE)
  if (!base::is.null(co$nonzero))
    base::cat("  *: the interval excludes 0 -- positive even allowing for the non-uniqueness\n")
  base::invisible(x)
}

#' @title Plot the identification region and intervals of an nmfkc.ident object
#' @description
#' Draws how far \eqn{X} and \eqn{C} (\eqn{B} for a fit without covariates)
#' can move without changing the fit. The ranges are the identification region
#' when it was computed (\code{Y} and \code{A} given, \code{nrow(A) == rank});
#' otherwise they are the ranges found by the algebraic check, which are inner
#' bounds (the region is at least this wide). Each line is the range of one
#' entry on its own: the entries do not reach their ends together, since one
#' remixing moves the whole matrix.
#' \itemize{
#'   \item \code{type = "C"}: one panel per entry of \eqn{C}, bases in rows and
#'     covariates in columns: the range (black; a dot when it is a single
#'     point), the profile-likelihood interval if computed (gray) and the
#'     estimate (red cross). Panels whose interval excludes 0 are shaded.
#'   \item \code{type = "X"}: one panel per basis, the range of every row of
#'     \eqn{X} (black) and the estimate (red cross). The row names are
#'     written once, left of the first panel, so that many bases fit side by
#'     side on a page; when they do not fit the height, every k-th row is
#'     named. Beyond 40 rows the row numbers are written at a few
#'     \code{pretty()} positions instead; \code{rows} draws a subset.
#'   \item \code{type = "B"}: one strip per basis, the columns of \eqn{C} along
#'     the horizontal axis (for a fit without covariates, \eqn{C = B} has one
#'     column per sample), the range of every entry (black line) and the
#'     estimate (small red point). Suited to many columns. When \code{at}
#'     is given, the horizontal axis is left undrawn so that
#'     \code{axis(1, ...)} can label it (e.g. with dates).
#'   \item \code{type = "C.prob"}, \code{"B.prob"}: the same with every
#'     column of \eqn{C} divided by its sum, i.e. the share of each basis in
#'     a covariate, or the membership probability of a sample (\code{B.prob}
#'     of the fit) without covariates. The column sums do not change under a
#'     remixing, so these ranges are exact rescalings and lie between 0 and 1; the
#'     profile-likelihood intervals, which are for \eqn{C} itself, are not
#'     drawn. No \code{"X.prob"} is offered: the row sums of \eqn{X} change
#'     under a remixing, so its row shares do not rescale.
#' }
#' The estimate is the repaired fit when the region was computed, otherwise
#' the fitted matrices. The colours named above are the defaults; \code{col}
#' changes them, and \code{cex} and \code{lwd} the sizes of the marks and the
#' widths of the lines.
#' @param x An object of class \code{"nmfkc.ident"}.
#' @param type \code{"C"}, \code{"X"}, \code{"B"}, \code{"C.prob"} or
#'   \code{"B.prob"}. Default: \code{"C"} when \eqn{C} has at most 12
#'   columns, otherwise \code{"B"}.
#' @param bases Bases to draw, by index or name (columns of \eqn{X}, rows of
#'   \eqn{C}). Default: all of them for \code{"C"} and \code{"X"} (a basis
#'   that cannot move shows as dots); for \code{"B"} the bases whose range
#'   exceeds \code{move.tol} of a column sum (all of them if none does).
#' @param at For \code{type = "B"}: horizontal positions of the columns
#'   (e.g. times). Default \code{1:ncol(C)}.
#' @param panel For \code{type = "B"}: optional \code{function(q)} called in
#'   each strip after drawing (e.g. to add \code{abline(v = ...)}); \code{q} is
#'   the index of the basis.
#' @param col Colours, as a named vector with any of \code{range} (the
#'   ranges; default \code{"black"}), \code{estimate} (the estimate marks;
#'   \code{"red3"}), \code{interval} (the profile-likelihood intervals of
#'   \code{type = "C"}; \code{"gray55"}) and \code{shade} (the panels whose
#'   interval excludes 0; \code{"gray90"}), e.g.
#'   \code{col = c(estimate = "blue")}. Elements not given keep their
#'   default, and \code{NA} leaves that part out. The legend below the plot
#'   names a colour given by name (e.g. "blue cross: estimate") and
#'   otherwise only the mark (e.g. "cross: estimate" for \code{"#1f77b4"}).
#' @param cex Sizes of the marks, as a named vector with any of
#'   \code{estimate} (the estimate marks) and \code{range} (the dot drawn
#'   when a range is a single point; types \code{"C"} and \code{"X"}), e.g.
#'   \code{cex = c(estimate = 0.5)}. Defaults: \code{estimate} 1.6
#'   (\code{"C"}), 0.9 (\code{"X"}), 0.5 (\code{"B"}); \code{range} 0.8
#'   (\code{"C"}), 0.45 (\code{"X"}).
#' @param lwd Line widths, as a named vector with any of \code{range} (the
#'   range lines), \code{interval} (the profile-likelihood intervals of
#'   \code{type = "C"}) and \code{estimate} (the strokes of the estimate
#'   crosses; types \code{"C"} and \code{"X"}), e.g.
#'   \code{lwd = c(range = 2)}. Defaults: \code{range} 6 (\code{"C"}), 4
#'   (\code{"X"}), \code{par("lwd")} (\code{"B"}); \code{interval} 2.5;
#'   \code{estimate} 2 (\code{"C"}), 1.5 (\code{"X"}). For \code{col},
#'   \code{cex} and \code{lwd}, elements not given keep their default and
#'   elements that do not apply to the type are ignored.
#' @param rows For \code{type = "X"}: rows of \eqn{X} to draw, by index or
#'   name, in the order given from the top (e.g. \code{rows = 1:30}).
#'   Default: all of them.
#' @param ... Not used.
#' @return \code{x}, invisibly. For \code{type = "B"} the graphical
#'   parameters are left as set, with the last strip active, so that
#'   \code{axis(1, ...)} can relabel the horizontal axis afterwards.
#' @seealso \code{\link{nmfkc.ident}}, \code{\link{summary.nmfkc.ident}}
#' @export
plot.nmfkc.ident <- function(x, type = NULL, bases = NULL, at = NULL, panel = NULL, col = NULL,
                             cex = NULL, lwd = NULL, rows = NULL, ...) {
  pc <- .nid.plot.col(col); cc <- pc$col; cw <- pc$word
  ucex <- .nid.plot.size(cex, "cex", base::c("range", "estimate"))
  ulwd <- .nid.plot.size(lwd, "lwd", base::c("range", "interval", "estimate"))
  ## legend words: the colour's name if it has one, otherwise the mark alone
  nm.rng <- function(noun) if (base::nzchar(cw[["range"]])) cw[["range"]] else noun
  nm.est <- function(noun) base::paste(base::c(if (base::nzchar(cw[["estimate"]])) cw[["estimate"]], noun), collapse = " ")
  legtxt <- function(...) base::paste(base::c(...), collapse = ";  ")
  reg <- !base::is.null(x$Xi.hat)
  if (reg) {
    Xl <- x$X.region.lower; Xu <- x$X.region.upper; Cl <- x$C.region.lower; Cu <- x$C.region.upper
  } else if (!base::is.null(x$algebraic)) {
    a <- x$algebraic
    Xl <- a$X.lower; Xu <- a$X.upper; Cl <- a$C.lower; Cu <- a$C.upper
  } else stop("no identification region and no algebraic check in x")
  what <- if (reg) "identification region" else "range found by the algebraic check (at least this wide)"
  Q <- base::nrow(x$C); R <- base::ncol(x$C)
  if (base::is.null(type)) type <- if (R <= 12) "C" else "B"
  type <- base::match.arg(type, c("C", "X", "B", "C.prob", "B.prob"))
  prob <- type %in% c("C.prob", "B.prob"); type <- base::sub(".prob", "", type, fixed = TRUE)
  ## mark sizes and line widths of this type (NA: not drawn in it), with the caller's values
  sz <- base::switch(type,
    C = base::list(cex = base::c(range = 0.8, estimate = 1.6), lwd = base::c(range = 6, interval = 2.5, estimate = 2)),
    X = base::list(cex = base::c(range = 0.45, estimate = 0.9), lwd = base::c(range = 4, interval = NA, estimate = 1.5)),
    B = base::list(cex = base::c(range = NA, estimate = 0.5),
                   lwd = base::c(range = graphics::par("lwd"), interval = NA, estimate = NA)))
  cx <- sz$cex; cx[base::names(ucex)] <- ucex
  lw <- sz$lwd; lw[base::names(ulwd)] <- ulwd
  bn <- base::rownames(x$C); if (base::is.null(bn)) bn <- base::paste0("Basis", base::seq_len(Q))
  sel <- if (base::is.null(bases)) {
    if (type == "B") {
      cs <- base::pmax(base::colSums(x$C), 1e-300)
      mt <- if (!base::is.null(x$algebraic$move.tol)) x$algebraic$move.tol else 1e-3
      mv <- base::which(base::apply(base::sweep(Cu - Cl, 2, cs, "/"), 1, base::max) > mt)
      if (base::length(mv)) mv else base::seq_len(Q)
    } else base::seq_len(Q)
  } else if (base::is.character(bases)) base::match(bases, bn) else bases
  if (!base::length(sel) || base::anyNA(sel) || base::any(sel < 1 | sel > Q))
    stop("bases must name or index bases of the fit")
  sel <- base::as.integer(sel)
  if (type == "X") {   # rows of X to draw (checked before the graphical parameters are touched)
    rn <- base::rownames(x$X)
    rr <- if (base::is.null(rows)) base::seq_len(base::nrow(x$X))
          else if (base::is.character(rows)) base::match(rows, rn) else rows
    if (!base::length(rr) || base::anyNA(rr) || base::any(rr < 1 | rr > base::nrow(x$X)))
      stop("rows must name or index rows of X")
    rr <- base::as.integer(rr)
  }
  ## shares of the column sums of C (membership probabilities for B): the column sums
  ## do not change under a remixing (1'T^-1 = 1'), so the ranges scale exactly
  Cv <- x$C; scl <- ""
  if (prob) {
    cs0 <- base::pmax(base::colSums(x$C), 1e-300)
    Cv <- base::sweep(x$C, 2, cs0, "/"); Cl <- base::sweep(Cl, 2, cs0, "/"); Cu <- base::sweep(Cu, 2, cs0, "/")
    scl <- ", as shares of the column sums"
  }
  op <- graphics::par(no.readonly = TRUE)
  if (type != "B") base::on.exit(graphics::par(op))
  if (type == "X") {
    P <- base::length(rr); yy <- base::rev(base::seq_len(P)); nq <- base::length(sel)
    rl <- if (base::is.null(rn)) base::as.character(rr) else rn[rr]
    ## up to 40 rows the row names are written once, left of the first panel (thinned if they
    ## do not fit its height); beyond, the row numbers at pretty() positions
    lab <- P <= 40
    graphics::layout(base::matrix(base::seq_len(nq), 1),
                     widths = base::c(if (lab) 1.6 else 1.25, base::rep(1, nq - 1)))
    graphics::par(oma = base::c(1.5, 0, 0, 0), mgp = base::c(1.5, 0.4, 0), tcl = -0.25)
    tk <- NULL   # rows with a tick (and a label in the first panel), set in the first panel
    for (i in base::seq_len(nq)) {
      q <- sel[i]
      graphics::par(mar = base::c(2.5, if (i == 1) (if (lab) 4 else 2.5) else 0.5, 2, 0.5))
      graphics::plot(NA, xlim = base::c(0, base::max(Xu[rr, sel]) * 1.05), ylim = base::c(0.5, P + 0.5),
                     yaxt = "n", xlab = "", ylab = "", main = base::colnames(x$X)[q], cex.main = 0.9,
                     cex.axis = 0.75)
      if (base::is.null(tk)) tk <- if (lab) {
        ## every k-th row when the names would overlap: rows are par("pin")[2] / P inches
        ## apart, a name needs about 1.4 times the height of "M"
        k <- base::max(1, base::ceiling(1.4 * graphics::strheight("M", units = "inches", cex = 0.6) /
                                          (graphics::par("pin")[2] / P)))
        base::seq(1, P, by = k)
      } else {
        pr <- base::pretty(base::c(1, P)); base::unique(base::c(1, pr[pr >= 1 & pr <= P]))
      }
      graphics::axis(2, at = yy[tk], labels = if (i > 1) FALSE else if (lab) rl[tk] else base::as.character(rr[tk]),
                     las = 1, cex.axis = 0.6)
      w <- Xu[rr, q] - Xl[rr, q]; pt <- w < 1e-6 * base::max(Xu)
      if (!base::is.na(cc[["range"]])) {
        graphics::segments(Xl[rr, q][!pt], yy[!pt], Xu[rr, q][!pt], yy[!pt], lwd = lw[["range"]], lend = 1,
                           col = cc[["range"]])
        graphics::points(Xl[rr, q][pt], yy[pt], pch = 19, cex = cx[["range"]], col = cc[["range"]])
      }
      if (!base::is.na(cc[["estimate"]]))
        graphics::points(x$X[rr, q], yy, pch = 4, col = cc[["estimate"]], cex = cx[["estimate"]], lwd = lw[["estimate"]])
    }
    graphics::mtext(legtxt(if (!base::is.na(cc[["range"]]))
                             base::paste0(nm.rng("lines"), ": ", what, " of X (a dot: a single point)"),
                           if (!base::is.na(cc[["estimate"]])) base::paste0(nm.est("cross"), ": estimate")),
                    side = 1, outer = TRUE, cex = 0.75)
    return(base::invisible(x))
  }
  if (type == "B") {
    draw.axis <- base::is.null(at)   # with `at` given, the caller labels the axis
    if (base::is.null(at)) at <- base::seq_len(R)
    if (base::length(at) != R) stop("at must have one position per column of C")
    graphics::par(mfrow = base::c(base::length(sel), 1), mar = base::c(0.5, 4, 0.5, 1), oma = base::c(4, 0, 0, 0),
                  cex.lab = 0.8)
    for (i in base::seq_along(sel)) {
      q <- sel[i]
      graphics::plot(at, Cv[q, ], type = "n", xaxt = "n", las = 1, ylab = bn[q], xlab = "",
                     ylim = if (prob) base::c(0, 1) else base::range(Cl[q, ], Cu[q, ], Cv[q, ]))
      if (!base::is.na(cc[["range"]])) graphics::segments(at, Cl[q, ], at, Cu[q, ], col = cc[["range"]], lwd = lw[["range"]])
      if (!base::is.na(cc[["estimate"]])) graphics::points(at, Cv[q, ], pch = 16, cex = cx[["estimate"]], col = cc[["estimate"]])
      if (!base::is.null(panel)) panel(q)
      if (i == base::length(sel) && draw.axis) graphics::axis(1)
    }
    graphics::mtext(legtxt(if (!base::is.na(cc[["range"]]))
                             base::paste0(nm.rng("lines"), ": ", what, " of each entry", scl),
                           if (!base::is.na(cc[["estimate"]])) base::paste0(nm.est("points"), ": estimate")),
                    side = 1, outer = TRUE, line = 2.5, cex = 0.7)
    return(base::invisible(x))
  }
  hasci <- !base::is.null(x$C.ci.lower) && !prob   # intervals are for C itself, not its shares
  hi <- if (hasci) x$C.ci.upper else Cu
  hs <- hi[sel, , drop = FALSE]
  xmax <- if (prob) 1 else base::max(hs[base::is.finite(hs)], Cv[sel, ]) * 1.05
  nq <- base::length(sel)
  graphics::layout(base::matrix(base::seq_len(nq * R), nq, R, byrow = TRUE), widths = base::c(1.35, base::rep(1, R - 1)))
  graphics::par(oma = base::c(2, 0, 0.5, 0), mgp = base::c(1.5, 0.4, 0), tcl = -0.25)
  for (i in base::seq_len(nq)) for (k in base::seq_len(R)) {
    q <- sel[i]
    graphics::par(mar = base::c(if (i == nq) 2 else 0.6, if (k == 1) 5 else 0.5, if (i == 1) 1.4 else 0.3, 0.5))
    graphics::plot(NA, xlim = base::c(0, xmax), ylim = base::c(0.5, 1.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    if (hasci && x$C.ci.lower[q, k] > 0 && !base::is.na(cc[["shade"]]))
      graphics::rect(-xmax, 0, 2 * xmax, 2, col = cc[["shade"]], border = NA)
    if (hasci && !base::is.na(cc[["interval"]]))
      graphics::segments(x$C.ci.lower[q, k], 1, base::min(x$C.ci.upper[q, k], 2 * xmax), 1,
                         lwd = lw[["interval"]], col = cc[["interval"]], lend = 1)
    if (!base::is.na(cc[["range"]])) {
      graphics::segments(Cl[q, k], 1, Cu[q, k], 1, lwd = lw[["range"]], lend = 1, col = cc[["range"]])
      if (Cu[q, k] - Cl[q, k] < 1e-9 * xmax)
        graphics::points(Cl[q, k], 1, pch = 19, cex = cx[["range"]], col = cc[["range"]])
    }
    if (!base::is.na(cc[["estimate"]]))
      graphics::points(Cv[q, k], 1, pch = 4, col = cc[["estimate"]], cex = cx[["estimate"]], lwd = lw[["estimate"]])
    graphics::axis(1, labels = i == nq, cex.axis = 0.75)
    if (k == 1) graphics::mtext(bn[q], side = 2, line = 0.5, las = 1, adj = 1, cex = 0.7)
    if (i == 1) graphics::mtext(base::colnames(x$C)[k], side = 3, line = 0.2, cex = 0.75)
    graphics::box()
  }
  leg <- legtxt(if (!base::is.na(cc[["range"]])) base::paste0(nm.rng("thick line"), ": ", what, scl),
                if (hasci && !base::is.na(cc[["interval"]]))
                  base::sprintf("%s: %g%% interval",
                                if (base::nzchar(cw[["interval"]])) cw[["interval"]] else "thin line", 100 * x$ci.level),
                if (!base::is.na(cc[["estimate"]])) base::paste0(nm.est("x"), ": estimate"),
                if (hasci && !base::is.na(cc[["shade"]])) "shaded: interval excludes 0")
  graphics::mtext(leg, side = 1, outer = TRUE, line = 0.6, cex = 0.7)
  base::invisible(x)
}

## colours of plot.nmfkc.ident(): the defaults, the caller's overrides, and the word
## naming each colour in the legend ("" for a colour without a name, e.g. "#1f77b4")
.nid.plot.col <- function(col) {
  def <- base::c(range = "black", estimate = "red3", interval = "gray55", shade = "gray90")
  word <- base::c(range = "black", estimate = "red", interval = "gray", shade = "")
  if (base::is.null(col)) return(base::list(col = def, word = word))
  nm <- base::names(col)
  if (!base::is.atomic(col) || base::is.null(nm) || base::anyNA(nm) || !base::all(nm %in% base::names(def)) ||
      base::anyDuplicated(nm))
    stop("col must be a named vector with names among \"range\", \"estimate\", \"interval\" and \"shade\"")
  col <- base::as.character(col); base::names(col) <- nm
  ok <- base::vapply(col, function(v) base::is.na(v) ||
                       !base::inherits(base::try(grDevices::col2rgb(v), silent = TRUE), "try-error"), TRUE)
  if (!base::all(ok)) stop("col: not a colour: ", base::paste(col[!ok], collapse = ", "))
  for (k in nm) if (!base::identical(col[[k]], def[[k]]))
    word[[k]] <- if (!base::is.na(col[[k]]) && col[[k]] %in% grDevices::colors()) base::sub("[0-9]+$", "", col[[k]]) else ""
  def[nm] <- col
  base::list(col = def, word = word)
}

## mark sizes (cex) or line widths (lwd) of plot.nmfkc.ident(): a named vector of
## positive numbers over `allowed`; the defaults depend on the type and are set there
.nid.plot.size <- function(v, arg, allowed) {
  if (base::is.null(v)) return(base::numeric(0))
  nm <- base::names(v)
  if (!base::is.numeric(v) || base::is.null(nm) || base::anyNA(nm) || !base::all(nm %in% allowed) ||
      base::anyDuplicated(nm))
    stop(base::sprintf("%s must be a named numeric vector with names among %s", arg,
                       base::paste0("\"", allowed, "\"", collapse = ", ")))
  if (base::any(!base::is.finite(v) | v <= 0)) stop(base::sprintf("%s must be positive", arg))
  v
}

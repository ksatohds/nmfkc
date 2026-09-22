#' Choose the number of factors for NMF-FFB by element-wise cross-validation
#'
#' Step 1 of the NMF-FFB workflow: select \eqn{Q} by element-wise
#' cross-validation of the \emph{feed-forward} fit.  Held-out entries of
#' \eqn{Y_1} are predicted from the remaining ones, so the criterion measures how
#' well a \eqn{Q}-factor non-negative basis describes the outcomes, which is what
#' \eqn{Q} is for.
#'
#' Feedback cannot be cross-validated, and that is not a limitation of the
#' implementation.  By the reduced-form equivalence, the equilibrium mapping
#' \eqn{M = (I - X\Theta_1)^{-1}X\Theta_2} is reproduced exactly by a
#' feed-forward model with the exogenous matrix \eqn{L_Q\Theta_2}: the two models
#' predict \eqn{Y_1} from \eqn{Y_2} identically, so no prediction criterion --
#' column-wise CV, element-wise CV, or a test-sample error -- can prefer one over
#' the other.  Feedback is identified by the \emph{conditional covariance}, not by
#' prediction, which is why it is addressed by \code{\link{nmf.ffb.test}} and not
#' here.
#'
#' @param Y1,Y2 Endogenous and exogenous matrices (variables in rows, units in
#'   columns).
#' @param rank Candidate values of \eqn{Q}.  Defaults to
#'   \code{seq_len(min(3, nrow(Y1)))}.
#' @param X.init Stage-1 initialization, passed to \code{\link{nmfkc}}.
#' @param X.L2.ortho Orthogonality penalty on the basis in stage 1.  Keep the
#'   value that will be used in \code{\link{nmf.ffb}}: \eqn{Q} and the penalty
#'   are not separable.
#' @param epsilon,maxit Convergence tolerance and iteration cap for stage 1.
#' @param ... Passed to \code{\link{nmfkc.ecv}}: \code{nfolds}, \code{seed},
#'   and \code{cores} (or \code{ncores}) to evaluate the rank x fold grid in
#'   parallel, defaulting to \code{getOption("mc.cores", 1L)} as elsewhere in the
#'   package.
#'
#' @return The object returned by \code{\link{nmfkc.ecv}}: the held-out error by
#'   rank, with the selected \code{rank}.
#'
#' @section Workflow:
#' \preformatted{
#' ecv <- nmf.ffb.ecv(Y1, Y2, rank = 1:5)   # 1. choose Q            <- this function
#' fit <- nmf.ffb(Y1, Y2, rank = ecv$rank)  # 2. estimate; BIC selects the support
#' tst <- nmf.ffb.test(fit, Y1, Y2)         # 3. test the feed-forward null
#' dgn <- nmf.ffb.diagnostics(fit)          # 4. cycles, spectral radius, identifiability
#' inf <- nmf.ffb.inference(fit, Y1, Y2)    # 5. intervals for the retained entries
#' }
#'
#' @examples
#' set.seed(1)
#' Y <- t(as.matrix(iris[, 1:4]))
#' Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
#' ecv <- nmf.ffb.ecv(Y1, Y2, rank = 1:2, nfolds = 3)
#' ecv$rank
#'
#' @seealso \code{\link{nmf.ffb}}, \code{\link{nmf.ffb.test}},
#'   \code{\link{nmfkc.ecv}}
#' @export
nmf.ffb.ecv <- function(Y1, Y2,
                        rank = NULL,
                        X.init = "nndsvd",
                        X.L2.ortho = 100.0,
                        epsilon = 1e-6,
                        maxit = 5000,
                        ...) {
  ## One implementation: nmf.ffb.cv(method = "fiml") has delegated to nmfkc.ecv() since 0.9.8, and this is
  ## the same call under the name that says what it does.  nmf.ffb.cv() is kept for the published
  ## multiplicative-update path (see NEWS 0.9.8).
  nmf.ffb.cv(Y1, Y2, rank = rank, X.init = X.init, X.L2.ortho = X.L2.ortho,
             epsilon = epsilon, maxit = maxit, method = "fiml", ...)
}

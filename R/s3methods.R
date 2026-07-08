# ============================================================
# S3 generic methods: coef, fitted, residuals, plot, summary
# ============================================================

## Internal: row order for a coefficients table by grouping factor.
## Every parameter (Theta / C) table has a row factor (the basis / X slot)
## and a column factor (the covariate / A slot). In the symmetric-network
## case A = t(X), so the two factors are stored as "Basis.row"/"Basis.col"
## instead of "Basis"/"Covariate", but the grouping rule is identical:
##   by = "covariate" (default): list all rows within each column factor
##     (column outer, row inner) -- the historical display order.
##   by = "basis": list all columns within each row factor (row outer).
## The original appearance order of each factor is preserved.
.coef.order.by <- function(cf, by = c("covariate", "basis")) {
  by <- match.arg(by)
  rowf <- if (!is.null(cf$Basis))     cf$Basis     else cf$Basis.row  # X slot
  colf <- if (!is.null(cf$Covariate)) cf$Covariate else cf$Basis.col  # A slot
  if (is.null(rowf) || is.null(colf)) return(seq_len(nrow(cf)))
  b  <- factor(rowf, levels = unique(rowf))
  cv <- factor(colf, levels = unique(colf))
  if (by == "covariate") order(cv, b) else order(b, cv)
}

# --- plot (convergence) ---

#' @title Plot convergence diagnostics for NMF models
#' @description
#' Plots the objective function value over iterations for \code{nmfre} and
#' \code{nmf.sem} objects. (For \code{nmfkc} and \code{nmfae}, plot methods
#' are defined in their respective source files.)
#'
#' @param x A fitted model object.
#' @param ... Additional graphical arguments passed to \code{\link{plot}}.
#' @return Invisible \code{NULL}.
#' @seealso \code{\link{nmfre}}, \code{\link{nmf.sem}}
#' @name plot.nmfre
#' @examples
#' \donttest{
#' set.seed(1)
#' Y <- matrix(runif(20), nrow = 4)
#' A <- diag(5)
#' res <- nmfre(Y, A, rank = 2)
#' plot(res)
#' }
#'
NULL

#' @rdname plot.nmfre
#' @export
plot.nmfre <- function(x, ...) {
  extra_args <- list(...)
  ## Default to the marginal negative log-likelihood, which the ECM decreases
  ## monotonically. The fixed-lambda penalized objective (objfunc.iter) is NOT
  ## monotone across outer iterations because lambda = sigma2/tau2 is updated
  ## (the penalty jumps by (lambda' - lambda) ||U||^2), so it is unsuitable for
  ## illustrating monotone convergence. Fall back to objfunc.iter for objects
  ## fitted before nll.trace was recorded.
  use_nll <- !is.null(x$nll.trace) && any(is.finite(x$nll.trace))
  args <- list(x = if (use_nll) x$nll.trace else x$objfunc.iter)
  if (is.null(extra_args$main)) {
    r2 <- if (!is.null(x$r.squared) && is.finite(x$r.squared)) round(x$r.squared, 3) else "NA"
    args$main <- paste0("R2 = ", r2)
  }
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- if (use_nll) "marginal NLL" else "objfunc"
  if (is.null(extra_args$type)) args$type <- "l"
  all_args <- c(args, extra_args)
  do.call("plot", all_args)
  invisible(NULL)
}

#' @rdname plot.nmfre
#' @param which For \code{plot.nmf.sem}: which objective to plot.
#'   One of \code{"full"} (default; \code{loss + penalties}, the actual
#'   monotonically-decreasing quantity that the multiplicative updates
#'   minimize), \code{"reconstruction"} (Frobenius distance only,
#'   \eqn{\| Y_1 - X B \|_F^2}), or \code{"both"} (overlay both with a
#'   legend).  \code{"both"} is useful for diagnosing whether
#'   regularization is actively shaping the solution: if the two curves
#'   diverge, the penalties are pulling the optimizer away from the
#'   pure least-squares minimum.
#' @export
plot.nmf.sem <- function(x, ..., which = c("full", "reconstruction", "both")) {
  which <- match.arg(which)
  extra_args <- list(...)

  ## Pick the iteration trace(s) to plot.  Older nmf.sem objects may
  ## carry only x$objfunc (reconstruction loss); fall back gracefully.
  has_full <- !is.null(x$objfunc.full)
  if (which == "full" && !has_full) which <- "reconstruction"

  if (which == "both" && has_full) {
    y_full <- x$objfunc.full
    y_rec  <- x$objfunc
    iter_idx <- seq_along(y_full)
    main_default <- sprintf("MAE = %s, SC.cov = %s",
                            round(x$MAE, 3), round(x$SC.cov, 3))
    if (is.null(extra_args$main)) extra_args$main <- main_default
    if (is.null(extra_args$xlab)) extra_args$xlab <- "iter"
    if (is.null(extra_args$ylab)) extra_args$ylab <- "objective"
    do.call(plot, c(list(iter_idx, y_full, type = "l", lwd = 2,
                          col = "black"), extra_args))
    graphics::lines(iter_idx, y_rec, col = "tomato", lwd = 1.5, lty = 2)
    graphics::legend("topright",
                     legend = c("loss + penalties (full)", "reconstruction only"),
                     col    = c("black", "tomato"),
                     lty    = c(1, 2), lwd = c(2, 1.5),
                     bty    = "n", cex = 0.85)
  } else {
    y <- if (which == "full") x$objfunc.full else x$objfunc
    args <- list(x = y)
    if (is.null(extra_args$main)) {
      label <- if (which == "full") "loss + penalties"
               else "reconstruction only"
      args$main <- sprintf("%s | MAE = %s, SC.cov = %s",
                           label, round(x$MAE, 3), round(x$SC.cov, 3))
    }
    if (is.null(extra_args$xlab)) args$xlab <- "iter"
    if (is.null(extra_args$ylab))
      args$ylab <- if (which == "full") "objfunc.full" else "objfunc"
    if (is.null(extra_args$type)) args$type <- "l"
    all_args <- c(args, extra_args)
    do.call("plot", all_args)
  }
  invisible(NULL)
}


# --- summary.nmf.sem ---

#' @title Summary method for nmf.sem objects
#' @description
#' Produces a formatted summary of a fitted NMF-FFB model, including
#' matrix dimensions, convergence, stability diagnostics, fit statistics,
#' and inference results (if available).
#'
#' @param object An object of class \code{"nmf.ffb"} (or legacy
#'   \code{"nmf.sem"}) returned by \code{\link{nmf.ffb}} /
#'   \code{\link{nmf.sem}}.
#' @param ... Not used.
#' @return An object of class \code{"summary.nmf.sem"} (the fitted model
#'   tagged for printing); printed by \code{\link{print.summary.nmf.sem}}.
#' @seealso \code{\link{nmf.ffb}}, \code{\link{nmf.ffb.inference}}
#' @export
#' @examples
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
#' result <- nmf.ffb(Y1, Y2, rank = 2, maxit = 500)
#' summary(result)
#'
summary.nmf.sem <- function(object, ...) {
  class(object) <- "summary.nmf.sem"
  object
}

#' @title Print method for summary.nmf.sem objects
#' @description
#' Prints the NMF-FFB model summary (dimensions, convergence, stability
#' diagnostics, fit statistics, and inference results if available).
#' @param x An object of class \code{"summary.nmf.sem"} returned by
#'   \code{\link{summary.nmf.sem}}.
#' @param ... Not used.
#' @return Invisible \code{x}.
#' @seealso \code{\link{summary.nmf.sem}}
#' @export
print.summary.nmf.sem <- function(x, ...) {
  object <- x
  P1 <- nrow(object$X)
  Q  <- ncol(object$X)
  P2 <- ncol(object$C2)

  cat(sprintf("NMF-FFB: Y1(%d,N) = X(%d,%d) [C1(%d,%d) Y1 + C2(%d,%d) Y2]\n",
              P1, P1, Q, Q, P1, Q, P2))
  cat(sprintf("Iterations: %d\n", object$iter))

  cat("\nStability diagnostics:\n")
  cat(sprintf("  Spectral radius(XC1): %.4f %s\n",
              object$XC1.radius,
              if (object$XC1.radius < 1) "(stable)" else "(UNSTABLE)"))
  cat(sprintf("  ||XC1||_1:            %.4f\n", object$XC1.norm1))
  cat(sprintf("  Amplification:        %.4f (bound: %.4f)\n",
              object$amplification, object$amplification.bound))

  cat("\nFit statistics:\n")
  if (!is.null(object$SC.map) && is.finite(object$SC.map))
    cat(sprintf("  SC.map (mapping correlation):    %.4f\n", object$SC.map))
  if (!is.null(object$SC.cov) && is.finite(object$SC.cov))
    cat(sprintf("  SC.cov (covariance correlation): %.4f\n", object$SC.cov))
  if (!is.null(object$MAE)   && is.finite(object$MAE))
    cat(sprintf("  MAE (mean absolute error):       %.4f\n", object$MAE))
  if (!is.null(object$effective.rank) && is.finite(object$effective.rank))
    cat(sprintf("  Effective Rank:                  %.2f / %d  (%.1f%%)\n",
                object$effective.rank, Q, 100 * object$effective.rank / Q))

  # Coefficients from inference
  if (!is.null(object$coefficients) && is.data.frame(object$coefficients)) {
    cf <- object$coefficients

    ## Bootstrap meta-info (new full-pair-bootstrap inference; v0.6.8+)
    has_boot <- !is.null(object$bootstrap.B)
    if (has_boot) {
      cat(sprintf("\nBootstrap inference (X-fixed full pair bootstrap):\n"))
      cat(sprintf("  B = %d, valid = %d, threshold = %g, ci.level = %g\n",
                  object$bootstrap.B,
                  if (!is.null(object$bootstrap.n.valid)) object$bootstrap.n.valid else NA_integer_,
                  if (!is.null(object$bootstrap.threshold)) object$bootstrap.threshold else NA_real_,
                  if (!is.null(object$bootstrap.ci.level)) object$bootstrap.ci.level else NA_real_))
    }

    ## Print one block per Type ("C1" feedback / "C2" exogenous) when
    ## present.  Falls back to a single block for legacy inference output
    ## that lacks the Type column.
    print_block <- function(block, title) {
      if (nrow(block) == 0L) return(invisible(NULL))
      cat(sprintf("\n%s\n", title))
      rnames <- paste0(block$Covariate, " -> ", block$Basis)
      est <- formatC(block$Estimate, format = "g", digits = 4, width = 10)
      cl  <- if ("CI_low"  %in% names(block))
               formatC(block$CI_low,  format = "g", digits = 3, width = 10) else NULL
      cu  <- if ("CI_high" %in% names(block))
               formatC(block$CI_high, format = "g", digits = 3, width = 10) else NULL
      sup <- if ("support_rate" %in% names(block))
               formatC(block$support_rate, format = "f", digits = 3, width = 8) else NULL
      pv  <- block$p_value
      pv_str <- ifelse(!is.finite(pv), "      NA",
                 ifelse(pv < 2.2e-16, "  <2e-16",
                   formatC(pv, format = "g", digits = 4, width = 8)))
      stars <- if ("sig" %in% names(block)) {
                 block$sig
               } else {
                 ifelse(!is.finite(pv), " ",
                   ifelse(pv < 0.001, "***",
                     ifelse(pv < 0.01, "**",
                       ifelse(pv < 0.05, "*",
                         ifelse(pv < 0.1, ".", " ")))))
               }
      max_lw <- max(nchar(rnames))
      hdr_parts <- c(formatC("Estimate", width = 10))
      if (!is.null(cl)) hdr_parts <- c(hdr_parts, formatC("CI_low",  width = 10))
      if (!is.null(cu)) hdr_parts <- c(hdr_parts, formatC("CI_high", width = 10))
      if (!is.null(sup)) hdr_parts <- c(hdr_parts, formatC("support", width = 8))
      hdr_parts <- c(hdr_parts, formatC("Pr(>0)", width = 8), "")
      cat(sprintf("%s %s\n", formatC("", width = max_lw),
                  paste(hdr_parts, collapse = " ")))
      for (i in seq_along(rnames)) {
        row_parts <- c(est[i])
        if (!is.null(cl))  row_parts <- c(row_parts, cl[i])
        if (!is.null(cu))  row_parts <- c(row_parts, cu[i])
        if (!is.null(sup)) row_parts <- c(row_parts, sup[i])
        row_parts <- c(row_parts, pv_str[i], stars[i])
        cat(sprintf("%s %s\n",
                    formatC(rnames[i], width = max_lw),
                    paste(row_parts, collapse = " ")))
      }
    }

    if (!is.null(cf$Type)) {
      ## v0.6.8+ inference: separate C1 / C2 blocks
      cf_C1 <- cf[cf$Type == "C1", , drop = FALSE]
      cf_C2 <- cf[cf$Type == "C2", , drop = FALSE]
      print_block(cf_C1, "C1 Coefficients (Y1 -> Factor; feedback Theta_1):")
      print_block(cf_C2, "C2 Coefficients (Y2 -> Factor; exogenous Theta_2):")
    } else {
      ## Legacy inference output (only C2)
      print_block(cf, "C2 Coefficients (Covariate -> Basis):")
    }
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  invisible(object)
}


# --- coef ---

#' @title Extract coefficients from NMF models
#' @description
#' Returns the \code{coefficients} data frame from a fitted NMF model
#' that has been passed through an inference function
#' (\code{\link{nmfkc.inference}}, \code{\link{nmfae.inference}},
#' \code{\link{nmfre.inference}}).
#'
#' If inference has not been run, returns the parameter matrix \eqn{C}
#' (\eqn{\Theta}) directly.
#'
#' For \code{nmf.sem} objects, returns \eqn{C_2} (exogenous block) as fallback.
#'
#' @param object A fitted model object of class \code{"nmf"}, \code{"nmfkc"},
#'   \code{"nmfae"}, \code{"nmfre"}, or \code{"nmf.sem"}.
#' @param ... Not used.
#' @return A data frame of coefficients (if inference was performed),
#'   or the parameter matrix \eqn{C}.
#' @seealso \code{\link{nmfkc.inference}}, \code{\link{nmfae.inference}},
#'   \code{\link{nmfre.inference}}, \code{\link{nmf.sem.inference}}
#' @name coef.nmf
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, rank = 1)
#' coef(result)  # returns C matrix
#'
#' result2 <- nmfkc.inference(result, Y, A)
#' coef(result2)  # returns coefficients data frame
#'
NULL

#' @rdname coef.nmf
#' @export
coef.nmf <- function(object, ...) {
  if (!is.null(object$coefficients)) object$coefficients else object$C
}

#' @rdname coef.nmf
#' @export
coef.nmf.sem <- function(object, ...) {
  ## With inference: return the unified coefficients data frame (rows for
  ## both C1 and C2, with CI / support_rate / sig columns).
  if (!is.null(object$coefficients)) return(object$coefficients)

  ## Without inference: return a long-format data frame with rows for
  ## every entry of C1 (feedback) and C2 (exogenous), so the column
  ## layout matches `coef(res_inf)`.  Users can post-hoc filter via
  ## `subset(coef(res), Type == "C1")` regardless of whether
  ## nmf.sem.inference() has been run.
  C1 <- object$C1
  C2 <- object$C2
  if (is.null(C1) || is.null(C2)) {
    ## Fallback (shouldn't normally happen for an nmf.sem result)
    return(if (!is.null(C2)) C2 else C1)
  }
  Q  <- nrow(C1)
  P1 <- ncol(C1)
  P2 <- ncol(C2)
  bas <- if (!is.null(rownames(C1))) rownames(C1) else paste0("Factor", seq_len(Q))
  y1  <- if (!is.null(colnames(C1))) colnames(C1) else paste0("Y1_", seq_len(P1))
  y2  <- if (!is.null(colnames(C2))) colnames(C2) else paste0("Y2_", seq_len(P2))
  C1_block <- data.frame(
    Type      = "C1",
    Basis     = rep(bas, times = P1),
    Covariate = rep(y1,  each  = Q),
    Estimate  = as.vector(C1),
    stringsAsFactors = FALSE
  )
  C2_block <- data.frame(
    Type      = "C2",
    Basis     = rep(bas, times = P2),
    Covariate = rep(y2,  each  = Q),
    Estimate  = as.vector(C2),
    stringsAsFactors = FALSE
  )
  out <- rbind(C1_block, C2_block)
  rownames(out) <- NULL
  out
}


# --- fitted ---

#' @title Extract fitted values from NMF models
#' @description
#' Returns the reconstructed matrix \eqn{\hat{Y} = X B} from a fitted
#' NMF model.
#'
#' For \code{nmf.sem} objects, returns the equilibrium prediction
#' \eqn{\hat{Y}_1 = M_{model} Y_2} if available. Supply \code{Y1} and
#' \code{Y2} to get the direct reconstruction
#' \eqn{X (C_1 Y_1 + C_2 Y_2)} instead.
#'
#' @param object A fitted model object of class \code{"nmf"}, \code{"nmfkc"},
#'   \code{"nmfae"}, \code{"nmfre"}, or \code{"nmf.sem"}.
#' @param ... For \code{nmf.sem}: optionally \code{Y1} and \code{Y2}.
#' @return The fitted matrix \eqn{X B}.
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfae}}, \code{\link{nmfre}},
#'   \code{\link{nmf.sem}}, \code{\link{residuals.nmf}}
#' @name fitted.nmf
#' @examples
#' result <- nmfkc(matrix(runif(50), 5, 10), rank = 2)
#' fitted(result)
#'
NULL

#' @rdname fitted.nmf
#' @export
fitted.nmf <- function(object, ...) {
  object$XB
}

#' @rdname fitted.nmf
#' @export
fitted.nmfae <- function(object, ...) {
  object$Y1hat
}

#' @rdname fitted.nmf
#' @param type For \code{nmfre} objects: \code{"blup"} (default, the fit
#'   \eqn{X(\Theta A + U)} including random effects) or \code{"fixed"}
#'   (\eqn{X\Theta A}, fixed effects only). Chosen so that
#'   \code{Y - fitted(object, type)} equals \code{residuals(object, Y, type)}.
#' @export
fitted.nmfre <- function(object, type = c("blup", "fixed"), ...) {
  type <- match.arg(type)
  if (type == "blup") object$XB.blup else object$XB
}

#' @rdname fitted.nmf
#' @export
fitted.nmf.sem <- function(object, ...) {
  extra <- list(...)
  Y1 <- extra$Y1
  Y2 <- extra$Y2
  if (!is.null(Y1) && !is.null(Y2)) {
    # Direct reconstruction: X(C1*Y1 + C2*Y2)
    object$X %*% (object$C1 %*% Y1 + object$C2 %*% Y2)
  } else {
    # Equilibrium prediction: M.model %*% Y2 (requires Y2)
    if (!is.null(Y2) && !anyNA(object$M.model)) {
      object$M.model %*% Y2
    } else {
      stop("Supply Y1 and Y2 for direct reconstruction, or Y2 for equilibrium prediction.")
    }
  }
}


# --- residuals ---

#' @title Extract residuals from NMF models
#' @description
#' Returns the residual matrix \eqn{Y - \hat{Y}} from a fitted NMF model.
#' Requires the original observation matrix \code{Y} to be supplied.
#'
#' For \code{nmfre} objects, residuals are computed from the BLUP
#' reconstruction (\eqn{Y - X(B_{blup})}) by default. Set
#' \code{type = "fixed"} to use fixed-effects only.
#'
#' @param object A fitted model object.
#' @param Y The original observation matrix used for fitting.
#' @param type For \code{nmfre} objects: \code{"blup"} (default) or
#'   \code{"fixed"}.
#' @param ... Not used.
#' @return The residual matrix.
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfae}}, \code{\link{nmfre}},
#'   \code{\link{nmf.sem}}, \code{\link{fitted.nmf}}
#' @name residuals.nmf
#' @examples
#' Y <- matrix(runif(50), 5, 10)
#' result <- nmfkc(Y, rank = 2)
#' residuals(result, Y)
#'
NULL

#' @rdname residuals.nmf
#' @export
residuals.nmf <- function(object, Y, ...) {
  Y - object$XB
}

#' @rdname residuals.nmf
#' @export
residuals.nmfae <- function(object, Y, ...) {
  Y - object$Y1hat
}

#' @rdname residuals.nmf
#' @export
residuals.nmfre <- function(object, Y, type = c("blup", "fixed"), ...) {
  type <- match.arg(type)
  if (type == "blup") {
    Y - object$XB.blup
  } else {
    Y - object$XB
  }
}

#' @rdname residuals.nmf
#' @export
residuals.nmf.sem <- function(object, Y, ...) {
  ## Delegate to fitted.nmf.sem so residuals use the SAME reconstruction as
  ## fitted (direct if Y1/Y2 given via ..., else the Y2-equilibrium form).
  ## Y is the observed response block (Y1). nmf.sem/nmf.ffb need Y2 (and,
  ## for the direct form, Y1) supplied via ... -- see fitted.nmf.sem.
  Y - fitted(object, ...)
}


# --- print.nmf.inference ---

#' @title Print method for NMF inference objects
#' @description
#' Prints a summary of any NMF inference result object
#' (\code{"nmfkc.inference"} or \code{"nmfae.inference"}).
#' @param x An object of class \code{"nmf.inference"}.
#' @param ... Additional arguments passed to the corresponding
#'   \code{print.summary.*} method.
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @seealso \code{\link{nmfkc.inference}}, \code{\link{nmfae.inference}}
#' @export
print.nmf.inference <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

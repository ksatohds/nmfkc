# ============================================================
# S3 generic methods: coef, fitted, residuals, plot, summary
# ============================================================

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
#' @name plot.nmfre
#' @examples
#' \donttest{
#' res <- nmfre(Y, A, rank = 2, wild.bootstrap = FALSE)
#' plot(res)
#' }
#'
NULL

#' @rdname plot.nmfre
#' @export
plot.nmfre <- function(x, ...) {
  extra_args <- list(...)
  args <- list(x = x$objfunc.iter)
  if (is.null(extra_args$main)) {
    r2 <- if (!is.null(x$r.squared) && is.finite(x$r.squared)) round(x$r.squared, 3) else "NA"
    args$main <- paste0("R2 = ", r2)
  }
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  if (is.null(extra_args$type)) args$type <- "l"
  all_args <- c(args, extra_args)
  do.call("plot", all_args)
  invisible(NULL)
}

#' @rdname plot.nmfre
#' @export
plot.nmf.sem <- function(x, ...) {
  extra_args <- list(...)
  args <- list(x = x$objfunc)
  if (is.null(extra_args$main)) {
    args$main <- paste0("MAE = ", round(x$MAE, 3), ", SC.cov = ", round(x$SC.cov, 3))
  }
  if (is.null(extra_args$xlab)) args$xlab <- "iter"
  if (is.null(extra_args$ylab)) args$ylab <- "objfunc"
  if (is.null(extra_args$type)) args$type <- "l"
  all_args <- c(args, extra_args)
  do.call("plot", all_args)
  invisible(NULL)
}


# --- summary.nmf.sem ---

#' @title Summary method for nmf.sem objects
#' @description
#' Produces a formatted summary of a fitted NMF-SEM model, including
#' matrix dimensions, convergence, stability diagnostics, fit statistics,
#' and inference results (if available).
#'
#' @param object An object of class \code{"nmf.sem"} returned by
#'   \code{\link{nmf.sem}}.
#' @param ... Not used.
#' @return Invisible \code{object}.
#' @seealso \code{\link{nmf.sem}}, \code{\link{nmf.sem.inference}}
#' @export
#' @examples
#' Y <- t(iris[, -5])
#' Y1 <- Y[1:2, ]; Y2 <- Y[3:4, ]
#' result <- nmf.sem(Y1, Y2, rank = 2, maxit = 500)
#' summary(result)
#'
summary.nmf.sem <- function(object, ...) {
  P1 <- nrow(object$X)
  Q  <- ncol(object$X)
  P2 <- ncol(object$C2)

  cat(sprintf("NMF-SEM: Y1(%d,N) = X(%d,%d) [C1(%d,%d) Y1 + C2(%d,%d) Y2]\n",
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
  if (!is.na(object$SC.cov))
    cat(sprintf("  SC.cov (covariance correlation): %.4f\n", object$SC.cov))
  if (!is.na(object$MAE))
    cat(sprintf("  MAE (mean absolute error):       %.4f\n", object$MAE))

  # Coefficients from inference
  if (!is.null(object$coefficients) && is.data.frame(object$coefficients)) {
    cat("\nC2 Coefficients (Exogenous -> Factor):\n")
    cf <- object$coefficients
    p_side <- if (!is.null(object$C2.p.side)) object$C2.p.side else "one.sided"
    p_header <- if (p_side == "one.sided") "Pr(>z)" else "Pr(>|z|)"

    rnames <- paste0(cf$Exogenous, ":", cf$Factor)
    est <- formatC(cf$Estimate, format = "f", digits = 4, width = 10)
    se  <- formatC(cf$SE, format = "f", digits = 4, width = 10)
    zv  <- formatC(cf$z_value, format = "f", digits = 2, width = 7)
    pv  <- cf$p_value
    pv_str <- ifelse(!is.finite(pv), "      NA",
               ifelse(pv < 2.2e-16, "  <2e-16",
                 formatC(pv, format = "g", digits = 4, width = 8)))
    stars <- ifelse(!is.finite(pv), " ",
               ifelse(pv < 0.001, "***",
                 ifelse(pv < 0.01, "**",
                   ifelse(pv < 0.05, "*",
                     ifelse(pv < 0.1, ".", " ")))))

    has_bse <- "BSE" %in% names(cf) && any(is.finite(cf$BSE))
    max_lw <- max(nchar(rnames))

    if (has_bse) {
      bse <- formatC(cf$BSE, format = "f", digits = 4, width = 8)
      hdr <- sprintf("%s %s %s %s %s %s",
                     formatC("Estimate", width = 10),
                     formatC("Std. Error", width = 10),
                     formatC("(Boot)", width = 8),
                     formatC("z value", width = 7),
                     formatC(p_header, width = 8), "")
      cat(sprintf("%s %s\n", formatC("", width = max_lw), hdr))
      for (i in seq_along(rnames)) {
        cat(sprintf("%s %s %s %s %s %s %s\n",
                    formatC(rnames[i], width = max_lw),
                    est[i], se[i], bse[i], zv[i], pv_str[i], stars[i]))
      }
    } else {
      hdr <- sprintf("%s %s %s %s %s",
                     formatC("Estimate", width = 10),
                     formatC("Std. Error", width = 10),
                     formatC("z value", width = 7),
                     formatC(p_header, width = 8), "")
      cat(sprintf("%s %s\n", formatC("", width = max_lw), hdr))
      for (i in seq_along(rnames)) {
        cat(sprintf("%s %s %s %s %s %s\n",
                    formatC(rnames[i], width = max_lw),
                    est[i], se[i], zv[i], pv_str[i], stars[i]))
      }
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
#' @param object A fitted model object.
#' @param ... Not used.
#' @return A data frame of coefficients (if inference was performed),
#'   or the parameter matrix \eqn{C}.
#' @seealso \code{\link{nmfkc.inference}}, \code{\link{nmfae.inference}},
#'   \code{\link{nmfre.inference}}
#' @name coef.nmfkc
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, Q = 1)
#' coef(result)  # returns C matrix
#'
#' result2 <- nmfkc.inference(result, Y, A)
#' coef(result2)  # returns coefficients data frame
#'
NULL

#' @rdname coef.nmfkc
#' @export
coef.nmfkc <- function(object, ...) {
  if (!is.null(object$coefficients)) object$coefficients else object$C
}

#' @rdname coef.nmfkc
#' @export
coef.nmfae <- function(object, ...) {
  if (!is.null(object$coefficients)) object$coefficients else object$C
}

#' @rdname coef.nmfkc
#' @export
coef.nmfre <- function(object, ...) {
  if (!is.null(object$coefficients)) object$coefficients else object$C
}


#' @rdname coef.nmfkc
#' @export
coef.nmf.sem <- function(object, ...) {
  if (!is.null(object$coefficients)) object$coefficients else object$C2
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
#' @param object A fitted model object.
#' @param ... For \code{nmf.sem}: optionally \code{Y1} and \code{Y2}.
#' @return The fitted matrix \eqn{X B}.
#' @name fitted.nmfkc
#' @examples
#' result <- nmfkc(matrix(runif(50), 5, 10), Q = 2)
#' fitted(result)
#'
NULL

#' @rdname fitted.nmfkc
#' @export
fitted.nmfkc <- function(object, ...) {
  object$XB
}

#' @rdname fitted.nmfkc
#' @export
fitted.nmfae <- function(object, ...) {
  # nmfae stores reconstruction in XB (= X1 %*% C %*% X2 %*% Y2)
  object$XB
}

#' @rdname fitted.nmfkc
#' @export
fitted.nmfre <- function(object, ...) {
  # XB.blup includes random effects; XB is fixed-effects only
  object$XB
}

#' @rdname fitted.nmfkc
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
#' @name residuals.nmfkc
#' @examples
#' Y <- matrix(runif(50), 5, 10)
#' result <- nmfkc(Y, Q = 2)
#' residuals(result, Y)
#'
NULL

#' @rdname residuals.nmfkc
#' @export
residuals.nmfkc <- function(object, Y, ...) {
  Y - object$XB
}

#' @rdname residuals.nmfkc
#' @export
residuals.nmfae <- function(object, Y, ...) {
  Y - object$XB
}

#' @rdname residuals.nmfkc
#' @export
residuals.nmfre <- function(object, Y, type = c("blup", "fixed"), ...) {
  type <- match.arg(type)
  if (type == "blup") {
    Y - object$XB.blup
  } else {
    Y - object$XB
  }
}

#' @rdname residuals.nmfkc
#' @export
residuals.nmf.sem <- function(object, Y, ...) {
  extra <- list(...)
  Y2 <- extra$Y2
  if (is.null(Y2)) stop("Supply Y2 to compute residuals.")
  if (!anyNA(object$M.model)) {
    Y - object$M.model %*% Y2
  } else {
    stop("Leontief inverse is NA; cannot compute residuals.")
  }
}

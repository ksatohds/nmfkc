# ============================================================
# S3 generic methods: coef, fitted, residuals
# ============================================================

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


# --- fitted ---

#' @title Extract fitted values from NMF models
#' @description
#' Returns the reconstructed matrix \eqn{\hat{Y} = X B} from a fitted
#' NMF model.
#'
#' @param object A fitted model object.
#' @param ... Not used.
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

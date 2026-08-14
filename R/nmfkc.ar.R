# nmfkc.ar.R — AR/VAR related functions
# nmfkc.ar, nmfkc.ar.predict, nmfkc.ar.degree.cv, nmfkc.ar.stationarity, nmfkc.ar.DOT

#' @title Construct observation and covariate matrices for a vector autoregressive model
#' @description
#' \code{nmfkc.ar} generates the observation matrix and covariate matrix
#' corresponding to a specified autoregressive lag order.
#'
#' If the input \code{Y} is a \code{ts} object, its time properties are preserved
#' in the \code{"tsp_info"} attribute, adjusted for the lag.
#' Additionally, the column names of \code{Y} and \code{A} are set to the corresponding time points.
#'
#' @param Y An observation matrix (P x N) or a \code{ts} object.
#'   If \code{Y} is a \code{ts} object (typically N x P), it is automatically transposed to match the (P x N) format.
#' @param degree The lag order of the autoregressive model. The default is 1.
#' @param intercept Logical. If TRUE (default), an intercept term is added to the covariate matrix.
#'
#' @return A list containing:
#' \item{Y}{Observation matrix (P x N_A) used for NMF. Includes adjusted \code{"tsp_info"} attribute and time-based column names.}
#' \item{A}{Covariate matrix (R x N_A) constructed according to the specified lag order. Includes adjusted \code{"tsp_info"} attribute and time-based column names.}
#' \item{A.columns}{Index matrix used to generate \code{A}.}
#' \item{degree.max}{Maximum lag order.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ar.degree.cv}}, \code{\link{nmfkc.ar.stationarity}}, \code{\link{nmfkc.ar.DOT}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to multivariate time series data as a vector autoregression model.
#'   \emph{Japanese Journal of Statistics and Data Science}. arXiv:2501.17446.
#'   \doi{10.1007/s42081-025-00314-0}
#' @examples
#' # Example using AirPassengers (ts object)
#' d <- AirPassengers
#' ar_data <- nmfkc.ar(d, degree = 2)
#' dim(ar_data$Y)
#' dim(ar_data$A)
#'
#' # Example using matrix input
#' Y <- matrix(1:20, nrow = 2)
#' ar_data <- nmfkc.ar(Y, degree = 1)
#' ar_data$degree.max
#'
#' @export
nmfkc.ar <- function(Y, degree=1, intercept=TRUE){

  # --- 1. Time Series Object Handling & Standardization ---
  # Initialize basic tsp info
  tsp_start <- 1
  tsp_freq <- 1

  # [FIX 1] Record if input was originally a ts object
  is_ts_object <- stats::is.ts(Y)

  if(is_ts_object){
    # Capture original time properties BEFORE transforming Y
    y_tsp_orig <- stats::tsp(Y)
    tsp_start <- y_tsp_orig[1]
    tsp_freq <- y_tsp_orig[3]

    # Standard 'ts' objects are Time(rows) x Variables(cols).
    # 'nmfkc' requires Variables(rows) x Time(cols).
    Y <- t(as.matrix(Y))

  } else {
    # Handle vector input (not ts) -> 1 x N matrix
    if(is.vector(Y)) Y <- matrix(Y, nrow=1)
    if(!is.matrix(Y)) Y <- as.matrix(Y)

    # If Y is a regular matrix, we assume it is already Variables x Time.
    # Default tsp (Start=1, Freq=1) is used.
  }

  # --- 2. NA Check ---
  if(any(is.na(Y))){
    stop("Y contains missing values (NA). In NMF-VAR, lagged Y is used to construct the covariate matrix A, which cannot contain missing values. Please impute Y before calling nmfkc.ar().")
  }

  P <- nrow(Y) # Features/Variables
  N <- ncol(Y) # Time points

  # --- 3. Degree Validity Check ---
  if (degree < 1) {
    stop("The 'degree' (lag order) must be 1 or greater.")
  }
  if (degree >= N) {
    stop(paste0("The 'degree' (", degree, ") must be strictly less than the number of time points in Y (", N, ")."))
  }

  # Length of the time series after lagging (N_A)
  N_A <- N - degree
  t_start_idx <- degree + 1

  # --- 4. Calculate Adjusted tsp for Y and A ---
  tsp_start_new <- tsp_start + (degree / tsp_freq)
  tsp_end_new <- tsp_start_new + (N_A - 1) / tsp_freq
  y_tsp_new <- c(tsp_start_new, tsp_end_new, tsp_freq)

  # --- 5. Generate Time Sequence for Column Names [FIXED] ---
  # Determine correct column names
  if (!is_ts_object && !is.null(colnames(Y))) {
    # Case A: Not a ts object AND has existing column names.
    # Preserve existing names by slicing them.
    time_names <- colnames(Y)[t_start_idx : N]
  } else {
    # Case B: Is a ts object OR no column names exist.
    # Generate numeric time sequence based on tsp info.
    time_seq <- seq(from = tsp_start_new, by = 1/tsp_freq, length.out = N_A)
    time_names <- as.character(round(time_seq, 5))
  }

  # --- 6. A.columns.index ---
  A.columns.index <- matrix(0, nrow = degree, ncol = N_A)
  for(i in 1:degree) {
    A.columns.index[i, ] <- (t_start_idx - i) : (N - i)
  }

  # --- 7. A: Construct covariate matrix A ---
  A_data <- Y[, as.vector(A.columns.index)]
  A <- matrix(A_data, nrow = P * degree, ncol = N_A, byrow = FALSE)

  # Set row labels
  if(is.null(rownames(Y))) rownames(Y) <- 1:P
  label <- unlist(lapply(1:degree, function(i) paste0(rownames(Y), "_", i)))
  rownames(A) <- label

  # Add Intercept
  if(intercept){
    A <- rbind(A, matrix(1, nrow=1, ncol=ncol(A)))
    rownames(A) <- c(label,"(Intercept)")
  }

  # --- 8. Ya (Non-lagged Observation Matrix) ---
  Ya <- Y[, t_start_idx : N, drop=FALSE]

  # --- 9. Set Column Names ---
  colnames(Ya) <- time_names
  colnames(A) <- time_names

  # --- 10. Metadata Calculation and Attribute Storage ---
  degree.max <- min(ncol(Ya), floor(10*log10(ncol(Ya))))

  # Store AR parameters
  attr(A, "degree") <- degree
  attr(A, "intercept") <- intercept
  attr(A, "function.name") <- "nmfkc.ar"

  # Store Adjusted Time Properties
  attr(A, "tsp_info") <- y_tsp_new
  attr(Ya, "tsp_info") <- y_tsp_new

  list(Y = Ya,
       A = A,
       A.columns = A.columns.index,
       degree.max = degree.max)
}




#' @title Forecast future values for NMF-VAR model
#' @description
#' \code{nmfkc.ar.predict} computes multi-step-ahead forecasts for a fitted NMF-VAR model
#' using recursive forecasting.
#'
#' If the fitted model contains time series property information (from \code{nmfkc.ar}),
#' the forecasted values will have appropriate time-based column names.
#'
#' @param x An object of class \code{nmfkc} (the fitted model).
#' @param Y The historical observation matrix used for fitting (or at least the last \code{degree} columns).
#' @param degree Optional integer. Lag order (D). If \code{NULL} (default), it is inferred
#'   from \code{x$A.attr} (when available) or from the dimensions of \code{x$C}.
#' @param n.ahead Integer (>=1). Number of steps ahead to forecast.
#'
#' @return A list with components:
#' \item{pred}{A \eqn{P \times n.ahead} matrix of predicted values. Column names are future time points if time information is available.}
#' \item{time}{A numeric vector of future time points corresponding to the columns of \code{pred}.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ar}}
#' @examples
#' # Forecast AirPassengers
#' d <- AirPassengers
#' ar_data <- nmfkc.ar(d, degree = 2)
#' result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#' pred <- nmfkc.ar.predict(result, Y = matrix(d, nrow = 1), degree = 2, n.ahead = 3)
#' pred$pred
#'
#' @export
nmfkc.ar.predict <- function(x, Y, degree = NULL, n.ahead = 1){
  # --- Basic checks ---
  if (!inherits(x, "nmfkc")) stop("Argument 'x' must be an object of class 'nmfkc'.")
  if (is.vector(Y)) Y <- matrix(Y, nrow = 1)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.numeric(Y)) stop("'Y' must be numeric.")
  if (length(n.ahead) != 1 || !is.finite(n.ahead) || n.ahead < 1) {
    stop("'n.ahead' must be a positive integer.")
  }
  n.ahead <- as.integer(n.ahead)

  # Fit objects & dims
  Xfit <- x$X
  Cfit <- x$C
  if (is.null(Xfit) || is.null(Cfit)) stop("The fitted object 'x' must contain 'X' and 'C'.")
  P_fit <- nrow(Xfit)
  Q_fit <- ncol(Xfit)
  R_fit <- ncol(Cfit)

  P <- nrow(Y)
  if (P != P_fit) {
    stop(sprintf("Row dimension mismatch: nrow(Y)=%d but nrow(x$X)=%d.", P, P_fit))
  }

  # --- Degree & intercept detection ---
  A.attributes <- x$A.attr
  has_meta <- !is.null(A.attributes) &&
    !is.null(A.attributes$`function.name`) &&
    A.attributes$`function.name` == "nmfkc.ar"

  # 1) user-specified degree has the highest priority
  if (!is.null(degree)) {
    if (length(degree) != 1 || !is.finite(degree) || degree < 1) {
      stop("'degree' must be a positive integer.")
    }
    degree_final <- as.integer(degree)
    # intercept: if metadata available, trust it; otherwise infer by dimension
    if (has_meta && !is.null(A.attributes$intercept)) {
      has_intercept_final <- isTRUE(A.attributes$intercept)
    } else {
      # Infer from C dimensions: R_fit == P*D (+1 if intercept)
      if ((R_fit - 1) %% P_fit == 0 && (R_fit - 1) / P_fit >= 1) {
        # ambiguous if user degree doesn't match dims; check and warn/stop
        if ((R_fit - 1) / P_fit == degree_final) {
          has_intercept_final <- TRUE
        } else if (R_fit %% P_fit == 0 && R_fit / P_fit == degree_final) {
          has_intercept_final <- FALSE
        } else {
          stop(sprintf("Dimension mismatch for supplied degree=%d with C (ncol=%d).", degree_final, R_fit))
        }
      } else if (R_fit %% P_fit == 0 && R_fit / P_fit >= 1) {
        if (R_fit / P_fit == degree_final) {
          has_intercept_final <- FALSE
        } else {
          stop(sprintf("Dimension mismatch for supplied degree=%d with C (ncol=%d).", degree_final, R_fit))
        }
      } else {
        stop("Could not reconcile 'degree' with the dimensions of 'C'.")
      }
    }
  } else if (has_meta) {
    degree_final <- as.integer(A.attributes$degree)
    has_intercept_final <- isTRUE(A.attributes$intercept)
  } else {
    # 3) fall back to dimension inference from C
    if ((R_fit - 1) %% P_fit == 0 && (R_fit - 1) / P_fit >= 1) {
      has_intercept_final <- TRUE
      degree_final <- as.integer((R_fit - 1) / P_fit)
    } else if (R_fit %% P_fit == 0 && R_fit / P_fit >= 1) {
      has_intercept_final <- FALSE
      degree_final <- as.integer(R_fit / P_fit)
    } else {
      stop("Could not infer lag order and intercept status. Please supply 'degree'.")
    }
  }

  # --- Data sufficiency check ---
  if (ncol(Y) < degree_final) {
    stop(sprintf("Not enough historical data in Y (%d columns) for degree=%d.", ncol(Y), degree_final))
  }

  # --- Forecast loop ---
  preds <- matrix(0, nrow = P, ncol = n.ahead)
  rownames(preds) <- rownames(Y)
  colnames(preds) <- paste0("t+", seq_len(n.ahead)) # Default names

  current_Y <- Y
  # Expected length of newA
  expected_R <- P_fit * degree_final + if (has_intercept_final) 1L else 0L
  if (expected_R != R_fit) {
    stop(sprintf("Internal check failed: expected R=%d (from degree/intercept) but ncol(C)=%d.", expected_R, R_fit))
  }

  for (step in seq_len(n.ahead)) {
    # Build lag vector: [Y_t-1; Y_t-2; ...; Y_t-degree] stacked by rows (P each)
    a_vec <- numeric(P * degree_final + if (has_intercept_final) 1L else 0L)
    pos <- 0L
    for (k in 1:degree_final) {
      col_idx <- ncol(current_Y) - k + 1L
      yk <- current_Y[, col_idx, drop = FALSE]
      a_vec[(pos + 1L):(pos + P)] <- as.numeric(yk)
      pos <- pos + P
    }
    if (has_intercept_final) a_vec[length(a_vec)] <- 1

    # Safety check for newA length
    if (length(a_vec) != R_fit) {
      stop(sprintf("Length of constructed lag vector (%d) != ncol(C) (%d).", length(a_vec), R_fit))
    }

    newA <- matrix(a_vec, ncol = 1)
    pred_step <- Xfit %*% (Cfit %*% newA)  # (P x Q) %*% (Q x R) %*% (R x 1) = (P x 1)
    preds[, step] <- as.numeric(pred_step)
    current_Y <- cbind(current_Y, pred_step)
  }

  # --- Time inference logic ---
  future_times <- NULL

  # Strategy 1: Use "tsp_info" from attributes (Robust method)
  if (has_meta && !is.null(A.attributes$tsp_info)) {
    # tsp_info = c(start, end, freq)
    t_info <- A.attributes$tsp_info
    last_time <- t_info[2]
    freq <- t_info[3]

    # If user passed a subset Y (not the full training Y), last_time in attribute might be old.
    # We should check if colnames(Y) can be parsed to sync with Y's actual end time.
    # However, if Y is just the tail of training data, we can try to infer dt from freq.

    if (freq <= 0) stop("'freq' in tsp_info must be positive.")
    dt <- 1 / freq

    # Try to find the last time point of current Y to start prediction from
    last_time_Y <- NULL
    if (!is.null(colnames(Y))) {
      num_cols <- suppressWarnings(as.numeric(colnames(Y)))
      if (!any(is.na(num_cols))) {
        last_time_Y <- utils::tail(num_cols, 1)
      }
    }

    # If Y colnames are not numeric, fallback to the model's last_time (assuming prediction continues from training)
    # OR if Y seems to be the full training set.
    if (is.null(last_time_Y)) {
      start_time_pred <- last_time + dt
    } else {
      start_time_pred <- last_time_Y + dt
    }

    future_times <- seq(from = start_time_pred, by = dt, length.out = n.ahead)

  }
  # Strategy 2: Infer from numeric colnames of Y (Fallback method)
  else if (!is.null(colnames(Y))) {
    times <- suppressWarnings(as.numeric(colnames(Y)))
    if (!any(is.na(times)) && length(times) >= 2) {
      n_check <- min(length(times), 10L)
      recent_times <- utils::tail(times, n_check)
      diffs <- diff(recent_times)
      if (length(diffs) > 0) {
        tol <- 0.01 * abs(mean(diffs)) + 1e-8  # ~1% tolerance
        if (max(diffs) - min(diffs) < tol) {
          dt <- mean(diffs)
          last_time <- utils::tail(times, 1)
          future_times <- last_time + seq_len(n.ahead) * dt
        }
      }
    }
  }

  if (!is.null(future_times)) {
    # Rounding to keep colnames clean (consistent with nmfkc.ar)
    time_names <- as.character(round(future_times, 5))
    colnames(preds) <- time_names
  }

  return(list(pred = preds, time = future_times))
}










#' @title Optimize lag order for the autoregressive model
#' @description
#' \code{nmfkc.ar.degree.cv} selects the optimal lag order for an autoregressive model
#' by applying cross-validation over candidate degrees.
#'
#' This function accepts both standard matrices (Variables x Time) and \code{ts} objects
#' (Time x Variables). \code{ts} objects are automatically transposed internally.
#'
#' @param Y Observation matrix \eqn{Y(P,N)} or a \code{ts} object.
#' @param rank Rank of the basis matrix. For backward compatibility,
#'   \code{Q} is accepted via \code{...}.
#' @param degree A vector of candidate lag orders to be evaluated.
#' @param intercept Logical. If TRUE (default), an intercept is added to the covariate matrix.
#' @param plot Logical. If TRUE (default), a plot of the objective function values is drawn.
#' @param ... Additional arguments passed to \code{nmfkc.cv}.  Also accepts
#'   \code{cores} (\code{getOption("mc.cores", 1L)}) to evaluate the candidate
#'   degrees in parallel; results are identical to the sequential loop for any
#'   \code{cores} (PSOCK cluster on Windows, forking elsewhere).
#'
#' @return A list with components:
#' \item{degree}{The lag order that minimizes the cross-validation objective function.}
#' \item{degree.max}{Maximum recommended lag order, computed as \eqn{10 \log_{10}(N)}
#'   following the \code{ar} function in the \pkg{stats} package.}
#' \item{objfunc}{Objective function values for each candidate lag order.}
#' @seealso \code{\link{nmfkc.ar}}, \code{\link{nmfkc.cv}}
#' @export
#' @examples
#' # Example using ts object directly
#' d <- AirPassengers
#'
#' # Selection of degree (using ts object)
#' # Note: Y is automatically transposed if it is a ts object
#' nmfkc.ar.degree.cv(Y=d, rank=1, degree=11:14)

nmfkc.ar.degree.cv <- function(Y, rank=1, degree=1:2, intercept=TRUE, plot=TRUE, ...){
  extra_cv <- list(...)
  if (!is.null(extra_cv$Q)) rank <- extra_cv$Q
  Q <- rank

  # --- 1. Time Series Object Handling & Standardization ---
  # Ensure Y is (P x N) matrix before processing
  if(stats::is.ts(Y)){
    # ts objects are Time x Var (N x P) -> Transpose to Var x Time (P x N)
    Y <- t(as.matrix(Y))
  } else {
    if(is.vector(Y)) Y <- matrix(Y, nrow=1)
    if(!is.matrix(Y)) Y <- as.matrix(Y)
  }

  # --- 2. Main Processing ---
  # Remove backward-compat aliases from extra_args to avoid duplicate argument errors
  cores <- if (!is.null(extra_cv$cores)) extra_cv$cores else getOption("mc.cores", 1L)
  extra_args <- extra_cv[!names(extra_cv) %in% c("Q", "cores")]

  # Each degree is an independent nmfkc.ar + (shuffle=FALSE, self-seeded)
  # nmfkc.cv fit; .nmfkc.parlapply preserves input order so the objfunc vector,
  # first-min selection and plot are identical to the sequential loop for any
  # `cores`. Timing / progress messages stay inside the per-degree task.
  run_degree <- function(i){
    start.time <- Sys.time()
    current_degree <- degree[i]
    message(paste0("degree=", current_degree, "..."), appendLF = FALSE)

    res <- tryCatch({
      # 1) build lagged Y & A
      # Y is already standardized to (P x N), so nmfkc.ar handles it as a matrix
      a <- nmfkc.ar(Y = Y, degree = current_degree, intercept = intercept)

      # 2) prepare args for CV (block CV for time series)
      main_args <- list(Y = a$Y, A = a$A, Q = Q)
      all_args  <- c(extra_args, main_args, list(shuffle = FALSE))

      # 3) run CV
      result.cv <- suppressMessages(do.call("nmfkc.cv", all_args))

      # 4) success payload
      list(ok = TRUE, obj = result.cv$objfunc, err = NULL)
    }, error = function(e){
      warning(paste0("Skipping degree=", current_degree, " due to error: ", e$message), call. = FALSE)
      list(ok = FALSE, obj = NA_real_, err = e)
    })

    end.time <- Sys.time()
    diff.time <- difftime(end.time, start.time, units = "sec")
    diff.time.st <- ifelse(diff.time <= 180, paste0(round(diff.time, 1), "sec"),
                           paste0(round(diff.time/60, 1), "min"))
    message(if (res$ok) diff.time.st else "Skipped (Error)")
    list(obj = res$obj, ok = res$ok)
  }

  results <- .nmfkc.parlapply(seq_along(degree), run_degree,
                              cores = if (length(degree) > 1L) cores else 1L)
  objfuncs <- vapply(results, function(z) z$obj, numeric(1))
  success_status <- vapply(results, function(z) z$ok, logical(1))

  valid_indices <- which(success_status)
  if(length(valid_indices) == 0) stop("Cross-validation failed for all candidate degrees.")

  i0 <- valid_indices[which.min(objfuncs[valid_indices])]
  best.degree <- degree[i0]

  # Calculate degree.max based on the standardized Y (ncol = Time points N)
  degree.max <- min(ncol(Y), floor(10 * log10(ncol(Y))))

  if(plot){
    plot(degree, objfuncs, type = "l", col = 2,
         xlab = paste0("degree (max=", degree.max, ")"),
         ylab = "objfunc")
    graphics::points(degree[valid_indices], objfuncs[valid_indices], cex = 1, col = 2)
    graphics::points(degree[i0], objfuncs[i0], cex = 3, col = 2)
    graphics::text(degree, objfuncs, degree, pos = 3)
  }

  names(objfuncs) <- degree
  list(degree = best.degree, degree.max = degree.max, objfunc = objfuncs)
}









#' Split a fitted NMF-VAR object into its VAR parts (Internal)
#'
#' Parser used by \code{\link{nmfkc.ar.stationarity}}.  The covariate matrix built by
#' \code{\link{nmfkc.ar}} stacks the lag-\eqn{d} block of \eqn{Y} in rows
#' \eqn{(d-1)P+1,\dots,dP} and (optionally) puts the intercept in the last
#' row, so the columns of \eqn{\Theta} (\code{x$C}) follow the same layout.
#'
#' \code{nmfkc.ar} records \code{degree}, \code{intercept} and
#' \code{function.name} as attributes of \eqn{A}, which \code{\link{nmfkc}}
#' keeps in \code{x$A.attr}.  Those are authoritative and are used whenever
#' present (the same policy as \code{\link{nmfkc.ar.predict}}); the layout is
#' only inferred from the dimensions of \eqn{\Theta} as a fallback.  That
#' fallback is genuinely ambiguous when \eqn{P = 1}, since \eqn{PD} and
#' \eqn{PD+1} are then both divisible by \eqn{P}, so it warns.
#'
#' @param x A fitted \code{"nmfkc"} object from a \code{\link{nmfkc.ar}} design.
#' @return A list with \code{X} (\eqn{P\times Q}), \code{Theta.lags}
#'   (\eqn{Q\times PD}), \code{theta0} (\eqn{Q\times 1} intercept or
#'   \code{NULL}), \code{P}, \code{Q}, \code{D}, \code{intercept}.
#' @keywords internal
#' @noRd
.nmfvar.parts <- function(x) {
  if (!inherits(x, "nmfkc")) stop("'x' must be an object of class 'nmfkc'.")
  if (is.null(x$X) || is.null(x$C)) stop("'x' must contain 'X' and 'C' components (from nmfkc.ar).")
  X <- x$X                                   # P x Q
  Theta <- x$C                               # Q x (P*D [+1])
  P <- nrow(X); Q <- ncol(X)
  if (P == 0 || Q == 0) stop("'X' must have positive dimensions.")
  total_cols <- ncol(Theta)

  att <- x$A.attr
  fname <- if (!is.null(att)) att$function.name else NULL
  if (!is.null(fname) && !identical(as.character(fname), "nmfkc.ar"))
    stop("'x' was not fitted on an nmfkc.ar() design (A came from ",
         as.character(fname)[1], "); the VAR interpretation does not apply.")

  if (!is.null(att) && !is.null(att$degree)) {
    ## Metadata recorded by nmfkc.ar(): authoritative, and unambiguous at P = 1.
    D <- as.integer(att$degree)
    has_intercept <- if (!is.null(att$intercept)) as.logical(att$intercept)
                     else (total_cols == P * D + 1L)
    if (total_cols != P * D + as.integer(has_intercept))
      stop("Dimensions of 'C' (", total_cols, " columns) are inconsistent with the ",
           "recorded design (P = ", P, ", degree = ", D,
           ", intercept = ", has_intercept, ").")
  } else {
    ## Fallback: infer from dimensions.  P * D and P * D + 1 are both multiples
    ## of P when P = 1, so the split cannot be recovered there -- assume an
    ## intercept (the nmfkc.ar() default) and say so.
    has_intercept <- (total_cols - 1) %% P == 0
    if (P == 1L && has_intercept)
      warning("P = 1 and no nmfkc.ar() metadata found: cannot tell an intercept ",
              "column from a further lag; assuming an intercept is present ",
              "(degree = ", total_cols - 1L, ").")
    D <- if (has_intercept) (total_cols - 1) %/% P else total_cols %/% P
  }
  if (D <= 0) stop("Cannot determine lag order D from dimensions of 'C'.")

  list(X = X,
       Theta.lags = if (has_intercept) Theta[, 1:(total_cols - 1), drop = FALSE] else Theta,
       theta0 = if (has_intercept) Theta[, total_cols, drop = FALSE] else NULL,
       P = P, Q = Q, D = D, intercept = has_intercept)
}


#' @title Latent transition matrices of an NMF-VAR model
#' @description
#' This function is \strong{experimental}. The interface may change in future
#' versions: argument names, defaults and the contents of the returned object
#' are not yet stable.
#'
#' Reduces a fitted NMF-VAR to the \eqn{Q\times Q} operators that drive its
#' coefficient vector.  Writing \eqn{\bm b_t=\Theta\bm a_t} for the fitted
#' scores, the residual-free recursion is
#' \deqn{\bm b_t=\sum_{d=1}^{D}G_d\,\bm b_{t-d}+\bm\theta,\qquad G_d=\Theta_d X
#'   \ \ (Q\times Q),}
#' and the \eqn{G_d} carry the model's whole deterministic structure ---
#' transition weights, propagation, roots, the long-run mean.  They stay
#' readable where the \eqn{P\times P} matrices \eqn{\Xi_d = X\Theta_d} do not:
#' \eqn{4\times4} instead of \eqn{47\times47}.
#'
#' Use \code{\link{plot.nmfkc.ar.latent}} to draw them,
#' \code{\link{nmfkc.ar.latent.inference}} to attach standard errors, and
#' \code{\link{nmfkc.ar.stationarity}} if the stationarity verdict is all that
#' is wanted.
#'
#' @details
#' \strong{The latent process is VARMA, not VAR.}  The display above drops the
#' observation residual.  With \eqn{\bm y_t = X\bm b_t + \bm e_t} the exact
#' identity is
#' \deqn{\bm b_t=\sum_{d=1}^{D}G_d\,\bm b_{t-d}+\bm\theta
#'   +\sum_{d=1}^{D}\Theta_d\,\bm e_{t-d},}
#' i.e. a VARMA\eqn{(D,D)} in general (verified to machine precision; dropping
#' the moving-average term leaves an error of order the residuals themselves).
#' \eqn{G_d} is therefore the \strong{autoregressive operator of the fitted VAR
#' reduced to the latent space}, not the coefficient matrix of a latent VAR that
#' the scores obey exactly.  In particular an entry of \eqn{G_d} should not be
#' read as Granger causality between latent conditions: the moving-average term
#' is common to all coordinates and is not conditioned on.
#'
#' \strong{Direction convention}: \eqn{(G_d)_{q q'}} is the effect of
#' \eqn{b_{q',\,t-d}} on \eqn{b_{q,\,t}} --- \strong{row = effect} (at \eqn{t}),
#' \strong{column = cause} (at \eqn{t-d}).
#'
#' \strong{Identification.}  The factorization is unique only up to
#' \eqn{(X,\Theta)\to(XT,T^{-1}\Theta)}, under which \eqn{G_d\to T^{-1}G_dT}.
#' The eigenvalues, \eqn{\rho}, \eqn{\Xi_d} and \eqn{\bm\mu_y} are therefore
#' invariant, while individual entries of \eqn{G_d} are not.  Column
#' normalization removes only the scale part of \eqn{T}, and a separable fitted
#' \eqn{X} does not certify the rest: with \eqn{X=I_2},
#' \eqn{\Theta=\left(\begin{smallmatrix}2&3\\1&1\end{smallmatrix}\right)} and
#' \eqn{T=\left(\begin{smallmatrix}0.9&0.1\\0.1&0.9\end{smallmatrix}\right)} the
#' basis is perfectly separable and column-normalized, yet \eqn{XT\ge0},
#' \eqn{T^{-1}\Theta\ge0}, the product is unchanged and \eqn{T} is no
#' permutation.  Read the entries of \eqn{G_d} as descriptive and the invariants
#' as inferential; see \code{\link{nmfkc.ar.latent.inference}}.
#'
#' \strong{Oscillation.}  A complex pair means a cycle, but a \emph{negative
#' real} root also oscillates --- it alternates in sign each period, i.e. with
#' period 2.  \code{cycle.period} covers the first case and \code{alternating}
#' the second; \code{non.oscillatory} is \code{TRUE} only when neither occurs.
#' Positive real roots rule out oscillation but \strong{not overshoot}: a
#' non-normal \eqn{G} can send a coordinate up before it comes down, which is
#' exactly what the off-diagonal impulse responses of the Canada fit do.  Note
#' also that for \eqn{Q=2}, \eqn{D=1} and \eqn{G\ge0} the roots are real by
#' construction (the discriminant is \eqn{(a-d)^2+4bc\ge0}), so in that setting
#' the absence of a cycle is a property of the design, not a finding.
#'
#' @param x A fitted \code{"nmfkc"} object obtained from an
#'   \code{\link{nmfkc.ar}} design (\code{nmfkc(ar$Y, ar$A, rank = Q)}).
#' @return An object of class \code{"nmfkc.ar.latent"}: a list with
#' \item{G}{List of \eqn{D} latent transition matrices \eqn{G_d}
#'   (\eqn{Q\times Q}), named \code{lag1}, ..., with \code{dimnames} = basis
#'   labels (row = effect, column = cause).}
#' \item{G.sum}{\eqn{\sum_d G_d}.}
#' \item{eigenvalues}{The \eqn{DQ} eigenvalues of the latent companion matrix.
#'   The \eqn{PD\times PD} companion of the \eqn{\Xi_d} has the same non-zero
#'   eigenvalues, so these are all of the model's roots.}
#' \item{spectral.radius}{\eqn{\rho}, the largest modulus among them.}
#' \item{spectral.radius.sum}{\eqn{\rho(\sum_d G_d)}; below 1 iff stationary.}
#' \item{stationary}{Logical; \code{TRUE} when \code{spectral.radius < 1}.}
#' \item{cycle.period}{Period \eqn{2\pi/|\arg\lambda|} implied by the dominant
#'   complex root, or \code{NA} when no root is complex.}
#' \item{alternating}{\code{TRUE} when some root is real and negative, i.e. the
#'   decay alternates in sign (period 2).}
#' \item{non.oscillatory}{\code{TRUE} only when every root is real and positive.
#'   Named after the roots, not the trajectory --- overshoot is still possible.}
#' \item{theta0}{Latent intercept \eqn{\bm\theta} (or \code{NULL}).}
#' \item{mu.b, mu.y}{Long-run means \eqn{(I-\sum_d G_d)^{-1}\bm\theta} and
#'   \eqn{X\bm\mu_b}; \code{NA} when the fit is not stationary or has no
#'   intercept.}
#' \item{p.star}{Composition of the long-run mean,
#'   \eqn{\bm\mu_b/(\bm 1'\bm\mu_b)} --- equivalently the fixed point of the
#'   \emph{residual-free} composition dynamics.  It is \strong{not} in general
#'   the limit of \eqn{\bm b_t/(\bm 1'\bm b_t)} for the stochastic process, nor
#'   its expectation: the ratio of expectations is not the expectation of the
#'   ratio.}
#' \item{colsum}{Column sums \eqn{c_j} of \eqn{\sum_d\Xi_d = X\sum_d\Theta_d},
#'   computed as \code{colSums(X) \%*\% Lambda} (length \eqn{P}).  Being column
#'   sums of a non-negative matrix whose radius is \eqn{\rho(\sum_d G_d)}, they
#'   bracket it.  The unweighted \code{colSums(Lambda)} would only be correct
#'   when \eqn{\bm 1'X=\bm 1'}.}
#' \item{colsum.max}{\eqn{\max_j c_j}; below 1 it certifies stationarity.}
#' \item{df}{Effective dimension \eqn{Q(P+m-Q)} of the fit, where \eqn{m} is the
#'   number of columns of \eqn{\Theta} (\eqn{PD+1} with an intercept, \eqn{PD}
#'   without).  These are the degrees of freedom to use when comparing
#'   \eqn{R^2} or information criteria across \eqn{Q}: the naive count
#'   \eqn{Q(P+m)} ignores the \eqn{Q^2} rotational indeterminacy.}
#' \item{separability}{Per-basis maximum of the row-normalised \eqn{X}; a value
#'   near 1 means that basis owns a variable loading almost only on it.  This is
#'   a \emph{necessary}-looking but not sufficient condition for uniqueness ---
#'   see the counterexample under Identification.}
#' \item{X, Theta}{The basis and the full coefficient matrix the quantities
#'   above were computed from, carried along so that
#'   \code{\link{nmfkc.ar.latent.inference}} does not need the fitted object a
#'   second time.}
#' \item{X.restriction}{The constraint \code{\link{nmfkc}} placed on the columns
#'   of \eqn{X}, which decides how much of \eqn{T} is pinned down.}
#' \item{dims}{Named vector of \code{P}, \code{Q}, \code{D}.}
#' @seealso \code{\link{plot.nmfkc.ar.latent}},
#'   \code{\link{nmfkc.ar.latent.inference}},
#'   \code{\link{nmfkc.ar.stationarity}}, \code{\link{nmfkc.ar}},
#'   \code{\link{nmfkc.ar.DOT}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to multivariate time series data as a vector autoregression model.
#'   \emph{Japanese Journal of Statistics and Data Science}. arXiv:2501.17446.
#'   \doi{10.1007/s42081-025-00314-0}
#' @examples
#' d <- AirPassengers
#' ar_data <- nmfkc.ar(d, degree = 2)
#' result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#' lat <- nmfkc.ar.latent(result)
#' lat$G$lag1          # latent transition matrix of lag 1
#' lat$eigenvalues     # all roots of the model
#' lat                 # full report
#'
#' @export
nmfkc.ar.latent <- function(x) {
  prt <- .nmfvar.parts(x)
  X <- prt$X; Theta_lags <- prt$Theta.lags
  P <- prt$P; Q <- prt$Q; D <- prt$D

  labels <- colnames(X)
  if (is.null(labels)) labels <- paste0("Basis", seq_len(Q))

  ## The rho(sum_d G_d) < 1 equivalence and the column-sum bracket rest on
  ## Perron-Frobenius, i.e. on Xi_d = X Theta_d >= 0.  That holds for an NMF
  ## fit, but a signed variant inherits class "nmfkc" too, so check.
  if (any(X < 0) || any(Theta_lags < 0))
    warning("X or Theta has negative entries: the eigenvalues are still exact, ",
            "but the rho(sum_d G_d) equivalence and the column-sum bound ",
            "assume non-negative coefficients and may not hold.")

  ## G_d = Theta_d X : row = effect at t, column = cause at t-d
  G <- lapply(seq_len(D), function(d) {
    g <- Theta_lags[, ((d - 1) * P + 1):(d * P), drop = FALSE] %*% X
    dimnames(g) <- list(labels, labels)
    g
  })
  names(G) <- paste0("lag", seq_len(D))

  Lam <- Theta_lags[, seq_len(P), drop = FALSE]                 # sum_d Theta_d
  if (D > 1) for (d in 2:D)
    Lam <- Lam + Theta_lags[, ((d - 1) * P + 1):(d * P), drop = FALSE]
  G.sum <- Lam %*% X
  dimnames(G.sum) <- list(labels, labels)

  ## All of the model's roots.  The DQ x DQ companion of the G_d and the
  ## PD x PD companion of the Xi_d share their non-zero eigenvalues (the latter
  ## has PD - DQ extra zeros), so the cheap route loses nothing.
  ev  <- .nmfvar.roots(G, Q, D)
  rho <- max(Mod(ev))
  rho.sum <- max(Mod(eigen(G.sum, only.values = TRUE)$values))
  stationary <- (rho < 1)

  ## A cycle needs a complex pair, but oscillation does not: a NEGATIVE real
  ## root alternates in sign every period, i.e. it oscillates with period 2.
  ## "all roots real" therefore does not mean monotone decay.
  cplx <- ev[abs(Im(ev)) > 1e-8]
  period <- if (length(cplx)) {
    k <- which.max(Mod(cplx))
    2 * pi / abs(Arg(cplx[k]))
  } else NA_real_
  alternating <- any(abs(Im(ev)) <= 1e-8 & Re(ev) < -1e-8)
  ## Positive real roots rule out OSCILLATION, not overshoot: a non-normal G
  ## can send a component up before it comes down (the off-diagonal IRF of the
  ## Canada fit peaks at h = 3 with both roots real and positive).  So this is
  ## deliberately named after the roots, not after the trajectory.
  non.oscillatory <- !length(cplx) && !alternating

  ## Perron-Frobenius bracket.  The bound applies to the column sums of the
  ## non-negative matrix whose radius we want, i.e. of sum_d Xi_d = X Lambda:
  ##   1_P' (X Lambda) = (1_P' X) Lambda = colSums(X) %*% Lambda,
  ## and rho(X Lambda) = rho(Lambda X) = rho(sum_d G_d).  Using colSums(Lambda)
  ## on its own is the special case colSums(X) = 1 and is WRONG otherwise:
  ## with X = (20, 20)' and Lambda = (0.04, 0.04) it reports max_j c_j = 0.04
  ## and would certify stationarity, while rho = 1.6.
  cs <- as.numeric(colSums(X) %*% Lam)
  names(cs) <- colnames(Lam)

  ## Long-run means: defined only for a stationary fit with an intercept.
  theta0 <- prt$theta0
  mu.b <- NA_real_; mu.y <- NA_real_; p.star <- NA_real_
  if (!is.null(theta0)) {
    if (stationary) {
      mu.b <- solve(diag(Q) - G.sum, as.numeric(theta0))
      names(mu.b) <- labels
      mu.y <- as.numeric(X %*% mu.b)
      names(mu.y) <- rownames(X)
      if (sum(mu.b) > 0) {
        p.star <- mu.b / sum(mu.b)
        names(p.star) <- labels
      }
    } else {
      warning("The fit is not stationary (spectral radius >= 1); long-run means are NA.")
    }
  }

  ## Separability diagnostic: rows of X normalised to sum 1; a column whose
  ## maximum is close to 1 owns an (almost) pure observation.
  rs <- rowSums(X)
  sep <- if (all(rs > 0)) apply(X / rs, 2, max) else rep(NA_real_, Q)
  names(sep) <- labels

  structure(list(
    G = G, G.sum = G.sum,
    eigenvalues = ev, spectral.radius = rho,
    spectral.radius.sum = rho.sum, stationary = stationary,
    cycle.period = period, alternating = alternating,
    non.oscillatory = non.oscillatory,
    theta0 = theta0, mu.b = mu.b, mu.y = mu.y, p.star = p.star,
    colsum = cs, colsum.max = max(cs),
    ## PQ (basis) + Q*ncol(Theta) (coefficients) - Q^2 (rotational
    ## indeterminacy); ncol(Theta) is PD + 1 with an intercept, PD without.
    df = Q * (P + ncol(Theta_lags) + as.integer(prt$intercept) - Q),
    separability = sep,
    ## The pieces the fit was read from, so downstream callers (notably
    ## nmfkc.ar.latent.inference()) do not need the fitted object again.
    ## Theta is the FULL coefficient matrix (lags and, when present, the
    ## intercept column), i.e. the x$C that was parsed.
    X = X, Theta = x$C,
    ## How much of the rotational freedom nmfkc() pinned down; entry-level
    ## inference on G_d needs the column scale fixed.
    X.restriction = if (is.null(x$X.restriction)) NA_character_ else x$X.restriction,
    ## Settings a bootstrap re-fit must reproduce, or it would resample a
    ## DIFFERENT estimator than the one being described.  Only what nmfkc()
    ## records can be recovered here; anything else the user passed (X.init,
    ## penalties, ...) has to be supplied through refit.args=.
    ## With X.restriction = "fixed" the basis is held at X.init throughout, so
    ## the returned X *is* that matrix (verified bit-identical).  Carry it, or a
    ## re-fit would hold X at whatever nmfkc() initialises to instead -- a
    ## different matrix entirely, and therefore a different estimator.
    fit.args = c(list(method = x$method, X.restriction = x$X.restriction),
                 if (isTRUE(x$X.restriction == "fixed")) list(X.init = x$X)),
    call = x$call,
    dims = c(P = P, Q = Q, D = D)),
    class = "nmfkc.ar.latent")
}


#' @title Is the factorization unique up to a permutation? (Internal)
#' @description
#' Two one-sided conditions together force \eqn{T} in
#' \eqn{(X,\Theta)\to(XT,T^{-1}\Theta)} to be a permutation:
#' \itemize{
#'   \item an \strong{anchor row} for every column of \eqn{X} (a variable loading
#'     on that basis alone) makes \eqn{XT\ge0} imply \eqn{T\ge0}, because the
#'     anchor row of \eqn{XT} is a row of \eqn{T};
#'   \item a \strong{pure column} for every row of \eqn{\Theta} (a covariate
#'     whose coefficient is non-zero for that basis alone) makes
#'     \eqn{T^{-1}\Theta\ge0} imply \eqn{T^{-1}\ge0}, because that column of
#'     \eqn{T^{-1}\Theta} is a multiple of a column of \eqn{T^{-1}}.
#' }
#' A non-negative matrix whose inverse is also non-negative is monomial, and the
#' column-sum constraint fixes the remaining scale, so \eqn{T} is a permutation.
#'
#' Either side alone is not enough.  With \eqn{X=I_2},
#' \eqn{\Theta=[2\ 3;\ 1\ 1]} and \eqn{T=[0.9\ 0.1;\ 0.1\ 0.9]} the anchor
#' condition holds and \eqn{T\ge0}, but \eqn{\Theta} has no pure column,
#' \eqn{T^{-1}} has negative entries, and both factorizations are admissible.
#' @param X Basis matrix (\eqn{P\times Q}).
#' @param Theta Full coefficient matrix (\eqn{Q\times m}).
#' @param sep.tol A column of \eqn{X} counts as anchored when the largest
#'   row-normalised loading reaches \code{sep.tol}.
#' @param zero.tol Entries of \eqn{\Theta} at or below this are treated as zero
#'   when looking for pure columns.
#' \strong{This is an approximate diagnostic, not a certificate.}  The argument
#' above needs an \emph{exact} anchor row (loading exactly 1 after row
#' normalisation) and an \emph{exact} pure column (the other entries exactly 0).
#' Multiplicative updates give small positive leakage instead, so a fit that
#' clears \code{sep.tol} / \code{zero.tol} is merely close to a configuration the
#' theorem covers.  The returned margins say how close.
#' @return List with \code{anchor} (logical, per column of \eqn{X}), \code{pure}
#'   (logical, per row of \eqn{\Theta}), \code{ok}, and the two margins
#'   \code{anchor.margin} (largest \eqn{1-} loading among the anchors, 0 when
#'   exact) and \code{pure.margin} (largest leaked entry in the pure columns,
#'   0 when exact).
#' @keywords internal
#' @noRd
.nmfvar.identified <- function(X, Theta, sep.tol = 0.999, zero.tol = 1e-3) {
  Q  <- ncol(X)
  rs <- rowSums(X)
  Xn <- if (all(rs > 0)) X / rs else X
  best <- vapply(seq_len(Q), function(q) {
    v <- Xn[, q]; v <- v[is.finite(v)]
    if (length(v)) max(v) else NA_real_
  }, numeric(1))
  anchor <- !is.na(best) & best >= sep.tol
  ## How far the best anchor is from an exact one (0 = exact).
  anchor.margin <- if (any(anchor)) max(1 - best[anchor]) else NA_real_

  nz   <- Theta > zero.tol
  pcol <- lapply(seq_len(nrow(Theta)), function(q)
    which(nz[q, ] & colSums(nz) == 1L))
  pure <- vapply(pcol, function(j) length(j) > 0L, logical(1))
  ## Largest entry that the pure-column test treated as zero (0 = exact).
  pure.margin <- if (any(pure)) {
    j <- unlist(pcol)
    off <- Theta[, j, drop = FALSE][!nz[, j, drop = FALSE]]
    if (length(off)) max(off) else 0
  } else NA_real_

  list(anchor = anchor, pure = pure, ok = all(anchor) && all(pure),
       anchor.margin = anchor.margin, pure.margin = pure.margin)
}


#' @title Roots of the latent companion matrix (Internal)
#' @param G List of \eqn{D} latent transition matrices \eqn{G_d}.
#' @param Q Number of bases.
#' @param D Number of lags.
#' @return The \eqn{DQ} eigenvalues of the companion matrix.
#' @keywords internal
#' @noRd
.nmfvar.roots <- function(G, Q, D) {
  comp <- matrix(0, nrow = D * Q, ncol = D * Q)
  for (d in seq_len(D)) comp[1:Q, ((d - 1) * Q + 1):(d * Q)] <- G[[d]]
  if (D > 1) comp[(Q + 1):(D * Q), 1:(Q * (D - 1))] <- diag(Q * (D - 1))
  eigen(comp, only.values = TRUE)$values
}


#' @title Print method for nmfkc.ar.latent objects
#' @param x An object of class \code{"nmfkc.ar.latent"}.
#' @param digits Number of digits used when printing the matrices.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.nmfkc.ar.latent <- function(x, digits = 3, ...) {
  cat("Latent transition matrices  G_d = Theta_d %*% X",
      "  (row = effect at t, column = cause at t-d)\n\n")
  for (nm in names(x$G)) {
    cat(nm, ":\n", sep = "")
    print(round(x$G[[nm]], digits))
    cat("\n")
  }
  D <- x$dims[["D"]]; Q <- x$dims[["Q"]]
  cat(sprintf("roots: %d (companion %dx%d)   rho = %.4f  ->  %s\n",
              length(x$eigenvalues), D * Q, D * Q, x$spectral.radius,
              if (x$stationary) "stationary" else "NOT stationary"))
  if (is.finite(x$cycle.period))
    cat(sprintf("dominant complex root: period %.1f periods\n", x$cycle.period))
  else if (isTRUE(x$alternating))
    cat("all roots real but one is negative: alternating decay (period 2)\n")
  else
    cat("all roots real and positive: non-oscillatory",
        "(overshoot is still possible)\n")
  if (Q == 2L && D == 1L)
    cat("  (note: with Q = 2, D = 1 and G >= 0 the roots are real by",
        "construction --\n   the absence of a cycle is not a finding)\n")
  if (length(x$mu.b) > 1 || !all(is.na(x$mu.b))) {
    cat("\nLong-run mean of the latent scores (mu.b):\n")
    print(round(x$mu.b, digits))
    if (!all(is.na(x$p.star))) {
      cat("Long-run composition (p*):\n")
      print(round(x$p.star, digits))
    }
  }
  cat(sprintf("\nEffective dimension df = Q(P + PD + 1 - Q) = %d\n", x$df))
  cat("For the stationarity verdict see nmfkc.ar.stationarity(),",
      "for standard errors nmfkc.ar.latent.inference().\n")
  invisible(x)
}


#' @title Stationarity of an NMF-VAR model
#' @description
#' Reports whether a fitted NMF-VAR is stationary, and how cheaply that could be
#' decided.  Three numbers are returned, in increasing order of sharpness:
#' \code{colsum.max} \eqn{=\max_j c_j}, the largest column sum of
#' \eqn{\sum_d\Xi_d = X\sum_d\Theta_d}, certifies stationarity when below 1 with
#' no eigenvalue computation at all (a sufficient condition only, so a value
#' above 1 is merely inconclusive); \code{spectral.radius.sum}
#' \eqn{=\rho(\sum_d G_d)} is below 1 exactly when the model is stationary (and
#' is bracketed by \eqn{\min_j c_j} and \eqn{\max_j c_j}); and
#' \code{spectral.radius} is the spectral radius of the companion matrix itself.
#'
#' Note that the \eqn{c_j} are weighted by \eqn{\bm 1'X}, i.e.
#' \code{colSums(X) \%*\% Lambda}.  Using the raw column sums of
#' \eqn{\sum_d\Theta_d} is correct only when \eqn{\bm 1'X=\bm 1'} and can
#' certify stationarity for a non-stationary fit otherwise.
#'
#' This function answers the yes/no question only.  For the dynamics behind the
#' number --- the transition matrices, the roots, the impulse responses, the
#' long-run mean --- use \code{\link{nmfkc.ar.latent}}, and for the standard
#' error of \eqn{\rho} use \code{\link{nmfkc.ar.latent.inference}}.
#'
#' @param x A fitted \code{"nmfkc"} object obtained from an
#'   \code{\link{nmfkc.ar}} design, or an object from
#'   \code{\link{nmfkc.ar.latent}}.
#' @param method Route used for \code{spectral.radius}.  \code{"latent"}
#'   (default) uses the \eqn{DQ\times DQ} companion matrix of the \eqn{G_d};
#'   \code{"companion"} uses the \eqn{PD\times PD} companion matrix of the
#'   \eqn{\Xi_d}.  The two companion matrices share the same non-zero
#'   eigenvalues, so the results agree up to floating-point rounding, but
#'   \code{"latent"} costs \eqn{O((DQ)^3)} instead of \eqn{O((PD)^3)} --- for a
#'   47-variable, 7-lag, rank-4 model that is a \eqn{28\times28} eigenproblem
#'   instead of a \eqn{329\times329} one.  Use \code{"companion"} to reproduce
#'   the value bit-for-bit as computed by earlier versions of the package.
#' @return An object of class \code{"nmfkc.ar.stationarity"}: a list with
#' \item{spectral.radius}{Spectral radius of the companion matrix (route
#'   selected by \code{method}).  A value below 1 indicates stationarity.}
#' \item{stationary}{Logical; \code{TRUE} when \code{spectral.radius < 1}.}
#' \item{method}{The route actually used.}
#' \item{spectral.radius.sum}{\eqn{\rho(\sum_d G_d)}; below 1 iff stationary.}
#' \item{colsum}{Column sums \eqn{c_j} of \eqn{\sum_d\Xi_d}, i.e.
#'   \code{colSums(X) \%*\% Lambda} (length \eqn{P}).}
#' \item{colsum.max}{\eqn{\max_j c_j}; below 1 it certifies stationarity.}
#' \item{dims}{Named vector of \code{P}, \code{Q}, \code{D}.}
#' @seealso \code{\link{nmfkc.ar.latent}},
#'   \code{\link{nmfkc.ar.latent.inference}}, \code{\link{nmfkc.ar}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to multivariate time series data as a vector autoregression model.
#'   \emph{Japanese Journal of Statistics and Data Science}. arXiv:2501.17446.
#'   \doi{10.1007/s42081-025-00314-0}
#' @examples
#' # Check stationarity of fitted AR model
#' d <- AirPassengers
#' ar_data <- nmfkc.ar(d, degree = 2)
#' result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#' st <- nmfkc.ar.stationarity(result)
#' st$spectral.radius
#' st$stationary
#' st                 # full report
#'
#' @export
nmfkc.ar.stationarity <- function(x, method = c("latent", "companion")) {
  method <- match.arg(method)
  lat <- if (inherits(x, "nmfkc.ar.latent")) x else nmfkc.ar.latent(x)
  P <- lat$dims[["P"]]; D <- lat$dims[["D"]]

  rho <- if (method == "latent") lat$spectral.radius else {
    ## The PD x PD route the package used before the latent representation
    ## existed, kept so that a published number reproduces bit-for-bit.
    X <- lat$X
    Theta_lags <- lat$Theta[, seq_len(P * D), drop = FALSE]
    Xi <- lapply(seq_len(D), function(d)
      X %*% Theta_lags[, ((d - 1) * P + 1):(d * P), drop = FALSE])
    cm <- matrix(0, nrow = P * D, ncol = P * D)
    for (d in seq_len(D)) cm[1:P, ((d - 1) * P + 1):(d * P)] <- Xi[[d]]
    if (D > 1) cm[(P + 1):(P * D), 1:(P * (D - 1))] <- diag(P * (D - 1))
    max(Mod(eigen(cm, only.values = TRUE)$values))
  }

  structure(list(
    spectral.radius = rho, stationary = (rho < 1), method = method,
    spectral.radius.sum = lat$spectral.radius.sum,
    colsum = lat$colsum, colsum.max = lat$colsum.max,
    dims = lat$dims),
    class = "nmfkc.ar.stationarity")
}


#' @title Print method for nmfkc.ar.stationarity objects
#' @param x An object of class \code{"nmfkc.ar.stationarity"}.
#' @param digits Ignored; kept for compatibility with the print generic.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.nmfkc.ar.stationarity <- function(x, digits = 3, ...) {
  nn <- if (x$method == "latent") x$dims[["D"]] * x$dims[["Q"]]
        else x$dims[["P"]] * x$dims[["D"]]
  cat(sprintf("rho(companion, %dx%d, method=\"%s\") = %.4f  ->  %s\n",
              nn, nn, x$method, x$spectral.radius,
              if (x$stationary) "stationary" else "NOT stationary"))
  cat(sprintf("rho(sum_d G_d)        = %.4f   (< 1 iff stationary)\n",
              x$spectral.radius.sum))
  cat(sprintf("bracket: min_j c_j = %.4f  <=  rho(sum_d G_d)  <=  max_j c_j = %.4f\n",
              min(x$colsum), x$colsum.max))
  cat(if (x$colsum.max < 1)
        "   max_j c_j < 1: stationarity certified without an eigenvalue\n"
      else "   max_j c_j >= 1: the cheap bound is inconclusive, rho decides\n")
  cat("\nFor the dynamics behind this number see nmfkc.ar.latent(),",
      "for a test of H0: rho >= 1 see nmfkc.ar.latent.inference().\n")
  invisible(x)
}


#' @title Plot the stationarity bracket of an NMF-VAR model
#' @description
#' Draws the Perron-Frobenius bracket \eqn{\min_j c_j\le\rho(\sum_d
#' G_d)\le\max_j c_j} against the unit boundary, where \eqn{c_j} are the column
#' sums of \eqn{\sum_d\Theta_d}.  It shows at a glance whether the cheap
#' sufficient condition settles stationarity or an eigenvalue computation is
#' needed.  For the model's dynamics see \code{\link{plot.nmfkc.ar.latent}}.
#' @param x An object from \code{\link{nmfkc.ar.stationarity}}.
#' @param ... Passed to the underlying plot call.
#' @return \code{invisible(NULL)}; called for the plot.
#' @seealso \code{\link{nmfkc.ar.stationarity}},
#'   \code{\link{plot.nmfkc.ar.latent}}
#' @examples
#' set.seed(1)
#' Y   <- matrix(abs(rnorm(4 * 60)) + 1, 4, 60)
#' ar  <- nmfkc.ar(Y, degree = 2)
#' fit <- nmfkc(ar$Y, ar$A, rank = 2, verbose = FALSE)
#' plot(nmfkc.ar.stationarity(fit))
#' @export
plot.nmfkc.ar.stationarity <- function(x, ...) {
  extra <- base::list(...)
  dflt <- function(...) {
    d <- base::list(...)
    c(d[base::setdiff(base::names(d), base::names(extra))], extra)
  }
  cj <- x$colsum
  base::do.call("plot", dflt(x = NA, y = NA,
                             xlim = c(0, base::max(1.25, base::max(cj) * 1.1)),
                             ylim = c(0.4, 1.6), yaxt = "n", ylab = "",
                             xlab = "spectral radius / column sums",
                             main = "Stationarity bracket:  min c_j <= rho <= max c_j"))
  graphics::abline(v = 1, col = "firebrick", lwd = 2, lty = 2)
  graphics::segments(base::min(cj), 1, base::max(cj), 1, lwd = 8,
                     col = "grey80", lend = 1)
  graphics::points(cj, base::rep(1, base::length(cj)), pch = "|", cex = 1.3,
                   col = "grey35")
  graphics::points(x$spectral.radius.sum, 1, pch = 19, col = "steelblue4", cex = 1.7)
  graphics::text(x$spectral.radius.sum, 1.18,
                 base::sprintf("rho = %.3f", x$spectral.radius.sum),
                 col = "steelblue4", cex = 0.9)
  graphics::title(sub = if (x$colsum.max < 1)
                    "max c_j < 1: stationarity certified without an eigenvalue"
                  else "max c_j >= 1: the cheap bound is inconclusive, rho decides",
                  cex.sub = 0.85)
  base::invisible(NULL)
}


#' @title Permutation number \code{i} of \code{seq_len(Q)} (Internal)
#' @param Q Size of the permuted set.
#' @param i Index in lexicographic order, \code{1 .. factorial(Q)}.
#' @return An integer permutation vector of length \code{Q}.
#' @keywords internal
#' @noRd
.nmfvar.perm <- function(Q, i) {
  el <- seq_len(Q); out <- integer(Q); i <- i - 1L
  for (k in seq_len(Q)) {
    f <- factorial(Q - k); j <- i %/% f + 1L; i <- i %% f
    out[k] <- el[j]; el <- el[-j]
  }
  out
}


#' @title Plot the latent dynamics of an NMF-VAR model
#' @description
#' Draws the model's own dynamics --- what the fit \emph{does} --- rather than
#' the realized series.  All of it lives in the \eqn{Q\times Q} latent
#' transition matrices \eqn{G_d=\Theta_d X} returned by
#' \code{\link{nmfkc.ar.latent}}, so every panel stays readable at a size the
#' observed space cannot: a 47-region, 7-lag, rank-4 fit has \eqn{47^2=2209}
#' observed impulse responses but only \eqn{4^2=16} latent ones.
#'
#' @details
#' The available panels are
#' \describe{
#'   \item{\code{"roots"}}{(default) The \eqn{DQ} eigenvalues of the latent
#'     companion matrix on the unit circle.  Complex roots are drawn in a
#'     different colour and their implied period is reported: a fit can sit near
#'     the boundary either because it cycles or because it is seasonal, and the
#'     spectral radius alone does not say which.}
#'   \item{\code{"graph"}}{The transition graph of \eqn{G_d}: self-loops are
#'     persistence, arrows \eqn{q'\to q} are spillover, widths are proportional
#'     to the entries.  Unlike a \eqn{y\to b\to y} bipartite picture this graph
#'     is closed, so an asymmetric feedback loop is visible as such.}
#'   \item{\code{"irf"}}{Latent impulse responses \eqn{\Psi_h}, a \eqn{Q\times Q}
#'     panel.  Since \eqn{G_d\ge 0} the responses never change sign, so the
#'     shape is a decay or a hump and its half-life can be read off.}
#'   \item{\code{"lag"}}{\eqn{(G_d)_{qq'}} as a function of \eqn{d}: how many
#'     periods a shock takes to arrive.  With \eqn{Q=1} this is a single bar
#'     plot of the lag profile.}
#'   \item{\code{"phase"}}{Phase portrait of \eqn{b_t\mapsto Gb_{t-1}+\theta}
#'     with the fixed point \eqn{\mu_b} and the eigen-directions.  Requires
#'     \eqn{D=1} and \eqn{Q=2}.  This is where a monotone relaxation and a true
#'     cycle look different; two coefficient series plotted against time do not
#'     distinguish them.}
#' }
#'
#' @param x An object from \code{\link{nmfkc.ar.latent}}.
#' @param type Which panel to draw; see Details.
#' @param lag Which lag to use for \code{"graph"} and \code{"phase"}.  Default
#'   \code{NULL} means \eqn{\sum_d G_d} for \code{"graph"} (the total one-step
#'   influence) and lag 1 for \code{"phase"}.
#' @param horizon Number of periods for \code{"irf"}.  Default 16.
#' @param ... Passed to the underlying plot call.  \code{"phase"} also accepts
#'   \code{start}, a matrix of initial \eqn{b_0} values (one per row).
#' @return \code{invisible(NULL)}; called for the plot.  \code{"irf"} and
#'   \code{"lag"} set up their own multi-panel layout and restore the previous
#'   one on exit.
#' @seealso \code{\link{nmfkc.ar.latent}},
#'   \code{\link{nmfkc.ar.latent.inference}},
#'   \code{\link{plot.nmfkc.ar.stationarity}}, \code{\link{nmfkc.ar.DOT}}
#' @examples
#' set.seed(1)
#' Y   <- matrix(abs(rnorm(4 * 60)) + 1, 4, 60)
#' ar  <- nmfkc.ar(Y, degree = 3)
#' fit <- nmfkc(ar$Y, ar$A, rank = 2, verbose = FALSE)
#' lat <- nmfkc.ar.latent(fit)
#' plot(lat)                  # roots on the unit circle
#' plot(lat, type = "graph")  # latent transition graph
#' plot(lat, type = "irf")    # latent impulse responses
#' @export
plot.nmfkc.ar.latent <- function(x,
                                 type = c("roots", "graph", "irf", "lag",
                                          "phase"),
                                 lag = NULL, horizon = 16L, ...) {
  type  <- base::match.arg(type)
  extra <- base::list(...)
  Q <- x$dims[["Q"]]; D <- x$dims[["D"]]
  lab <- base::colnames(x$G.sum)
  ## Named defaults the caller may override through `...`.
  dflt <- function(...) {
    d <- base::list(...)
    c(d[base::setdiff(base::names(d), base::names(extra))], extra)
  }

  if (type == "roots") {
    ev <- x$eigenvalues
    cx <- base::abs(base::Im(ev)) > 1e-8
    a  <- base::seq(0, 2 * base::pi, length.out = 400)
    base::do.call("plot", dflt(x = base::cos(a), y = base::sin(a), type = "l",
                               asp = 1, col = "grey55", lty = 2, xlab = "Re",
                               ylab = "Im", xlim = c(-1.15, 1.15),
                               ylim = c(-1.15, 1.15),
                               main = "Roots of the latent companion matrix"))
    graphics::abline(h = 0, v = 0, col = "grey88")
    graphics::points(base::Re(ev), base::Im(ev), pch = 19, cex = 1.2,
                     col = base::ifelse(cx, "darkorange3", "firebrick"))
    sub <- base::sprintf("rho = %.4f   (%d roots)", base::max(Mod(ev)), base::length(ev))
    if (base::any(cx)) {
      k <- base::which(cx)[base::which.max(Mod(ev[cx]))]
      sub <- base::paste0(sub, base::sprintf("   dominant cycle: period %.1f",
                                             2 * base::pi / base::abs(base::Arg(ev[k]))))
    }
    graphics::title(sub = sub, cex.sub = 0.85)
    return(base::invisible(NULL))
  }

  if (type == "graph") {
    G <- if (base::is.null(lag)) x$G.sum else x$G[[lag]]
    ttl <- if (base::is.null(lag)) "Latent transition graph (sum of lags)"
           else base::sprintf("Latent transition graph, lag %d", lag)
    r <- 0.20
    pos <- if (Q == 1) base::cbind(0, 0) else {
      ang <- if (Q == 2) c(base::pi, 0)
             else base::seq(base::pi / 2, -3 * base::pi / 2, length.out = Q + 1)[base::seq_len(Q)]
      base::cbind(base::cos(ang), base::sin(ang))
    }
    W <- base::max(G, na.rm = TRUE); if (!base::is.finite(W) || W <= 0) W <- 1
    base::do.call("plot", dflt(x = NA, y = NA, xlim = c(-1.7, 1.7),
                               ylim = c(-1.35, 1.5), asp = 1, axes = FALSE,
                               xlab = "", ylab = "", main = ttl))
    ## Curve so that q -> q' and q' -> q stay apart: an asymmetric pair is the
    ## whole point of drawing the closed graph.
    bez <- function(p0, p1, curv = 0.28, n = 80) {
      m <- (p0 + p1) / 2; d <- p1 - p0
      c2 <- m + curv * c(-d[2], d[1])
      t <- base::seq(0, 1, length.out = n)
      base::cbind((1 - t)^2 * p0[1] + 2 * (1 - t) * t * c2[1] + t^2 * p1[1],
                  (1 - t)^2 * p0[2] + 2 * (1 - t) * t * c2[2] + t^2 * p1[2])
    }
    for (i in base::seq_len(Q)) for (j in base::seq_len(Q)) {
      w <- G[i, j]; if (!base::is.finite(w) || w <= 0) next
      lwd <- 0.8 + 7 * w / W
      if (i == j) {
        ctr <- if (Q == 1) c(0, 0.37) else pos[i, ] * (1 + r + 0.17)
        aa <- base::seq(0, 2 * base::pi, length.out = 90)
        graphics::lines(ctr[1] + 0.17 * base::cos(aa), ctr[2] + 0.17 * base::sin(aa),
                        lwd = lwd, col = "grey30")
        graphics::text(ctr[1], ctr[2] + 0.30, base::sprintf("%.2f", w),
                       cex = 0.9, col = "grey20")
      } else {
        b  <- bez(pos[j, ], pos[i, ])
        d0 <- base::sqrt(base::rowSums((b - base::matrix(pos[j, ], base::nrow(b), 2, byrow = TRUE))^2))
        d1 <- base::sqrt(base::rowSums((b - base::matrix(pos[i, ], base::nrow(b), 2, byrow = TRUE))^2))
        k  <- base::which(d0 > r & d1 > r)
        if (base::length(k) < 3) next
        graphics::lines(b[k, ], lwd = lwd, col = "steelblue4")
        n <- base::length(k)
        graphics::arrows(b[k[n - 1], 1], b[k[n - 1], 2], b[k[n], 1], b[k[n], 2],
                         length = 0.12, lwd = lwd, col = "steelblue4")
        mid <- b[k[base::round(n / 2)], ]
        graphics::text(mid[1], mid[2], base::sprintf("%.2f", w), cex = 0.9,
                       col = "steelblue4", pos = if (mid[2] > 0) 3 else 1)
      }
    }
    aa <- base::seq(0, 2 * base::pi, length.out = 100)
    for (i in base::seq_len(Q)) {
      graphics::polygon(pos[i, 1] + r * base::cos(aa), pos[i, 2] + r * base::sin(aa),
                        col = "white", border = "grey20")
      graphics::text(pos[i, 1], pos[i, 2], lab[i], cex = 1.0)
    }
    graphics::mtext("edge = G[effect, cause];  self-loop = persistence",
                    side = 1, cex = 0.85)
    return(base::invisible(NULL))
  }

  if (type == "irf") {
    H <- base::as.integer(horizon)
    Psi <- base::vector("list", H + 1L)
    Psi[[1L]] <- base::diag(Q)
    for (h in base::seq_len(H)) {
      M <- base::matrix(0, Q, Q)
      for (d in base::seq_len(base::min(h, D))) M <- M + x$G[[d]] %*% Psi[[h - d + 1L]]
      Psi[[h + 1L]] <- M
    }
    ymax <- base::max(base::vapply(Psi, base::max, base::numeric(1)), na.rm = TRUE)
    op <- graphics::par(mfrow = c(Q, Q), mar = c(4, 4, 2.4, 1))
    base::on.exit(graphics::par(op), add = TRUE)
    for (i in base::seq_len(Q)) for (j in base::seq_len(Q)) {
      y <- base::vapply(Psi, function(M) M[i, j], base::numeric(1))
      base::do.call("plot", dflt(x = 0:H, y = y, type = "h", lwd = 3,
                                 col = "steelblue4", ylim = c(0, ymax),
                                 xlab = "horizon h", ylab = "response",
                                 main = base::sprintf("shock in %s  ->  %s", lab[j], lab[i])))
      graphics::points(0:H, y, pch = 16, cex = 0.6, col = "steelblue4")
      graphics::abline(h = 0, col = "grey70")
    }
    return(base::invisible(NULL))
  }

  if (type == "lag") {
    prof <- function(i, j) base::vapply(x$G, function(M) M[i, j], base::numeric(1))
    ymax <- base::max(base::vapply(x$G, base::max, base::numeric(1)), na.rm = TRUE)
    one <- function(i, j, main) {
      v <- prof(i, j)
      bp <- graphics::barplot(v, names.arg = base::seq_len(D), col = "steelblue4",
                              border = NA, ylim = c(0, ymax * 1.15),
                              xlab = "lag d", ylab = expression(G[d]), main = main)
      k <- v >= 0.005                       # anything smaller just prints "0.00"
      if (base::any(k)) graphics::text(bp[k], v[k], base::sprintf("%.2f", v[k]),
                                       pos = 3, cex = 0.85)
      graphics::abline(h = 0, col = "grey60")
    }
    if (Q == 1L) {
      one(1L, 1L, "Lag profile of the latent transition")
    } else {
      op <- graphics::par(mfrow = c(Q, Q), mar = c(4, 4, 2.4, 1))
      base::on.exit(graphics::par(op), add = TRUE)
      for (i in base::seq_len(Q)) for (j in base::seq_len(Q))
        one(i, j, base::sprintf("%s  ->  %s", lab[j], lab[i]))
    }
    return(base::invisible(NULL))
  }

  ## type == "phase"
  if (D != 1L || Q != 2L)
    stop("type = \"phase\" needs a rank-2 model with a single lag (D = 1, Q = 2).")
  if (base::is.null(x$theta0))
    stop("type = \"phase\" needs an intercept; refit nmfkc.ar() with intercept = TRUE.")
  G <- x$G[[1L]]; th <- x$theta0; mu <- x$mu.b
  if (!base::all(base::is.finite(mu)))
    stop("type = \"phase\" needs a stationary fit (the fixed point is undefined).")
  start <- if (!base::is.null(extra$start)) base::as.matrix(extra$start) else {
    s <- base::max(mu) * 2
    base::rbind(c(0.05, 0.9), c(0.9, 0.05), c(0.95, 0.95), c(0.05, 0.05)) * s
  }
  extra$start <- NULL
  trace <- base::lapply(base::seq_len(base::nrow(start)), function(i) {
    o <- base::matrix(NA_real_, 41L, 2L); o[1L, ] <- start[i, ]
    for (t in base::seq_len(40L)) o[t + 1L, ] <- base::drop(G %*% o[t, ] + th)
    o
  })
  rng <- base::range(c(base::unlist(trace), mu, 0), na.rm = TRUE)
  base::do.call("plot", dflt(x = NA, y = NA, xlim = rng, ylim = rng,
                             xlab = base::sprintf("b[%s]", lab[1]),
                             ylab = base::sprintf("b[%s]", lab[2]),
                             main = "Phase portrait of the latent VAR"))
  ee <- base::eigen(G)
  for (k in 1:2) {
    v <- base::Re(ee$vectors[, k]); v <- v / base::sqrt(base::sum(v^2))
    d <- base::diff(rng)
    graphics::segments(mu[1] - d * v[1], mu[2] - d * v[2],
                       mu[1] + d * v[1], mu[2] + d * v[2], col = "grey75", lty = 2)
    graphics::text(mu[1] + 0.42 * d * v[1], mu[2] + 0.42 * d * v[2],
                   base::sprintf("lambda=%.2f", base::Re(ee$values[k])),
                   col = "grey45", cex = 0.8)
  }
  cols <- grDevices::hcl.colors(base::max(base::length(trace), 2L), "Dark 3")
  for (i in base::seq_along(trace)) {
    graphics::lines(trace[[i]], col = cols[i], lwd = 1.6)
    graphics::points(trace[[i]][1L, 1L], trace[[i]][1L, 2L], pch = 1,
                     col = cols[i], cex = 1.1)
    graphics::points(trace[[i]], pch = 16, cex = 0.35, col = cols[i])
  }
  graphics::points(mu[1], mu[2], pch = 8, cex = 2, lwd = 2, col = "firebrick")
  ## Real roots do not imply monotone decay: a negative real root alternates.
  ## And at Q = 2, D = 1 with G >= 0 the roots CANNOT be complex, so reporting
  ## "no cycle" as a finding there would overclaim.
  lam <- ee$values
  graphics::title(sub = if (base::any(base::abs(base::Im(lam)) > 1e-8))
                    "complex eigenvalues: the fit cycles"
                  else if (base::any(base::Re(lam) < -1e-8))
                    "a negative real root: decay alternates in sign (period 2)"
                  else
                    "real positive roots: non-oscillatory (forced here by Q = 2, D = 1, G >= 0)",
                  cex.sub = 0.85)
  base::invisible(NULL)
}


#' @title Bootstrap inference for the latent VAR of an NMF-VAR model
#' @description
#' This function is \strong{experimental}. The interface may change in future
#' versions: argument names, defaults and the contents of the returned object
#' are not yet stable.
#'
#' \code{nmfkc.ar.latent.inference} attaches a wild bootstrap to the quantities
#' reported by \code{\link{nmfkc.ar.latent}}: standard errors and intervals for
#' the spectral radius, the long-run mean, the long-run composition and the
#' entries of the latent transition matrices \eqn{G_d}.
#'
#' @details
#' \strong{What is identified, and why that decides what can be inferred.}  The
#' factorization admits \eqn{(X,\Theta)\to(XT,T^{-1}\Theta)} for invertible
#' \eqn{T}, so an entry of \eqn{\Theta} means nothing across replicates;
#' \code{\link{nmfkc.inference}} deals with this by conditioning on the
#' estimated \eqn{X}.  This function instead re-estimates the basis in every
#' replicate, propagating its uncertainty.  That splits the output in two:
#' \itemize{
#'   \item \strong{Inferential.}  \eqn{\rho}, the eigenvalues and the observed
#'     long-run mean \eqn{\bm\mu_y=X\bm\mu_b} are invariant under the whole
#'     family (\eqn{G_d\to T^{-1}G_dT} leaves the spectrum alone), so their
#'     replicates need no alignment and their intervals mean what they say.
#'   \item \strong{Conditional on identification.}  The entries
#'     \eqn{(G_d)_{qq'}} and the composition \eqn{\bm p^\ast} need \eqn{T} to be
#'     a permutation, and three conditions together deliver that: the column
#'     scale fixed (\code{X.restriction = "colSums"}), an \strong{anchor row} for
#'     every column of \eqn{X}, and a \strong{pure column} for every row of
#'     \eqn{\Theta}.  The anchor rows force \eqn{T\ge0}, the pure columns force
#'     \eqn{T^{-1}\ge0}, and a non-negative matrix with a non-negative inverse is
#'     monomial.  \strong{Either one-sided condition alone is not enough}: with
#'     \eqn{X=I_2},
#'     \eqn{\Theta=\left(\begin{smallmatrix}2&3\\1&1\end{smallmatrix}\right)} and
#'     \eqn{T=\left(\begin{smallmatrix}0.9&0.1\\0.1&0.9\end{smallmatrix}\right)}
#'     the basis is perfectly separable and column-normalized, yet \eqn{XT\ge0},
#'     \eqn{T^{-1}\Theta\ge0}, \eqn{X\Theta} is unchanged, \eqn{T} is no
#'     permutation, and \eqn{G} differs between the two representatives ---
#'     \eqn{\Theta} has no pure column there.  \code{identified.approx} reports
#'     the verdict and \code{identification} the two sides plus how far each is
#'     from exact; the function warns either way, naming the side that fails or
#'     the margins when it only passes approximately.  \strong{The thresholds
#'     make this a diagnostic, not a certificate}: the argument needs exact
#'     anchors and exact zeros, which multiplicative updates do not deliver.
#'     Replicates are permutation-aligned to the original basis regardless.
#'     \code{X.restriction = "fixed"} is the one clean case for \eqn{Q>1} --- the
#'     basis is not estimated, so \eqn{T=I}.
#' }
#'
#' \strong{\code{prob.nonstationary} is not a unit-root test.}  It is the
#' bootstrap tail probability \eqn{P^\ast(\rho^\ast\ge1)}, with replicates drawn
#' around the \emph{estimated} model.  Nothing imposes \eqn{H_0:\rho=1}, so it
#' does not have the interpretation of a p-value for that hypothesis; a proper
#' test would have to resample under the null.
#'
#' \strong{The zeros of \eqn{\Theta} are only meaningful under
#' misspecification.}  The strict-complementarity diagnostic reports the dual
#' margin \eqn{\delta_{\mathrm{dual}}=\min|\Lambda^\ast|} over the zero
#' coefficients and the primal margin
#' \eqn{\delta_{\mathrm{prim}}=\min\Theta} over the positive ones.  Both must be
#' bounded away from zero for the support to be recovered exactly and for the
#' active coefficients to attain the distribution that knowing the support would
#' give.  If the model is correctly specified, \eqn{\Xi_0=X\Theta} makes
#' \eqn{\Lambda^\ast=0} exactly, so strict complementarity fails at
#' \emph{every} zero: \eqn{P(\hat\theta_k=0)} then tends to a constant in
#' \eqn{(0,1)} (one half for an isolated zero) instead of to one, and no
#' resampling scheme is consistent there.  A small \code{kkt$ratio} is therefore
#' not a numerical defect --- it says the zeros are not estimating a sign
#' restriction.  What makes them estimate something is misspecification: they
#' recover the coordinates whose \emph{unconstrained} population coefficient is
#' strictly negative.
#'
#' \strong{Bootstrap validity at the boundary.}  A naive nonparametric bootstrap
#' is inconsistent for a parameter on the boundary.  The scheme used here
#' resamples the residuals and re-imposes the non-negativity constraint in every
#' re-fit, which is the practical remedy; a moving-block scheme would in
#' addition preserve the serial dependence that the fixed design discards.
#'
#' \strong{Bootstrap p-values invert the centred distribution.}  Every replicate
#' of a non-negative coefficient is \eqn{\ge 0}, so comparing the raw replicates
#' with zero would return \eqn{p = 0} for every entry.  The reported
#' \code{G.p.value} is the basic-bootstrap two-sided \eqn{p} obtained from
#' \eqn{P^\ast(\hat G^\ast \ge 2\hat G)}.
#'
#' \strong{Resampling scheme.}  With \eqn{\hat Y = X\Theta A} and residuals
#' \eqn{E=Y-\hat Y}, each replicate forms \eqn{Y^\ast=\hat Y+E\odot W} and
#' re-fits \code{\link{nmfkc}} at the same rank.  \eqn{W} carries one multiplier
#' per time point (\code{wild.unit = "column"}, the default, which preserves the
#' contemporaneous covariance between series) or one per cell.  Negative entries
#' of \eqn{Y^\ast} are truncated at 0 to keep the response non-negative; the
#' truncation rate is reported, and a large one is a warning sign for the whole
#' procedure.
#'
#' \strong{This is a fixed-design bootstrap}: the covariate matrix \eqn{A} is
#' held at the observed lags rather than rebuilt from \eqn{Y^\ast}, so the
#' intervals are conditional on the observed past and do not propagate the
#' randomness of the lagged design itself.  They may therefore be somewhat
#' optimistic.
#'
#' @param object An object from \code{\link{nmfkc.ar.latent}}.
#' @param Y,A The observation and covariate matrices the model was fitted on,
#'   i.e. the \code{Y} and \code{A} returned by \code{\link{nmfkc.ar}}.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{\code{wild.B}}{Number of bootstrap replicates.  Default 500.}
#'     \item{\code{wild.dist}}{Multiplier distribution, \code{"rademacher"}
#'       (default) or \code{"exp"} (mean-centred).}
#'     \item{\code{wild.unit}}{\code{"column"} (default, one multiplier per time
#'       point) or \code{"element"} (one per cell).}
#'     \item{\code{wild.level}}{Confidence level.  Default 0.95.}
#'     \item{\code{wild.seed}}{Seed for the bootstrap.  Default 123.}
#'     \item{\code{fit}}{The fitted \code{"nmfkc"} object.  Only needed when
#'       \code{object} does not carry \code{X} / \code{Theta}.}
#'     \item{\code{truncate}}{Logical; clip \eqn{Y^\ast} at 0 (default
#'       \code{TRUE}).}
#'     \item{\code{epsilon}, \code{maxit}}{Convergence control of the re-fits,
#'       default \code{1e-8} and \code{20000} --- tight, because a coefficient
#'       heading for the non-negativity boundary moves slowly and a loose
#'       tolerance biases anything that is compared with zero.}
#'     \item{\code{kkt.tol}, \code{kkt.ratio.min}}{Thresholds of the strict
#'       complementarity diagnostic (defaults \code{1e-3} and \code{5}).}
#'     \item{\code{sep.tol}}{Row-normalised loading a column of \eqn{X} must
#'       reach for that basis to count as anchored.  Default \code{0.999}.  The
#'       zero threshold for the pure-column test on \eqn{\Theta} is
#'       \code{kkt.tol}.}
#'     \item{\code{refit.args}}{Named list of extra arguments for the bootstrap
#'       re-fits.  \code{method} and \code{X.restriction} are taken from the
#'       original fit automatically, so the replicates use the same estimator;
#'       anything else the original call set (\code{X.init}, penalties, ...)
#'       cannot be recovered from the fitted object and must be passed here.
#'       The function warns when it sees such arguments in the recorded call.}
#'     \item{\code{cores}}{Number of workers for the re-fits, default
#'       \code{getOption("mc.cores", 1L)} (sequential).  The bootstrap
#'       multipliers are drawn before the loop and each re-fit is deterministic
#'       given its data, so the result does not depend on \code{cores}.}
#'   }
#' @return An object of class \code{"nmfkc.ar.latent.inference"}: a list with
#' \item{coefficients}{Data frame of the \eqn{DQ^2} entries of the \eqn{G_d},
#'   with columns \code{Lag}, \code{Basis} (effect at \eqn{t}),
#'   \code{Covariate} (cause at \eqn{t-d}), \code{Estimate}, \code{BSE},
#'   \code{CI_low}, \code{CI_high}, \code{p_value}.  Descriptive unless
#'   \code{identified.approx}.}
#' \item{G, G.se, G.ci.lower, G.ci.upper, G.p.value}{The same quantities as
#'   lists of \eqn{Q\times Q} matrices, one per lag.}
#' \item{spectral.radius}{The radius of the \eqn{DQ\times DQ} latent companion
#'   matrix --- the same quantity \code{\link{nmfkc.ar.latent}} reports, not
#'   \eqn{\rho(\sum_d G_d)}, which differs when \eqn{D>1}.  With
#'   \code{spectral.radius.se}, \code{spectral.radius.ci.lower} / \code{.upper}
#'   and \code{spectral.radius.boot.mean}.}
#' \item{prob.nonstationary}{Bootstrap tail probability
#'   \eqn{P^\ast(\rho^\ast\ge1)}.  \strong{Not} a p-value for
#'   \eqn{H_0:\rho\ge1} --- see Details.}
#' \item{complex, complex.frac}{Whether the point estimate has complex roots,
#'   and the fraction of replicates that do.}
#' \item{mu.y}{With \code{.ci.lower} / \code{.ci.upper}; invariant, so
#'   inferential.}
#' \item{p.star}{With \code{.ci.lower} / \code{.ci.upper}; descriptive unless
#'   \code{identified.approx}.}
#' \item{kkt}{Strict-complementarity diagnostic: \code{n.boundary},
#'   \code{all.positive}, \code{ratio}, and the two margins \code{delta.dual}
#'   and \code{delta.prim}.}
#' \item{identified.approx}{Logical, an \strong{approximate diagnostic}:
#'   \code{TRUE} for \eqn{Q=1}, for \code{X.restriction = "fixed"} (where
#'   \eqn{X} is not estimated, so \eqn{T=I}), or when the per-column scale is
#'   fixed (\code{"colSums"} / \code{"colSqSums"}) and the anchor-row and
#'   pure-column tests both pass at \code{sep.tol} / \code{kkt.tol}.  It is not a
#'   certificate --- the uniqueness argument needs \emph{exact} anchors and
#'   \emph{exact} zeros, and multiplicative updates leak a little.  When
#'   \code{FALSE} the entry-level summaries and \code{p.star} are descriptive;
#'   \eqn{\rho}, the eigenvalues and \eqn{\bm\mu_y} remain inferential either
#'   way.}
#' \item{identification}{The two sides separately: \code{anchor} (per column of
#'   \eqn{X}), \code{pure} (per row of \eqn{\Theta}), \code{ok}, and the margins
#'   \code{anchor.margin} / \code{pure.margin} saying how far the fit is from an
#'   exact configuration (0 = exact).}
#' \item{refit.args}{The arguments the bootstrap re-fits actually used, so the
#'   replicates can be checked against the original estimator.}
#'   Plus the bookkeeping entries \code{wild.B}, \code{wild.B.requested},
#'   \code{n.fail}, \code{truncate.rate}, \code{wild.unit}, \code{wild.dist},
#'   \code{wild.level}, \code{dims}.
#' @seealso \code{\link{nmfkc.ar.latent}}, \code{\link{nmfkc.ar.stationarity}},
#'   \code{\link{nmfkc.ar}}, \code{\link{nmfkc.inference}}
#' @references
#' Satoh, K. (2025). Applying non-negative matrix factorization with covariates
#'   to multivariate time series data as a vector autoregression model.
#'   \emph{Japanese Journal of Statistics and Data Science}. arXiv:2501.17446.
#'   \doi{10.1007/s42081-025-00314-0}
#' @examples
#' \donttest{
#' set.seed(1)
#' Y   <- matrix(abs(rnorm(4 * 40)) + 1, 4, 40)
#' ar  <- nmfkc.ar(Y, degree = 1)
#' fit <- nmfkc(ar$Y, ar$A, rank = 2, verbose = FALSE)
#' lat <- nmfkc.ar.latent(fit)
#' nmfkc.ar.latent.inference(lat, ar$Y, ar$A, wild.B = 50)
#' }
#' @export
nmfkc.ar.latent.inference <- function(object, Y, A, ...) {
  if (!inherits(object, "nmfkc.ar.latent"))
    stop("'object' must come from nmfkc.ar.latent().")
  ex <- base::list(...)
  wild.B     <- if (!is.null(ex$wild.B))     ex$wild.B     else 500L
  wild.dist  <- if (!is.null(ex$wild.dist))
                  match.arg(ex$wild.dist, c("rademacher", "exp")) else "rademacher"
  wild.unit  <- if (!is.null(ex$wild.unit))
                  match.arg(ex$wild.unit, c("column", "element")) else "column"
  wild.level <- if (!is.null(ex$wild.level)) ex$wild.level else 0.95
  wild.seed  <- if (!is.null(ex$wild.seed))  ex$wild.seed  else 123
  truncate   <- if (!is.null(ex$truncate))   ex$truncate   else TRUE
  epsilon    <- if (!is.null(ex$epsilon))    ex$epsilon    else 1e-8
  maxit      <- if (!is.null(ex$maxit))      ex$maxit      else 20000L
  kkt.tol    <- if (!is.null(ex$kkt.tol))    ex$kkt.tol    else 1e-3
  kkt.ratio.min <- if (!is.null(ex$kkt.ratio.min)) ex$kkt.ratio.min else 5
  sep.tol    <- if (!is.null(ex$sep.tol))    ex$sep.tol    else 0.999
  ## The re-fits must use the SAME estimator as the original fit, otherwise the
  ## bootstrap describes a different one.  What nmfkc() records (method,
  ## X.restriction) is carried on the object and reused; anything else the
  ## original call set must be supplied here.
  refit.args <- object$fit.args
  refit.args <- refit.args[!vapply(refit.args, is.null, logical(1))]
  if (!is.null(ex$refit.args)) refit.args[names(ex$refit.args)] <- ex$refit.args
  ## Warn about non-default settings in the original call that we cannot see.
  if (!is.null(object$call)) {
    ## Arguments that cannot change the estimator (labels and reporting) plus the
    ## ones we set or inherit ourselves.
    seen <- c("Y", "A", "rank", "data", "epsilon", "maxit", "verbose",
              "prefix", "print.trace", "detail", "", names(refit.args))
    cl <- as.list(object$call)
    ## A literal seed can be inherited outright; anything else has to be passed.
    if (!is.null(cl$seed) && is.numeric(cl$seed) && is.null(refit.args$seed)) {
      refit.args$seed <- cl$seed
      seen <- c(seen, "seed")
    }
    unseen <- setdiff(names(cl)[-1], seen)
    if (length(unseen))
      warning("The original fit was called with ", paste(unseen, collapse = ", "),
              ", which the bootstrap re-fits cannot reproduce automatically. ",
              "Pass them through refit.args= so the replicates use the same ",
              "estimator.")
  }
  ## Keep our own seeding out of the caller's random stream.
  .rng <- .nmfkc.rng.save(wild.seed)
  on.exit(.nmfkc.rng.restore(.rng), add = TRUE)
  if (!is.null(wild.seed)) set.seed(wild.seed)

  X <- object$X; Theta <- object$Theta
  if (is.null(X) || is.null(Theta)) {
    if (is.null(ex$fit))
      stop("'object' does not carry X/Theta; pass the fitted nmfkc object via fit=.")
    X <- ex$fit$X; Theta <- ex$fit$C
  }
  Y <- as.matrix(Y); A <- as.matrix(A)
  P <- object$dims[["P"]]; Q <- object$dims[["Q"]]; D <- object$dims[["D"]]
  n <- ncol(Y)
  intercept <- (ncol(Theta) == P * D + 1L)

  Gsum <- function(Xm, Cm) {
    L <- Cm[, seq_len(P), drop = FALSE]
    if (D > 1) for (d in 2:D) L <- L + Cm[, ((d - 1) * P + 1):(d * P), drop = FALSE]
    L %*% Xm
  }
  ## The roots must be the SAME quantity nmfkc.ar.latent() reports, i.e. the
  ## eigenvalues of the DQ x DQ companion matrix of the G_d -- not those of
  ## sum_d G_d.  For D > 1 the two differ: with G_1 = 0.2, G_2 = 0.1 the
  ## companion radius is 0.4317 (roots 0.4317, -0.2317) while
  ## rho(G_1 + G_2) = 0.3.  They agree on the stationary/non-stationary verdict
  ## under non-negativity, but the number, its SE and the implied period do not.
  Gd.list <- function(Xm, Cm) base::lapply(base::seq_len(D), function(d)
    Cm[, ((d - 1) * P + 1):(d * P), drop = FALSE] %*% Xm)
  ev.of <- function(Xm, Cm) .nmfvar.roots(Gd.list(Xm, Cm), Q, D)
  ## G_d as one stacked Q x (QD) matrix, which is the convenient shape for
  ## permutation-aligning a replicate: G_d -> P' G_d P acts on both margins.
  Gd.of <- function(Xm, Cm) do.call(cbind, Gd.list(Xm, Cm))

  ## Can an ENTRY of G_d be given a standard error?  Only if the factorization
  ## is unique up to a permutation, which takes three things: the column scale
  ## fixed, an anchor row for every column of X (forcing T >= 0) and a pure
  ## column for every row of Theta (forcing T^-1 >= 0).  See
  ## .nmfvar.identified() for why either one-sided condition alone is not
  ## enough -- separability of X on its own admits T = [.9 .1; .1 .9].
  xr <- object$X.restriction
  ## Which constraints actually fix the PER-COLUMN scale of X?  Only the two
  ## that normalise each column: "colSums" and "colSqSums".  "totalSum"
  ## constrains the grand total, so columns can still trade scale between them,
  ## and "none" constrains nothing.  "fixed" is the opposite extreme -- X is not
  ## estimated at all, so T = I and there is no rotational freedom to certify.
  scale.fixed <- isTRUE(xr %in% c("colSums", "colSqSums"))
  X.held      <- isTRUE(xr == "fixed")
  ident <- .nmfvar.identified(X, Theta, sep.tol = sep.tol, zero.tol = kkt.tol)
  ## APPROXIMATE, not a theorem.  The uniqueness argument needs an EXACT anchor
  ## row (loading exactly 1 after row normalisation) and an EXACT pure column
  ## (the other entries exactly 0).  A fit that merely comes within sep.tol and
  ## kkt.tol is close to such a configuration but does not satisfy its
  ## hypotheses, so this is evidence about the fit, not a guarantee.
  identified.approx <- (Q == 1L) || X.held || (scale.fixed && ident$ok)
  if (Q > 1L && !X.held && !identified.approx) {
    why <- c(if (!scale.fixed) paste0("X.restriction = \"", xr,
               "\" does not fix the per-column scale of X"),
             if (!all(ident$anchor)) paste0("basis ",
               paste(which(!ident$anchor), collapse = ", "),
               " has no anchor row in X (so T >= 0 is not forced)"),
             if (!all(ident$pure)) paste0("basis ",
               paste(which(!ident$pure), collapse = ", "),
               " has no pure column in Theta (so T^-1 >= 0 is not forced)"))
    warning("The factorization is not even approximately unique up to a ",
            "permutation: ", paste(why, collapse = "; "), ". The entries of ",
            "G_d and p.star are therefore DESCRIPTIVE; rho, the eigenvalues ",
            "and mu.y are invariant under the whole family and remain the ",
            "inferential targets.")
  } else if (Q > 1L && !X.held) {
    warning("The anchor-row / pure-column diagnostic passes only APPROXIMATELY ",
            "(anchor margin ", signif(ident$anchor.margin, 3),
            ", pure margin ", signif(ident$pure.margin, 3), "). The uniqueness ",
            "argument needs exact anchors and exact zeros, so treat the ",
            "entry-level summaries as well-behaved but not certified.")
  }

  ## Strict-complementarity (oracle) diagnostic.  In the population,
  ##   Lambda* = X'(Xi_0 - X Theta) Sigma_A,
  ## and the oracle property (exact support recovery, plus the asymptotic
  ## distribution of the active coefficients being the one that knowing the
  ## support would give) needs BOTH margins away from zero:
  ##   delta.dual = min |Lambda*| over the zero coefficients,
  ##   delta.prim = min Theta over the positive ones.
  ## Note what a vanishing dual margin means.  If the model is correctly
  ## specified, Xi_0 = X Theta and Lambda* = 0 exactly, so strict
  ## complementarity fails at EVERY zero -- P(theta_k = 0) then tends to a
  ## constant in (0, 1) (1/2 for an isolated zero) rather than to 1.  A small
  ## ratio is therefore not a numerical defect; it says the zeros are not
  ## estimating a sign restriction, and no resampling scheme is consistent for
  ## them.  What makes the zeros meaningful is misspecification: they estimate
  ## the set of coordinates whose UNCONSTRAINED population coefficient is
  ## strictly negative.
  Lam <- t(X) %*% (X %*% Theta %*% A - Y) %*% t(A) / n   # = -Lambda*
  S   <- Theta > kkt.tol
  kkt <- list(n.boundary = sum(!S), all.positive = NA, ratio = NA_real_,
              delta.dual = NA_real_, delta.prim = NA_real_)
  if (any(!S) && any(S)) {
    kkt$all.positive <- all(Lam[!S] > 0)
    kkt$delta.dual   <- min(Lam[!S])
    kkt$delta.prim   <- min(Theta[S])
    kkt$ratio        <- kkt$delta.dual / max(abs(Lam[S]))
    if (!isTRUE(kkt$all.positive) || kkt$ratio < kkt.ratio.min)
      warning("Strict complementarity does not hold (KKT ratio = ",
              signif(kkt$ratio, 3), "). Either the model is correctly ",
              "specified there -- in which case the population multipliers are ",
              "exactly zero and the boundary coefficients are not consistently ",
              "estimated as zeros -- or the margin is too small to tell. ",
              "No resampling scheme is consistent at such coordinates.")
  }

  ## Point estimates
  ev0   <- ev.of(X, Theta)
  rho0  <- max(Mod(ev0))
  cplx0 <- any(abs(Im(ev0)) > 1e-8)
  th0   <- if (intercept) Theta[, P * D + 1L] else rep(0, Q)
  mu.b0 <- if (rho0 < 1) solve(diag(Q) - Gsum(X, Theta), th0) else rep(NA_real_, Q)
  p0    <- if (all(is.finite(mu.b0)) && sum(mu.b0) > 0) mu.b0 / sum(mu.b0)
           else rep(NA_real_, Q)

  Gd0 <- Gd.of(X, Theta)

  ## G_d and p* are invariant only up to the ordering of the bases, so each
  ## replicate is matched to the original X first (rho and mu.y need no such
  ## alignment).
  perms <- if (Q <= 6)
    do.call(rbind, lapply(seq_len(factorial(Q)), function(i) .nmfvar.perm(Q, i)))
  else NULL
  align <- function(Xs) {
    if (!is.null(perms)) {
      d <- apply(perms, 1, function(p) sum((Xs[, p, drop = FALSE] - X)^2))
      perms[which.min(d), ]
    } else {
      p <- integer(Q); avail <- seq_len(Q)
      for (q in seq_len(Q)) {
        d <- vapply(avail, function(j) sum((Xs[, j] - X[, q])^2), numeric(1))
        p[q] <- avail[which.min(d)]; avail <- setdiff(avail, p[q])
      }
      p
    }
  }

  Yhat <- X %*% Theta %*% A
  E    <- Y - Yhat

  ## The multipliers are drawn up front, in the order the sequential loop drew
  ## them, so that the replicate datasets do not depend on `cores`.  What is
  ## left per replicate is a re-fit that is deterministic given its data (nmfkc
  ## seeds itself), which is what makes the parallel run result-identical.
  m  <- if (wild.unit == "column") n else P * n
  Ys.all <- base::lapply(base::seq_len(wild.B), function(b) {
    w <- base::switch(wild.dist,
                      rademacher = base::sample(c(-1, 1), m, TRUE),
                      exp        = stats::rexp(m) - 1)
    Ws <- if (wild.unit == "column") base::matrix(w, P, n, byrow = TRUE)
          else base::matrix(w, P, n)
    Yhat + E * Ws
  })
  trunc <- base::sum(base::vapply(Ys.all, function(Ys) base::mean(Ys < 0),
                                  base::numeric(1)))

  one_boot <- function(Ys) {
    if (truncate) Ys[Ys < 0] <- 0
    ## detail = "fast": only X and C are read from the re-fit, so the O(N^2)
    ## sample-clustering criteria would be computed and discarded wild.B times.
    rr <- base::try(base::do.call(nmfkc, c(base::list(
                      Y = Ys, A = A, rank = Q, epsilon = epsilon,
                      maxit = maxit, detail = "fast", verbose = FALSE),
                      refit.args)),
                    silent = TRUE)
    if (base::inherits(rr, "try-error"))
      return(base::list(rho = 0, cplx = FALSE, muy = NULL, p = NULL, G = NULL,
                        fail = TRUE))

    Xs <- rr$X; Cs <- rr$C
    ev <- ev.of(Xs, Cs)
    pm <- align(Xs)
    ## G_d -> P' G_d P: the permutation acts on the effect margin (rows) and on
    ## the cause margin (columns) of every lag block alike.
    Gs <- Gd.of(Xs, Cs)[pm, base::rep(pm, D) + base::rep((base::seq_len(D) - 1L) * Q,
                                                        each = Q), drop = FALSE]
    out <- base::list(rho = base::max(Mod(ev)),
                      cplx = base::any(base::abs(Im(ev)) > 1e-8),
                      muy = NULL, p = NULL, G = Gs, fail = FALSE)
    if (out$rho < 1) {
      ths <- if (intercept) Cs[, P * D + 1L] else base::rep(0, Q)
      mb  <- base::try(base::solve(base::diag(Q) - Gsum(Xs, Cs), ths), silent = TRUE)
      if (!base::inherits(mb, "try-error")) {
        out$muy <- base::as.numeric(Xs %*% mb)
        if (base::sum(mb) > 0) out$p <- (mb / base::sum(mb))[pm]
      }
    }
    out
  }

  cores <- if (!base::is.null(ex$cores)) ex$cores else base::getOption("mc.cores", 1L)
  res <- .nmfkc.parlapply(Ys.all, one_boot, cores = cores,
                          envir = base::environment())

  rho.b <- base::vapply(res, function(r) r$rho, base::numeric(1))
  cplx.b <- base::vapply(res, function(r) r$cplx, base::logical(1))
  fail.b <- base::vapply(res, function(r) r$fail, base::logical(1))
  nfail  <- base::sum(fail.b)
  muy.b <- base::matrix(NA_real_, wild.B, P); p.b <- base::matrix(NA_real_, wild.B, Q)
  G.b   <- base::matrix(NA_real_, wild.B, Q * Q * D)
  for (b in base::seq_len(wild.B)) {
    if (!base::is.null(res[[b]]$muy)) muy.b[b, ] <- res[[b]]$muy
    if (!base::is.null(res[[b]]$p))   p.b[b, ]   <- res[[b]]$p
    if (!base::is.null(res[[b]]$G))   G.b[b, ]   <- base::as.numeric(res[[b]]$G)
  }

  ## Use the failure flag, not rho > 0: a replicate whose Theta collapses to
  ## zero has rho = 0 legitimately (and is stationary), so the old test threw
  ## away a valid draw and inflated the failure count.
  ok <- !fail.b
  a2 <- (1 - wild.level) / 2
  qci <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) c(NA_real_, NA_real_)
    else stats::quantile(v, c(a2, 1 - a2), names = FALSE)
  }
  rho.ci <- qci(rho.b[ok])
  muy.ci <- apply(muy.b, 2, qci)
  pst.ci <- apply(p.b,   2, qci)

  ## Entries of G_d.  The p-value inverts the CENTRED replicate distribution:
  ## every replicate of a non-negative coefficient is >= 0, so comparing the raw
  ## replicates with zero would return p = 0 everywhere.
  g0     <- base::as.numeric(Gd0)
  G.se   <- base::apply(G.b, 2, function(v) stats::sd(v[base::is.finite(v)]))
  G.ci   <- base::apply(G.b, 2, qci)
  G.p    <- base::vapply(base::seq_along(g0), function(i) {
    v <- G.b[, i]; v <- v[base::is.finite(v)]
    if (!base::length(v)) return(NA_real_)
    base::min(1, 2 * base::min(base::mean(v >= 2 * g0[i]),
                               1 - base::mean(v >= 2 * g0[i])))
  }, base::numeric(1))

  lab  <- base::colnames(Gd0)[base::seq_len(Q)]
  if (base::is.null(lab)) lab <- base::paste0("Basis", base::seq_len(Q))
  ## as.numeric() unrolled column-major over a Q x (QD) matrix, so the index
  ## runs effect-fastest within (cause, lag).
  idx  <- base::expand.grid(Basis = lab, Covariate = lab, Lag = base::seq_len(D),
                            stringsAsFactors = FALSE)
  spl  <- function(v) {
    out <- base::lapply(base::seq_len(D), function(d) {
      m <- base::matrix(v[((d - 1) * Q * Q + 1):(d * Q * Q)], Q, Q,
                        dimnames = base::list(lab, lab))
      m
    })
    base::names(out) <- base::paste0("lag", base::seq_len(D))
    out
  }
  coefficients <- base::data.frame(
    Lag = idx$Lag, Basis = idx$Basis, Covariate = idx$Covariate,
    Estimate = g0, BSE = G.se,
    CI_low = G.ci[1L, ], CI_high = G.ci[2L, ], p_value = G.p,
    stringsAsFactors = FALSE)

  structure(list(
    coefficients = coefficients,
    G = spl(g0), G.se = spl(G.se),
    G.ci.lower = spl(G.ci[1L, ]), G.ci.upper = spl(G.ci[2L, ]),
    G.p.value = spl(G.p),
    identified.approx = identified.approx, identification = ident,
    refit.args = refit.args,
    spectral.radius           = rho0,
    spectral.radius.se        = stats::sd(rho.b[ok]),
    spectral.radius.ci.lower  = rho.ci[1L],
    spectral.radius.ci.upper  = rho.ci[2L],
    spectral.radius.boot.mean = mean(rho.b[ok]),
    ## NOT a p-value: the replicates are drawn around the ESTIMATED model, not
    ## under an imposed H0: rho = 1, so this is the bootstrap tail probability
    ## P*(rho* >= 1) and nothing stronger.  A real unit-root test would have to
    ## resample under the null.
    prob.nonstationary = mean(rho.b[ok] >= 1),
    complex = cplx0, complex.frac = mean(cplx.b[ok]),
    mu.y = as.numeric(X %*% mu.b0),
    mu.y.ci.lower = muy.ci[1L, ], mu.y.ci.upper = muy.ci[2L, ],
    p.star = p0,
    p.star.ci.lower = pst.ci[1L, ], p.star.ci.upper = pst.ci[2L, ],
    kkt = kkt,
    wild.B = sum(ok), wild.B.requested = wild.B, n.fail = nfail,
    truncate.rate = trunc / wild.B,
    wild.unit = wild.unit, wild.dist = wild.dist, wild.level = wild.level,
    dims = object$dims),
    class = "nmfkc.ar.latent.inference")
}


#' @title Print method for nmfkc.ar.latent.inference objects
#' @param x An object of class \code{"nmfkc.ar.latent.inference"}.
#' @param digits Number of digits used when printing.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.nmfkc.ar.latent.inference <- function(x, digits = 4, ...) {
  cat("Bootstrap inference for the latent VAR\n")
  cat(sprintf("  B = %d (%d failed), unit = \"%s\", dist = \"%s\", truncation rate = %.2f%%\n",
              x$wild.B, x$n.fail, x$wild.unit, x$wild.dist, 100 * x$truncate.rate))
  cat("\nLatent transition matrices  G_d = Theta_d %*% X",
      "  (row = effect at t, column = cause at t-d)\n")
  if (!isTRUE(x$identified.approx))
    cat("  NOTE: the anchor-row / pure-column diagnostic FAILS, so these\n",
        "       entries (and p.star) are DESCRIPTIVE. rho, the eigenvalues\n",
        "       and mu.y are invariant and inferential.\n", sep = "")
  else if (x$dims[["Q"]] > 1L)
    cat("  NOTE: the diagnostic passes only APPROXIMATELY (anchor margin ",
        signif(x$identification$anchor.margin, 3), ", pure margin ",
        signif(x$identification$pure.margin, 3), ");\n",
        "       the uniqueness argument needs exact anchors and exact zeros.\n",
        sep = "")
  for (nm in names(x$G)) {
    cat("\n", nm, ":\n", sep = "")
    print(round(x$G[[nm]], 3))
    cat("SE:\n"); print(round(x$G.se[[nm]], 3))
    cat("p-value:\n"); print(round(x$G.p.value[[nm]], 3))
  }
  cat(sprintf("\nrho(G) = %.*f   SE = %.*f   %g%% CI = [%.*f, %.*f]\n",
              digits, x$spectral.radius, digits, x$spectral.radius.se,
              100 * x$wild.level, digits, x$spectral.radius.ci.lower,
              digits, x$spectral.radius.ci.upper))
  cat(sprintf("  bootstrap mean = %.*f  (bias %+.*f)\n",
              digits, x$spectral.radius.boot.mean, digits,
              x$spectral.radius.boot.mean - x$spectral.radius))
  cat(sprintf("  bootstrap tail P*(rho* >= 1) = %.3f\n", x$prob.nonstationary))
  cat("  (a tail probability around the fitted model, NOT a unit-root test:",
      "the\n   replicates are not drawn under an imposed H0: rho = 1)\n")
  cat(sprintf("\ncomplex eigenvalues: point estimate %s, in %.1f%% of replicates\n",
              if (x$complex) "yes" else "no", 100 * x$complex.frac))
  if (any(is.finite(x$p.star.ci.lower))) {
    cat("\nlong-run composition p*:\n")
    print(round(cbind(estimate = x$p.star,
                      lower = x$p.star.ci.lower, upper = x$p.star.ci.upper), digits))
  }
  if (is.finite(x$kkt$ratio))
    cat(sprintf("\nstrict complementarity: %d boundary coords, all multipliers positive = %s, ratio = %.1f\n",
                x$kkt$n.boundary, x$kkt$all.positive, x$kkt$ratio))
  invisible(x)
}
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------




############################################################
## nmfkc.ar.DOT  (Graphviz visualization for NMF-AR / VAR models)
############################################################

#' Generate a Graphviz DOT Diagram for NMF-AR / NMF-VAR Models
#'
#' @description
#' Produces a Graphviz DOT script for visualizing autoregressive
#' NMF-with-covariates models constructed via \code{nmfkc.ar} + \code{nmfkc}.
#'
#' The diagram displays three types of directed relationships:
#' \itemize{
#'   \item Lagged predictors: \eqn{T_{t-k} \rightarrow X},
#'   \item Current latent factors: \eqn{X \rightarrow T_t},
#'   \item Optional intercept effects: \code{Const -> X}.
#' }
#'
#' Importantly, *no direct edges from lagged variables to current outputs*
#' (\eqn{T_{t-k} \rightarrow T_t}) are drawn, in accordance with the NMF-AR
#' formulation.
#'
#' Each block of lagged variables is displayed in its own DOT subgraph
#' (e.g., “T-1”, “T-2”, ...), while latent factor nodes and current-time
#' outputs are arranged in separate clusters.
#'
#' @param result A fitted \code{nmfkc} object representing the AR model.
#'   Must contain matrices \code{X} and \code{C}.
#' @param degree Maximum AR lag to visualize.
#' @param intercept Logical; if \code{TRUE}, draws intercept nodes for
#'   columns named "(Intercept)" in matrix \code{C}.
#'   Default is \code{TRUE} when an intercept column is detected in \code{C},
#'   \code{FALSE} otherwise (auto-detected).
#' @param threshold Minimum coefficient magnitude required to draw an edge.
#' @param rankdir Graphviz rank direction (e.g., \code{"RL"}, \code{"LR"}, \code{"TB"}).
#' @param fill Logical; whether nodes are filled with color.
#' @param weight_scale_xy Scaling factor for edges \eqn{X \rightarrow T}.
#' @param weight_scale_lag Scaling factor for lagged edges \eqn{T-k \rightarrow X}.
#' @param weight_scale_int Scaling factor for intercept edges.
#' @param hide.isolated Logical. If \code{TRUE} (default), Y nodes that have no
#'   edges at or above \code{threshold} are excluded from the graph.
#'
#' @return A character string representing a Graphviz DOT file.
#' @seealso \code{\link{nmfkc.ar}}, \code{\link{nmfkc}}, \code{\link{plot.nmfkc.DOT}}
#' @examples
#' d <- AirPassengers
#' ar_data <- nmfkc.ar(d, degree = 2)
#' result <- nmfkc(ar_data$Y, ar_data$A, rank = 1)
#' dot <- nmfkc.ar.DOT(result, degree = 2)
#' cat(dot)
#'
#' @export
nmfkc.ar.DOT <- function(result,
                         degree    = 1,
                         intercept = any(colnames(result$C) == "(Intercept)"),
                         threshold = 0.1,
                         rankdir   = "RL",
                         fill      = TRUE,
                         weight_scale_xy  = 5,
                         weight_scale_lag = 5,
                         weight_scale_int = 3,
                         hide.isolated   = TRUE) {

  ## -------------------------------------------------------------
  ## Extract required AR components
  ## -------------------------------------------------------------
  X     <- result$X
  C_raw <- result$C

  if (is.null(X) || is.null(C_raw)) {
    stop("result must contain matrices X and C.")
  }

  C <- as.matrix(C_raw)

  ## Ensure node labels exist
  if (is.null(rownames(X))) rownames(X) <- paste0("Y", seq_len(nrow(X)))
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  if (is.null(colnames(C))) {
    colnames(C)     <- paste0("A_", seq_len(ncol(C)))
    colnames(C_raw) <- colnames(C)
  }

  X_lab <- colnames(X)
  Y_lab <- rownames(X)
  C_lab <- colnames(C_raw)

  X_ids <- .nmfkc_dot_sanitize_id(X_lab)
  Y_ids <- .nmfkc_dot_sanitize_id(Y_lab)
  C_ids <- .nmfkc_dot_sanitize_id(colnames(C))

  ## Helper: strip numeric suffix _k for lag grouping
  base_from_name <- function(nm) sub("_([0-9]+)$", "", nm)

  has_intercept_col <- any(C_lab == "(Intercept)")

  ## -------------------------------------------------------------
  ## Determine lag structure
  ## -------------------------------------------------------------
  A_labels_raw <- unique(base_from_name(C_lab))
  A_labels_raw <- setdiff(A_labels_raw, "(Intercept)")

  A_per_lag <- length(A_labels_raw)
  if (A_per_lag == 0L) {
    stop("Unable to determine covariate base names for lag separation.")
  }

  total_cols        <- ncol(C)
  total_cols_no_int <- if (has_intercept_col) total_cols - 1L else total_cols
  Dmax              <- floor(total_cols_no_int / A_per_lag)

  if (Dmax < 1L) {
    stop("Insufficient columns in C to construct at least one lag block.")
  }

  D <- min(Dmax, as.integer(degree))

  ## -------------------------------------------------------------
  ## Compute edge width scaling maxima
  ## -------------------------------------------------------------
  max_X   <- suppressWarnings(max(X,   na.rm = TRUE))
  max_C   <- suppressWarnings(max(C[, C_lab != "(Intercept)"], na.rm = TRUE))
  max_int <- if (has_intercept_col) {
    suppressWarnings(max(C[, C_lab == "(Intercept)"], na.rm = TRUE))
  } else {
    NA_real_
  }

  pw     <- .nmfkc_dot_penwidth
  digits <- .nmfkc_dot_digits_from_threshold(threshold)
  fmtc   <- function(v) .nmfkc_dot_format_coef(v, digits)

  ## -------------------------------------------------------------
  ## Filter isolated Y nodes (hide.isolated)
  ## -------------------------------------------------------------
  idx_Y <- seq_along(Y_lab)
  if (isTRUE(hide.isolated)) {
    used_Y <- apply(X, 1L, function(row) any(row >= threshold, na.rm = TRUE))
    idx_Y <- which(used_Y)
  }

  ## -------------------------------------------------------------
  ## DOT header
  ## -------------------------------------------------------------
  scr <- .nmfkc_dot_header(graph_name = "NMF_AR", rankdir = rankdir)

  ## -------------------------------------------------------------
  ## Current-time output cluster (T)
  ## -------------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  // Current-time outputs (T)\n',
    .nmfkc_dot_cluster_nodes(
      cluster_id  = "Y",
      title       = "T",
      node_ids    = Y_ids[idx_Y],
      node_labels = Y_lab[idx_Y],
      shape       = "box",
      fill        = fill,
      fillcolor   = "lightblue",
      line_width  = 1.5
    )
  )

  ## -------------------------------------------------------------
  ## Latent factors cluster (X)
  ## -------------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  // Latent variables (X)\n',
    .nmfkc_dot_cluster_nodes(
      cluster_id  = "X",
      title       = "Latent Variables",
      node_ids    = X_ids,
      node_labels = X_lab,
      shape       = "ellipse",
      fill        = fill,
      fillcolor   = "wheat",
      line_width  = 1.0
    )
  )

  ## -------------------------------------------------------------
  ## Edge defaults
  ## -------------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  edge [fontname="Arial", fontsize=8, arrowhead=open];\n'
  )

  ## -------------------------------------------------------------
  ## 1. X → T edges (factor loadings)
  ## -------------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  // X -> T edges (factor loadings)\n',
    '  edge [color="gray0", fontcolor="gray0", style=solid];\n'
  )

  if (is.finite(max_X) && max_X > 0) {
    for (i in idx_Y) {
      for (j in seq_len(ncol(X))) {
        val <- X[i, j]
        if (is.finite(val) && val >= threshold) {
          pen <- pw(val, max_X, weight_scale_xy)
          lab <- fmtc(val)
          scr <- paste0(
            scr,
            sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                    X_ids[j], Y_ids[i], lab, pen)
          )
        }
      }
    }
  }

  ## -------------------------------------------------------------
  ## 2. Optional intercept nodes
  ## -------------------------------------------------------------
  if (isTRUE(intercept) && has_intercept_col && is.finite(max_int) && max_int > 0) {

    scr <- paste0(scr, '\n  // Intercept nodes\n')
    int_col <- which(C_lab == "(Intercept)")

    for (j in seq_len(ncol(X))) {
      val <- C[j, int_col]
      if (is.finite(val) && val >= threshold) {
        pen <- pw(val, max_int, weight_scale_int)
        lab <- fmtc(val)
        node_id <- paste0("Const", j)

        scr <- paste0(
          scr,
          sprintf('  %s [shape=circle, label="%s"];\n', node_id, lab),
          sprintf('  %s -> %s [penwidth=%.2f];\n', node_id, X_ids[j], pen)
        )
      }
    }
  }

  ## -------------------------------------------------------------
  ## 3. Lag blocks and T-k → X edges
  ## -------------------------------------------------------------
  for (k in seq_len(D)) {

    start <- (k - 1L) * A_per_lag
    cols  <- (start + 1L):(start + A_per_lag)

    Ck     <- C[, cols, drop = FALSE]
    Ck_ids <- C_ids[cols]
    Ck_lab <- base_from_name(C_lab[cols])

    ## Skip lag block if no coefficient exceeds threshold
    if (suppressWarnings(max(Ck, na.rm = TRUE)) < threshold) {
      next
    }

    ## ---- Cluster for lag k ----
    st <- sprintf('  subgraph cluster_C%d {label="T-%d" style="rounded";\n', k, k)

    if (fill) {
      st <- paste0(
        st,
        '    node [shape=box, style="filled,rounded", fillcolor="lightcoral", color=black, penwidth=1.5];\n'
      )
    } else {
      st <- paste0(
        st,
        '    node [shape=box, style="rounded", color=black, penwidth=1.5];\n'
      )
    }

    for (j in seq_len(ncol(Ck))) {
      if (suppressWarnings(max(Ck[, j], na.rm = TRUE)) >= threshold) {
        st <- paste0(
          st,
          sprintf('    %s [label="%s"];\n', Ck_ids[j], Ck_lab[j])
        )
      }
    }
    st  <- paste0(st, "  }\n")

    scr <- paste0(scr, "\n  // Lag block T-", k, "\n", st)

    ## ---- Lag → X edges ----
    scr <- paste0(
      scr,
      '  // T-', k, ' -> X edges\n',
      '  edge [color=black, fontcolor=black, style=solid];\n'
    )

    if (is.finite(max_C) && max_C > 0) {
      for (q in seq_len(nrow(Ck))) {
        for (j in seq_len(ncol(Ck))) {
          val <- Ck[q, j]
          if (is.finite(val) && val >= threshold) {
            pen <- pw(val, max_C, weight_scale_lag)
            lab <- fmtc(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                      Ck_ids[j], X_ids[q], lab, pen)
            )
          }
        }
      }
    }
  }

  result <- paste0(scr, "}\n")
  class(result) <- c("nmfkc.ar.DOT", "nmfkc.DOT")
  result
}

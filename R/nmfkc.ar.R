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
#' result <- nmfkc(ar_data$Y, ar_data$A, Q = 1)
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
#' @param Q Rank of the basis matrix.
#' @param degree A vector of candidate lag orders to be evaluated.
#' @param intercept Logical. If TRUE (default), an intercept is added to the covariate matrix.
#' @param plot Logical. If TRUE (default), a plot of the objective function values is drawn.
#' @param ... Additional arguments passed to \code{nmfkc.cv}.
#'
#' @return A list with components:
#' \item{degree}{The lag order that minimizes the cross-validation objective function.}
#' \item{degree.max}{Maximum recommended lag order, computed as \eqn{10 \log_{10}(N)}
#'   following the \code{ar} function in the \pkg{stats} package.}
#' \item{objfunc}{Objective function values for each candidate lag order.}
#' @seealso \code{\link{nmfkc.ar}}, \code{\link{nmfkc.cv}}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example using ts object directly
#' d <- AirPassengers
#'
#' # Selection of degree (using ts object)
#' # Note: Y is automatically transposed if it is a ts object
#' nmfkc.ar.degree.cv(Y=d, Q=1, degree=11:14)

nmfkc.ar.degree.cv <- function(Y, Q=1, degree=1:2, intercept=TRUE, plot=TRUE, ...){

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
  extra_args <- list(...)
  objfuncs <- numeric(length(degree))
  success_status <- logical(length(degree))

  for(i in seq_along(degree)){
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

    objfuncs[i] <- res$obj
    success_status[i] <- res$ok

    end.time <- Sys.time()
    diff.time <- difftime(end.time, start.time, units = "sec")
    diff.time.st <- ifelse(diff.time <= 180, paste0(round(diff.time, 1), "sec"),
                           paste0(round(diff.time/60, 1), "min"))
    message(if (res$ok) diff.time.st else "Skipped (Error)")
  }

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









#' @title Check stationarity of an NMF-VAR model
#' @description
#' \code{nmfkc.ar.stationarity} assesses the dynamic stability of a VAR model
#' by computing the spectral radius of its companion matrix.
#' It returns both the spectral radius and a logical indicator of stationarity.
#'
#' @param x The return value of \code{nmfkc} for a VAR model.
#'
#' @return A list with components:
#' \item{spectral.radius}{Numeric. The spectral radius of the companion matrix. A value less than 1 indicates stationarity.}
#' \item{stationary}{Logical. \code{TRUE} if the spectral radius is less than 1 (i.e., the system is stationary), \code{FALSE} otherwise.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ar}}
#' @examples
#' # Check stationarity of fitted AR model
#' d <- AirPassengers
#' ar_data <- nmfkc.ar(d, degree = 2)
#' result <- nmfkc(ar_data$Y, ar_data$A, Q = 1)
#' nmfkc.ar.stationarity(result)
#'
#' @export

nmfkc.ar.stationarity <- function(x){
  if (!inherits(x, "nmfkc")) stop("'x' must be an object of class 'nmfkc'.")
  if (is.null(x$X) || is.null(x$C)) stop("'x' must contain 'X' and 'C' components (from nmfkc.ar).")
  X <- x$X  # P × Q
  Theta <- x$C  # Q × (P * D [+1] )
  P <- nrow(X)
  Q <- ncol(X)
  if (P == 0 || Q == 0) stop("'X' must have positive dimensions.")
  total_cols <- ncol(Theta)
  has_intercept <- (total_cols - 1) %% P == 0
  D <- if (has_intercept) (total_cols - 1) %/% P else total_cols %/% P
  if (D <= 0) stop("Cannot determine lag order D from dimensions of 'C'.")
  Theta_lags <- if (has_intercept) Theta[, 1:(total_cols - 1), drop = FALSE] else Theta
  Xi_list <- lapply(1:D, function(d){
    cols <- ((d - 1) * P + 1):(d * P)
    X %*% Theta_lags[, cols]})
  companion_matrix <- matrix(0, nrow = P * D, ncol = P * D)
  for (d in 1:D)companion_matrix[1:P, ((d - 1) * P + 1):(d * P)] <- Xi_list[[d]]
  if (D > 1)companion_matrix[(P + 1):(P * D), 1:(P * (D - 1))] <- diag(P * (D - 1))
  rho <- max(Mod(eigen(companion_matrix)$values))
  return(list(spectral.radius = rho, stationary = (rho < 1)))
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
#' @param x A fitted \code{nmfkc} object representing the AR model.
#'   Must contain matrices \code{X} and \code{C}.
#' @param degree Maximum AR lag to visualize.
#' @param intercept Logical; if \code{TRUE}, draws intercept nodes for
#'   columns named "(Intercept)" in matrix \code{C}.
#' @param threshold Minimum coefficient magnitude required to draw an edge.
#' @param rankdir Graphviz rank direction (e.g., \code{"RL"}, \code{"LR"}, \code{"TB"}).
#' @param fill Logical; whether nodes are filled with color.
#' @param weight_scale_xy Scaling factor for edges \eqn{X \rightarrow T}.
#' @param weight_scale_lag Scaling factor for lagged edges \eqn{T-k \rightarrow X}.
#' @param weight_scale_int Scaling factor for intercept edges.
#'
#' @return A character string representing a Graphviz DOT file.
#' @examples
#' d <- AirPassengers
#' ar_data <- nmfkc.ar(d, degree = 2)
#' result <- nmfkc(ar_data$Y, ar_data$A, Q = 1)
#' dot <- nmfkc.ar.DOT(result, degree = 2)
#' cat(dot)
#'
#' @export
nmfkc.ar.DOT <- function(x,
                         degree    = 1,
                         intercept = FALSE,
                         threshold = 0.1,
                         rankdir   = "RL",
                         fill      = TRUE,
                         weight_scale_xy  = 5,
                         weight_scale_lag = 5,
                         weight_scale_int = 3) {

  ## -------------------------------------------------------------
  ## Extract required AR components
  ## -------------------------------------------------------------
  X     <- x$X
  C_raw <- x$C

  if (is.null(X) || is.null(C_raw)) {
    stop("x must contain matrices X and C.")
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
      node_ids    = Y_ids,
      node_labels = Y_lab,
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
    for (i in seq_len(nrow(X))) {
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

  paste0(scr, "}\n")
}

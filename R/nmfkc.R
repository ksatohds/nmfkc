.onAttach <- function(libname, pkgname) {
  pkg_date <- utils::packageDescription("nmfkc", fields = "Date")
  formatted_date <- format(as.Date(pkg_date), "%d %b %Y")
  packageStartupMessage(
    paste("Package: nmfkc (Version", as.character(utils::packageVersion("nmfkc")),
          ", released on", formatted_date, ")")
  )
  packageStartupMessage("https://ksatohds.github.io/nmfkc/")
}




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
#' @param Q Rank of the basis matrix. Must satisfy \eqn{Q \le \min(P,N)}.
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
#' @export

nmfkc.ar.stationarity <- function(x){
  X <- x$X  # P × Q
  Theta <- x$C  # Q × (P * D [+1] )
  P <- nrow(X)
  Q <- ncol(X)
  total_cols <- ncol(Theta)
  has_intercept <- (total_cols - 1) %% P == 0
  D <- if (has_intercept) (total_cols - 1) %/% P else total_cols %/% P
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





#' @title Create a kernel matrix from covariates
#' @description
#' \code{nmfkc.kernel} constructs a kernel matrix from covariate matrices.
#' It supports Gaussian, Exponential, Periodic, Linear, Normalized Linear, and Polynomial kernels.
#'
#' @param U Covariate matrix \eqn{U(K,N) = (u_1, \dots, u_N)}. Each row may be normalized in advance.
#' @param V Covariate matrix \eqn{V(K,M) = (v_1, \dots, v_M)}, typically used for prediction. If \code{NULL}, the default is \code{U}.
#' @param kernel Kernel function to use. Default is \code{"Gaussian"}. Options are \code{"Gaussian"}, \code{"Exponential"}, \code{"Periodic"}, \code{"Linear"}, \code{"NormalizedLinear"}, and \code{"Polynomial"}.
#' @param ... Additional arguments passed to the specific kernel function (e.g., \code{beta}, \code{degree}).
#'
#' @return Kernel matrix \eqn{A(N,M)}.
#' @seealso \code{\link{nmfkc.kernel}}, \code{\link{nmfkc.cv}}
#' @export
#' @source Satoh, K. (2024). Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model.
#'   arXiv preprint arXiv:2403.05359. \url{https://arxiv.org/abs/2403.05359}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' A <- nmfkc.kernel(U,V,beta=28/1000)
#' dim(A)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.kernel <- function(U, V = NULL,
                         kernel = c("Gaussian","Exponential","Periodic",
                                    "Linear","NormalizedLinear","Polynomial"),...){
  k_params <- list(...)
  U <- as.matrix(U); storage.mode(U) <- "double"
  if (is.null(V)) V <- U else V <- as.matrix(V)
  storage.mode(V) <- "double"
  kernel <- match.arg(kernel)
  G <- crossprod(U, V)

  # Determine the specific parameter used for the kernel (e.g., beta or degree)
  # This section extracts the effective parameter value to store in attributes.
  effective_param <- NULL
  if (kernel %in% c("Gaussian", "Exponential")) {
    effective_param <- if (!is.null(k_params$beta)) k_params$beta else 0.5
  } else if (kernel == "Periodic") {
    effective_param <- if (!is.null(k_params$beta)) k_params$beta else c(1,1)
  } else if (kernel == "Polynomial") {
    # For Polynomial, the effective param is (beta, degree) pair
    effective_param <- list(beta = if (!is.null(k_params$beta)) k_params$beta else 0,
                            degree = if (!is.null(k_params$degree)) k_params$degree else 2)
  }

  if (kernel %in% c("Gaussian","Exponential","Periodic")) {
    u2 <- colSums(U * U)
    v2 <- colSums(V * V)
    D2 <- outer(u2, v2, "+") - 2 * G
    D2 <- pmax(D2, 0)
  }

  K <- switch(kernel,
              Gaussian ={
                beta <- if (!is.null(k_params$beta)) k_params$beta else 0.5
                exp(-beta * D2)
              },
              Exponential = {
                beta <- if (!is.null(k_params$beta)) k_params$beta else 0.5
                d <- sqrt(D2)
                exp(-beta * d)
              },
              Periodic = {
                beta <- if (!is.null(k_params$beta)) k_params$beta else c(1,1)
                if (length(beta) < 2L)
                  stop("For 'Periodic', set beta as c(beta1, beta2).")
                d <- sqrt(D2)
                exp(-beta[1] * (sin(beta[2] * d)^2))
              },
              Linear = G,
              NormalizedLinear = {
                nu <- sqrt(colSums(U * U))
                nv <- sqrt(colSums(V * V))
                nu[nu == 0] <- .Machine$double.eps
                nv[nv == 0] <- .Machine$double.eps
                G / outer(nu, nv)
              },
              Polynomial = {
                beta <- if (!is.null(k_params$beta)) k_params$beta else 0
                degree <- if (!is.null(k_params$degree)) k_params$degree else 2
                (G + beta)^degree
              }
  )
  K[K < 0] <- 0
  dimnames(K) <- list(colnames(U), colnames(V))

  # --- Store Kernel Metadata as Attributes (NEW BLOCK) ---
  # These attributes allow nmfkc to identify the kernel type and parameters used.
  attr(K, "params") <- effective_param
  attr(K, "kernel") <- kernel
  attr(K, "function.name") <- "nmfkc.kernel"
  # ------------------------------------------------------

  return(K)
}



#' @title Estimate Gaussian/RBF kernel parameter beta from covariates (supports landmarks)
#'
#' @description
#' Computes a data-driven reference scale for the Gaussian/RBF kernel from covariates
#' using a robust "median nearest-neighbor (or nearest-landmark) distance" heuristic,
#' and returns the corresponding kernel parameter \eqn{\beta}.
#'
#' The Gaussian/RBF kernel is assumed to be written in the form
#' \deqn{k(u,v) = \exp\{-\beta \|u-v\|^2\} = \exp\{-\|u-v\|^2/(2\sigma^2)\},}
#' hence \eqn{\beta = 1/(2\sigma^2)}. This function first estimates a typical distance
#' scale \eqn{\sigma_0} by the median of distances, then sets \eqn{\beta_0 = 1/(2\sigma_0^2)}.
#'
#' If \code{Uk} is \code{NULL}, \eqn{\sigma_0} is estimated as the median of
#' nearest-neighbor distances within \code{U} (excluding self-distance).
#' If \code{Uk} is provided, \eqn{\sigma_0} is estimated as the median of
#' nearest-landmark distances from each sample in \code{U} to its closest landmark in \code{Uk}.
#'
#' To control memory usage for large \code{N} (and \code{M}), distances are computed in blocks.
#' Optionally, columns of \code{U} can be randomly subsampled via \code{sample_size} to reduce cost.
#'
#' @details
#' \strong{Candidate grid:}
#' Along with \code{beta}, the function returns \code{beta_candidates}, a small logarithmic grid
#' suitable for cross-validation.
#'
#' In the landmark case (\code{Uk} provided), the grid is designed to be symmetric on the
#' bandwidth scale \eqn{\sigma} around \eqn{\sigma_0} over one decade:
#' \deqn{\sigma = \sigma_0 \times 10^{t}, \quad t \in \{-1,-2/3,-1/3,0,1/3,2/3,1\}.}
#' Using \eqn{\beta = 1/(2\sigma^2)}, this corresponds to
#' \deqn{\beta = \beta_0 \times 10^{-2t}.}
#'
#' When \code{Uk} is \code{NULL}, a shorter coarse grid may be returned (see \code{Value}).
#'
#' \strong{Notes:}
#' \itemize{
#'   \item When \code{Uk} is identical to \code{U}, the function detects this case and excludes
#'         self-distances (distance 0) to avoid \eqn{\sigma_0=0}.
#'   \item \code{sample_size} performs random subsampling without setting a seed. For reproducible
#'         results, set \code{set.seed()} before calling this function.
#' }
#'
#' @param U A numeric matrix of covariates (\code{K x N}); columns are samples.
#' @param Uk An optional numeric matrix of landmarks (\code{K x M}); columns are landmark points.
#'   If provided, distances are computed from samples in \code{U} to landmarks in \code{Uk}.
#' @param block_size Integer. Number of columns of \code{U} processed per block when computing
#'   distances (controls memory usage). If \code{N <= 1000}, it is automatically set to \code{N}.
#' @param block_size_Uk Integer. Number of columns of \code{Uk} processed per block when \code{Uk}
#'   is not \code{NULL} (controls memory usage). If \code{M <= 2000}, it is automatically set to \code{M}.
#' @param sample_size Integer or \code{NULL}. If not \code{NULL}, randomly subsamples this many columns
#'   of \code{U} (without replacement) before computing distances, to reduce computational cost.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{beta}: Estimated kernel parameter \eqn{\beta_0 = 1/(2\sigma_0^2)}.
#'   \item \code{beta_candidates}: Numeric vector of candidate \eqn{\beta} values (logarithmic grid)
#'         intended for cross-validation.
#'   \item \code{dist_median}: The estimated distance scale \eqn{\sigma_0} (median of nearest-neighbor
#'         or nearest-landmark distances).
#'   \item \code{block_size_used}: The effective block size(s) used. Either a scalar (no \code{Uk}) or
#'         a named vector \code{c(U=..., Uk=...)} when \code{Uk} is provided.
#'   \item \code{sample_size_used}: The number of columns of \code{U} actually used (after subsampling).
#'   \item \code{uk_is_u}: Logical flag indicating whether \code{Uk} was detected as identical to \code{U}
#'         (only returned when \code{Uk} is provided).
#' }
#'
#' @examples
#' # Basic (nearest-neighbor within U)
#' # beta_info <- nmfkc.kernel.beta.nearest.med(U)
#' # beta0 <- beta_info$beta
#' # betas <- beta_info$beta_candidates
#'
#' # With landmarks (nearest-landmark distances)
#' # beta_info <- nmfkc.kernel.beta.nearest.med(U, Uk)
#' # beta0 <- beta_info$beta
#' # betas <- beta_info$beta_candidates
#'
#' @export
nmfkc.kernel.beta.nearest.med <- function(
    U,
    Uk = NULL,
    block_size = 1000,
    block_size_Uk = 2000,
    sample_size = NULL
){
  U <- as.matrix(U)
  storage.mode(U) <- "double"
  N <- ncol(U)
  if (N < 2) stop("U must have at least 2 columns.")

  # ---- optional subsampling over U columns (for speed) ----
  if (!is.null(sample_size)) {
    sample_size <- as.integer(sample_size)
    if (sample_size <= 1) stop("sample_size must be >= 2")
    if (sample_size < N) {
      idxU <- sample.int(N, sample_size)
      U <- U[, idxU, drop = FALSE]
      N <- ncol(U)
    }
  }
  sample_size_used <- N

  # If Uk is NULL, behave like the original function: NN within U (exclude self).
  if (is.null(Uk)) {
    X  <- t(U)                    # N x K
    if (N <= 1000) block_size <- N
    XX <- rowSums(X * X)          # length N
    min_d2 <- rep(Inf, N)

    for (i in seq(1, N, by = block_size)) {
      i2 <- min(i + block_size - 1, N)
      Xi <- X[i:i2, , drop = FALSE]
      Xi_norm <- rowSums(Xi * Xi)

      dist2 <- outer(Xi_norm, XX, "+") - 2 * Xi %*% t(X)
      idx <- i:i2
      dist2[cbind(seq_along(idx), idx)] <- Inf  # exclude self
      dist2[dist2 < 0] <- 0

      nn_local <- apply(dist2, 1, min)
      min_d2[idx] <- pmin(min_d2[idx], nn_local)

      rm(Xi, Xi_norm, dist2); gc(FALSE)
    }

    d_med <- stats::median(sqrt(min_d2))
    beta  <- 1 / (2 * d_med^2)
    return(list(
      beta = beta,
      beta_candidates = beta * 10^c(-2:1),
      dist_median = d_med,
      block_size_used = block_size,
      sample_size_used = sample_size_used
    ))
  }

  # ---- Uk provided: nearest landmark distance (U vs Uk) ----
  Uk <- as.matrix(Uk)
  storage.mode(Uk) <- "double"
  if (nrow(Uk) != nrow(U)) stop("nrow(Uk) must equal nrow(U).")

  M <- ncol(Uk)
  if (M < 1) stop("Uk must have at least 1 column.")

  # We'll compute, for each column of U, min_j ||u_i - uk_j||^2.
  # Block over U (columns) and optionally over Uk to control memory.
  if (N <= 1000) block_size <- N
  if (M <= 2000) block_size_Uk <- M

  U2  <- colSums(U * U)           # length N
  Uk2 <- colSums(Uk * Uk)         # length M
  min_d2 <- rep(Inf, N)

  # Detect "Uk is U" case (same object / identical values) to exclude self-distance.
  uk_is_u <- (M == N) && isTRUE(all.equal(Uk, U, check.attributes = FALSE))

  for (i in seq(1, N, by = block_size)) {
    i2 <- min(i + block_size - 1, N)
    Ui <- U[, i:i2, drop = FALSE]
    Ui2 <- U2[i:i2]

    # we may need to accumulate min over Uk blocks
    cur_min <- rep(Inf, ncol(Ui))

    for (j in seq(1, M, by = block_size_Uk)) {
      j2 <- min(j + block_size_Uk - 1, M)
      Ukj <- Uk[, j:j2, drop = FALSE]
      Ukj2 <- Uk2[j:j2]

      # dist2: (block_U) x (block_Uk)
      # dist2 = ||u||^2 + ||uk||^2 - 2 uk^T u
      G <- crossprod(Ukj, Ui)  # (block_Uk) x (block_U)
      dist2 <- outer(Ui2, Ukj2, "+") - 2 * t(G)
      dist2[dist2 < 0] <- 0

      # exclude self-distances only when Uk == U
      if (uk_is_u) {
        # global indices in U: col(Ui) = (i..i2), col(Ukj) = (j..j2)
        # where they overlap, set that cell to Inf
        gi <- i:i2
        gj <- j:j2
        common <- intersect(gi, gj)
        if (length(common) > 0) {
          # map global index -> local positions
          li <- match(common, gi)
          lj <- match(common, gj)
          dist2[cbind(li, lj)] <- Inf
        }
      }

      cur_min <- pmin(cur_min, apply(dist2, 1, min))

      rm(Ukj, Ukj2, G, dist2); gc(FALSE)
    }

    min_d2[i:i2] <- pmin(min_d2[i:i2], cur_min)
    rm(Ui, Ui2, cur_min); gc(FALSE)
  }

  d_med <- stats::median(sqrt(min_d2))
  beta  <- 1 / (2 * d_med^2)

  t <- c(-1, -2/3, -1/3, 0, 1/3, 2/3, 1)
  betas <- beta * 10^(-2*t)

  list(
    beta = beta,
    beta_candidates = betas,
    dist_median = d_med,
    block_size_used = c(U = block_size, Uk = block_size_Uk),
    sample_size_used = sample_size_used,
    uk_is_u = uk_is_u
  )
}





#' @title Optimize beta of the Gaussian kernel function by cross-validation
#' @description
#' \code{nmfkc.kernel.beta.cv} selects the optimal beta parameter of the kernel function by applying cross-validation over a set of candidate values.
#'
#' @param Y Observation matrix \eqn{Y(P,N)}.
#' @param Q Rank of the basis matrix. Must satisfy \eqn{Q \le \min(P,N)}.
#' @param U Covariate matrix \eqn{U(K,N) = (u_1, \dots, u_N)}. Each row may be normalized in advance.
#' @param V Covariate matrix \eqn{V(K,M) = (v_1, \dots, v_M)}, typically used for prediction. If \code{NULL}, the default is \code{U}.
#' @param beta A numeric vector of candidate kernel parameters to evaluate via cross-validation.
#' @param plot Logical. If TRUE (default), plots the objective function values for each candidate \code{beta}.
#' @param ... Additional arguments passed to \code{nmfkc.cv}.
#'
#' @return A list with components:
#' \item{beta}{The beta value that minimizes the cross-validation objective function.}
#' \item{objfunc}{Objective function values for each candidate \code{beta}.}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' nmfkc.kernel.beta.cv(Y,Q=1,U,V,beta=25:30/1000)
#' A <- nmfkc.kernel(U,V,beta=28/1000)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.kernel.beta.cv <- function(Y,Q=2,U,V=NULL,beta=NULL,plot=TRUE,...){
  extra_args <- list(...)
  kernel_arg_names <- names(formals(nmfkc.kernel))
  cv_arg_names <- names(formals(nmfkc.cv))
  kernel_args_for_call <- extra_args[names(extra_args) %in% kernel_arg_names]
  cv_args_for_call <- extra_args[names(extra_args) %in% cv_arg_names]

  if(is.null(beta)){
    if(is.null(V)) V <- U
    med_args <- c(list(U = V), extra_args[names(extra_args) %in% names(formals(nmfkc.kernel.beta.nearest.med))])
    result.beta <- do.call("nmfkc.kernel.beta.nearest.med", med_args)
    beta <- result.beta$beta_candidates
  }
  objfuncs <- numeric(length(beta))
  for(i in 1:length(beta)){
    start.time <- Sys.time()
    message(paste0("beta=",beta[i],"..."),appendLF=FALSE)

    kernel_args <- c(
      list(U = U, V = V, beta = beta[i]),
      kernel_args_for_call
    )
    A <- do.call("nmfkc.kernel", kernel_args)

    cv_args <- c(
      list(Y = Y, A = A, Q = Q),
      cv_args_for_call
    )
    result <- do.call("nmfkc.cv", cv_args)

    objfuncs[i] <- result$objfunc
    end.time <- Sys.time()
    diff.time <- difftime(end.time,start.time,units="sec")
    diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                           paste0(round(diff.time/60,1),"min"))
    message(diff.time.st)
  }
  i0 <- which.min(objfuncs)
  beta.best <- beta[i0]
  if(plot){
    plot(beta,objfuncs,type="l",col=2,xlab="beta",ylab="objfunc",log="x")
    graphics::points(beta[i0],objfuncs[i0],cex=3,col=2)
    graphics::text(beta,objfuncs,format(beta,scientific=TRUE,digits=5))
  }
  names(objfuncs) <- beta
  result <- list(beta=beta.best,objfunc=objfuncs)
  return(result)
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




#' @title Initialize W (X) matrix using NNDSVDar (Vectorized)
#' @description
#'   Internal function to compute the NNDSVDar (Nonnegative Double SVD and Random)
#'   initialization for the basis matrix X. This method fills the zero entries
#'   generated by basic NNDSVD with small random values scaled by the matrix average,
#'   improving stability and convergence.
#' @param Y Input matrix (P x N)
#' @param Q Rank (number of components)
#' @return X (P x Q) non-negative initial basis matrix
#' @keywords internal
#' @noRd
.nndsvdar <- function(Y, Q) {
  P <- nrow(Y)
  N <- ncol(Y)
  s <- svd(Y, nu = Q, nv = Q)
  W <- matrix(0, P, Q)
  W[, 1] <- sqrt(s$d[1]) * abs(s$u[, 1])
  if (Q > 1) {
    idx <- 2:Q
    U <- s$u[, idx, drop = FALSE]
    V <- s$v[, idx, drop = FALSE]
    D_sqrt <- sqrt(s$d[idx])
    U_p <- pmax(U, 0)
    U_n <- pmax(-U, 0)
    V_p <- pmax(V, 0)
    V_n <- pmax(-V, 0)
    norm_Up <- sqrt(colSums(U_p^2))
    norm_Un <- sqrt(colSums(U_n^2))
    norm_Vp <- sqrt(colSums(V_p^2))
    norm_Vn <- sqrt(colSums(V_n^2))
    Mp <- norm_Up * norm_Vp
    Mn <- norm_Un * norm_Vn
    eps <- .Machine$double.eps
    norm_Up_safe <- norm_Up; norm_Up_safe[norm_Up == 0] <- eps
    norm_Un_safe <- norm_Un; norm_Un_safe[norm_Un == 0] <- eps
    W_pos <- sweep(U_p, 2, (D_sqrt * Mp / norm_Up_safe), FUN = "*")
    W_neg <- sweep(U_n, 2, (D_sqrt * Mn / norm_Un_safe), FUN = "*")
    use_pos <- (Mp > Mn)
    W_combined <- W_pos
    W_combined[, !use_pos] <- W_neg[, !use_pos]
    W[, idx] <- W_combined
  }
  W[is.nan(W)] <- 0
  W[W < 0] <- 0
  avg_Y <- mean(Y)
  idx_zero <- which(W == 0)
  n_zero <- length(idx_zero)
  fill_values <- stats::runif(n_zero) * avg_Y / 100
  W[idx_zero] <- fill_values
  return(W)
}





#' @title Compute a simplified/approximate silhouette coefficient (Internal)
#' @description
#' This internal function computes an approximate version of the silhouette coefficient
#' for clustering results. Unlike the standard definition (e.g., \code{cluster::silhouette}),
#' it determines the nearest neighboring cluster (\code{qn}) by comparing the
#' **Euclidean distances between each sample and the cluster mean vectors (centroids)**,
#' rather than comparing the average distance to all members of other clusters.
#'
#' However, the calculation of \code{a(i)} (average distance to members of the
#' assigned cluster) and \code{b(i)} (average distance to members of the nearest
#' neighboring cluster) uses the pairwise distances of all samples, consistent
#' with the standard definition. This combined approach reduces the computational
#' complexity of the \code{qn} determination step while maintaining consistency
#' for \code{a(i)} and \code{b(i)} using pre-computed distance matrix \code{D}.
#'
#' @param B.prob The coefficient probability matrix (Q x N).
#' @param B.cluster The hard clustering vector derived from \code{B.prob}.
#'
#' @return A list containing:
#' \item{cluster}{Sorted cluster labels.}
#' \item{silhouette}{Sorted silhouette values.}
#' \item{silhouette.means}{Mean of the cluster indices.}
#' \item{silhouette.mean}{Overall mean silhouette score.}
#' @keywords internal
#' @noRd
.silhouette.simple <- function(B.prob,B.cluster){
  if(is.matrix(B.prob)){Q <- nrow(B.prob)}else{Q <- 1}
  if(Q==1){
    return(list(cluster=NA,silhouette=NA,
                silhouette.means=NA,silhouette.mean=NA))
  }else{
    index <-!is.na(B.cluster)
    B.prob <- B.prob[,index, drop=FALSE]
    B.cluster <- B.cluster[index]
    N_samples <- ncol(B.prob)
    D <- as.matrix(stats::dist(t(B.prob)))
    cluster.means <- matrix(0,nrow=Q,ncol=Q)
    ns <- NULL
    cluster.list <- NULL
    for(q in 1:Q){
      ns <- c(ns,sum(B.cluster==q))
      cluster.list <- c(cluster.list,list(which(B.cluster==q)))
      cluster.means[,q] <- rowMeans(B.prob[,B.cluster==q,drop=F])
    }
    B_prob_sq_sum <- colSums(B.prob^2)
    cluster_means_sq_sum <- colSums(cluster.means^2)
    d2_matrix <- outer(B_prob_sq_sum, rep(1, Q)) +
      outer(rep(1, N_samples), cluster_means_sq_sum) -
      2 * t(B.prob) %*% cluster.means
    d2_matrix[d2_matrix < 0] <- 0
    si <- 0*B.cluster
    neighbor.cluster <- 0*B.cluster
    for(q in 1:Q){
      q_indices <- cluster.list[[q]]
      for(i_idx in 1:length(q_indices)){
        i <- q_indices[i_idx]
        di <- d2_matrix[i, ]
        qn <- order(di)[ifelse(order(di)[1]==q, 2, 1)]
        neighbor.cluster[i] <- qn
        if(ns[q]==1){
          si[i] <- 0
        }else{
          ai_distances <- D[i, q_indices]
          ai <- sum(ai_distances) / (ns[q]-1)
          qn_indices <- cluster.list[[qn]]
          bi_distances <- D[i, qn_indices]
          bi <- sum(bi_distances) / ns[qn]
          si[i] <- (bi-ai)/max(ai,bi)
        }
      }
    }
    si.mean <- mean(si)
    si.sort.cluster.means <- 0*ns
    for(q in 1:Q){
      si.sort.cluster.means[q] <- mean(cluster.list[[q]])
    }
    si.sort <- NULL
    si.sort.cluster <- NULL
    for(q in 1:Q){
      si.sort <- c(si.sort,sort(si[cluster.list[[q]]],decreasing=T))
      si.sort.cluster <- c(si.sort.cluster,rep(q,length(cluster.list[[q]])))
    }
    return(list(cluster=si.sort.cluster,silhouette=si.sort,
                silhouette.means=si.sort.cluster.means,silhouette.mean=si.mean))
  }
}





#' @title Parse formula and prepare Y and A matrices
#' @description
#' Internal function to handle formula input, parse variables, and generate
#' the transposed Y and A matrices (features/covariates in rows, samples in columns)
#' as required by the core \code{nmfkc} function. Supports 'data' mode
#' (variable names and dot notation) and direct expression evaluation mode.
#' @param formula A formula object.
#' @param data Optional data frame or environment.
#' @return A list containing the prepared matrices: \code{Y} (P x N) and \code{A} (R x N or NULL).
#' @keywords internal
#' @noRd
.nmfkc_parse_formula <- function(formula, data) {

  # Helper function to abbreviate column messages (shows max_show columns)
  .abbreviate_msg <- function(cols, prefix, max_show = 5) {
    total_len <- base::length(cols)
    if (total_len > max_show) {
      msg <- base::paste0(base::paste(cols[1:max_show], collapse = ", "),
                          base::sprintf("... (Total: %d)", total_len))
    } else {
      msg <- base::paste(cols, collapse = ", ")
    }
    base::message(base::paste(prefix, "created using columns:", msg))
  }

  # Capture the formula object and set the environment
  f <- formula
  f_env <- if (base::missing(data)) base::environment(f) else base::environment()

  # Extract the expressions for Y and A
  Y_expr <- f[[2]]
  A_expr <- f[[3]]

  # --- Input Mode Branching ---

  if (!base::missing(data)) {
    # MODE 1: data provided (Variable names and dot notation)
    data <- base::as.data.frame(data)
    all_cols <- base::names(data)

    # Check if A_expr is structurally present in the formula object.
    A_is_structurally_present <- base::length(f) >= 3

    # Extract Y and A expressions as clean strings using deparse
    Y_expr_str <- base::paste(base::deparse(Y_expr), collapse = " ")

    A_expr_str <- if (A_is_structurally_present) {
      base::paste(base::deparse(A_expr), collapse = " ")
    } else {
      NA # Use simple NA for missing structural part
    }

    # Determine dot notation and missing status based on the safe string A_expr_str
    Y_is_dot <- (Y_expr_str == ".")
    A_is_missing <- base::is.na(A_expr_str)
    A_is_dot <- (!A_is_missing && A_expr_str == ".")

    # --- Check for explicit A omission symbols (0 or -1) ---
    A_is_explicitly_omitted <- (!A_is_missing &&
                                  (A_expr_str == "0" || A_expr_str == "-1" ||
                                     A_expr_str == "0 + ." || A_expr_str == "-1 + ." ||
                                     A_expr_str == ". + 0" || A_expr_str == ". + -1"))

    # Update A_is_missing status if explicitly omitted
    A_is_missing <- A_is_missing || A_is_explicitly_omitted

    Y_cols <- NULL
    A_cols <- NULL

    # Helper to clean and split variable names from expression string
    .extract_cols <- function(expr_str) {
      # Split by '+', '~', ' ' and remove empty elements. This handles "A1 + A2" structure.
      cols <- base::unlist(base::strsplit(expr_str, "[ +~]", fixed = FALSE))
      cols <- cols[cols != ""]
      # Remove explicit omission symbols if they were part of the split result
      cols <- cols[!(cols %in% base::c("0", "-1"))]
      return(cols)
    }

    if (!Y_is_dot) { Y_cols <- .extract_cols(Y_expr_str) }

    # A_cols is only extracted if A is NOT missing and NOT explicitly omitted (and not dot)
    if (!A_is_missing && !A_is_dot) {
      A_cols <- .extract_cols(A_expr_str)
    }

    # Error Check: Both sides cannot be '.'
    if (Y_is_dot && A_is_dot) {
      base::stop("Formula error: '.' cannot be used on both the left (Y) and right (A) sides simultaneously when 'data' is provided.")
    }

    used_cols <- base::unique(base::c(Y_cols, A_cols))
    remaining_cols <- base::setdiff(all_cols, used_cols)

    # Assign dot notation variables
    if (Y_is_dot) {
      Y_cols <- remaining_cols
      if (base::length(Y_cols) == 0) { base::stop("Formula error: '.' for Y resulted in no remaining variables.") }
      .abbreviate_msg(Y_cols, "Y")
    } else if (A_is_dot) {
      A_cols <- remaining_cols
      if (base::length(A_cols) == 0) { base::stop("Formula error: '.' for A resulted in no remaining variables.") }
      # NOTE: Message output for A is suppressed here to avoid duplication
      #       and is handled in the A_mat creation block below.
    }

    # Matrix Creation
    Y_mat <- data[, Y_cols, drop = FALSE]

    # --- A Matrix Finalization and Message Output ---
    if (A_is_missing) { # Catches structural omission, NA, or explicit 0/-1
      A_mat <- NULL
      base::message("A (covariate matrix) is omitted. Performing standard NMF (Y ~ X B).")
    } else {
      # This block now runs if A is explicitly defined with variables or is dot notation
      A_mat <- data[, A_cols, drop = FALSE]
      .abbreviate_msg(A_cols, "A") # Message output for A is performed here ONCE
    }

  } else {
    # MODE 2: data omitted (Direct matrix expression evaluation)

    if (base::is.symbol(Y_expr) && base::as.character(Y_expr) == ".") {
      base::stop("Formula error: '.' is not supported for direct matrix evaluation mode (without 'data' argument).")
    }

    Y_mat <- base::tryCatch({ base::as.matrix(base::eval(Y_expr, envir = f_env)) },
                            error = function(e) { base::stop(base::paste("Error evaluating Y expression:", base::conditionMessage(e))) })

    A_mat <- base::tryCatch({
      # Check if A_expr exists (length >= 3) or is explicitly 0/-1
      A_is_structurally_present <- base::length(f) >= 3
      if (!A_is_structurally_present || (base::is.numeric(A_expr) && A_expr == 0) || (base::is.numeric(A_expr) && A_expr == -1)) { NULL
      } else { base::as.matrix(base::eval(A_expr, envir = f_env)) }
    }, error = function(e) { base::stop(base::paste("Error evaluating A expression:", base::conditionMessage(e))) })

    if (base::is.null(A_mat)) {
      base::message("A (covariate matrix) is omitted. Performing standard NMF (Y ~ X B).")
    }
    # NOTE: No column abbreviation message needed for Mode 2 as matrices are evaluated directly.
  }

  # Transpose: R rows=samples -> nmfkc columns=samples (P x N or R x N)
  # NOTE: R data frames/matrices are typically samples x variables. nmfkc requires variables x samples.
  Y <- base::t(Y_mat)

  if (!base::is.null(A_mat)) {
    A <- base::t(A_mat)
    if (base::ncol(Y) != base::ncol(A)) {
      base::stop(base::paste0("Dimension error: Number of columns (samples) in Y (", base::ncol(Y), ") must match number of columns (samples) in A (", base::ncol(A), ")."))
    }
  } else {
    A <- NULL
  }

  return(base::list(Y = Y, A = A))
}












#' @title Optimize NMF with kernel covariates (Full Support for Missing Values)
#' @description
#' \code{nmfkc} fits a nonnegative matrix factorization with kernel covariates
#' under the tri-factorization model \eqn{Y \approx X C A = X B}.
#'
#' This function supports two major input modes:
#' 1. **Matrix Mode (Existing)**: \code{nmfkc(Y=matrix, A=matrix, ...)}
#' 2. **Formula Mode (New)**: \code{nmfkc(formula=Y_vars ~ A_vars, data=df, rank=Q, ...)}
#'
#' The rank of the basis matrix can be specified using either the \code{rank} argument
#' (preferred for formula mode) or the hidden \code{Q} argument (for backward compatibility).
#'
#' @param Y Observation matrix, OR a formula object if \code{data} is supplied.
#' @param A Covariate matrix. Default is \code{NULL} (no covariates).
#' @param rank Integer. The rank of the basis matrix \eqn{X} (Q). Preferred over \code{Q}.
#' @param data Optional. A data frame from which variables in the formula should be taken.
#' @param epsilon Positive convergence tolerance.
#' @param maxit Maximum number of iterations.
#' @param ... Additional arguments passed for fine-tuning regularization, initialization, constraints,
#'   and output control. This includes the backward-compatible arguments \code{Q} and \code{method}.
#'   \itemize{
#'     \item \code{Y.weights}: Optional numeric matrix (P x N) or vector (length N).
#'       0 indicates missing/ignored values. If NULL (default), weights are automatically
#'       set to 0 for NAs in Y, and 1 otherwise.
#'     \item \code{X.L2.ortho}: Nonnegative penalty parameter for the orthogonality of \eqn{X} (default: 0).
#'       It minimizes the off-diagonal elements of the Gram matrix \eqn{X^\top X}, reducing the correlation
#'       between basis vectors (conceptually minimizing \eqn{\| X^\top X - \mathrm{diag}(X^\top X) \|_F^2}).
#'       (Formerly \code{lambda.ortho}).
#'     \item \code{B.L1}: Nonnegative penalty parameter for L1 regularization on \eqn{B = C A} (default: 0).
#'       Promotes **sparsity in the coefficients**. (Formerly \code{gamma}).
#'     \item \code{C.L1}: Nonnegative penalty parameter for L1 regularization on \eqn{C} (default: 0).
#'       Promotes **sparsity in the parameter matrix**. (Formerly \code{lambda}).
#'     \item \code{Q}: Backward-compatible name for the rank of the basis matrix (Q).
#'     \item \code{method}: Objective function: Euclidean distance \code{"EU"} (default) or Kullback–Leibler divergence \code{"KL"}.
#'     \item \code{X.restriction}: Constraint for columns of \eqn{X}. Options: \code{"colSums"} (default), \code{"colSqSums"}, \code{"totalSum"}, or \code{"fixed"}.
#'     \item \code{X.init}: Method for initializing the basis matrix \eqn{X}. Options: \code{"kmeans"} (default), \code{"runif"}, \code{"nndsvd"}, or a user-specified matrix.
#'     \item \code{nstart}: Number of random starts for \code{kmeans} when initializing \eqn{X} (default: 1).
#'     \item \code{seed}: Integer seed for reproducibility (default: 123).
#'     \item \code{C.init}: Optional numeric matrix giving the initial value of the parameter matrix \eqn{C}
#'       (i.e., \eqn{\Theta}). If \code{A} is \code{NULL}, \code{C} has dimension \eqn{Q \times N} (equivalently \eqn{B});
#'       otherwise, \code{C} has dimension \eqn{Q \times K} where \eqn{K = nrow(A)}. Default initializes all entries to 1.
#'     \item \code{prefix}: Prefix for column names of \eqn{X} and row names of \eqn{B} (default: "Basis").
#'     \item \code{print.trace}: Logical. If \code{TRUE}, prints progress every 10 iterations (default: \code{FALSE}).
#'     \item \code{print.dims}: Logical. If \code{TRUE} (default), prints matrix dimensions and elapsed time.
#'     \item \code{save.time}: Logical. If \code{TRUE} (default), skips some post-computations (e.g., CPCC, silhouette) to save time.
#'     \item \code{save.memory}: Logical. If \code{TRUE}, performs only essential computations (implies \code{save.time = TRUE}) to reduce memory usage (default: \code{FALSE}).
#'   }
#' @return A list with components:
#' \item{call}{The matched call, as captured by `match.call()`.}
#' \item{dims}{A character string summarizing the matrix dimensions of the model.}
#' \item{runtime}{A character string summarizing the computation time.}
#' \item{X}{Basis matrix. Column normalization depends on \code{X.restriction}.}
#' \item{B}{Coefficient matrix \eqn{B = C A}.}
#' \item{XB}{Fitted values for \eqn{Y}.}
#' \item{C}{Parameter matrix.}
#' \item{B.prob}{Soft-clustering probabilities derived from columns of \eqn{B}.}
#' \item{B.cluster}{Hard-clustering labels (argmax over \eqn{B.prob} for each column).}
#' \item{X.prob}{Row-wise soft-clustering probabilities derived from \eqn{X}.}
#' \item{X.cluster}{Hard-clustering labels (argmax over \eqn{X.prob} for each row).}
#' \item{A.attr}{List of attributes of the input covariate matrix \code{A}, containing metadata like lag order and intercept status if created by \code{nmfkc.ar} or \code{nmfkc.kernel}.}
#' \item{objfunc}{Final objective value.}
#' \item{objfunc.iter}{Objective values by iteration.}
#' \item{r.squared}{Coefficient of determination \eqn{R^2} between \eqn{Y} and \eqn{X B}.}
#' \item{sigma}{The residual standard error, representing the typical deviation of the observed values \eqn{Y} from the fitted values \eqn{X B}.}
#' \item{criterion}{A list of selection criteria, including \code{ICp}, \code{CPCC}, \code{silhouette}, \code{AIC}, and \code{BIC}.}
#' @seealso \code{\link{nmfkc.cv}}, \code{\link{nmfkc.rank}}, \code{\link{nmfkc.kernel}}, \code{\link{nmfkc.ar}}, \code{\link{predict.nmfkc}}
#' @export
#' @source Satoh, K. (2024). Applying Non-negative Matrix Factorization with Covariates
#'   to the Longitudinal Data as Growth Curve Model. arXiv:2403.05359.
#'   \url{https://arxiv.org/abs/2403.05359}
#' @references
#' Ding, C., Li, T., Peng, W., & Park, H. (2006). Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering.
#'   In \emph{Proceedings of the 12th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining} (pp. 126–135).
#'   \doi{10.1145/1150402.1150420}
#' Potthoff, R. F., & Roy, S. N. (1964). A generalized multivariate analysis of variance model useful especially for growth curve problems.
#'   \emph{Biometrika}, 51, 313–326. \doi{10.2307/2334137}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example 1. Matrix Mode (Existing)
#' library(nmfkc)
#' X <- cbind(c(1,0,1),c(0,1,0))
#' B <- cbind(c(1,0),c(0,1),c(1,1))
#' Y <- X %*% B
#' rownames(Y) <- paste0("P",1:nrow(Y))
#' colnames(Y) <- paste0("N",1:ncol(Y))
#' print(X); print(B); print(Y)
#' library(nmfkc)
#' res <- nmfkc(Y,Q=2,epsilon=1e-6)
#' res$X
#' res$B
#'
#' # Example 2. Formula Mode (New)
#' # dummy_data <- data.frame(Y1=rpois(10,5), Y2=rpois(10,10), A1=1:10, A2=rnorm(10,5))
#' # res_f <- nmfkc(Y1 + Y2 ~ A1 + A2, data=dummy_data, rank=2)
#'
nmfkc <- function(Y, A=NULL, rank=NULL, data, epsilon=1e-4, maxit=5000, ...){
  # A small constant for numerical stability to prevent division by zero and log(0).
  .eps <- 1e-10

  extra_args <- base::list(...)

  if(!base::is.null(A)) {
    if(any(is.na(A))) base::stop("Covariate matrix A contains NAs. Please impute or remove them.")
    if(base::min(A, na.rm=TRUE)<0) base::stop("The matrix A should be non-negative.")
  }
  if(base::min(Y, na.rm=TRUE)<0) base::stop("The matrix Y should be non-negative.")

  # --- 1. Parameter Extraction ---
  Q_hidden <- if (!base::is.null(extra_args$Q)) extra_args$Q else NULL
  Q_val <- if (!base::is.null(rank)) rank else if (!base::is.null(Q_hidden)) Q_hidden else 2
  Q <- Q_val

  C.L1 <- if (!base::is.null(extra_args$C.L1)) extra_args$C.L1 else 0
  B.L1 <- if (!base::is.null(extra_args$B.L1)) extra_args$B.L1 else 0
  X.L2.ortho <- if (!base::is.null(extra_args$X.L2.ortho)) extra_args$X.L2.ortho else 0

  if (C.L1 == 0 && !base::is.null(extra_args$lambda)) C.L1 <- extra_args$lambda
  if (B.L1 == 0 && !base::is.null(extra_args$gamma)) B.L1 <- extra_args$gamma
  if (X.L2.ortho == 0 && !base::is.null(extra_args$lambda.ortho)) X.L2.ortho <- extra_args$lambda.ortho

  method <- if (!base::is.null(extra_args$method)) extra_args$method else "EU"
  X.restriction <- if (!base::is.null(extra_args$X.restriction)) extra_args$X.restriction else "colSums"
  X.init <- if (!base::is.null(extra_args$X.init)) extra_args$X.init else "kmeans"
  nstart <- if (!base::is.null(extra_args$nstart)) extra_args$nstart else 1
  seed <- if (!base::is.null(extra_args$seed)) extra_args$seed else 123
  C.init <- if (!is.null(extra_args$C.init)) extra_args$C.init else NULL

  prefix <- if (!base::is.null(extra_args$prefix)) extra_args$prefix else "Basis"
  print.trace <- if (!base::is.null(extra_args$print.trace)) extra_args$print.trace else FALSE
  print.dims <- if (!base::is.null(extra_args$print.dims)) extra_args$print.dims else TRUE
  save.time <- if (!base::is.null(extra_args$save.time)) extra_args$save.time else TRUE
  save.memory <- if (!base::is.null(extra_args$save.memory)) extra_args$save.memory else FALSE

  Y.weights <- if (!base::is.null(extra_args$Y.weights)) extra_args$Y.weights else NULL

  # --- 2. Input Data Preparation ---
  if (base::inherits(Y, "formula")) {
    data_list <- .nmfkc_parse_formula(formula = Y, data = data)
    Y <- data_list$Y
    A <- data_list$A
  } else {
    if(base::is.vector(Y)) Y <- base::matrix(Y,nrow=1)
    if(!base::is.matrix(Y)) Y <- base::as.matrix(Y)
  }

  # === Weights Handling ===
  # 1. Vector Expansion (Column-wise weights)
  if (!is.null(Y.weights) && is.vector(Y.weights)) {
    if (length(Y.weights) == ncol(Y)) {
      # Expand vector to matrix: each row gets the same weight vector
      Y.weights <- matrix(Y.weights, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
    } else if (length(Y.weights) == 1) {
      Y.weights <- matrix(Y.weights, nrow = nrow(Y), ncol = ncol(Y))
    } else {
      stop("Length of Y.weights vector must match ncol(Y) (or be 1).")
    }
  }

  # 2. Check Dimensions & Handle NAs
  if (is.null(Y.weights)) {
    if (any(is.na(Y))) {
      Y.weights <- matrix(1, nrow=nrow(Y), ncol=ncol(Y))
      Y.weights[is.na(Y)] <- 0
      Y[is.na(Y)] <- 0
      if(print.dims) message("Notice: Missing values (NA) in Y were treated as weights=0.")
    } else {
      Y.weights <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
    }
  } else {
    if (!is.matrix(Y.weights)) Y.weights <- as.matrix(Y.weights)
    if (!all(dim(Y.weights) == dim(Y))) stop("Dimension mismatch between Y and Y.weights.")
    Y.weights[is.na(Y.weights)] <- 0
    Y[is.na(Y)] <- 0
    Y <- Y * Y.weights
  }

  # --- 3. Algorithm Setup ---
  X.restriction <- base::match.arg(X.restriction, base::c("colSums", "colSqSums", "totalSum","fixed"))
  xnorm <- base::switch(X.restriction,
                        colSums   = function(X) base::sweep(X, 2, base::colSums(X), "/"),
                        colSqSums = function(X) base::sweep(X, 2, base::sqrt(base::colSums(X^2)), "/"),
                        totalSum  = function(X) X / base::sum(X),
                        fixed = function(X) X
  )

  if(base::is.null(A)){
    dims <- base::sprintf("Y(%d,%d)~X(%d,%d)B(%d,%d)",
                          base::nrow(Y),base::ncol(Y),base::nrow(Y),Q,Q,base::ncol(Y))
  }else{
    dims <- base::sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)=XB(%d,%d)",
                          base::nrow(Y),base::ncol(Y),base::nrow(Y),Q,Q,base::nrow(A),base::nrow(A),base::ncol(Y),Q,base::ncol(Y))
  }
  if(print.dims) base::message(base::paste0(dims,"..."),appendLF=FALSE)
  start.time <- base::Sys.time()

  if(!base::is.null(A)) if(base::min(A)<0) base::stop("The matrix A should be non-negative.")

  # [Fix 1] Added na.rm=TRUE to min(Y) check
  if(base::min(Y, na.rm=TRUE)<0) base::stop("The matrix Y should be non-negative.")

  # Initialize X
  is.X.scalar <- FALSE
  if(nrow(Y)>=2){
    if(ncol(Y)>=Q){
      if(ncol(Y)==Q){
        X <- Y
      }else{
        # [Fix 2] Create Y_init for initialization (impute NAs with row means)
        Y_init <- Y
        if (is.matrix(Y.weights) && any(Y.weights == 0)) {
          # Calculate weighted mean (Since Y is zero-filled, sum(Y) is equivalent to sum(Y*W))
          row_means <- rowSums(Y) / (rowSums(Y.weights) + .eps)
          mask_binary <- (Y.weights > 0)
          idx_missing <- which(!mask_binary, arr.ind = TRUE)
          if (nrow(idx_missing) > 0) {
            Y_init[idx_missing] <- row_means[idx_missing[,1]]
          }
        }

        if (is.matrix(X.init)) {
          X <- X.init
        } else if (is.character(X.init)) {
          if (!is.null(seed)) set.seed(seed)
          if (X.init == "kmeans") {
            # Use Y_init
            res.kmeans <- tryCatch(stats::kmeans(t(Y_init), centers = Q, iter.max = maxit, nstart = nstart),
                                   error = function(e) NULL)
            if(!is.null(res.kmeans)) X <- t(res.kmeans$centers) else X <- matrix(stats::runif(nrow(Y) * Q), nrow = nrow(Y), ncol = Q)
          } else if (X.init == "runif") {
            X <- matrix(stats::runif(nrow(Y) * Q), nrow = nrow(Y), ncol = Q)
          } else {
            if (Q > min(nrow(Y), ncol(Y))) {
              X <- matrix(stats::runif(nrow(Y) * Q), nrow = nrow(Y), ncol = Q)
            } else {
              X <- .nndsvdar(Y_init, Q)
            }
          }
        }
      }
    }else{
      stop("It does not hold Q<=N (ncol(Y)).")
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
    is.X.scalar <- TRUE
  }
  X <- xnorm(X)

  # [FIX: Initialization of tX]
  # Initialize tX here so it exists even if the X update loop is skipped (e.g., scalar X)
  tX <- t(X)

  if(is.null(A)){
    if(is.null(C.init)) C <- matrix(1, nrow=Q, ncol=ncol(Y)) else C <- C.init
  }else{
    if(is.null(C.init)) C <- matrix(1, nrow=Q, ncol=nrow(A)) else C <- C.init
  }
  hasA <- !is.null(A)

  ones_QN <- matrix(1, nrow=Q, ncol=ncol(Y))
  if(hasA) {
    At <- t(A)
    ones_QN_At <- ones_QN %*% At
  }

  epsilon.iter <- Inf
  objfunc.iter <- 0*(1:maxit)
  i_end <- NULL

  # --- 4. Main Loop (Weighted) ---
  for(i in 1:maxit){
    if(is.null(A)) B <- C else B <- C %*% A
    XB <- X %*% B
    if(print.trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))

    if(method=="EU"){
      if(!is.X.scalar & X.restriction!="fixed"){
        num_X <- (Y.weights * Y) %*% t(B)
        den_X <- (Y.weights * XB) %*% t(B)
        if (X.L2.ortho > 0) {
          XtX <- crossprod(X); diag(XtX) <- 0
          den_X <- den_X + X.L2.ortho * (X %*% XtX)
        }
        X <- X * (num_X / (den_X + .eps))
        X <- xnorm(X)
        tX <- t(X)
      }
      if(is.null(A)) {
        num_C <- tX %*% (Y.weights * Y)
        den_C <- tX %*% (Y.weights * XB)
        if (C.L1 != 0) den_C <- den_C + (C.L1/2) * ones_QN
        if (B.L1 != 0) den_C <- den_C + (B.L1/2) * ones_QN
        C <- C * (num_C / (den_C + .eps))
      } else {
        num_C <- tX %*% (Y.weights * Y) %*% At
        den_C <- tX %*% (Y.weights * XB) %*% At
        if (C.L1 != 0) den_C <- den_C + (C.L1/2) * matrix(1, nrow=Q, ncol=nrow(A))
        if (B.L1 != 0) den_C <- den_C + (B.L1/2) * ones_QN_At
        C <- C * (num_C / (den_C + .eps))
      }
      resid <- Y.weights * (Y - XB)
      obj <- sum(resid^2)

    }else{ # KL
      if(!is.X.scalar & X.restriction!="fixed"){
        ratio <- Y.weights * (Y / (XB + .eps))
        num_X <- ratio %*% t(B)
        den_X <- Y.weights %*% t(B)
        if (X.L2.ortho > 0) {
          XtX <- crossprod(X); diag(XtX) <- 0
          den_X <- den_X + X.L2.ortho * (X %*% XtX)
        }
        X <- X * (num_X / (den_X + .eps))
        X <- xnorm(X)
        tX <- t(X)
      }
      if(is.null(A)) {
        ratio <- Y.weights * (Y / (XB + .eps))
        num_C <- tX %*% ratio
        den_C <- tX %*% Y.weights
        if (C.L1 != 0) den_C <- den_C + C.L1 * ones_QN
        if (B.L1 != 0) den_C <- den_C + B.L1 * ones_QN
        C <- C * (num_C / (den_C + .eps))
      } else {
        ratio <- Y.weights * (Y / (XB + .eps))
        num_C <- tX %*% ratio %*% At
        den_C <- tX %*% Y.weights %*% At
        if (C.L1 != 0) den_C <- den_C + C.L1 * matrix(1, nrow=Q, ncol=nrow(A))
        if (B.L1 != 0) den_C <- den_C + B.L1 * ones_QN_At
        C <- C * (num_C / (den_C + .eps))
      }
      term1 <- - (Y.weights * Y) * log(XB + .eps)
      term2 <- Y.weights * XB
      obj <- sum(term1 + term2)
    }

    if (C.L1 != 0) obj <- obj + C.L1 * sum(C)
    if (B.L1 != 0) obj <- obj + if (hasA) B.L1 * sum(C %*% A) else B.L1 * sum(C)
    if (X.L2.ortho != 0) {
      XtX <- crossprod(X); diag(XtX) <- 0
      obj <- obj + (X.L2.ortho / 2) * sum(XtX^2)
    }
    objfunc.iter[i] <- obj

    if(i>=10){
      #epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1]) / pmax(abs(objfunc.iter[i]), 1)
      if(epsilon.iter <= abs(epsilon)){ i_end <- i; break }
    }
  }

  if(is.null(A)) B <- C else B <- C %*% A
  XB <- X %*% B

  if(method=="EU"){
    resid <- Y.weights * (Y - XB)
    objfunc <- sum(resid^2)
  } else {
    term1 <- - (Y.weights * Y) * log(XB + .eps)
    term2 <- Y.weights * XB
    objfunc <- sum(term1 + term2)
  }

  if(!is.null(i_end)){ objfunc.iter <- objfunc.iter[10:i_end]
  } else if (i >= 10){ objfunc.iter <- objfunc.iter[10:i]
  } else { objfunc.iter <- objfunc.iter[1:i] }

  if(ncol(X) > 1 & X.restriction != "fixed"){
    index <- order(matrix(1:nrow(X)/nrow(X),nrow=1) %*% X)
    X <- X[,index,drop=FALSE]; B <- B[index,,drop=FALSE]; C <- C[index,,drop=FALSE]
  }
  rownames(C) <- paste0(prefix,1:nrow(C))
  rownames(X) <- rownames(Y); colnames(X) <- paste0(prefix,1:ncol(X))
  rownames(B) <- paste0(prefix,1:nrow(B)); colnames(B) <- colnames(Y)

  N_obs <- sum(Y.weights > 0)

  if(X.restriction=="fixed") nparam.X <- 0 else if(X.restriction=="totalSum") nparam.X <- prod(dim(X)) - 1 else nparam.X <- (nrow(X)-1)*ncol(X)
  if(is.null(A)) nparam <- nparam.X + prod(dim(B)) else nparam <- nparam.X + prod(dim(C))

  # Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.
  if (method == "EU") {
    sigma2 <- objfunc / N_obs
    P <- nrow(Y)
    T <- ncol(Y)
    g_icp1 <- (P + T) / (P * T) * log(P * T)
    g_icp2 <- (P + T) / (P * T) * log(min(P, T))
    g_icp3 <- log(min(P, T)) / min(P, T)
    ICp1 <- log(sigma2) + nparam * g_icp1
    ICp2 <- log(sigma2) + nparam * g_icp2
    ICp3 <- log(sigma2) + nparam * g_icp3
    ICp <- min(ICp1, ICp2, ICp3)
    AIC <- N_obs * log(sigma2) + 2 * nparam
    BIC <- N_obs * log(sigma2) + nparam * log(N_obs)
  } else {
    ICp <- NA;ICp1 <- NA;ICp2 <-NA;ICp3 <- NA;
    AIC <- NA
    BIC <- NA
  }

  if(save.memory==FALSE){
    # [Fix 3] Calculate statistics using valid data only
    valid_idx <- (Y.weights > 0)
    if(any(valid_idx)) {
      r2 <- stats::cor(XB[valid_idx], Y[valid_idx])^2
      sigma <- stats::sd(Y[valid_idx] - XB[valid_idx])
      mae <- mean(abs(Y[valid_idx] - XB[valid_idx]))
    } else { r2 <- NA; sigma <- NA }

    B.prob <- t( t(B) / (colSums(B) + .eps) )
    if(Q > 1){
      B.prob.sd.min <- min(apply(B.prob,1,stats::sd))
      B.prob.max.mean <- mean(apply(B.prob,2,max))
      p <- B.prob + .eps # avoid log(0)
      B.prob.entropy.mean <- -mean(colSums(p * log(p))) / log(Q)
    } else {
      B.prob.sd.min <- 0
      B.prob.entropy.mean <- 0
      B.prob.max.mean <- 1
    }
    B.cluster <- apply(B.prob,2,which.max)
    B.cluster[colSums(B.prob)==0] <- NA
    X.prob <- X / (rowSums(X) + .eps)
    X.cluster <- apply(X.prob,1,which.max)

    # [Fix 4] Skip dist.cor if there are missing values
    if(save.time || any(Y.weights == 0)){
      silhouette <- NA; CPCC <- NA; dist.cor <- NA
    }else{
      silhouette <- .silhouette.simple(B.prob,B.cluster)
      dist.cor <- stats::cor(as.vector(stats::dist(t(Y))),as.vector(stats::dist(t(B))))
      if(Q>=2){
        M <- t(B.prob) %*% B.prob
        h.dist <- as.matrix(stats::cophenetic(stats::hclust(stats::as.dist(1-M))))
        up <- upper.tri(M)
        CPCC <- stats::cor(h.dist[up],(1-M)[up])
      }else{ CPCC <- NA }
    }
  }else{
    r2 <- NA; sigma <- NA; mae <- NA
    B.prob <- NA; B.cluster <- NA
    B.prob.sd.min <- NA; B.prob.max.mean <- NA; B.prob.entropy.mean <- NA
    XB <- NA; X.prob <- NA; X.cluster <- NA;
    silhouette <- NA; CPCC <- NA; dist.cor <- NA
  }
  if(epsilon.iter > abs(epsilon)) warning(paste0("maximum iterations (",maxit,") reached..."))
  end.time <- Sys.time()
  diff.time.st <- paste0(round(difftime(end.time,start.time,units="sec"),1),"sec")
  if(print.dims) message(diff.time.st)

  n.missing <- sum(Y.weights == 0)
  n.total <- prod(dim(Y))
  A.attr <- NULL
  if (!is.null(A)) A.attr <- attributes(A)

  result <- list(
    call      = match.call(),
    dims      = dims,
    runtime   = diff.time.st,
    method    = method,
    X         = X,
    B         = B,
    XB        = XB,
    C         = C,
    B.prob    = B.prob,
    B.cluster = B.cluster,
    X.prob    = X.prob,
    X.cluster = X.cluster,
    A.attr    = A.attr,
    n.missing = n.missing,
    n.total   = n.total,
    rank      = Q,
    objfunc   = objfunc,
    objfunc.iter = objfunc.iter,
    r.squared = r2,
    sigma     = sigma,
    mae =mae,
    criterion = list(
      B.prob.sd.min = B.prob.sd.min,
      B.prob.max.mean = B.prob.max.mean,
      B.prob.entropy.mean= B.prob.entropy.mean,
      ICp1 = ICp1, ICp2 = ICp2, ICp3 = ICp3, ICp = ICp,
      AIC  = AIC,  BIC  = BIC,
      silhouette = silhouette,
      CPCC       = CPCC,
      dist.cor   = dist.cor
    )
  )
  class(result) <- "nmfkc"
  return(result)
}















#' @title Plot method for objects of class \code{nmfkc}
#' @description
#' \code{plot.nmfkc} produces a diagnostic plot for the return value of
#' \code{nmfkc}, showing the objective function across iterations.
#'
#' @param x An object of class \code{nmfkc}, i.e., the return value of \code{nmfkc}.
#' @param ... Additional arguments passed to the base \code{\link{plot}} function.
#' @export
plot.nmfkc <- function(x,...){
  extra_args <- list(...)
  args <- list(x = x$objfunc.iter)
  if(is.null(extra_args$main)) args$main <- paste0("r.squared=",round(x$r.squared,3))
  if(is.null(extra_args$xlab)) args$xlab <- "iter"
  if(is.null(extra_args$ylab)) args$ylab <- "objfunc"
  all_args <- c(args, extra_args)
  do.call("plot",all_args)
}




#' @export
summary.nmfkc <- function(object, ...) {
  ans <- list()
  ans$call <- object$call
  ans$dims <- object$dims
  ans$rank <- object$rank
  ans$runtime <- object$runtime

  # Missing values
  if(!is.null(object$n.missing) && !is.null(object$n.total)){
    ans$n.missing <- object$n.missing
    ans$prop.missing <- object$n.missing / object$n.total * 100
  } else {
    ans$n.missing <- NULL
  }

  ans$method <- object$method
  ans$iter <- length(object$objfunc.iter)
  ans$objfunc <- object$objfunc
  ans$r.squared <- object$r.squared
  ans$sigma <- object$sigma
  ans$mae <- object$mae
  if(!is.null(object$criterion)){
    ans$ICp <- object$criterion$ICp
  } else {
    ans$ICp <- NULL
  }

  # --- Diagnostics (Sparsity & Clustering Quality) ---

  # 1. Basis (X)
  if (!is.null(object$X) && is.matrix(object$X)) {
    # Sparsity: Proportion of elements close to zero (< 1e-4)
    ans$X.sparsity <- mean(object$X < 1e-4)
  }

  # 2. Probabilities (B.prob)
  if (!is.null(object$B.prob)){
    # Sparsity
    ans$B.prob.sparsity <- mean(object$B.prob < 1e-4)
    ans$B.prob.entropy.mean <- object$criterion$B.prob.entropy.mean
    ans$B.prob.max.mean <- object$criterion$B.prob.max.mean
  }

  class(ans) <- "summary.nmfkc"
  return(ans)
}

#' @export
print.summary.nmfkc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Dimensions:", x$dims, "\n")
  if(!is.null(x$rank)) cat("Rank (Q):   ", x$rank, "\n")
  cat("Runtime:    ", x$runtime, "\n")
  if (!is.null(x$method)) cat("Method:     ", x$method, "\n")
  cat("Iterations: ", x$iter, "\n")

  if (!is.null(x$n.missing)) {
    cat("Missing:    ", x$n.missing,
        sprintf("(%.1f%%)", x$prop.missing), "\n")
  }

  cat("\nStatistics:\n")
  cat("  Objective function:  ", format(x$objfunc, digits = digits), "\n")
  cat("  Multiple R-squared:  ", format(x$r.squared, digits = digits), "\n")
  cat("  Residual Std Error:  ", format(x$sigma, digits = digits), "\n")
  cat("  Mean Absolute Error: ", format(x$mae, digits = digits), "\n")

  if (!is.null(x$ICp)) {
    cat("  ICp:                 ", format(x$ICp, digits = digits), "\n")
  }

  cat("\nStructure Diagnostics:\n")
  if (!is.null(x$X.sparsity)) {
    cat("  Basis (X) Sparsity:   ", sprintf("%.1f%%", x$X.sparsity * 100), "(< 1e-4)\n")
  }
  if (!is.null(x$B.prob.sparsity)) {
    cat("  Coef (B) Sparsity:    ", sprintf("%.1f%%", x$B.prob.sparsity * 100), "(< 1e-4)\n")
  }
  if(!is.null(x$B.prob.entropy.mean)){
    cat("  Clustering Entropy:   ", format(x$B.prob.entropy.mean, digits = digits),
        "(range: 0-1, closer to 0 is better)\n")
    cat("  Clustering Crispness: ", format(x$B.prob.max.mean, digits = digits),
        "(range: 0-1, closer to 1 is better)\n")
  }
  cat("\n")
  invisible(x)
}














#' @title Normalize a matrix to the range \eqn{[0,1]}
#' @description
#' \code{nmfkc.normalize} rescales the values of a matrix to lie between 0 and 1
#' using the column-wise minimum and maximum values of a reference matrix.
#'
#' @param x A numeric matrix (or vector) to be normalized.
#' @param ref A reference matrix from which the column-wise minima and maxima are taken.
#'   Default is \code{x}.
#'
#' @return A matrix of the same dimensions as \code{x}, with each column rescaled to the \eqn{[0,1]} range.
#' @seealso \code{\link{nmfkc.denormalize}}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' x <- nmfkc.normalize(iris[,-5])
#' apply(x,2,range)
nmfkc.normalize <- function(x,ref=x){
  if(is.vector(x)==TRUE){
    x <- matrix(x,ncol=1)
    ref <- matrix(ref,ncol=1)
  }
  r <- apply(ref,2,range)
  denom <- r[2, ] - r[1, ]
  denom[denom == 0] <- 1   # 0幅列はそのままにする（全0返しが嫌なら別方針に）
  y <- sweep(x, 2, r[1, ], FUN = "-")
  y <- sweep(y, 2, denom,   FUN = "/")
  return(y)
}




#' @title Denormalize a matrix from \eqn{[0,1]} back to its original scale
#' @description
#' \code{nmfkc.denormalize} rescales a matrix with values in \eqn{[0,1]} back to its
#' original scale using the column-wise minima and maxima of a reference matrix.
#'
#' @param x A numeric matrix (or vector) with values in \eqn{[0,1]} to be denormalized.
#' @param ref A reference matrix used to obtain the original column-wise minima
#'   and maxima. Must have the same number of columns as \code{x}.
#'
#' @return A numeric matrix with values transformed back to the original scale.
#' @seealso \code{\link{nmfkc.normalize}}
#' @export
#' @examples
#' x <- nmfkc.normalize(iris[, -5])
#' x_recovered <- nmfkc.denormalize(x, iris[, -5])
#' apply(x_recovered - iris[, -5], 2, max)
nmfkc.denormalize <- function(x, ref=x) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    ref <- matrix(ref, ncol = 1)
  }
  r <- apply(ref, 2, range)
  y <- sweep(x, 2, r[2, ] - r[1, ], FUN = "*")
  y <- sweep(y, 2, r[1, ], FUN = "+")
  return(y)
}




#' @title Create a class (one-hot) matrix from a categorical vector
#' @description
#' \code{nmfkc.class} converts a categorical or factor vector into a class matrix
#' (one-hot encoded representation), where each row corresponds to a category
#' and each column corresponds to an observation.
#'
#' @param x A categorical vector or a factor.
#'
#' @return A binary matrix with one row per unique category and one column per observation. Each column has exactly one entry equal to 1, indicating the category of the observation.
#' @seealso \code{\link{nmfkc}}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- nmfkc.class(iris$Species)
#' Y[,1:6]
nmfkc.class <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  lev <- levels(x)
  X <- outer(lev, x, "==")
  mode(X) <- "numeric"
  rownames(X) <- lev
  if(!is.null(names(x))) colnames(X) <- names(x) else colnames(X) <- 1:length(x)
  X
}





#' @title Prediction method for objects of class \code{nmfkc}
#' @description
#' \code{predict.nmfkc} generates predictions from an object of class \code{nmfkc},
#' either using the fitted covariates or a new covariate matrix.
#'
#' @param object An object of class \code{nmfkc}, i.e., the return value of \code{nmfkc}.
#' @param newA Optional. A new covariate matrix to be used for prediction.
#' @param type Type of prediction to return. Options are "response", "prob", "class".
#' @param ... Further arguments passed to or from other methods.
#' @export
predict.nmfkc <- function(object, newA = NULL, type = "response", ...) {
  x <- object
  .eps <- 1e-10

  if(is.null(newA)){
    if(type=="response"){
      result <- x$X %*% x$B
    }else{
      XB.prob <- x$X %*% x$B.prob
      if(type=="prob"){
        result <- XB.prob
      }else{
        result <- rownames(x$X)[apply(XB.prob,2,which.max)]
      }
    }
  }else{
    B <- x$C %*% newA
    if(type=="response"){
      result <- x$X %*% B
    }else{
      B.prob <- sweep(B, 2, colSums(B) + .eps, "/")
      XB.prob <- x$X %*% B.prob
      if(type=="prob"){
        result <- XB.prob
      }else{
        result <- rownames(x$X)[apply(XB.prob,2,which.max)]
      }
    }
  }
  return(result)
}








#' @title Perform k-fold cross-validation for NMF with kernel covariates
#'
#' @description
#' \code{nmfkc.cv} performs k-fold cross-validation for the tri-factorization model
#' \eqn{Y \approx X C A = X B}, where
#' \itemize{
#'   \item \eqn{Y(P,N)} is the observation matrix,
#'   \item \eqn{A(R,N)} is the covariate (or kernel) matrix,
#'   \item \eqn{X(P,Q)} is the basis matrix (with \eqn{Q \le \min(P,N)}),
#'   \item \eqn{C(Q,R)} is the parameter matrix, and
#'   \item \eqn{B(Q,N)} is the coefficient matrix (\eqn{B = C A}).
#' }
#' Given \eqn{Y} (and optionally \eqn{A}), \eqn{X} and \eqn{C} are fitted on each
#' training split and predictive performance is evaluated on the held-out split.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#' @param Q Rank of the basis matrix \eqn{X}; must satisfy \eqn{Q \le \min(P,N)}.
#' @param ... Additional arguments controlling CV and the internal \code{\link{nmfkc}} call:
#'   \describe{
#'     \item{\code{Y.weights}}{Optional numeric matrix or vector; 0 indicates missing/ignored values.}
#'     \item{\code{div}}{Number of folds (\eqn{k}); default: \code{5}.}
#'     \item{\code{seed}}{Integer seed for reproducible partitioning; default: \code{123}.}
#'     \item{\code{shuffle}}{Logical. If \code{TRUE} (default), randomly shuffles samples (standard CV);
#'       if \code{FALSE}, splits sequentially (block CV; recommended for time series).}
#'     \item{\emph{Arguments passed to} \code{\link{nmfkc}}}{e.g., \code{gamma} (\code{B.L1}), \code{epsilon},
#'       \code{maxit}, \code{method} (\code{"EU"} or \code{"KL"}), \code{X.restriction}, \code{X.init}, etc.}
#'   }
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{objfunc}}{Mean loss per valid entry over all folds (MSE for \code{method="EU"}).}
#'     \item{\code{sigma}}{Residual standard error (RMSE). Available only if \code{method="EU"}; on the same scale as \code{Y}.}
#'     \item{\code{objfunc.block}}{Loss for each fold.}
#'     \item{\code{block}}{Vector of fold indices (1, …, \code{div}) assigned to each column of \eqn{Y}.}
#'   }
#'
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.kernel.beta.cv}}, \code{\link{nmfkc.ar.degree.cv}}
#'
#' @examples
#' # Example 1 (with explicit covariates):
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' res <- nmfkc.cv(Y, A, Q = 1)
#' res$objfunc
#'
#' # Example 2 (kernel A and beta sweep):
#' Y <- matrix(cars$dist, nrow = 1)
#' U <- matrix(c(5, 10, 15, 20, 25), nrow = 1)
#' V <- matrix(cars$speed, nrow = 1)
#' betas <- 25:35/1000
#' obj <- numeric(length(betas))
#' for (i in seq_along(betas)) {
#'   A <- nmfkc.kernel(U, V, beta = betas[i])
#'   obj[i] <- nmfkc.cv(Y, A, Q = 1, div = 10)$objfunc
#' }
#' betas[which.min(obj)]
#'
#' @export

nmfkc.cv <- function(Y, A=NULL, Q=2, ...){
  # A small constant for numerical stability to prevent division by zero and log(0).
  .eps <- 1e-10

  extra_args <- list(...)

  div <- if (!is.null(extra_args$div)) extra_args$div else 5
  seed <- if (!is.null(extra_args$seed)) extra_args$seed else 123
  shuffle <- if (!is.null(extra_args$shuffle)) extra_args$shuffle else TRUE

  gamma <- if (!is.null(extra_args$B.L1)) extra_args$B.L1 else if (!is.null(extra_args$gamma)) extra_args$gamma else 0
  epsilon <- if (!is.null(extra_args$epsilon)) extra_args$epsilon else 1e-4
  maxit   <- if (!is.null(extra_args$maxit))   extra_args$maxit   else 5000
  method  <- if (!is.null(extra_args$method))  extra_args$method  else "EU"
  Y.weights <- if (!is.null(extra_args$Y.weights)) extra_args$Y.weights else NULL

  if(!is.matrix(Y)) Y <- as.matrix(Y)
  P <- nrow(Y)
  N <- ncol(Y)

  # === Weights Preparation (Same as nmfkc) ===
  # Expand vector to matrix for easier slicing
  if (!is.null(Y.weights) && is.vector(Y.weights)) {
    if (length(Y.weights) == N) {
      Y.weights <- matrix(Y.weights, nrow = P, ncol = N, byrow = TRUE)
    } else if (length(Y.weights) == 1) {
      Y.weights <- matrix(Y.weights, nrow = P, ncol = N)
    } else {
      stop("Length of Y.weights vector must match ncol(Y) (or be 1).")
    }
  }
  # Handle NAs
  if (is.null(Y.weights)) {
    if (any(is.na(Y))) {
      Y.weights <- matrix(1, nrow=P, ncol=N)
      Y.weights[is.na(Y)] <- 0
      Y[is.na(Y)] <- 0
    } else {
      Y.weights <- matrix(1, nrow=P, ncol=N) # Full 1s for easier slicing logic
    }
  } else {
    if (!is.matrix(Y.weights)) Y.weights <- as.matrix(Y.weights)
    Y.weights[is.na(Y.weights)] <- 0
    Y[is.na(Y)] <- 0
    Y <- Y * Y.weights
  }

  # --- Helper: Weighted Optimization of B given fixed X and Y_test ---
  optimize.B.from.Y <- function(result, Y_test, W_test, gamma, epsilon, maxit, method){
    X <- result$X
    # Initialize C (which acts as B here)
    C <- matrix(1, nrow=ncol(X), ncol=ncol(Y_test))

    # Precompute ones for penalty (must match test dimension)
    ones_QN <- matrix(1, nrow=ncol(X), ncol=ncol(Y_test))

    oldSum <- 0
    epsilon.iter <- Inf

    for(l in 1:maxit){
      B <- C
      XB <- X %*% B

      if(method=="EU"){
        # Weighted Update for B (EU)
        # Num: X^T (W * Y)
        num <- t(X) %*% (W_test * Y_test)
        # Den: X^T (W * XB)
        den <- t(X) %*% (W_test * XB)

        if(gamma != 0) den <- den + (gamma/2) * ones_QN
        C <- C * ( num / (den + .eps) )

      }else{ # KL
        # Weighted Update for B (KL)
        # Num: X^T (W * (Y/XB))
        ratio <- W_test * (Y_test / (XB + .eps))
        num <- t(X) %*% ratio
        # Den: X^T W
        den <- t(X) %*% W_test

        if(gamma != 0) den <- den + gamma * ones_QN
        C <- C * ( num / (den + .eps) )
      }

      newSum <- sum(C)
      if(l>=10){
        #epsilon.iter <- abs(newSum-oldSum)/(abs(newSum)+0.1)
        epsilon.iter <- abs(newSum-oldSum) / pmax(abs(newSum), 1)
        if(epsilon.iter <= abs(epsilon)) break
      }
      oldSum <- sum(C)
    }

    B <- C
    XB <- X %*% B
    # Note: colnames handling omitted for speed in CV
    return(list(B=B, XB=XB))
  }

  is.identity.matrix <- function(A, tol = .Machine$double.eps) {
    if (nrow(A) != ncol(A)) return(FALSE)
    isTRUE(all.equal(A, diag(nrow(A)), tolerance = tol))
  }

  if(is.null(A)){
    is_identity <- TRUE
    is_symmetric.matrix <- FALSE
  }else{
    is_identity <- is.identity.matrix(A)
    is_symmetric.matrix <- isSymmetric(A, tol=.Machine$double.eps)
    A.function <- attr(A, "function.name")
    is_kernel_matrix <- !is.null(A.function) && A.function == "nmfkc.kernel"
  }

  # Create Folds
  remainder <- N %% div
  division <- N %/% div
  block <- 0*(1:N)

  if(shuffle){
    set.seed(seed)
    perm_index <- sample(1:N, N, replace=FALSE)
  } else {
    perm_index <- 1:N
  }

  processed_count <- 0
  for(i in 1:(div-1)){
    plus <- ifelse(i <= remainder, 1, 0)
    chunk_size <- division + plus
    target_indices <- perm_index[(processed_count + 1):(processed_count + chunk_size)]
    block[target_indices] <- i
    processed_count <- processed_count + chunk_size
  }
  target_indices <- perm_index[(processed_count + 1):N]
  block[target_indices] <- div

  objfunc.block <- numeric(div)
  total_valid_obs <- 0 # Denominator for weighted mean error

  for(j in 1:div){
    # Slice Y and Weights
    # Train
    Y_train <- Y[,block!=j, drop=FALSE]
    W_train <- Y.weights[,block!=j, drop=FALSE]

    # Test
    Y_test <- Y[,block==j, drop=FALSE]
    W_test <- Y.weights[,block==j, drop=FALSE]

    # Slice A
    if(is_identity){
      A_train <- NULL
    }else{
      if(is_symmetric.matrix && is_kernel_matrix){
        A_train <- A[block!=j,block!=j, drop=FALSE]
        A_test  <- A[block!=j,block==j, drop=FALSE]
      }else{
        A_train <- A[,block!=j, drop=FALSE]
        A_test  <- A[,block==j, drop=FALSE]
      }
    }

    # Run NMF on Training set (passing W_train)
    nmfkc_args <- c(
      list(...),
      list(Y = Y_train, A = A_train, Q = Q, Y.weights = W_train, # Pass weights!
           seed=NULL, print.trace = FALSE, print.dims = FALSE,
           save.time = TRUE, save.memory = TRUE)
    )
    nmfkc_args$shuffle <- NULL

    # Suppress messages from inner nmfkc calls
    res_j <- suppressMessages(do.call("nmfkc", nmfkc_args))

    # Predict on Test set
    if(is_identity){
      # Standard NMF: Optimize B for test set using weights
      resj <- optimize.B.from.Y(res_j, Y_test, W_test, gamma, epsilon, maxit, method)
      XB_test <- resj$XB
    }else{
      # Covariate NMF: Predict using A_test
      # XB = X * C * A_test
      XB_test <- res_j$X %*% res_j$C %*% A_test
    }

    # Evaluate Error (Weighted)
    if(method=="EU"){
      resid <- W_test * (Y_test - XB_test)
      objfunc.block[j] <- sum(resid^2)
    }else{
      term1 <- - (W_test * Y_test) * log(XB_test + .eps)
      term2 <- W_test * XB_test
      objfunc.block[j] <- sum(term1 + term2)
    }

    total_valid_obs <- total_valid_obs + sum(W_test > 0)
  }

  # Mean error per valid observation
  # (Avoid division by zero if all weights are 0, though unlikely)
  objfunc <- sum(objfunc.block) / max(total_valid_obs, 1)

  # Calculate RMSE for EU (This corresponds to 'sigma')
  sigma <- if(method == "EU") sqrt(objfunc) else NA

  return(list(objfunc=objfunc, sigma=sigma, objfunc.block=objfunc.block, block=block))
}





#' @title Perform Element-wise Cross-Validation (Wold's CV)
#' @description
#' \code{nmfkc.ecv} performs k-fold cross-validation by randomly holding out
#' individual elements of the data matrix (element-wise), assigning them a
#' weight of 0 via `Y.weights`, and evaluating the reconstruction error on
#' those held-out elements.
#'
#' This method (also known as Wold's CV) is theoretically robust for determining
#' the optimal rank (Q) in NMF. This function supports vector input for `Q`,
#' allowing simultaneous evaluation of multiple ranks on the same folds.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix.
#' @param Q Vector of ranks to evaluate (e.g., 1:5).
#' @param div Number of folds (default: 5).
#' @param seed Integer seed for reproducibility.
#' @param ... Additional arguments passed to \code{\link{nmfkc}} (e.g., \code{method="EU"}).
#'
#' @return A list with components:
#' \item{objfunc}{Numeric vector containing the Mean Squared Error (MSE) for each Q.}
#' \item{sigma}{Numeric vector containing the Residual Standard Error (RMSE) for each Q. Only available if method="EU".}
#' \item{objfunc.fold}{List of length equal to Q vector. Each element contains the MSE values for the k folds.}
#' \item{folds}{A list of length \code{div}, containing the linear indices of held-out elements for each fold (shared across all Q).}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.cv}}
#' @export
nmfkc.ecv <- function(Y, A=NULL, Q=1:3, div=5, seed=123, ...){

  if(!is.matrix(Y)) Y <- as.matrix(Y)
  P <- nrow(Y)
  N <- ncol(Y)

  # --- Argument Handling ---
  extra_args <- list(...)

  # If user mistakenly passed 'rank' in ..., treat it as Q
  if (!is.null(extra_args$rank)) {
    Q <- extra_args$rank
  }

  # 1. Create Folds
  if (!is.null(seed)) set.seed(seed)
  valid_indices <- which(!is.na(Y))
  n_valid <- length(valid_indices)

  perm_indices <- sample(valid_indices)
  folds <- vector("list", div)
  chunk_size <- n_valid %/% div
  remainder <- n_valid %% div

  start_idx <- 1
  for(k in 1:div){
    current_size <- chunk_size + ifelse(k <= remainder, 1, 0)
    end_idx <- start_idx + current_size - 1
    folds[[k]] <- perm_indices[start_idx:end_idx]
    start_idx <- end_idx + 1
  }

  method <- if(!is.null(extra_args$method)) extra_args$method else "EU"

  # 2. Loop over Q vectors
  num_Q <- length(Q)
  result_objfunc <- numeric(num_Q)
  result_sigma   <- numeric(num_Q)
  result_fold    <- vector("list", num_Q)

  names(result_objfunc) <- paste0("Q=", Q)
  names(result_sigma)   <- paste0("Q=", Q)
  names(result_fold)    <- paste0("Q=", Q)

  message(paste0("Performing Element-wise CV for Q = ", paste(Q, collapse=","), " (", div, "-fold)..."))

  for(i in 1:num_Q){
    q_curr <- Q[i]
    objfunc.fold <- numeric(div)

    for(k in 1:div){
      test_idx <- folds[[k]]
      weights_train <- matrix(1, nrow=P, ncol=N)
      if(any(is.na(Y))) weights_train[is.na(Y)] <- 0
      weights_train[test_idx] <- 0

      # Prepare arguments for nmfkc
      # Remove 'Q' and 'rank' from extra_args to avoid conflicts
      nmfkc_clean_args <- extra_args
      nmfkc_clean_args$Q <- NULL
      nmfkc_clean_args$rank <- NULL

      nmfkc_args <- c(list(Y=Y, A=A, Q=q_curr, Y.weights=weights_train), nmfkc_clean_args)

      fit <- do.call("nmfkc", nmfkc_args)

      pred <- fit$XB
      residuals <- Y[test_idx] - pred[test_idx]
      objfunc.fold[k] <- mean(residuals^2)
    }

    result_fold[[i]]  <- objfunc.fold
    result_objfunc[i] <- mean(objfunc.fold)
    result_sigma[i]   <- if(method == "EU") sqrt(result_objfunc[i]) else NA
  }

  return(list(objfunc = result_objfunc,
              sigma = result_sigma,
              objfunc.fold = result_fold,
              folds = folds))
}



#' @title Rank selection diagnostics with graphical output
#' @description
#' \code{nmfkc.rank} provides diagnostic criteria for selecting the rank (\eqn{Q})
#' in NMF with kernel covariates. Several model selection measures are computed
#' (e.g., R-squared, silhouette, CPCC, ARI), and results can be visualized in a plot.
#'
#' By default (\code{save.time = FALSE}), this function also computes the
#' Element-wise Cross-Validation error (Wold's CV Sigma) using \code{\link{nmfkc.ecv}}.
#'
#' The plot explicitly marks the "BEST" rank based on two criteria:
#' \enumerate{
#'   \item **Elbow Method (Red)**: Based on the curvature of the R-squared values (always computed if Q > 2).
#'   \item **Min RMSE (Blue)**: Based on the minimum Element-wise CV Sigma (only if \code{save.time=FALSE}).
#' }
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#' @param rank A vector of candidate ranks to be evaluated.
#' @param save.time Logical. If \code{TRUE}, skips heavy computations like Element-wise CV.
#'   Default is \code{FALSE} (computes everything).
#' @param plot Logical. If \code{TRUE} (default), draws a plot of the diagnostic criteria.
#' @param ... Additional arguments passed to \code{\link{nmfkc}} and \code{\link{nmfkc.ecv}}.
#'   \itemize{
#'     \item \code{Q}: (Deprecated) Alias for \code{rank}.
#'   }
#'
#' @return A list containing:
#' \item{rank.best}{The estimated optimal rank. Prioritizes ECV minimum if available, otherwise R-squared Elbow.}
#' \item{criteria}{A data frame containing diagnostic metrics for each rank.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ecv}}
#' @export
#' @references
#' Brunet, J.P., Tamayo, P., Golub, T.R., Mesirov, J.P. (2004).
#' Metagenes and molecular pattern discovery using matrix factorization.
#' \emph{Proc. Natl. Acad. Sci. USA}, 101, 4164–4169.
#' \doi{10.1073/pnas.0308531101}
#' Punera, K., & Ghosh, J. (2008).
#' Consensus-based ensembles of soft clusterings.
#' \emph{Applied Artificial Intelligence}, 22(7–8), 780–810.
#' \doi{10.1080/08839510802170546}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' library(nmfkc)
#' Y <- t(iris[,-5])
#' # Full run (default)
#' nmfkc.rank(Y, rank=1:4)
#' # Fast run (skip ECV)
#' nmfkc.rank(Y, rank=1:4, save.time=TRUE)

nmfkc.rank <- function(Y, A=NULL, rank=1:2, save.time=FALSE, plot=TRUE, ...){

  extra_args <- list(...)

  # Backward Compatibility for Q
  if (!is.null(extra_args$Q)) rank <- extra_args$Q
  Q <- rank
  max_rank <- min(nrow(Y), ncol(Y))
  if (any(Q > max_rank)) {
    invalid_Q <- Q[Q > max_rank]
    warning(paste0("Rank(s) ", paste(invalid_Q, collapse=", "),
                   " exceed min(nrow(Y), ncol(Y)) = ", max_rank,
                   ". They have been removed."))
    Q <- Q[Q <= max_rank]
    if(length(Q) == 0) stop("No valid ranks specified (all exceed matrix dimensions).")
  }
  # ---------------------------------------------
  AdjustedRandIndex <- function(x){
    choose2 <- function(n) choose(n,2)
    a <- sum(apply(x,c(1,2),choose2))
    ab <- sum(sapply(rowSums(x),choose2))
    b <- ab-a
    ac <- sum(sapply(colSums(x),choose2))
    c <- ac-a
    total <- choose2(sum(x))
    d <- total-a-b-c
    (ri <- (a+d)/total)
    e <- ab*ac/total+(total-ab)*(total-ac)/total
    (ari <- (a+d-e)/(total-e))
    return(list(RI=ri,ARI=ari))
  }

  num_q <- length(Q)
  results_df <- data.frame(
    rank = Q,
    r.squared = numeric(num_q),
    ICp = numeric(num_q),
    AIC = numeric(num_q),
    BIC = numeric(num_q),
    B.prob.sd.min = numeric(num_q),
    B.prob.entropy.mean = numeric(num_q),
    B.prob.max.mean = numeric(num_q),
    ARI = numeric(num_q),
    silhouette = numeric(num_q),
    CPCC = numeric(num_q),
    dist.cor = numeric(num_q),
    sigma.ecv = numeric(num_q)
  )

  cluster.old <- NULL

  # --- Main Loop for Standard Metrics ---
  for(q_idx in 1:num_q){
    current_Q <- Q[q_idx]

    extra_args_nmfkc <- extra_args
    extra_args_nmfkc$save.memory <- FALSE
    extra_args_nmfkc$save.time <- save.time
    extra_args_nmfkc$Q <- NULL

    nmfkc_args <- c(list(Y = Y, A = A, rank = current_Q), extra_args_nmfkc)
    result <- do.call("nmfkc", nmfkc_args)

    results_df$r.squared[q_idx] <- result$r.squared
    results_df$ICp[q_idx] <- result$criterion$ICp
    results_df$AIC[q_idx] <- result$criterion$AIC
    results_df$BIC[q_idx] <- result$criterion$BIC
    results_df$B.prob.sd.min[q_idx] <- result$criterion$B.prob.sd.min
    results_df$B.prob.max.mean[q_idx] = result$criterion$B.prob.max.mean
    results_df$B.prob.entropy.mean[q_idx] = result$criterion$B.prob.entropy.mean

    if(save.time){
      results_df$CPCC[q_idx] <- NA
      results_df$dist.cor[q_idx] <- NA
      results_df$silhouette[q_idx] <- NA
    } else {
      results_df$CPCC[q_idx] <- result$criterion$CPCC
      results_df$dist.cor[q_idx] <- result$criterion$dist.cor
      results_df$silhouette[q_idx] <- result$criterion$silhouette$silhouette.mean
    }

    if(is.null(cluster.old)){
      results_df$ARI[q_idx] <- NA
    } else {
      df <- data.frame(old=cluster.old, new=result$B.cluster)
      df <- df[stats::complete.cases(df),]
      if (nrow(df) > 0 && length(unique(df$old)) > 1 && length(unique(df$new)) > 1) {
        f <- table(df$old, df$new)
        results_df$ARI[q_idx] <- AdjustedRandIndex(f)$ARI
      } else {
        results_df$ARI[q_idx] <- NA
      }
    }
    cluster.old <- result$B.cluster
  }

  # --- Element-wise CV (Wold's CV) ---
  rank.best.ecv <- NA
  if(!save.time){
    ecv_args <- list(Y = Y, A = A, rank = Q)
    extra_args_ecv <- extra_args
    extra_args_ecv$save.time <- NULL
    extra_args_ecv$save.memory <- NULL
    extra_args_ecv$Q <- NULL

    ecv_full_args <- c(ecv_args, extra_args_ecv)

    message("Running Element-wise CV (this may take time)...")
    ecv_res <- do.call("nmfkc.ecv", ecv_full_args)
    results_df$sigma.ecv <- ecv_res$sigma

    # Determine best rank by ECV
    idx_best_ecv <- which.min(results_df$sigma.ecv)
    rank.best.ecv <- results_df$rank[idx_best_ecv]
  } else {
    results_df$sigma.ecv <- NA
  }

  # --- Determine R-squared Best Rank (Elbow) ---
  rank.best.r2 <- NA
  if(num_q > 2){
    x <- 1:num_q
    y <- results_df$r.squared
    y_norm <- (y - min(y)) / (max(y) - min(y))
    x_norm <- (x - min(x)) / (max(x) - min(x))
    x1 <- x_norm[1]; y1 <- y_norm[1]
    x2 <- x_norm[num_q]; y2 <- y_norm[num_q]

    distances <- numeric(num_q)
    for(i in 1:num_q){
      distances[i] <- abs((y2-y1)*x_norm[i] - (x2-x1)*y_norm[i] + x2*y1 - y2*x1) / sqrt((y2-y1)^2 + (x2-x1)^2)
    }
    idx_best_r2 <- which.max(distances)
    rank.best.r2 <- results_df$rank[idx_best_r2]
  }

  # Decide final recommended rank (Prioritize ECV)
  rank.final <- if(!is.na(rank.best.ecv)) rank.best.ecv else rank.best.r2

  # --- Plotting ---
  if(plot){
    old_par <- graphics::par(mar = c(5, 4, 4, 5) + 0.1)
    on.exit(graphics::par(old_par))

    # 1. Left Axis Plot: 0-1 Metrics
    plot(results_df$rank, results_df$r.squared, type="l", col=2, lwd=3,
         xlab="Rank (Q)", ylab="Fit / Stability (0-1)", ylim=c(0,1),
         main="Rank Selection Diagnostics")
    for(q in results_df$rank) graphics::abline(v=q,col="gray90",lwd=0.5)
    graphics::points(results_df$rank, results_df$r.squared, pch=16, col=2, cex=0.8)
    graphics::text(results_df$rank, results_df$r.squared, results_df$rank, pos=3, col=2, cex=0.8)

    # Always mark R2 Elbow BEST if found
    if(!is.na(rank.best.r2)){
      idx_r2 <- which(results_df$rank == rank.best.r2)
      graphics::points(rank.best.r2, results_df$r.squared[idx_r2], pch=16, col="red", cex=1.5)
      graphics::text(rank.best.r2, results_df$r.squared[idx_r2], "Best (Elbow)", pos=1, col="red", cex=0.8)
    }

    legend_txt <- c("r.squared")
    legend_col <- c(2)
    legend_lty <- c(1)

    graphics::lines(results_df$rank, results_df$B.prob.sd.min, col="green1", lwd=2)
    legend_txt <- c(legend_txt, "B.prob.sd.min"); legend_col <- c(legend_col, "green1"); legend_lty <- c(legend_lty, 1)

    graphics::lines(results_df$rank, results_df$B.prob.max.mean, col="green3", lwd=2)
    legend_txt <- c(legend_txt, "B.prob.max.mean"); legend_col <- c(legend_col, "green3"); legend_lty <- c(legend_lty, 1)

        graphics::lines(results_df$rank, results_df$B.prob.entropy.mean, col="green4", lwd=2)
    legend_txt <- c(legend_txt, "B.prob.entropy.mean"); legend_col <- c(legend_col, "green4"); legend_lty <- c(legend_lty, 1)


    if (any(!is.na(results_df$ARI))) {
      graphics::lines(results_df$rank, results_df$ARI, col=4, lwd=2)
      legend_txt <- c(legend_txt, "ARI"); legend_col <- c(legend_col, 4); legend_lty <- c(legend_lty, 1)
    }
    if (any(!is.na(results_df$silhouette))) {
      graphics::lines(results_df$rank, results_df$silhouette, col=7, lwd=2)
      legend_txt <- c(legend_txt, "Silhouette"); legend_col <- c(legend_col, 7); legend_lty <- c(legend_lty, 1)
    }
    if (any(!is.na(results_df$CPCC))) {
      graphics::lines(results_df$rank, results_df$CPCC, col=6, lwd=2)
      legend_txt <- c(legend_txt, "CPCC"); legend_col <- c(legend_col, 6); legend_lty <- c(legend_lty, 1)
    }
    if (any(!is.na(results_df$dist.cor))) {
      graphics::lines(results_df$rank, results_df$dist.cor, col=5, lwd=2)
      legend_txt <- c(legend_txt, "dist.cor"); legend_col <- c(legend_col, 5); legend_lty <- c(legend_lty, 1)
    }

    # 2. Right Axis Plot: Sigma ECV (RMSE)
    if (any(!is.na(results_df$sigma.ecv))) {
      graphics::par(new = TRUE)
      graphics::plot(results_df$rank, results_df$sigma.ecv, type="l", col="blue", lwd=3,
                     axes=FALSE, xlab="", ylab="")

      graphics::points(results_df$rank, results_df$sigma.ecv, pch=16, col="blue", cex=0.8)
      graphics::text(results_df$rank, results_df$sigma.ecv, results_df$rank, pos=3, col="blue", cex=0.8)

      # Always mark ECV Min BEST
      if(!is.na(rank.best.ecv)){
        idx_ecv <- which(results_df$rank == rank.best.ecv)
        graphics::points(results_df$rank, results_df$sigma.ecv, pch=1, col="blue")
        graphics::points(rank.best.ecv, results_df$sigma.ecv[idx_ecv], pch=16, col="blue", cex=1.5)
        graphics::text(rank.best.ecv, results_df$sigma.ecv[idx_ecv], "Best (Min)", pos=1, col="blue", cex=0.8)
      }

      graphics::axis(side=4, col="blue", col.axis="blue")
      graphics::mtext("ECV Sigma (RMSE)", side=4, line=3, col="blue")

      legend_txt <- c(legend_txt, "ECV Sigma")
      legend_col <- c(legend_col, "blue")
      legend_lty <- c(legend_lty, 1)
    }

    graphics::legend("right", legend=legend_txt, col=legend_col, lty=legend_lty, lwd=2, bg="white", cex=0.7)
  }

  return(list(rank.best = rank.final, criteria = results_df))
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



#' @title Plot Diagnostics: Original, Fitted, and Residual Matrices as Heatmaps
#' @description
#' This function generates a side-by-side plot of three heatmaps: the original
#' observation matrix Y, the fitted matrix XB (from NMF), and the residual matrix E (Y - XB).
#' This visualization aids in diagnosing whether the chosen rank Q is adequate
#' by assessing if the residual matrix E appears to be random noise.
#'
#' The axis labels (X-axis: Samples, Y-axis: Features) are integrated into the main title of each plot
#' to maximize the plot area, reflecting the compact layout settings.
#'
#' @param Y The original observation matrix (P x N).
#' @param result The result object returned by the nmfkc function.
#' @param Y_XB_palette A vector of colors used for Y and XB heatmaps. Defaults to a white-orange-red gradient.
#' @param E_palette A vector of colors used for the residuals (E) heatmap. Defaults to a blue-white-red gradient.
#' @param ... Additional graphical parameters passed to the internal image calls.
#'
#' @return NULL. The function generates a plot.
#' @export
nmfkc.residual.plot <- function(Y, result,
                                Y_XB_palette = grDevices::colorRampPalette(c("white", "orange", "red"))(256),
                                E_palette = grDevices::colorRampPalette(c("blue", "white", "red"))(256), ...){
  if (!inherits(result, "nmfkc")) {
    stop("The 'result' argument must be an object of class 'nmfkc'.")
  }
  XB <- result$XB
  E <- Y - XB
  if(nrow(Y) != nrow(E) || ncol(Y) != ncol(E)){
    stop("Dimension mismatch between Y and result$XB. Check input matrices.")
  }
  old_par <- graphics::par(mfrow = c(1, 3), mar = c(1,0.5,5,0.5) + 0.1)
  on.exit(graphics::par(old_par))
  min_YX <- min(Y, XB, na.rm = TRUE)
  max_YX <- max(Y, XB, na.rm = TRUE)
  max_abs_E <- max(abs(E), na.rm = TRUE)
  min_E <- -max_abs_E
  max_E <- max_abs_E
  graphics::image(t(Y)[, nrow(Y):1],
                  col = Y_XB_palette,
                  zlim = c(min_YX, max_YX),
                  main = "1. Original Matrix Y\n X-axis: Samples (N), Y-axis: Features (P)",
                  xlab = "",
                  ylab = "",
                  axes = FALSE, ...)
  graphics::box()
  graphics::image(t(XB)[, nrow(XB):1],
                  col = Y_XB_palette,
                  zlim = c(min_YX, max_YX),
                  main = paste0("2. Fitted Matrix XB (Q=",ncol(result$X),")\n X-axis: Samples (N), Y-axis: Features (P)"),
                  xlab = "",
                  ylab = "",
                  axes = FALSE, ...)
  graphics::box()
  graphics::image(t(E)[, nrow(E):1],
                  col = E_palette,
                  zlim = c(min_E, max_E),
                  main = "3. Residual Matrix E (Y - XB)\n X-axis: Samples (N), Y-axis: Features (P)",
                  xlab = "",
                  ylab = "",
                  axes = FALSE, ...)
  graphics::box()
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
#' @title NMF-SEM Main Estimation Algorithm
#'
#' @description
#' Fits the NMF-SEM model
#' \deqn{
#'   Y_1 \approx X \bigl( \Theta_1 Y_1 + \Theta_2 Y_2 \bigr)
#' }
#' under non-negativity constraints with orthogonality and sparsity regularization.
#' The function returns the estimated latent factors, structural coefficient matrices,
#' and the implied equilibrium (input–output) mapping.
#'
#' At equilibrium, the model can be written as
#' \deqn{
#'   Y_1 \approx (I - X \Theta_1)^{-1} X \Theta_2 Y_2
#'   \equiv M_{\mathrm{model}} Y_2,
#' }
#' where \eqn{M_{\mathrm{model}} = (I - X \Theta_1)^{-1} X \Theta_2} is a
#' Leontief-type cumulative-effect operator in latent space.
#'
#' Internally, the latent feedback and exogenous loading matrices are stored as
#' \code{C1} and \code{C2}, corresponding to \eqn{\Theta_1} and \eqn{\Theta_2},
#' respectively.
#'
#' @param Y1 A non-negative numeric matrix of endogenous variables with
#'   \strong{rows = variables (P1), columns = samples (N)}.
#' @param Y2 A non-negative numeric matrix of exogenous variables with
#'   \strong{rows = variables (P2), columns = samples (N)}.
#'   Must satisfy \code{ncol(Y1) == ncol(Y2)}.
#' @param rank Integer; number of latent factors \eqn{Q}. If \code{NULL},
#'   \eqn{Q} is taken from a hidden argument in \code{...} or defaults to
#'   \code{nrow(Y2)}.
#' @param X.init Optional non-negative initialization for the basis matrix
#'   \code{X} (\eqn{P_1 \times Q}). If supplied, it is projected to be
#'   non-negative and column-normalized.
#' @param X.L2.ortho L2 orthogonality penalty for \code{X}. This controls
#'   the penalty term \eqn{\lambda_X \lVert X^\top X - \mathrm{diag}(X^\top X)
#'   \rVert_F^2}. Default: \code{100}.
#' @param C1.L1 L1 sparsity penalty for \code{C1} (i.e., \eqn{\Theta_1}).
#'   Default: \code{1.0}.
#' @param C2.L1 L1 sparsity penalty for \code{C2} (i.e., \eqn{\Theta_2}).
#'   Default: \code{0.1}.
#' @param epsilon Relative convergence threshold for the objective function.
#'   Iterations stop when the relative change in reconstruction loss falls
#'   below this value. Default: \code{1e-6}.
#' @param maxit Maximum number of iterations for the multiplicative updates.
#'   Default: \code{20000}.
#' @param seed Random seed used to initialize \code{X}, \code{C1}, and \code{C2}.
#'   Default: \code{123}.
#' @param ... Additional arguments. Currently used to pass a hidden rank
#'   \code{Q} (e.g., via \code{Q = 3}) if \code{rank} is \code{NULL}.
#'
#' @return A list with components:
#'   \item{X}{Estimated basis matrix (\eqn{P_1 \times Q}).}
#'   \item{C1}{Estimated latent feedback matrix (\eqn{\Theta_1}, \eqn{Q \times P_1}).}
#'   \item{C2}{Estimated exogenous loading matrix (\eqn{\Theta_2}, \eqn{Q \times P_2}).}
#'   \item{XC1}{Feedback matrix \eqn{X \Theta_1}.}
#'   \item{XC2}{Direct-effect matrix \eqn{X \Theta_2}.}
#'   \item{XC1.radius}{Spectral radius \eqn{\rho(X \Theta_1)}.}
#'   \item{XC1.norm1}{Induced 1-norm \eqn{\lVert X \Theta_1 \rVert_{1,\mathrm{op}}}.}
#'   \item{Leontief.inv}{Leontief-type inverse \eqn{(I - X \Theta_1)^{-1}.}}
#'   \item{M.model}{Equilibrium mapping
#'     \eqn{M_{\mathrm{model}} = (I - X \Theta_1)^{-1} X \Theta_2}.}
#'   \item{amplification}{Latent amplification factor
#'     \eqn{\lVert M_{\mathrm{model}} \rVert_{1,\mathrm{op}} /
#'          \bigl\lVert X \Theta_2 \bigr\rVert_{1,\mathrm{op}}}.}
#'   \item{amplification.bound}{Geometric-series upper bound
#'     \eqn{1 / (1 - \lVert X \Theta_1 \rVert_{1,\mathrm{op}})} if
#'     \eqn{\lVert X \Theta_1 \rVert_{1,\mathrm{op}} < 1}, otherwise \code{Inf}.}
#'   \item{Q}{Effective latent dimension used in the fit.}
#'   \item{SC.cov}{Correlation between sample and model-implied covariance
#'     (flattened) of \eqn{Y_1}.}
#'   \item{MAE}{Mean absolute error between \eqn{Y_1} and its equilibrium
#'     prediction \eqn{\hat Y_1 = M_{\mathrm{model}} Y_2}.}
#'   \item{objfunc}{Vector of reconstruction losses per iteration.}
#'   \item{objfunc.full}{Vector of penalized objective values per iteration.}
#'   \item{iter}{Number of iterations actually performed.}
#'
#' @export
nmf.sem <- function(
    Y1, Y2,
    rank = NULL,
    X.init = NULL,
    X.L2.ortho = 100.0,
    C1.L1 = 1.0,
    C2.L1 = 0.1,
    epsilon = 1e-6,
    maxit = 20000,
    seed  = 123,
    ...
) {
  # ------------------------------ checks ------------------------------
  if (!is.matrix(Y1)) Y1 <- as.matrix(Y1)
  if (!is.matrix(Y2)) Y2 <- as.matrix(Y2)

  if (any(!is.finite(Y1)) || any(!is.finite(Y2)))
    stop("Y1 and Y2 must not contain NA/NaN/Inf.")
  if (min(Y1) < 0 || min(Y2) < 0)
    stop("Y1 and Y2 must be non-negative.")
  if (ncol(Y1) != ncol(Y2))
    stop("ncol(Y1) must be equal to ncol(Y2).")

  extra_args <- list(...)
  Q_hidden <- if (!is.null(extra_args$Q)) extra_args$Q else NULL
  Q0 <- if (!is.null(rank)) rank else if (!is.null(Q_hidden)) Q_hidden else nrow(Y2)

  P1 <- nrow(Y1); P2 <- nrow(Y2); N <- ncol(Y1)
  Q  <- min(Q0, P1, N)
  if (Q < 1) stop("Effective rank Q must be >= 1.")

  # -------------------------- labels (output) -------------------------
  Y1_labels    <- if (!is.null(rownames(Y1))) rownames(Y1) else paste0("Y1_", 1:P1)
  Y2_labels    <- if (!is.null(rownames(Y2))) rownames(Y2) else paste0("Y2_", 1:P2)
  Basis_labels <- paste0("Factor", 1:Q)

  set.seed(seed)
  .eps <- 1e-10
  .xnorm  <- function(X) sweep(X, 2, pmax(colSums(X), .eps), "/")
  mat1norm <- function(A) max(colSums(abs(A)))

  # ---------------------------- init X,C1,C2 --------------------------
  if (is.null(X.init)) {
    X <- .nndsvdar(Y1, Q)   # exists in nmfkc
  } else {
    X <- as.matrix(X.init)
    if (!all(dim(X) == c(P1, Q))) {
      stop("X.init must have dimension (nrow(Y1) x rank).")
    }
    X[X < 0] <- 0
  }
  X <- .xnorm(X)

  C1 <- matrix(stats::runif(Q * P1, 0.01, 0.1), nrow = Q, ncol = P1)
  C2 <- matrix(stats::runif(Q * P2, 0.001, 0.01), nrow = Q, ncol = P2)

  min_dim <- min(Q, P2)
  for (i in 1:min_dim) {
    C2[i, i] <- stats::runif(1, 0.1, 0.2)
  }

  objfunc      <- numeric(maxit)
  objfunc.full <- numeric(maxit)

  # ----------------------------- main loop ----------------------------
  for (it in 1:maxit) {
    M  <- C1 %*% Y1 + C2 %*% Y2
    Mt <- t(M)

    # 2.1 update X
    Numerator_X       <- Y1 %*% Mt
    Denominator_X_rec <- X %*% M %*% Mt

    if (X.L2.ortho > 0) {
      XtX <- t(X) %*% X
      XtX_offdiag <- XtX
      diag(XtX_offdiag) <- 0
      Denominator_X_ortho <- X.L2.ortho * X %*% XtX_offdiag
    } else {
      Denominator_X_ortho <- 0
    }

    X <- X * (Numerator_X / (Denominator_X_rec + Denominator_X_ortho + .eps))
    X <- .xnorm(X)

    Xt  <- t(X)
    XtX <- Xt %*% X

    # 2.2 update C1
    Numerator_C1   <- Xt %*% Y1 %*% t(Y1)
    Denominator_C1 <- XtX %*% (C1 %*% Y1 + C2 %*% Y2) %*% t(Y1) + C1.L1 + .eps
    C1 <- C1 * (Numerator_C1 / Denominator_C1)

    # 2.3 update C2
    Numerator_C2   <- Xt %*% Y1 %*% t(Y2)
    Denominator_C2 <- XtX %*% (C1 %*% Y1 + C2 %*% Y2) %*% t(Y2) + C2.L1 + .eps
    C2 <- C2 * (Numerator_C2 / Denominator_C2)

    # loss + penalties
    XB <- X %*% (C1 %*% Y1 + C2 %*% Y2)
    loss_rec <- sum((Y1 - XB)^2)
    objfunc[it] <- loss_rec

    if (X.L2.ortho > 0) {
      XtX_off <- XtX
      diag(XtX_off) <- 0
      pen_X_ortho <- 0.5 * X.L2.ortho * sum(XtX_off^2)
    } else {
      pen_X_ortho <- 0
    }
    pen_C1_L1 <- C1.L1 * sum(C1)
    pen_C2_L1 <- C2.L1 * sum(C2)
    objfunc.full[it] <- loss_rec + pen_X_ortho + pen_C1_L1 + pen_C2_L1

    if (it >= 10) {
      epsilon_iter <- abs(objfunc[it] - objfunc[it - 1]) / pmax(abs(objfunc[it]), 1)
      if (epsilon_iter <= epsilon) break
    }
  }

  # ------------------ reorder factors (nmfkc centroid order) ----------
  centroid <- as.numeric((1:nrow(X)) / nrow(X)) %*% X
  index <- order(centroid)
  X  <- X[, index, drop = FALSE]
  C1 <- C1[index, , drop = FALSE]
  C2 <- C2[index, , drop = FALSE]

  # ------------------------------ names -------------------------------
  colnames(X)  <- Basis_labels
  rownames(C1) <- Basis_labels
  rownames(C2) <- Basis_labels
  rownames(X)  <- Y1_labels
  colnames(C1) <- Y1_labels
  colnames(C2) <- Y2_labels

  # -------------------- feedback + stability diagnostics --------------
  XC1  <- X %*% C1
  eigs <- eigen(XC1, only.values = TRUE)$values
  rho  <- max(abs(eigs))
  if (rho >= 1)
    warning("Leontief.inv may be unstable; spectral radius >= 1.")

  XC1_norm1 <- mat1norm(XC1)

  # -------------------- Leontief inverse + equilibrium mapping --------
  I_mat <- diag(nrow(XC1))
  XC2   <- X %*% C2

  Leontief.inv <- tryCatch(
    base::solve(I_mat - XC1),
    error = function(e) {
      warning("Failed to compute Leontief.inv via solve(I - XC1). Returning NA matrices.")
      matrix(NA_real_, nrow = nrow(XC1), ncol = ncol(XC1))
    }
  )
  M.model <- Leontief.inv %*% XC2

  amplification <- mat1norm(M.model) / (mat1norm(XC2) + .eps)
  amplification.bound <- if (XC1_norm1 < 1) 1 / (1 - XC1_norm1) else Inf

  # -------------------- fit indices (equilibrium prediction) ----------
  Y1_hat   <- M.model %*% Y2
  S.sample <- Y1 %*% t(Y1)
  S.model  <- Y1_hat %*% t(Y1_hat)
  SC.cov   <- stats::cor(as.numeric(S.sample), as.numeric(S.model))
  MAE      <- mean(abs(Y1 - Y1_hat))

  list(
    X                   = X,
    C1                  = C1,
    C2                  = C2,
    XC1                 = XC1,
    XC2                 = XC2,
    XC1.radius          = rho,
    XC1.norm1           = XC1_norm1,
    Leontief.inv        = Leontief.inv,
    M.model             = M.model,
    amplification       = amplification,
    amplification.bound = amplification.bound,
    Q                   = Q,
    SC.cov              = SC.cov,
    MAE                 = MAE,
    objfunc             = objfunc[1:it],
    objfunc.full        = objfunc.full[1:it],
    iter                = it
  )
}



#' @title Cross-Validation for NMF-SEM
#' @description
#' Performs K-fold cross-validation to evaluate the equilibrium mapping of
#' the NMF-SEM model.
#'
#' For each fold, \code{nmf.sem} is fitted on the training samples,
#' yielding an equilibrium mapping \eqn{\hat Y_1 = M_{\mathrm{model}} Y_2}.
#' The held-out endogenous variables \eqn{Y_1} are then predicted from \eqn{Y_2}
#' using this mapping, and the mean absolute error (MAE) over all entries in the
#' test block is computed. The returned value is the average MAE across folds.
#'
#' This implements the hyperparameter selection strategy described in the paper:
#' hyperparameters are chosen by predictive cross-validation rather than direct
#' inspection of the internal structural matrices.
#'
#' @param Y1 A non-negative numeric matrix of endogenous variables with
#'   \strong{rows = variables (P1), columns = samples (N)}.
#' @param Y2 A non-negative numeric matrix of exogenous variables with
#'   \strong{rows = variables (P2), columns = samples (N)}.
#'   Must satisfy \code{ncol(Y1) == ncol(Y2)}.
#' @param rank Integer; rank (number of latent factors) passed to \code{nmf.sem}.
#'   If \code{NULL}, \code{nmf.sem} decides the effective rank (via \code{...} or \code{nrow(Y2)}).
#' @param X.init Optional initialization for \code{X} (as in \code{nmf.sem}).
#' @param X.L2.ortho L2 orthogonality penalty for \code{X}.
#' @param C1.L1 L1 sparsity penalty for \code{C1} (\eqn{\Theta_1}).
#' @param C2.L1 L1 sparsity penalty for \code{C2} (\eqn{\Theta_2}).
#' @param epsilon Convergence threshold for \code{nmf.sem}.
#' @param maxit Maximum number of iterations for \code{nmf.sem}.
#' @param seed Master random seed for CV splitting and fold-specific calls to \code{nmf.sem}.
#'   If \code{NULL}, RNG is not controlled within folds.
#' @param div Number of CV folds. (Default: \code{5})
#' @param shuffle Logical; if \code{TRUE}, samples are randomly permuted
#'   before assigning to folds. (Default: \code{TRUE})
#' @param ... Additional arguments passed to \code{nmf.sem} (except for
#'   \code{rank}, \code{seed}, \code{div}, \code{shuffle}, which are handled here).
#'
#' @return A numeric scalar: mean MAE across CV folds.
#'
#' @export
nmf.sem.cv <- function(
    Y1, Y2,
    rank = NULL,
    X.init = NULL,
    X.L2.ortho = 100.0,
    C1.L1 = 0.5,        # L1 sparsity for C1 (Theta1)
    C2.L1 = 0.0,        # L1 sparsity for C2 (Theta2)
    epsilon = 1e-4,     # Convergence tolerance passed to nmf.sem
    maxit = 50000,
    seed = NULL,        # Master seed for CV (partition + fold seeds)
    div = 5,            # Number of CV folds
    shuffle = TRUE,     # Shuffle samples before assigning folds
    ...
){
  # ------------------------------------------------------------------
  # 1. Basic input checks
  #
  # NMF-SEM requires non-negative matrices. We also require that Y1 and Y2
  # share the same number of samples (columns) to allow paired CV splits.
  # ------------------------------------------------------------------
  if (!is.matrix(Y1)) Y1 <- as.matrix(Y1)
  if (!is.matrix(Y2)) Y2 <- as.matrix(Y2)
  if (any(!is.finite(Y1)) || any(!is.finite(Y2))) stop("Y1 and Y2 must not contain NA/NaN/Inf.")
  div <- as.integer(div)
  if (min(Y1) < 0 || min(Y2) < 0) stop("Y1 and Y2 must be non-negative.")
  if (ncol(Y1) != ncol(Y2)) {
    stop("Y1 and Y2 must have the same number of columns (samples).")
  }

  P1 <- nrow(Y1)
  P2 <- nrow(Y2)
  N  <- ncol(Y1)

  if (div < 2L) {
    stop("div (number of CV folds) must be >= 2.")
  }
  if (div > N) {
    stop("div (number of CV folds) must be <= number of samples.")
  }

  # ------------------------------------------------------------------
  # 2. Handle extra arguments for nmf.sem
  #
  # We collect additional arguments in 'extra_args' but explicitly remove
  # those that are managed at the CV level:
  #   - div, shuffle : used only here, not passed to nmf.sem.
  #   - rank        : passed explicitly from nmf.sem.cv.
  #   - seed        : fold-specific seeds are generated here.
  # ------------------------------------------------------------------
  extra_args <- list(...)

  extra_args$div     <- NULL
  extra_args$shuffle <- NULL
  extra_args$rank    <- NULL
  extra_args$seed    <- NULL

  # ------------------------------------------------------------------
  # 3. Set RNG for CV partition and per-fold seeds
  #
  # If a master 'seed' is given:
  #   - it is used to define the CV partition (sample permutation),
  #   - independent seeds for each nmf.sem run are drawn.
  # If 'seed' is NULL:
  #   - CV partition uses the current RNG state,
  #   - nmf.sem runs use whatever the global RNG state is at call time.
  # ------------------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Sample indices for fold division
  if (shuffle) {
    perm_index <- sample.int(N)
  } else {
    perm_index <- seq_len(N)
  }

  # Per-fold seeds for nmf.sem (optional, only if master seed specified)
  if (!is.null(seed)) {
    seeds_fold <- sample.int(.Machine$integer.max, div)
  } else {
    seeds_fold <- rep(NA_integer_, div)
  }

  # ------------------------------------------------------------------
  # 4. Create CV folds
  #
  # We assign approximately N/div samples to each fold, distributing
  # any remainder one-by-one to the earliest folds.
  # 'block[i] = k' means sample i belongs to fold k (as test set).
  # ------------------------------------------------------------------
  remainder <- N %% div
  division  <- N %/% div
  block     <- integer(N)

  processed_count <- 0L
  for (i in 1:(div - 1L)) {
    plus       <- ifelse(i <= remainder, 1L, 0L)
    chunk_size <- division + plus
    idx_range  <- (processed_count + 1L):(processed_count + chunk_size)
    target_idx <- perm_index[idx_range]
    block[target_idx] <- i
    processed_count   <- processed_count + chunk_size
  }
  # Last fold gets all remaining samples
  target_idx <- perm_index[(processed_count + 1L):N]
  block[target_idx] <- div

  # Per-fold CV loss (MAE) will be stored here
  objfunc.block <- numeric(div)

  # ------------------------------------------------------------------
  # 5. Cross-validation loop
  #
  # For each fold j:
  #   - train on all samples not in fold j,
  #   - test on samples in fold j,
  #   - fit nmf.sem on training data,
  #   - compute MAE on test block from equilibrium mapping M.model.
  # ------------------------------------------------------------------
  for (j in 1:div) {
    # Train / Test split
    train_idx <- block != j
    test_idx  <- block == j

    Y1_train <- Y1[, train_idx, drop = FALSE]
    Y1_test  <- Y1[, test_idx,  drop = FALSE]
    Y2_train <- Y2[, train_idx, drop = FALSE]
    Y2_test  <- Y2[, test_idx,  drop = FALSE]

    # Fold-specific seed for nmf.sem
    seed_j <- if (!is.null(seed)) seeds_fold[j] else NULL

    # Assemble arguments for nmf.sem
    nmf.sem.args <- c(
      extra_args,   # User-specified additional arguments (e.g., Q)
      list(
        Y1         = Y1_train,
        Y2         = Y2_train,
        rank       = rank,
        X.init     = X.init,
        X.L2.ortho = X.L2.ortho,
        C1.L1      = C1.L1,
        C2.L1      = C2.L1,
        epsilon    = epsilon,
        maxit      = maxit
      )
    )
    # Attach seed only when it is defined
    if (!is.null(seed_j)) {
      nmf.sem.args$seed <- seed_j
    }

    # Call nmf.sem on the training data (suppress messages for cleaner CV output)
    res_j <- suppressMessages(do.call("nmf.sem", nmf.sem.args))

    # If mapping is not usable, penalize this fold (do not crash CV)
    if (is.null(res_j$M.model) || any(!is.finite(res_j$M.model))) {
      objfunc.block[j] <- Inf
      next
    }
    Pre_test <- res_j$M.model %*% Y2_test
    if (any(!is.finite(Pre_test))) {
      objfunc.block[j] <- Inf
      next
    }
    objfunc.block[j] <- mean(abs(Y1_test - Pre_test))
  }

  # ------------------------------------------------------------------
  # 6. Aggregate CV score
  #
  # The overall CV criterion is the average MAE across all folds.
  # This is typically minimized over hyperparameter grids
  # (e.g., rank, X.L2.ortho, C1.L1, C2.L1) when tuning NMF-SEM.
  # ------------------------------------------------------------------
  objfunc <- mean(objfunc.block)
  return(objfunc)
}



#' @title Heuristic Variable Splitting for NMF-SEM
#'
#' @description
#' Infers a heuristic partition of observed variables into exogenous (\eqn{Y_2})
#' and endogenous (\eqn{Y_1}) blocks for use in NMF-SEM.
#' The method is based on positive-SEM logic, causal ordering, and optional
#' sign alignment using the first principal component (PC1).
#'
#' The procedure:
#' \itemize{
#'   \item internally standardizes variables (mean 0, sd 1),
#'   \item optionally flips signs so that most variables align positively with PC1,
#'   \item infers a causal ordering by repeatedly regressing each variable on the
#'         remaining ones and selecting the variable with the largest minimum
#'         standardized coefficient,
#'   \item determines an exogenous block by scanning the ordering from upstream
#'         and stopping at the first variable whose strongest parent coefficient
#'         exceeds \code{threshold}.
#' }
#'
#' If \code{n.exogenous} is supplied, it overrides the automatic threshold rule.
#'
#' @param x A numeric matrix or data frame with
#'   \strong{rows = samples} and \strong{columns = observed variables}.
#' @param n.exogenous Optional integer specifying the number of exogenous variables
#'   (\eqn{Y_2}). If \code{NULL}, the number is inferred automatically by the
#'   coefficient cut-off rule.
#' @param threshold Standardized regression-coefficient threshold used in the
#'   automatic exogenous–endogenous split. A variable is treated as endogenous
#'   once its maximum standardized parent coefficient exceeds this value.
#'   (Default: \code{0.1})
#' @param auto.flipped Logical; if \code{TRUE}, applies PC1-based automatic
#'   sign flipping after standardization to ensure consistent orientation.
#'   (Default: \code{TRUE})
#' @param verbose Logical; if \code{TRUE}, prints progress messages and the
#'   resulting variable split. (Default: \code{TRUE})
#'
#' @return A list with:
#'   \item{endogenous.variables}{
#'     Character vector of variables selected as endogenous (\eqn{Y_1}).}
#'   \item{exogenous.variables}{
#'     Character vector of variables selected as exogenous (\eqn{Y_2}).}
#'   \item{ordered.variables}{
#'     Variables in inferred causal order (from exogenous to endogenous).}
#'   \item{is.flipped}{
#'     Logical vector indicating which variables were sign-flipped during processing.}
#'   \item{n.exogenous}{
#'     Integer giving the number of exogenous variables.}
#'
#' @export
nmf.sem.split <- function(x, n.exogenous = NULL, threshold = 0.1,
                          auto.flipped = TRUE, verbose = TRUE) {

  if (!is.matrix(x) && !is.data.frame(x))
    stop("x must be a numeric matrix or data frame.")

  X_raw <- as.matrix(x)
  P <- ncol(X_raw)
  col_names <- colnames(X_raw)
  if (is.null(col_names)) {
    col_names <- paste0("V", 1:P)
    colnames(X_raw) <- col_names
  }

  # --------------------------------------------------------------------
  # Preprocessing Step 1: Standardize all variables
  #
  # Variables are centered and scaled (mean 0, sd 1). NMF-SEM requires
  # non-negative matrices, but the purpose of this function is only to
  # infer variable roles (Y1/Y2), so standardized values are allowed.
  #
  # Missing or NaN values resulting from constant columns are set to 0.
  # --------------------------------------------------------------------
  X_calc <- scale(X_raw, center = TRUE, scale = TRUE)
  X_calc[is.na(X_calc)] <- 0
  X_calc[is.nan(X_calc)] <- 0

  all_indices <- 1:P

  # --------------------------------------------------------------------
  # Preprocessing Step 2: Optional sign flipping based on PC1 alignment
  #
  # In positive SEM (and NMF-SEM), variables should ideally have
  # consistent sign orientation. To enforce this heuristic, variables
  # negatively correlated with the first principal component are flipped.
  #
  # This stabilizes the causal-ordering heuristic by avoiding mixtures
  # of arbitrary sign conventions in the raw data.
  # --------------------------------------------------------------------
  is.flipped <- rep(FALSE, P)
  names(is.flipped) <- col_names

  if (auto.flipped) {
    if (verbose) cat("Step 0: Checking correlations with PC1 (on standardized data)...\n")

    svd_res <- svd(X_calc)
    pc1 <- svd_res$u[, 1]

    cors <- stats::cor(X_calc, pc1)
    # Ensure majority alignment with PC1
    if (stats::median(cors, na.rm = TRUE) < 0) cors <- -cors

    flip_idx <- which(cors < 0)

    if (length(flip_idx) > 0) {
      is.flipped[flip_idx] <- TRUE
      X_calc[, flip_idx] <- -X_calc[, flip_idx]

      if (verbose) {
        cat(sprintf("   -> Detected %d flipped variables: %s\n",
                    length(flip_idx), paste(col_names[flip_idx], collapse=", ")))
      }
    }
  }

  # --------------------------------------------------------------------
  # Step 1: Causal ordering heuristic
  #
  # We infer an ordering of variables consistent with positive-SEM logic:
  # repeatedly select the variable that has the *largest minimum* coefficient
  # when regressed on the remaining variables. This favors variables that
  # are least explained by others → likely exogenous.
  #
  # The resulting order approximates a causal topological order in which
  # exogenous variables appear early and endogenous variables later.
  # --------------------------------------------------------------------
  if (verbose) cat("Step 1: Inferring Causal Ordering...\n")

  active_set <- all_indices
  ordering_reversed <- integer(P)

  for (t in 1:(P - 1)) {
    scores <- numeric(length(active_set))

    for (i in seq_along(active_set)) {
      target_col <- active_set[i]
      pred_cols <- active_set[-i]

      y_vec <- X_calc[, target_col]
      X_mat <- X_calc[, pred_cols, drop = FALSE]

      coefs <- tryCatch({
        stats::coef(stats::lm(y_vec ~ X_mat - 1))
      }, error = function(e) rep(NA_real_, length(pred_cols)))

      # Positive SEM → keep the smallest coefficient (weakest positive predictor)
      if (all(is.na(coefs))) {
        scores[i] <- -Inf
      } else {
        scores[i] <- min(coefs, na.rm = TRUE)
      }
    }

    best_idx <- which.max(scores)
    ordering_reversed[t] <- active_set[best_idx]
    active_set <- active_set[-best_idx]
  }
  ordering_reversed[P] <- active_set[1]

  # Causal order: exogenous → endogenous
  ordering_indices <- rev(ordering_reversed)

  # --------------------------------------------------------------------
  # Step 2: Automatic identification of exogenous variables
  #
  # Sweep through the causal ordering. For each variable, regress it on all
  # earlier variables. If its strongest parent coefficient exceeds the
  # threshold, the variable is considered endogenous.
  #
  # Variables before this point → exogenous (Y2)
  # Variables after this point → endogenous (Y1)
  #
  # If n.exogenous is given, it overrides this automatic rule.
  # --------------------------------------------------------------------
  if (is.null(n.exogenous)) {
    if (verbose) cat("Step 2: Detecting optimal cut-off for exogenous variables...\n")
    cutoff <- 1

    for (k in 2:(P - 1)) {
      curr_idx <- ordering_indices[k]
      parent_indices <- ordering_indices[1:(k - 1)]

      y_vec <- X_calc[, curr_idx]
      X_parents <- X_calc[, parent_indices, drop = FALSE]

      coefs <- stats::coef(stats::lm(y_vec ~ X_parents - 1))
      max_influence <- max(coefs, na.rm = TRUE)

      if (max_influence > threshold) {
        if (verbose)
          cat(sprintf("   [%d] %s : Max std.coef=%.3f -> Endogenous (Stop)\n",
                      k, col_names[curr_idx], max_influence))
        break
      } else {
        cutoff <- k
        if (verbose)
          cat(sprintf("   [%d] %s : Max std.coef=%.3f -> Exogenous (Continue)\n",
                      k, col_names[curr_idx], max_influence))
      }
    }
    n.exogenous <- cutoff
  }

  # --------------------------------------------------------------------
  # Step 3: Final classification into Y1 and Y2
  #
  # Variables appearing early in the ordering (determined by cut-off) are
  # treated as exogenous (Y2). The remainder are endogenous (Y1).
  #
  # Ordered list shows the full inferred causal sequence.
  # --------------------------------------------------------------------
  idx_exo <- ordering_indices[1:n.exogenous]
  idx_endo <- ordering_indices[(n.exogenous + 1):P]

  exogenous.variables <- col_names[idx_exo]
  endogenous.variables <- col_names[idx_endo]
  ordered.variables <- col_names[ordering_indices]

  if (verbose) {
    cat("\n--- Auto Split Result ---\n")
    cat(sprintf("Exogenous (Y2, n=%d): %s\n",
                n.exogenous, paste(exogenous.variables, collapse=", ")))
    cat(sprintf("Endogenous (Y1, n=%d): %s ...\n",
                length(endogenous.variables),
                paste(utils::head(endogenous.variables, 3), collapse=", ")))
  }

  return(list(
    endogenous.variables = endogenous.variables,
    exogenous.variables = exogenous.variables,
    ordered.variables = ordered.variables,
    is.flipped = is.flipped,
    n.exogenous = as.integer(n.exogenous)
  ))
}




############################################################
## Common DOT Helpers
############################################################

#' Determine the decimal digits based on a threshold
#'
#' This helper computes the number of decimal places that should be
#' used for formatting coefficient labels, based on the magnitude
#' of the threshold used for edge visualization.
#'
#' @param threshold Numeric scalar (>0).
#' @return Integer specifying the number of decimal places.
#' @keywords internal
#' @noRd
.nmfkc_dot_digits_from_threshold <- function(threshold) {
  if (!is.finite(threshold) || threshold <= 0) {
    return(2L)
  }
  if (threshold >= 1) {
    return(0L)
  }
  s <- format(threshold, scientific = FALSE, trim = TRUE)
  # Examples: "0.01" -> "01", "0.005" -> "005"
  s_dec <- sub("^0\\.", "", s)
  nchar(s_dec)
}


#' Format coefficient values for DOT edge labels
#'
#' If \code{digits} is \code{NULL}, the function uses the traditional
#' magnitude-based formatting rules.
#' If \code{digits} is provided, the coefficient is formatted with the
#' specified number of decimal places (typically derived from the threshold).
#'
#' @param x Numeric scalar.
#' @param digits Integer or \code{NULL}; number of decimal places.
#' @return Character string suitable for use inside DOT \code{label="..."}.
#' @keywords internal
#' @noRd
.nmfkc_dot_format_coef <- function(x, digits = NULL) {
  if (!is.finite(x)) return("")

  # Default behavior: magnitude-dependent auto formatting
  if (is.null(digits)) {
    ax <- abs(x)
    if (ax == 0) return("0")

    if (ax >= 1) {
      sprintf("%.2f", x)
    } else if (ax >= 1e-1) {
      sprintf("%.3f", x)
    } else if (ax >= 1e-2) {
      sprintf("%.4f", x)
    } else {
      sprintf("%.3e", x)
    }

    # Threshold-based fixed-digit formatting
  } else {
    if (digits <= 0) {
      sprintf("%.0f", x)
    } else {
      fmt <- paste0("%.", digits, "f")
      sprintf(fmt, x)
    }
  }
}


#' Sanitize character vectors for DOT node / cluster IDs
#'
#' Replace characters that are not alphanumeric, underscore, or dot
#' with underscores to generate safe DOT identifiers.
#'
#' @param x Character vector.
#' @return Character vector with sanitized identifiers.
#' @keywords internal
#' @noRd
.nmfkc_dot_sanitize_id <- function(x) {
  gsub("[^[:alnum:]_.]", "_", x, perl = TRUE)
}


#' Generate a standard DOT graph header for NMF-related diagrams
#'
#' Provides a consistent DOT header used by NMF visualization functions.
#'
#' @param graph_name Character; name of the DOT graph.
#' @param rankdir Character; Graphviz rank direction (e.g., "LR", "RL", "TB", "BT").
#' @param fontname Character; default font name used in the graph.
#'
#' @return Character scalar containing DOT header lines.
#' @keywords internal
#' @noRd
.nmfkc_dot_header <- function(graph_name = "NMF_GRAPH",
                              rankdir    = "LR",
                              fontname   = "Meiryo") {
  paste0(
    "digraph ", graph_name, " {\n",
    "  graph [rankdir=", rankdir, " compound=true];\n",
    "  splines=true; nodesep=0.4; ranksep=0.7; fontname=\"", fontname, "\";\n"
  )
}


#' Generate a DOT subgraph cluster with nodes
#'
#' Internal helper to define a cluster (subgraph) of nodes with shared
#' visual properties.
#' If \code{cluster_style = "none"}, nodes are listed individually without
#' creating a DOT subgraph.
#'
#' @param cluster_id Character; cluster identifier suffix (e.g. "Y", "X", "F").
#' @param title Character; cluster label (ignored if \code{cluster_style = "none"}).
#' @param node_ids Character vector; internal DOT node IDs.
#' @param node_labels Character vector; display labels for nodes.
#' @param shape Character; Graphviz node shape.
#' @param fill Logical; whether to use filled node shapes.
#' @param fillcolor Character; fill color (if \code{fill = TRUE}).
#' @param line_width Numeric; node border width.
#' @param indent Character; indentation prefix for formatting.
#' @param cluster_style Character; boundary style ("rounded", "dashed", "invis", "none").
#' @param cluster_color Character; cluster boundary color.
#' @param cluster_penwidth Numeric; cluster boundary width.
#'
#' @return Character scalar containing DOT code for the cluster.
#' @keywords internal
#' @noRd
.nmfkc_dot_cluster_nodes <- function(cluster_id,
                                     title,
                                     node_ids,
                                     node_labels,
                                     shape      = "box",
                                     fill       = TRUE,
                                     fillcolor  = "white",
                                     line_width = 1.5,
                                     indent     = "  ",
                                     cluster_style    = "rounded",
                                     cluster_color    = "black",
                                     cluster_penwidth = 1.0) {

  if (length(node_ids) == 0L) return("")
  if (length(node_ids) != length(node_labels)) {
    stop("node_ids and node_labels must have the same length.")
  }

  # Node style definition
  if (fill) {
    node_style <- sprintf(
      '%snode [shape=%s, style="filled,rounded", fillcolor="%s", color=black, penwidth=%.1f];\n',
      indent, shape, fillcolor, line_width
    )
  } else {
    node_style <- sprintf(
      '%snode [shape=%s, style="rounded", color=black, penwidth=%.1f];\n',
      indent, shape, line_width
    )
  }

  # Case 1: no cluster → output nodes directly
  if (identical(cluster_style, "none")) {
    scr <- node_style
    for (i in seq_along(node_ids)) {
      scr <- paste0(
        scr,
        sprintf('%s%s [label="%s"];\n', indent, node_ids[i], node_labels[i])
      )
    }
    return(scr)
  }

  # Case 2: cluster with boundary
  scr <- paste0(
    indent, "subgraph cluster_", cluster_id,
    '{label="', title, '"',
    ' style="', cluster_style, '"',
    ' color="', cluster_color, '"',
    ' penwidth=', sprintf("%.1f", cluster_penwidth), ";\n",
    node_style
  )

  for (i in seq_along(node_ids)) {
    scr <- paste0(
      scr,
      sprintf('%s  %s [label="%s"];\n', indent, node_ids[i], node_labels[i])
    )
  }
  paste0(scr, indent, "}\n")
}


#' Scale DOT edge width based on coefficient magnitude
#'
#' Computes an appropriate \code{penwidth} value for Graphviz edges
#' by scaling the coefficient relative to the maximum coefficient in
#' its category.
#'
#' @param value Numeric; coefficient value.
#' @param max_value Numeric; maximum coefficient within its group.
#' @param weight_scale Numeric; global scale factor for edge width.
#' @param min_pw Numeric; minimum penwidth.
#'
#' @return Numeric penwidth value.
#' @keywords internal
#' @noRd
.nmfkc_dot_penwidth <- function(value,
                                max_value,
                                weight_scale = 5,
                                min_pw       = 0.5) {
  if (!is.finite(max_value) || max_value <= 0) return(min_pw)
  if (!is.finite(value)     || value <= 0)     return(min_pw)
  max(min_pw, value * weight_scale / max_value)
}



############################################################
## 1. nmf.sem.DOT  (for NMF-SEM visualization)
############################################################

#' Generate a Graphviz DOT Diagram for an NMF-SEM Model
#'
#' @description
#' Creates a Graphviz DOT script that visualizes the structural network
#' estimated by \code{nmf.sem}.
#' The resulting diagram displays:
#' \itemize{
#'   \item endogenous observed variables (\eqn{Y_1}),
#'   \item exogenous observed variables (\eqn{Y_2}),
#'   \item latent factors (\eqn{F_1}, \dots, \eqn{F_Q}),
#' }
#' together with the non-negative path coefficients whose magnitudes
#' exceed a user-specified threshold.
#'
#' Directed edges represent estimated relationships:
#' \itemize{
#'   \item \eqn{Y_2 \rightarrow F_q}: entries of \code{C2} (exogenous loadings),
#'   \item \eqn{F_q \rightarrow Y_1}: rows of \code{X} (factor-to-endogenous mappings),
#'   \item \eqn{Y_1 \rightarrow F_q}: entries of \code{C1} (feedback paths).
#' }
#'
#' Edge widths are scaled by coefficient magnitude, and nodes are placed
#' in optional visual clusters. Only variables participating in
#' edges above the threshold are displayed, while latent factors are always shown.
#'
#' @param result A list returned by \code{nmf.sem}, containing matrices
#'   \code{X}, \code{C1}, and \code{C2}.
#' @param weight_scale Base scaling factor for edge widths.
#' @param weight_scale_y2f Optional override for scaling edges
#'   \eqn{Y_2 \rightarrow F_q}. Defaults to \code{weight_scale}.
#' @param weight_scale_fy1 Optional override for scaling edges
#'   \eqn{F_q \rightarrow Y_1}. Defaults to \code{weight_scale}.
#' @param weight_scale_feedback Optional override for scaling feedback edges
#'   \eqn{Y_1 \rightarrow F_q}. Defaults to \code{weight_scale}.
#' @param threshold Minimum coefficient value needed for an edge to be drawn.
#' @param rankdir Graphviz rank direction (e.g., \code{"LR"}, \code{"TB"}).
#' @param fill Logical; whether to use filled node shapes.
#' @param cluster.box Character string controlling the visibility and style
#'   of cluster frames around Y2, factors, and Y1 blocks.
#'   One of \code{"normal"}, \code{"faint"}, \code{"invisible"}, \code{"none"}.
#' @param cluster.labels Optional character vector of length 3 giving custom
#'   labels for the Y2, factor, and Y1 clusters.
#'
#' @return A character string representing a valid Graphviz DOT script.
#'
#' @export
nmf.sem.DOT <- function(result,
                        weight_scale          = 5,
                        weight_scale_y2f      = weight_scale,
                        weight_scale_fy1      = weight_scale,
                        weight_scale_feedback = weight_scale,
                        threshold             = 0.01,
                        rankdir               = "LR",
                        fill                  = TRUE,
                        cluster.box           = c("normal", "faint", "invisible", "none"),
                        cluster.labels        = NULL) {

  ## ---------------------------------------------------------------
  ## Cluster style selection
  ## ---------------------------------------------------------------
  cluster.box <- match.arg(cluster.box)

  cluster_style <- switch(cluster.box,
                          normal    = "rounded",
                          faint     = "rounded,dashed",
                          invisible = "rounded",
                          none      = "none")

  cluster_color <- switch(cluster.box,
                          normal    = "black",
                          faint     = "gray80",
                          invisible = "none",
                          none      = "none")

  cluster_penwidth <- switch(cluster.box,
                             normal    = 1.0,
                             faint     = 0.7,
                             invisible = 0.0,
                             none      = 0.0)

  ## ---------------------------------------------------------------
  ## Cluster titles
  ## ---------------------------------------------------------------
  default.labels <- c("Exogenous (Y2)", "Latent Factors", "Endogenous (Y1)")
  if (is.null(cluster.labels)) {
    titles <- default.labels
  } else {
    if (!is.character(cluster.labels) || length(cluster.labels) != 3L) {
      stop("cluster.labels must be a character vector of length 3: c(label_Y2, label_F, label_Y1).")
    }
    titles <- cluster.labels
  }

  ## ---------------------------------------------------------------
  ## Extract matrices and labels
  ## ---------------------------------------------------------------
  X  <- result$X
  C1 <- result$C1
  C2 <- result$C2

  if (is.null(X) || is.null(C1) || is.null(C2)) {
    stop("result must contain elements X, C1, and C2.")
  }

  Y1_labels <- rownames(X)
  Y2_labels <- colnames(C2)

  if (is.null(Y1_labels)) Y1_labels <- paste0("Y1_", seq_len(nrow(X)))
  if (is.null(Y2_labels)) Y2_labels <- paste0("Y2_", seq_len(ncol(C2)))

  P1 <- length(Y1_labels)
  P2 <- length(Y2_labels)
  Q  <- ncol(X)

  ## ---------------------------------------------------------------
  ## Identify nodes involved in edges >= threshold
  ## ---------------------------------------------------------------
  used_y2 <- apply(C2, 2L, function(col) any(col >= threshold, na.rm = TRUE))
  used_y1_from_X  <- apply(X,  1L, function(row) any(row >= threshold, na.rm = TRUE))
  used_y1_from_C1 <- apply(C1, 2L, function(col) any(col >= threshold, na.rm = TRUE))
  used_y1 <- used_y1_from_X | used_y1_from_C1

  idx_y1 <- which(used_y1)
  idx_y2 <- which(used_y2)

  ## ---------------------------------------------------------------
  ## Assign internal DOT node IDs
  ## ---------------------------------------------------------------
  Y1_ids <- .nmfkc_dot_sanitize_id(paste0("Y1_", seq_len(P1)))
  Y2_ids <- .nmfkc_dot_sanitize_id(paste0("Y2_", seq_len(P2)))
  F_ids  <- .nmfkc_dot_sanitize_id(paste0("F_",  seq_len(Q)))

  ## Node colors
  COLOR_Y2_NODE     <- "lightcoral"
  COLOR_Y1_NODE     <- "lightblue"
  COLOR_FACTOR_NODE <- "wheat"

  ## ---------------------------------------------------------------
  ## Header
  ## ---------------------------------------------------------------
  dot_script <- .nmfkc_dot_header(
    graph_name = "NMF_SEM_Full_Mechanism",
    rankdir    = rankdir
  )

  ## ---------------------------------------------------------------
  ## Y2 cluster
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    "\n  // Exogenous variables (Y2)\n",
    .nmfkc_dot_cluster_nodes(
      cluster_id        = "Y2",
      title             = titles[1],
      node_ids          = Y2_ids[idx_y2],
      node_labels       = Y2_labels[idx_y2],
      shape             = "box",
      fill              = fill,
      fillcolor         = COLOR_Y2_NODE,
      line_width        = 1.5,
      cluster_style     = cluster_style,
      cluster_color     = cluster_color,
      cluster_penwidth  = cluster_penwidth
    )
  )

  ## ---------------------------------------------------------------
  ## Y1 cluster
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    "\n  // Endogenous variables (Y1)\n",
    .nmfkc_dot_cluster_nodes(
      cluster_id        = "Y1",
      title             = titles[3],
      node_ids          = Y1_ids[idx_y1],
      node_labels       = Y1_labels[idx_y1],
      shape             = "box",
      fill              = fill,
      fillcolor         = COLOR_Y1_NODE,
      line_width        = 1.5,
      cluster_style     = cluster_style,
      cluster_color     = cluster_color,
      cluster_penwidth  = cluster_penwidth
    )
  )

  ## ---------------------------------------------------------------
  ## Latent factor cluster
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    "\n  // Latent Factors (F)\n",
    .nmfkc_dot_cluster_nodes(
      cluster_id        = "F",
      title             = titles[2],
      node_ids          = F_ids,
      node_labels       = paste0("Factor ", seq_len(Q)),
      shape             = "ellipse",
      fill              = fill,
      fillcolor         = COLOR_FACTOR_NODE,
      line_width        = 1.0,
      cluster_style     = cluster_style,
      cluster_color     = cluster_color,
      cluster_penwidth  = cluster_penwidth
    )
  )

  ## ---------------------------------------------------------------
  ## Edge defaults
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  edge [fontname="Meiryo", fontsize=8, arrowhead=open];\n'
  )

  pw     <- .nmfkc_dot_penwidth
  digits <- .nmfkc_dot_digits_from_threshold(threshold)
  fmtc   <- function(x) .nmfkc_dot_format_coef(x, digits)

  ## ---------------------------------------------------------------
  ## 1. Y2 → F edges (C2)
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  // 1. External Driving (Y2 -> Factor) [C2]\n',
    '  edge [color=black, fontcolor=black, style=solid];\n'
  )

  max_C2 <- suppressWarnings(max(C2, na.rm = TRUE))
  if (is.finite(max_C2) && max_C2 > 0) {
    for (q in seq_len(Q)) {
      for (p2 in idx_y2) {
        weight <- C2[q, p2]
        if (is.finite(weight) && weight >= threshold) {
          pen <- pw(weight, max_C2, weight_scale_y2f)
          lab <- fmtc(weight)
          path <- sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                          Y2_ids[p2], F_ids[q], lab, pen)
          dot_script <- paste0(dot_script, path)
        }
      }
    }
  }

  ## ---------------------------------------------------------------
  ## 2. F → Y1 edges (X)
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  // 2. Generation (Factor -> Y1) [X]\n',
    '  edge [color="gray0", fontcolor="gray0", style=solid];\n'
  )

  max_X <- suppressWarnings(max(X, na.rm = TRUE))
  if (is.finite(max_X) && max_X > 0) {
    for (q in seq_len(Q)) {
      for (p1 in idx_y1) {
        weight <- X[p1, q]
        if (is.finite(weight) && weight >= threshold) {
          pen <- pw(weight, max_X, weight_scale_fy1)
          lab <- fmtc(weight)
          path <- sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                          F_ids[q], Y1_ids[p1], lab, pen)
          dot_script <- paste0(dot_script, path)
        }
      }
    }
  }

  ## ---------------------------------------------------------------
  ## 3. Y1 → F edges (C1 feedback)
  ## ---------------------------------------------------------------
  dot_script <- paste0(
    dot_script,
    '\n  // 3. Internal Feedback (Y1 -> Factor) [C1]\n',
    '  edge [style=dashed, color="gray0", fontcolor="gray0"];\n'
  )

  max_C1 <- suppressWarnings(max(C1, na.rm = TRUE))
  if (is.finite(max_C1) && max_C1 > 0) {
    for (q in seq_len(Q)) {
      for (p1 in idx_y1) {
        weight <- C1[q, p1]
        if (is.finite(weight) && weight >= threshold) {
          pen <- pw(weight, max_C1, weight_scale_feedback)
          lab <- fmtc(weight)
          path <- sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                          Y1_ids[p1], F_ids[q], lab, pen)
          dot_script <- paste0(dot_script, path)
        }
      }
    }
  }

  paste0(dot_script, "}\n")
}





############################################################
## 2. nmfkc.DOT  (Static NMF / NMF-with-covariates visualization)
############################################################

#' Generate Graphviz DOT Scripts for NMF or NMF-with-Covariates Models
#'
#' @description
#' Produces a Graphviz DOT script visualizing the structure of an NMF model
#' (\eqn{Y \approx X C A}) or its simplified forms.
#'
#' Supported visualization types:
#' \itemize{
#'   \item \code{"YX"} — Standard NMF view: latent factors \eqn{X} map to observations \eqn{Y}.
#'   \item \code{"YA"} — Direct regression view: covariates \eqn{A} map directly to \eqn{Y}
#'         using the combined coefficient matrix \eqn{X C}.
#'   \item \code{"YXA"} — Full tri-factorization: \eqn{A \rightarrow C \rightarrow X \rightarrow Y}.
#' }
#'
#' Edge widths are scaled by coefficient magnitude, and nodes with no edges
#' above the threshold are omitted from the visualization.
#'
#' @param x The return value from \code{nmfkc}, containing matrices
#'   \code{X}, \code{B}, and optionally \code{C}.
#' @param type Character string specifying the visualization style:
#'   one of \code{"YX"}, \code{"YA"}, \code{"YXA"}.
#' @param threshold Minimum coefficient magnitude to display an edge.
#' @param rankdir Graphviz rank direction (e.g., \code{"LR"}, \code{"TB"}).
#' @param fill Logical; whether nodes should be drawn with filled shapes.
#' @param weight_scale Base scaling factor for edge widths.
#' @param weight_scale_ax Scaling factor for edges \eqn{A \rightarrow X} (type \code{"YXA"}).
#' @param weight_scale_xy Scaling factor for edges \eqn{X \rightarrow Y}.
#' @param weight_scale_ay Scaling factor for edges \eqn{A \rightarrow Y} (type \code{"YA"}).
#' @param Y.label Optional character vector for labels of Y nodes.
#' @param X.label Optional character vector for labels of X (latent factor) nodes.
#' @param A.label Optional character vector for labels of A (covariate) nodes.
#' @param Y.title Cluster title for Y nodes.
#' @param X.title Cluster title for X nodes.
#' @param A.title Cluster title for A nodes.
#'
#' @return A character string representing a Graphviz DOT script.
#'
#' @seealso \code{nmfkc}
#' @export
nmfkc.DOT <- function(
    x,
    type = c("YX","YA","YXA"),
    threshold = 0.01,
    rankdir   = "LR",
    fill      = TRUE,
    weight_scale    = 5,
    weight_scale_ax = weight_scale,
    weight_scale_xy = weight_scale,
    weight_scale_ay = weight_scale,
    Y.label = NULL, X.label = NULL, A.label = NULL,
    Y.title = "Observation (Y)",
    X.title = "Basis (X)",
    A.title = "Covariates (A)"
) {

  type <- match.arg(type)

  ## ---------------------------------------------------------
  ## Required matrices
  ## ---------------------------------------------------------
  X <- x$X
  B <- x$B
  if (is.null(X) || is.null(B)) {
    stop("x must contain X and B.")
  }

  ## If C exists and is a proper NMF-with-covariates factor:
  hasA <- !is.null(x$C) && ncol(x$C) != ncol(B)

  P <- nrow(X)
  Q <- ncol(X)

  ## ---------------------------------------------------------
  ## Labels
  ## ---------------------------------------------------------
  Y_labels <- if (is.null(Y.label)) rownames(X) else Y.label
  X_labels <- if (is.null(X.label)) colnames(X) else X.label

  if (is.null(Y_labels)) Y_labels <- paste0("Y", seq_len(P))
  if (is.null(X_labels)) X_labels <- paste0("Factor", seq_len(Q))

  Y_ids <- .nmfkc_dot_sanitize_id(paste0("Y_", seq_len(P)))
  X_ids <- .nmfkc_dot_sanitize_id(paste0("X_", seq_len(Q)))

  ## ---------------------------------------------------------
  ## Covariates and tri-factorization
  ## ---------------------------------------------------------
  if (hasA) {
    C <- as.matrix(x$C)          # Q x R
    A_cols <- ncol(C)
    A_labels <- if (is.null(A.label)) colnames(C) else A.label
    if (is.null(A_labels)) A_labels <- paste0("A", seq_len(A_cols))
    A_ids <- .nmfkc_dot_sanitize_id(paste0("A_", seq_len(A_cols)))

    ## Combined mapping for type = "YA"
    XC_mat <- X %*% C   # P x R

  } else if (type == "YX") {
    ## No A block needed
    C <- NULL
    A_cols <- 0L
    A_labels <- NULL
    A_ids <- NULL
    XC_mat <- NULL

  } else {
    stop("The model structure (matrix C) is incompatible with type 'YXA' or 'YA'.")
  }

  ## ---------------------------------------------------------
  ## DOT header
  ## ---------------------------------------------------------
  scr <- .nmfkc_dot_header(graph_name = "NMF", rankdir = rankdir)

  ## ---------------------------------------------------------
  ## Y node cluster
  ## ---------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  // Output variables (Y)\n',
    .nmfkc_dot_cluster_nodes(
      cluster_id  = "Y",
      title       = Y.title,
      node_ids    = Y_ids,
      node_labels = Y_labels,
      shape       = "box",
      fill        = fill,
      fillcolor   = "lightblue",
      line_width  = 1.5
    )
  )

  ## ---------------------------------------------------------
  ## X node cluster
  ## (Hidden when type = "YA")
  ## ---------------------------------------------------------
  if (type != "YA") {
    scr <- paste0(
      scr,
      '\n  // Latent factors (X)\n',
      .nmfkc_dot_cluster_nodes(
        cluster_id  = "X",
        title       = X.title,
        node_ids    = X_ids,
        node_labels = X_labels,
        shape       = "ellipse",
        fill        = fill,
        fillcolor   = "wheat",
        line_width  = 1.0
      )
    )
  }

  ## ---------------------------------------------------------
  ## A node cluster
  ## (Only for "YA" and "YXA")
  ## ---------------------------------------------------------
  if (type != "YX" && hasA) {
    scr <- paste0(
      scr,
      '\n  // Covariates (A)\n',
      .nmfkc_dot_cluster_nodes(
        cluster_id  = "A",
        title       = A.title,
        node_ids    = A_ids,
        node_labels = A_labels,
        shape       = "box",
        fill        = fill,
        fillcolor   = "lightcoral",
        line_width  = 1.5
      )
    )
  }

  ## ---------------------------------------------------------
  ## Edge defaults
  ## ---------------------------------------------------------
  scr <- paste0(
    scr,
    '\n  edge [fontname="Meiryo", fontsize=8, arrowhead=open];\n'
  )

  pw     <- .nmfkc_dot_penwidth
  digits <- .nmfkc_dot_digits_from_threshold(threshold)
  fmtc   <- function(x) .nmfkc_dot_format_coef(x, digits)


  ## =========================================================
  ## Case 1: Full tri-factorization (A → X → Y)
  ## =========================================================
  if (type == "YXA") {

    ## ---- A → X (C) edges ----
    scr <- paste0(
      scr,
      '\n  // A -> X edges (C)\n',
      '  edge [color=black, fontcolor=black, style=solid];\n'
    )

    max_C <- suppressWarnings(max(C, na.rm = TRUE))
    if (is.finite(max_C) && max_C > 0) {
      for (q in seq_len(Q)) {
        for (k in seq_len(A_cols)) {
          val <- C[q, k]
          if (is.finite(val) && val >= threshold) {
            pen <- pw(val, max_C, weight_scale_ax)
            lab <- fmtc(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                      A_ids[k], X_ids[q], lab, pen)
            )
          }
        }
      }
    }

    ## ---- X → Y edges ----
    scr <- paste0(
      scr,
      '\n  // X -> Y edges (X)\n',
      '  edge [color="gray0", fontcolor="gray0", style=solid];\n'
    )

    max_X <- suppressWarnings(max(X, na.rm = TRUE))
    if (is.finite(max_X) && max_X > 0) {
      for (i in seq_len(P)) {
        for (j in seq_len(Q)) {
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


    ## =========================================================
    ## Case 2: Direct regression view (A → Y)
    ## =========================================================
  } else if (type == "YA") {

    scr <- paste0(
      scr,
      '\n  // A -> Y edges (X %*% C)\n',
      '  edge [color=black, fontcolor=black, style=solid];\n'
    )

    max_XC <- suppressWarnings(max(XC_mat, na.rm = TRUE))
    if (is.finite(max_XC) && max_XC > 0) {
      for (i in seq_len(P)) {
        for (k in seq_len(A_cols)) {
          val <- XC_mat[i, k]
          if (is.finite(val) && val >= threshold) {
            pen <- pw(val, max_XC, weight_scale_ay)
            lab <- fmtc(val)
            scr <- paste0(
              scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                      A_ids[k], Y_ids[i], lab, pen)
            )
          }
        }
      }
    }


    ## =========================================================
    ## Case 3: Standard NMF (X → Y)
    ## =========================================================
  } else if (type == "YX") {

    scr <- paste0(
      scr,
      '\n  // X -> Y edges (X)\n',
      '  edge [color="gray0", fontcolor="gray0", style=solid];\n'
    )

    max_X <- suppressWarnings(max(X, na.rm = TRUE))
    if (is.finite(max_X) && max_X > 0) {
      for (i in seq_len(P)) {
        for (j in seq_len(Q)) {
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
  }

  paste0(scr, "}\n")
}





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
    '\n  edge [fontname="Meiryo", fontsize=8, arrowhead=open];\n'
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
    if (max(Ck, na.rm = TRUE) < threshold) {
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
      if (max(Ck[, j], na.rm = TRUE) >= threshold) {
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

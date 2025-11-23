.onAttach <- function(...) {
  packageStartupMessage("Last update on 22 NOV 2025")
  packageStartupMessage("https://github.com/ksatohds/nmfkc")
}


#' @title Construct observation and covariate matrices for a vector autoregressive model
#' @description
#' \code{nmfkc.ar} generates the observation matrix and covariate matrix
#' corresponding to a specified autoregressive lag order.
#'
#' @param Y An observation matrix, where columns are ordered by measurement time points.
#' @param degree The lag order of the autoregressive model. The default is 1.
#' @param intercept Logical. If TRUE (default), an intercept term is added to the covariate matrix.
#'
#' @return A list containing:
#' \item{Y}{Observation matrix constructed according to the specified lag order.}
#' \item{A}{Covariate matrix constructed according to the specified lag order.}
#' \item{A.columns}{Index matrix used to generate \code{A}.}
#' \item{degree.max}{Maximum lag order, set to \eqn{10 \log_{10}(N)} following the \code{ar} function in the \pkg{stats} package.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ar.degree.cv}}, \code{\link{nmfkc.ar.stationarity}}, \code{\link{nmfkc.ar.DOT}}
#' @export
#' @source Satoh, K. (2025). Applying Non-negative Matrix Factorization with Covariates
#'   to Multivariate Time Series Data as a Vector Autoregression Model.
#'    Japanese Journal of Statistics and Data Science. \url{https://doi.org/10.1007/s42081-025-00314-0}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' d <- AirPassengers
#' time <- time(ts(1:length(d),start=c(1949,1),frequency=12))
#' time.vec <- round(as.vector(t(time)),2)
#' Y0 <- matrix(as.vector(d),nrow=1)
#' colnames(Y0) <- time.vec
#' rownames(Y0) <- "t"
#' # nmf with covariates
#' Q <- 1
#' D <- 12
#' a <- nmfkc.ar(Y0,degree=D,intercept=TRUE); Y <- a$Y; A <- a$A
#' res <- nmfkc(Y=Y,A=A,Q=Q,prefix="Factor",epsilon=1e-6)
#' res$r.squared
#' # coefficients
#' print.table(round(res$C,2),zero.print="")
#' # fitted curve
#' plot(as.numeric(colnames(Y)),as.vector(Y),type="l",col=1,xlab="",ylab="AirPassengers")
#' lines(as.numeric(colnames(Y)),as.vector(res$XB),col=2)

nmfkc.ar <- function(Y,degree=1,intercept=T){
  if(is.vector(Y)) Y <- matrix(Y,nrow=1)
  if(!is.matrix(Y)) Y <- as.matrix(Y)

  # --- NA Check ---
  # Y must not contain missing values before generating the lagged matrix A.
  if(any(is.na(Y))){
    stop("Y contains missing values (NA). In NMF-VAR, lagged Y is used to construct the covariate matrix A, which cannot contain missing values. Please impute Y before calling nmfkc.ar().")
  }

  P <- nrow(Y) # Features/Variables
  N <- ncol(Y) # Time points

  # --- 1. Degree Validity Check (Guardrail) ---
  # Check if the lag order is valid with respect to the time series length N.
  if (degree < 1) {
    stop("The 'degree' (lag order) must be 1 or greater.")
  }
  if (degree >= N) {
    stop(paste0("The 'degree' (", degree, ") must be strictly less than the number of time points in Y (", N, ")."))
  }
  # -------------------------------------------------

  # Length of the time series after lagging (N_A)
  N_A <- N - degree
  t_start <- degree + 1

  # --- 2. A.columns.index: Construct index matrix for lagged Y (Fast) ---
  # A.columns.index: Index matrix for Y columns (degree rows x N_A columns)
  A.columns.index <- matrix(0, nrow = degree, ncol = N_A)
  for(i in 1:degree) {
    # Column indices corresponding to lag t-i for times t_start to N
    A.columns.index[i, ] <- (t_start - i) : (N - i)
  }

  # --- 3. A: Construct covariate matrix A in one operation (Fast) ---
  # 1. Slice the necessary data from Y using the vectorized indices
  A_data <- Y[, as.vector(A.columns.index)]

  # 2. Reshape A_data into the final A shape (P*degree rows x N_A columns)
  A <- matrix(A_data, nrow = P * degree, ncol = N_A, byrow = FALSE)

  # --- 4. Set row labels ---
  if(is.null(rownames(Y))) rownames(Y) <- 1:P
  label <- unlist(lapply(1:degree, function(i) paste0(rownames(Y), "_", i)))
  rownames(A) <- label

  # --- 5. Add Intercept (Safe Implementation) ---
  if(intercept){
    # Use matrix(1, ...) to safely add the intercept row, avoiding recycling dependence
    A <- rbind(A, matrix(1, nrow=1, ncol=ncol(A)))
    rownames(A) <- c(label,"(Intercept)")
  }

  # --- 6. Ya (Non-lagged Observation Matrix) ---
  # Y matrix for time points t_start to N
  Ya <- Y[, t_start : N, drop=FALSE]

  # --- 7. Metadata Calculation and Attribute Storage (Changed key to 'function') ---
  # Calculate degree.max based on R's stats::ar() empirical rule (10 * log10(N))
  degree.max <- min(ncol(Ya),floor(10*log10(ncol(Ya))))

  # Store AR parameters as attributes of A for simple model passing
  # These attributes ensure robust parameter transfer to nmfkc() and nmfkc.ar.predict().
  attr(A, "degree") <- degree
  attr(A, "intercept") <- intercept
  attr(A, "function.name") <- "nmfkc.ar" # Identifier for the model type (MODIFIED)

  # A.columns is returned for backward compatibility/diagnostics
  list(Y = Ya,
       A = A,
       A.columns = A.columns.index,
       degree.max = degree.max)
}


#' @title Generate DOT language scripts for vector autoregressive models
#' @description
#' \code{nmfkc.ar.DOT} generates scripts in the DOT language for visualizing
#' vector autoregressive models fitted using \code{nmfkc}.
#'
#' @param x The return value of \code{nmfkc} for a vector autoregressive model.
#' @param degree The maximum lag order to visualize. Default is 1.
#' @param intercept Logical. If TRUE, an intercept node is added. Default is FALSE.
#' @param digits Integer. Number of decimal places to display in edge labels.
#' @param threshold Numeric. Parameters greater than or equal to this threshold are displayed. Default is \eqn{10^{-\code{digits}}}.
#' @param rankdir Graph layout direction in DOT language. Default is "RL". Other options include "LR", "TB", and "BT".
#'
#' @return A character string containing a DOT script, suitable for use with the \pkg{DOT} package or Graphviz tools.
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ar}}
#' @export

#' @title Generate DOT language scripts for vector autoregressive models
#' @description
#' \code{nmfkc.ar.DOT} generates scripts in the DOT language for visualizing
#' vector autoregressive models fitted using \code{nmfkc}.
#'
#' @param x The return value of \code{nmfkc} for a vector autoregressive model.
#' @param degree The maximum lag order to visualize. Default is 1.
#' @param intercept Logical. If TRUE, an intercept node is added. Default is FALSE.
#' @param digits Integer. Number of decimal places to display in edge labels.
#' @param threshold Numeric. Parameters greater than or equal to this threshold are displayed. Default is \eqn{10^{-\code{digits}}}.
#' @param rankdir Graph layout direction in DOT language. Default is "RL". Other options include "LR", "TB", and "BT".
#'
#' @return A character string containing a DOT script, suitable for use with the \pkg{DOT} package or Graphviz tools.
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ar}}
#' @export

nmfkc.ar.DOT <- function(x,degree=1,intercept=FALSE,digits=1,threshold=10^(-digits),rankdir="RL"){
  X <- x$X; C <- round(x$C,digits); D <- min(ncol(C),degree)

  # --- MODIFIED: Sanitize Names for DOT Node IDs ---
  # Replace spaces, hyphens, and other non-alphanumeric chars (except .) with '_'
  # The original gsub('.', '', fixed=T) for dot removal is kept, but expanded.
  sanitize_dot_id <- function(names) {
    names <- gsub("[^[:alnum:]_.]", "_", names, perl=TRUE)
    return(names)
  }

  # Sanitize names for use as internal DOT Node IDs
  X_names <- sanitize_dot_id(colnames(X))
  Y_names <- sanitize_dot_id(rownames(X))
  C_col_names <- sanitize_dot_id(colnames(C))

  # Original label processing for display purposes
  rownames(X) <- gsub(".","",rownames(X),fixed=T)
  colnames(X) <- gsub(".","",colnames(X),fixed=T)
  rownames(C) <- gsub(".","",rownames(C),fixed=T)
  colnames(C) <- gsub(".","",colnames(C),fixed=T)
  # ---------------------------------------------------

  Alabels <- unique(gsub("_([0-9]+)","",colnames(C),fixed=F))
  index <-match("(Intercept)", Alabels)
  if(!is.na(index)) Alabels <- Alabels[-index]

  # --- 1. Graph Initialization ---
  scr <- paste0('digraph XCA {graph [rankdir=',rankdir,' compound=true]; \n')

  # --- 2. Y Nodes (Observation/Output) ---
  st <- 'subgraph cluster_Y{label="T" style="rounded"; \n'
  for(j in 1:nrow(X)){
    # Use sanitized Y_names as Node ID, and original rownames(X) as display label
    st <- paste0(st,sprintf('  %s [label="%s", shape=box]; \n',
                            Y_names[j],            # Safe Node ID
                            rownames(X)[j]))        # Display Label
  }
  st <- paste0(st,'}; \n'); scr <- paste0(scr,st)

  # --- 3. X Nodes (Latent Variables) ---
  st <- 'subgraph cluster_X{label="Latent Variables" style="rounded"; \n'
  for(j in 1:ncol(X)){
    # Use sanitized X_names as Node ID, and original colnames(X) as display label
    st <- paste0(st,sprintf('  %s [label="%s", shape=ellipse]; \n',
                            X_names[j],            # Safe Node ID
                            colnames(X)[j]))        # Display Label
  }
  st <- paste0(st,'}; \n'); scr <- paste0(scr,st)

  # --- 4. Edge: X to Y ---
  for(i in 1:nrow(X))for(j in 1:ncol(X)){
    if(X[i,j]>=threshold){
      # Use sanitized IDs for edges
      st <- sprintf(paste0('%s -> %s [label="%.',digits,'f"]; \n'),
                    X_names[j],Y_names[i],X[i,j]); scr <- paste0(scr,st)}
  }

  # --- 5. Intercept ---
  if(intercept==TRUE){
    for(i in 1:ncol(X))
      if(C[i,ncol(C)]>=threshold){
        # Use sanitized X_names for target node
        st <- sprintf(paste0('Const%d [shape=circle label="%.',digits,'f"]; '),i,C[i,ncol(C)])
        st <- paste0(st,sprintf(paste0('Const%d -> %s; \n'),i,X_names[i])) # Safe Node ID
        scr <- paste0(scr,st)
      }
  }

  # --- 6. Edge: T-k to X (Covariates) ---
  klist <- NULL; ktoplist <- NULL
  for(k in 1:D){
    Ck <- C[,(k-1)*length(Alabels)+1:length(Alabels)]
    if(is.matrix(Ck)==FALSE){
      Ck <- matrix(Ck,nrow=1)
      colnames(Ck) <- C_col_names[k] # Use sanitized name
      rownames(Ck) <- X_names[1]      # Use sanitized name
    }
    if(max(Ck)>=threshold){
      klist <- c(klist,k)
      st <- sprintf('subgraph cluster_C%d{label="T-%d" style="rounded"; \n',k,k)
      ktop <- NULL
      for(j in 1:ncol(Ck)){
        if(max(Ck[,j])>=threshold){
          if(is.null(ktop)==TRUE){
            ktop <- j
            ktoplist <- c(ktoplist,ktop)
          }
          # Extract display label and sanitize Node ID
          alabel <- gsub("_([0-9]+)","",colnames(C)[j],fixed=F)
          node_id <- C_col_names[(k-1)*length(Alabels)+j] # Safe Node ID

          # QUOTE LABEL: Protect A label in DOT
          st <- paste0(st,sprintf('  %s [label="%s",shape=box]; \n',
                                  node_id,         # Safe Node ID
                                  alabel))          # Display Label
        }
      }
      st <- paste0(st,"}; \n");scr <- paste0(scr,st)
    }
    for(i in 1:nrow(Ck))for(j in 1:ncol(Ck)){
      if(Ck[i,j]>=threshold){
        # Use sanitized IDs for edges
        st <- sprintf(paste0('%s -> %s [label="%.',digits,'f"]; \n'),
                      C_col_names[(k-1)*length(Alabels)+j], X_names[i], Ck[i,j]) # Safe IDs
        scr <- paste0(scr,st)
      }
    }
  }

  # --- 7. Invisible Edges for Cluster Ordering ---
  if(length(klist)>=2){
    for(k in 2:length(klist)){
      # Note: This block uses the original, unsanitized Alabels for label consistency,
      # but should use sanitized IDs for the nodes themselves.
      # Since Alabels are derived from sanitized names previously, we rely on the node IDs being safe.

      # We use the sanitized C_col_names for the actual node IDs:
      start_node_id <- C_col_names[(klist[k]-1)*length(Alabels)+ktoplist[k]]
      end_node_id <- C_col_names[(klist[k-1]-1)*length(Alabels)+ktoplist[k-1]]

      st <- sprintf('%s -> %s [ltail=cluster_C%d lhead=cluster_C%d style=invis]; \n',
                    start_node_id, end_node_id, # Safe Node IDs
                    klist[k], klist[k-1])
      scr <- paste0(scr,st)
    }
  }

  scr <- paste0(scr,"} \n")
  return(scr)
}



#' @title Forecast future values for NMF-VAR model
#' @description
#' \code{nmfkc.ar.predict} computes multi-step ahead forecasts for a fitted NMF-VAR model
#' using recursive forecasting strategies.
#'
#' @param x An object of class \code{nmfkc} (the fitted model).
#' @param Y The historical observation matrix used for fitting (or at least the last \code{degree} columns).
#' @param degree The lag order (D) used in the model. If \code{NULL} (default), it is automatically inferred from the dimensions of \code{x$X} and \code{x$C}.
#' @param n.ahead Integer. The number of steps ahead to forecast.
#'
#' @return A list containing:
#' \item{pred}{A matrix of predicted values with \code{n.ahead} columns.}
#' \item{time}{A numeric vector of future time points, if \code{colnames(Y)} are numeric time values. Otherwise \code{NULL}.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ar}}
#' @export

nmfkc.ar.predict <- function(x, Y, degree=NULL, n.ahead=1){
  if(!is.matrix(Y)) Y <- as.matrix(Y)

  # --- Initial Class Check ---
  if (!inherits(x, "nmfkc")) {
    stop("Argument 'x' must be an object of class 'nmfkc'.")
  }

  C <- x$C
  P <- nrow(Y) # The number of observed variables/features

  # --- 1. Degree & Intercept Detection (Robust, Accessing x$A.attr) ---

  # 1.1. Read parameters from x$A.attr (Prioritized check)
  A.attributes <- x$A.attr

  # Check if A.attr exists and contains VAR metadata
  is_metadata_available <- !is.null(A.attributes) &&
    !is.null(A.attributes$function.name) &&
    A.attributes$function.name == "nmfkc.ar"

  if (is_metadata_available) {
    # CASE A: VAR Metadata found in x$A.attr (Robust)
    degree_from_attr <- A.attributes$degree
    intercept_from_attr <- A.attributes$intercept

    if (is.null(degree_from_attr) || is.null(intercept_from_attr)) {
      stop("VAR model metadata (degree or intercept status) is missing in x$A.attr.")
    }
    degree_final <- degree_from_attr
    has_intercept_final <- intercept_from_attr

    # Check for user-provided 'degree' argument which overrides A's metadata.
    if (!is.null(degree) && degree != degree_final) {
      warning(paste0("The 'degree' argument (", degree, ") provided by the user contradicts the degree stored in the model (", degree_final, "). Using user-provided degree for prediction."))
      degree_final <- degree

      # Fallback to dimension check for the intercept status based on user's degree
      R <- ncol(C)
      if (R == P * degree_final + 1) {
        has_intercept_final <- TRUE
      } else if (R == P * degree_final) {
        has_intercept_final <- FALSE
      } else {
        stop(paste0("Dimension mismatch: ncol(C)=", R,
                    " does not match P*degree+1 (", P*degree_final+1,
                    ") or P*degree (", P*degree_final, ") for the user-provided 'degree'."))
      }
    }

  } else if (!is.null(degree)) {
    # CASE B: Metadata not found or not VAR, but user explicitly provided 'degree'.
    # Fallback to dimension check based on user input.
    degree_final <- degree
    R <- ncol(C)

    if (R == P * degree_final + 1) {
      has_intercept_final <- TRUE
    } else if (R == P * degree_final) {
      has_intercept_final <- FALSE
    } else {
      stop(paste0("Dimension mismatch: ncol(C)=", R,
                  " does not match P*degree+1 (", P*degree_final+1,
                  ") or P*degree (", P*degree_final, "). Check 'degree' argument."))
    }
  } else {
    # CASE C: No metadata and no user input. Fallback to the original dimension check.

    R <- ncol(C)
    P_fit <- nrow(x$X)

    is_intercept_model <- ((R - 1) %% P_fit == 0)
    is_no_intercept_model <- (R %% P_fit == 0)

    if (P_fit == 1) {
      has_intercept_final <- TRUE
      degree_final <- (R - 1) / P_fit
      if(R == 1) { has_intercept_final <- FALSE; degree_final <- 1 }
    } else {
      if (is_intercept_model && !is_no_intercept_model) {
        has_intercept_final <- TRUE
        degree_final <- (R - 1) / P_fit
      } else if (!is_intercept_model && is_no_intercept_model) {
        has_intercept_final <- FALSE
        degree_final <- R / P_fit
      } else {
        # Ambiguous or non-VAR case
        stop("The model dimensions lead to an ambiguous lag order and intercept status. Please ensure the model 'x' was created using nmfkc.ar() and nmfkc(), or manually specify the 'degree' argument.")
      }
    }
  }

  degree <- degree_final
  has_intercept <- has_intercept_final

  # --- 2. Historical Data Check ---
  if(ncol(Y) < degree){
    stop(paste0("Not enough historical data in Y (", ncol(Y), " columns) for the model's lag degree (", degree, "). Y must contain at least ", degree, " columns for the first forecast step."))
  }

  # --- 3. Forecasting Loop (Recursive Prediction) ---
  preds <- matrix(0, nrow=P, ncol=n.ahead)
  colnames(preds) <- paste0("t+", 1:n.ahead)
  rownames(preds) <- rownames(Y)

  current_Y <- Y

  for(step in 1:n.ahead){

    a_vec <- numeric(0)
    # Collect lagged values (A matrix rows)
    for(k in 1:degree){
      col_idx <- ncol(current_Y) - k + 1
      a_vec <- c(a_vec, current_Y[, col_idx])
    }

    if(has_intercept) a_vec <- c(a_vec, 1)
    newA <- matrix(a_vec, ncol=1) # Covariate vector for prediction time t+step

    # Predict next step: X * C * A(t+step)
    pred_step <- x$X %*% x$C %*% newA
    preds[, step] <- pred_step

    # Add predicted value to the history for the next recursive step
    current_Y <- cbind(current_Y, pred_step)
  }

  # --- 4. Time Inference ---
  future_times <- NULL
  if(!is.null(colnames(Y))){
    times <- suppressWarnings(as.numeric(colnames(Y)))

    if(!any(is.na(times)) && length(times) >= 2){
      n_check <- min(length(times), 10)
      recent_times <- utils::tail(times, n_check)

      diffs <- diff(recent_times)

      # Check if time points are equally spaced (tolerance: 1e-5 relative error)
      if(all(abs(diffs - mean(diffs)) < 1e-5 * abs(mean(diffs)) + 1e-8)){
        dt <- mean(diffs)

        last_time <- utils::tail(times, 1)

        future_times <- last_time + (1:n.ahead) * dt
        colnames(preds) <- future_times
      }
    }
  }

  return(list(pred = preds, time = future_times))
}





#' @title Optimize lag order for the autoregressive model
#' @description
#' \code{nmfkc.ar.degree.cv} selects the optimal lag order for an autoregressive model
#' by applying cross-validation over candidate degrees.
#'
#' @param Y Observation matrix \eqn{Y(P,N)}.
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
#' # Example.
#' d <- AirPassengers
#' time <- time(ts(1:length(d),start=c(1949,1),frequency=12))
#' time.vec <- round(as.vector(t(time)),2)
#' Y0 <- matrix(as.vector(d),nrow=1)
#' colnames(Y0) <- time.vec
#' rownames(Y0) <- "t"
#' # selection of degree
#' nmfkc.ar.degree.cv(Y=Y0,Q=1,degree=11:14)

nmfkc.ar.degree.cv <- function(Y,Q=1,degree=1:2,intercept=T,plot=TRUE,...){
  extra_args <- list(...)
  objfuncs <- numeric(length(degree)) # Use numeric(length(degree)) for clean initialization

  # --- FIX 1: Initialize status variable outside the loop ---
  success_status <- logical(length(degree))

  for(i in 1:length(degree)){
    start.time <- Sys.time()
    # --- FIX 1: Define current_degree inside the loop ---
    current_degree <- degree[i]

    message(paste0("degree=",current_degree,"..."),appendLF=FALSE)

    tryCatch({
      # 1. Execute nmfkc.ar (Construct Y and A)
      a <- nmfkc.ar(Y=Y,degree=current_degree,intercept=intercept)

      # 2. Setup arguments for nmfkc.cv
      main_args <- list(Y = a$Y,A = a$A,Q = Q)
      # Use shuffle=FALSE for block CV, suitable for time series data
      all_args <- c(extra_args, main_args, list(shuffle = FALSE))

      # 3. Execute nmfkc.cv
      # Note: We suppress messages from inner calls to keep the log clean
      result.cv <- suppressMessages(do.call("nmfkc.cv", all_args))

      # 4. Record successful results
      objfuncs[i] <- result.cv$objfunc/ncol(a$Y)
      success_status[i] <- TRUE

      # 5. Log processing time
      end.time <- Sys.time()
      diff.time <- difftime(end.time,start.time,units="sec")
      diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                             paste0(round(diff.time/60,1),"min"))
      message(diff.time.st)

    }, error = function(e) {
      # Handle error: Record NA, issue a warning, and continue the loop
      warning(paste0("Skipping degree=", current_degree, " due to error: ", e$message), call. = FALSE)
      objfuncs[i] <- NA
      success_status[i] <- FALSE
      message("Skipped (Error)")
    })
  }

  # --- FIX 2: Post-Processing executed ONCE outside the loop ---
  valid_indices <- which(success_status)
  if (length(valid_indices) == 0) {
    stop("Cross-validation failed for all candidate degrees.")
  }

  # Find the minimum objective function among successful runs
  i0 <- valid_indices[which.min(objfuncs[valid_indices])]
  best.degree <- degree[i0]

  # Calculate degree.max (empirical rule)
  degree.max <- min(ncol(Y),floor(10*log10(ncol(Y))))

  # Plotting (Only successful runs are plotted/labeled)
  if(plot){
    # Plot objective function values
    plot(degree,objfuncs,type="l",col=2,xlab=paste0("degree (max=",degree.max,")"),ylab="objfunc")
    graphics::points(degree[valid_indices],objfuncs[valid_indices],cex=1,col=2)

    # Highlight the best degree
    graphics::points(degree[i0],objfuncs[i0],cex=3,col=2)
    graphics::text(degree,objfuncs,degree, pos=3)
  }

  names(objfuncs) <- degree
  result <- list(degree=best.degree,degree.max=degree.max,objfunc=objfuncs)
  return(result)
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
  message(paste0("P=",P,",Q=",Q,",D=",D,",intercept=",ifelse(has_intercept,"T","F")))
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
  attr(K, "function") <- "nmfkc.kernel"
  # ------------------------------------------------------

  return(K)
}



#' @title Estimate kernel parameter beta from covariates
#' @description
#' \code{nmfkc.kernel.beta.nearest.med} estimates the Gaussian kernel
#' parameter \eqn{\beta} by computing the median of nearest-neighbor
#' distances among covariates. This is useful for setting the scale
#' parameter in kernel-based NMF with covariates.
#'
#' @param U covariate matrix \eqn{U(K,N)=(u_1,\dots,u_N)},
#'   where each column corresponds to an individual. Each row may be
#'   normalized in advance.
#' @param block_size number of samples to process at once.
#'   If \eqn{N \le 1000}, it is automatically set to \eqn{N}.
#'
#' @return A list with components:
#' \item{beta}{estimated kernel parameter \eqn{\beta=1/(2 d_{med}^2)}}
#' \item{beta_candidates}{a numeric vector of candidate values obtained by
#'   multiplying the estimate \eqn{\beta} by powers of 10,
#'   i.e.\ \eqn{\{\beta \cdot 10^{-2},\,\beta \cdot 10^{-1},\,\beta,\,\beta \cdot 10^{1}\}}}
#' \item{dist_median}{the median nearest-neighbor distance}
#' \item{block_size_used}{actual block size used in computation}
#' @details
#' The function computes all pairwise squared distances between columns of
#' \eqn{U}, excludes self-distances, and takes the median of the nearest-neighbor
#' distances (after square root). This median is then used to set \eqn{\beta}.
#' @seealso \code{\link{nmfkc.kernel}}, \code{\link{nmfkc.kernel.beta.cv}}
#' @export

nmfkc.kernel.beta.nearest.med <- function(U, block_size=1000){
  U <- as.matrix(U)
  N <- ncol(U)
  X <- t(U)
  if (N <= 1000) block_size <- N
  XX <- rowSums(X * X)
  min_d2 <- rep(Inf, N)
  for (i in seq(1, N, by = block_size)) {
    i2 <- min(i + block_size - 1, N)
    Xi <- X[i:i2, , drop = FALSE]
    Xi_norm <- rowSums(Xi * Xi)
    dist2 <- outer(Xi_norm, rep(1, N)) +
      outer(rep(1, nrow(Xi)), XX) -
      2 * Xi %*% t(X)
    idx <- i:i2
    dist2[cbind(seq_along(idx), idx)] <- Inf
    dist2[dist2 < 0] <- 0
    nn_local <- apply(dist2, 1, min)
    min_d2[idx] <- pmin(min_d2[idx], nn_local)
    rm(Xi, Xi_norm, dist2); gc(FALSE)
  }
  d_med <- stats::median(sqrt(min_d2))
  beta  <- 1 / (2 * d_med^2)
  beta_candidates <- beta*10^c(-2:1)
  list(beta = beta, beta_candidates=beta_candidates, dist_median = d_med, block_size_used = block_size)
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
#'       Enforces columns of \eqn{X} to be orthogonal (conceptually \eqn{X^\top X \to I}). (Formerly \code{lambda.ortho}).
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
      Y.weights <- 1
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
    if(min(nrow(Y),ncol(Y))>=Q){
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
            # Use Y_init
            X <- .nndsvdar(Y_init, Q)
          }
        }
      }
    }else{
      stop("It does not hold Q<=min(P,N).")
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
    is.X.scalar <- TRUE
  }
  X <- xnorm(X)

  if(is.null(A)) C <- matrix(1,nrow=ncol(X),ncol=ncol(Y)) else C <- matrix(1,nrow=ncol(X),ncol=nrow(A))
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
      }
      if(is.null(A)) {
        num_C <- t(X) %*% (Y.weights * Y)
        den_C <- t(X) %*% (Y.weights * XB)
        if (C.L1 != 0) den_C <- den_C + (C.L1/2) * ones_QN
        if (B.L1 != 0) den_C <- den_C + (B.L1/2) * ones_QN
        C <- C * (num_C / (den_C + .eps))
      } else {
        num_C <- t(X) %*% (Y.weights * Y) %*% At
        den_C <- t(X) %*% (Y.weights * XB) %*% At
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
      }
      if(is.null(A)) {
        ratio <- Y.weights * (Y / (XB + .eps))
        num_C <- t(X) %*% ratio
        den_C <- t(X) %*% Y.weights
        if (C.L1 != 0) den_C <- den_C + C.L1 * ones_QN
        if (B.L1 != 0) den_C <- den_C + B.L1 * ones_QN
        C <- C * (num_C / (den_C + .eps))
      } else {
        ratio <- Y.weights * (Y / (XB + .eps))
        num_C <- t(X) %*% ratio %*% At
        den_C <- t(X) %*% Y.weights %*% At
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
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
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

  if(method=="EU"){
    sigma2 <- objfunc / N_obs
    ICp <- log(objfunc/N_obs) + nparam * log(N_obs)
    AIC <- N_obs * log(sigma2) + 2*nparam
    BIC <- N_obs * log(sigma2) + nparam * log(N_obs)
  }else{
    ICp <- NA; AIC <- 2*objfunc + 2*nparam; BIC <- 2*objfunc + nparam*log(N_obs)
  }

  if(save.memory==FALSE){
    # [Fix 3] Calculate statistics using valid data only
    valid_idx <- (Y.weights > 0)
    if(any(valid_idx)) {
      r2 <- stats::cor(XB[valid_idx], Y[valid_idx])^2
      sigma <- stats::sd(Y[valid_idx] - XB[valid_idx])
    } else { r2 <- NA; sigma <- NA }

    B.prob <- t( t(B) / (colSums(B) + .eps) )
    B.prob.sd.min <- min(apply(B.prob,1,stats::sd))
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
    r2 <- NA; sigma <- NA; B.prob <- NA; B.prob.sd.min <- NA; B.cluster <- NA
    XB <- NA; X.prob <- NA; X.cluster <- NA; silhouette <- NA; CPCC <- NA; dist.cor <- NA
  }

  if(epsilon.iter > abs(epsilon)) warning(paste0("maximum iterations (",maxit,") reached..."))
  end.time <- Sys.time()
  diff.time.st <- paste0(round(difftime(end.time,start.time,units="sec"),1),"sec")
  if(print.dims) message(diff.time.st)

  n.missing <- sum(Y.weights == 0)
  n.total <- prod(dim(Y))
  A.attr <- NULL
  if (!is.null(A)) A.attr <- attributes(A)

  result <- list(call=match.call(), dims=dims, runtime=diff.time.st,
                 X=X, B=B, XB=XB, C=C,
                 B.prob=B.prob, B.cluster=B.cluster,
                 X.prob=X.prob, X.cluster=X.cluster,
                 A.attr=A.attr,
                 n.missing = n.missing,n.total = n.total,rank = Q,
                 objfunc=objfunc, objfunc.iter=objfunc.iter, r.squared=r2, sigma=sigma,
                 criterion=list(B.prob.sd.min=B.prob.sd.min, ICp=ICp, AIC=AIC, BIC=BIC,
                                silhouette=silhouette, CPCC=CPCC, dist.cor=dist.cor))
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

  ans$iter <- length(object$objfunc.iter)
  ans$objfunc <- object$objfunc
  ans$r.squared <- object$r.squared
  ans$sigma <- object$sigma

  if(!is.null(object$criterion)){
    ans$aic <- object$criterion$AIC
    ans$bic <- object$criterion$BIC
  } else {
    ans$aic <- NULL; ans$bic <- NULL
  }

  # --- Diagnostics (Sparsity & Clustering Quality) ---

  # 1. Basis (X)
  if (!is.null(object$X) && is.matrix(object$X)) {
    # Sparsity: Proportion of elements close to zero (< 1e-4)
    ans$X.sparsity <- mean(object$X < 1e-4)
  }

  # 2. Probabilities (B.prob)
  if (!is.null(object$B.prob) && is.matrix(object$B.prob)) {
    # Sparsity
    ans$B.prob.sparsity <- mean(object$B.prob < 1e-4)

    # Entropy (Normalized): 0=Crisp, 1=Uniform
    Q <- nrow(object$B.prob)
    if(Q > 1){
      p <- object$B.prob + 1e-10 # avoid log(0)
      entropies <- -colSums(p * log(p)) / log(Q)
      ans$B.prob.entropy <- mean(entropies)
      # Crispness: Mean Maximum Probability (Closer to 1 is better)
      ans$B.prob.crispness <- mean(apply(object$B.prob, 2, max))
    } else {
      ans$B.prob.entropy <- 0
      ans$B.prob.crispness <- 1
    }
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
  cat("Iterations: ", x$iter, "\n")

  if (!is.null(x$n.missing)) {
    cat("Missing:    ", x$n.missing,
        sprintf("(%.1f%%)", x$prop.missing), "\n")
  }

  cat("\nStatistics:\n")
  cat("  Objective function: ", format(x$objfunc, digits = digits), "\n")
  cat("  Multiple R-squared: ", format(x$r.squared, digits = digits), "\n")
  cat("  Residual Std Error: ", format(x$sigma, digits = digits), "\n")

  if (!is.null(x$aic)) {
    cat("  AIC:                ", format(x$aic, digits = digits), "\n")
    cat("  BIC:                ", format(x$bic, digits = digits), "\n")
  }

  cat("\nStructure Diagnostics:\n")
  if (!is.null(x$X.sparsity)) {
    cat("  Basis (X) Sparsity:   ", sprintf("%.1f%%", x$X.sparsity * 100), "(< 1e-4)\n")
  }
  if (!is.null(x$B.prob.sparsity)) {
    cat("  Coef (B) Sparsity:    ", sprintf("%.1f%%", x$B.prob.sparsity * 100), "(< 1e-4)\n")
    if(!is.null(x$B.prob.entropy)){
      cat("  Clustering Entropy:   ", format(x$B.prob.entropy, digits = 3), "(0=Crisp, 1=Ambiguous)\n")
      cat("  Clustering Crispness: ", format(x$B.prob.crispness, digits = 3), "(Mean Max Prob, >0.8 is good)\n")
    }
  }
  cat("\n")
  invisible(x)
}





#' @title Prediction method for objects of class \code{nmfkc}
#' @description
#' \code{predict.nmfkc} generates predictions from an object of class \code{nmfkc},
#' either using the fitted covariates or a new covariate matrix.
#'
#' @param x An object of class \code{nmfkc}, i.e., the return value of \code{nmfkc}.
#' @param newA Optional. A new covariate matrix to be used for prediction.
#'   If \code{NULL} (default), the fitted covariates are used.
#' @param type Type of prediction to return. Options are:
#'   \itemize{
#'     \item \code{"response"} (default): Returns the reconstructed matrix \eqn{X B}.
#'     \item \code{"prob"}: Returns probabilities using \code{B.prob} instead of \code{B}.
#'     \item \code{"class"}: Returns the class label corresponding to the maximum in each column of \code{B.prob}.
#'   }
#' @seealso \code{\link{nmfkc}}
#' @export
predict.nmfkc <- function(x,newA=NULL,type="response"){
  # A small constant for numerical stability to prevent division by zero and log(0).
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





#' @title Generate DOT language scripts for NMF models
#' @description
#' \code{nmfkc.DOT} generates scripts in the DOT language for visualizing
#' the NMF model structure (\eqn{Y \approx X C A}) or its simplified forms.
#'
#' @param x The return value of \code{nmfkc}.
#' @param type Character string specifying the visualization type.
#'   Options are:
#'   \itemize{
#'     \item \code{"YX"} (Default): Standard NMF view: \eqn{B \to X \to Y}.
#'     \item \code{"YXA"} : Full visualization of the tri-factorization \eqn{A \to C \to X \to Y}.
#'     \item \code{"YA"}: Direct regression view: \eqn{A \to Y}, where coefficients are from \eqn{X C}.
#'   }
#' @param digits Integer. Number of decimal places to display in edge labels. Default is 2.
#' @param threshold Numeric. Parameters greater than or equal to this threshold are displayed. Default is \eqn{10^{-\code{digits}}}.
#' @param rankdir Graph layout direction in DOT language. Default is "LR".
#' @param Y.label Character vector for row names of Y/X (features). If NULL, uses \code{rownames(x$X)}.
#' @param X.label Character vector for column names of X/rows of B (latent factors). If NULL, uses \code{colnames(x$X)}.
#' @param A.label Character vector for row names of A/columns of C (covariates). If NULL, uses \code{colnames(x$C)}.
#' @param Y.title Title for the Y node cluster. Default is "Observation (Y)".
#' @param X.title Title for the X node cluster. Default is "Basis (X)".
#' @param A.title Title for the A node cluster. Default is "Covariates (A)".
#' @param min.penwidth Numeric. Minimum line thickness for the path (default: 1.0).
#' @param max.penwidth Numeric. Maximum line thickness for the path (default: 5.0).
#'
#' @return A character string containing a DOT script, suitable for use with the \pkg{DOT} package or Graphviz tools.
#' @seealso \code{\link{nmfkc}}
#' @export
nmfkc.DOT <- function(x, type = c("YX","YA","YXA"), digits = 2, threshold = 10^(-digits), rankdir = "LR",
                      Y.label = NULL, X.label = NULL, A.label = NULL,
                      Y.title = "Observation (Y)", X.title = "Basis (X)", A.title = "Covariates (A)",
                      min.penwidth = 1.0, max.penwidth = 5.0) {
  type <- match.arg(type)
  X <- x$X
  B <- x$B
  hasA <- !is.null(x$C) && ncol(x$C) != ncol(B)
  P <- nrow(X)
  Q <- ncol(X)

  calculate_penwidth <- function(coeff_matrix, value) {
    if (is.null(coeff_matrix)) return(min.penwidth)
    max_coeff <- max(coeff_matrix, na.rm = TRUE)
    if (max_coeff <= threshold) return(min.penwidth)
    penwidth <- min.penwidth +
      (max.penwidth - min.penwidth) * (value - threshold) / (max_coeff - threshold)
    return(max(min.penwidth, penwidth))
  }

  # --- NEW: Sanitize Names for DOT Node IDs (Problem 7 Fix) ---
  # Replace non-alphanumeric chars (except underscore and dot) with '_'
  sanitize_dot_id <- function(names) {
    names <- gsub("[^[:alnum:]_.]", "_", names, perl=TRUE)
    return(names)
  }

  # Sanitize names for use as internal DOT Node IDs
  Y_names_id <- sanitize_dot_id(if (is.null(Y.label)) rownames(X) else Y.label)
  X_names_id <- sanitize_dot_id(if (is.null(X.label)) colnames(X) else X.label)

  # Use user-provided labels or matrix names for *display*
  Y_labels_display <- if (is.null(Y.label)) rownames(X) else Y.label
  X_labels_display <- if (is.null(X.label)) colnames(X) else X.label
  # ---------------------------------------------------

  # --- 1. Label Assignment & Coefficient Setup ---
  if (hasA) {
    C <- round(x$C, digits)
    A_cols_NMF <- ncol(C)
    A_labels_display <- if (is.null(A.label)) colnames(C) else A.label
    A_names_id <- sanitize_dot_id(A_labels_display) # Sanitize A names too
    XC_mat <- X %*% C
  } else if (type == "YX") {
    C <- B
    A_cols_NMF <- ncol(C) # Set A_cols_NMF to prevent error later
    A_labels_display <- NULL
    XC_mat <- NULL
    A_names_id <- NULL
  } else {
    stop("The model structure (A is missing) is incompatible with the selected type ('YXA' or 'YA').")
  }

  # --- 2. Graph Initialization ---
  scr <- paste0('digraph NMF {graph [rankdir=', rankdir, ' compound=true]; \n')

  # --- 3. Define Y Nodes (Observation/Output) ---
  # Y.title
  st <- paste0('subgraph cluster_Y{label="', Y.title, '" style="rounded"; \n')
  for (j in 1:P) {
    # Use Y%d as Node ID and Y_labels_display[j] as display label
    st <- paste0(st, sprintf('  Y%d [label="%s", shape=box]; \n',
                             j,                   # Numerical Node ID
                             Y_labels_display[j])) # Display Label
  }
  st <- paste0(st, '}; \n'); scr <- paste0(scr, st)

  # --- 3. Define X Nodes (Latent Variables/Basis) ---
  st <- paste0('subgraph cluster_X{label="', X.title, '" style="rounded"; \n')
  for (j in 1:Q) {
    st <- paste0(st, sprintf('  X%d [label="%s", shape=ellipse]; \n',
                             j,                   # Numerical Node ID
                             X_labels_display[j])) # Display Label
  }
  st <- paste0(st, '}; \n'); scr <- paste0(scr, st)

  # --- 4. Draw Paths based on 'type' ---
  if (type == "YXA") {
    # 4.1. Define A Nodes (Covariates/Input)
    st <- paste0('subgraph cluster_A{label="', A.title, '" style="rounded"; \n')
    for (j in 1:A_cols_NMF) {
      st <- paste0(st, sprintf('  A%d [label="%s", shape=box]; \n',
                               j,                   # Numerical Node ID
                               A_labels_display[j])) # Display Label
    }
    st <- paste0(st, '}; \n'); scr <- paste0(scr, st)

    # 4.2. Edges from A to X (via C) - Using numerical IDs
    for (i in 1:Q) {
      for (j in 1:A_cols_NMF) {
        coeff_val <- C[i, j]
        if (coeff_val >= threshold) {
          penwidth <- calculate_penwidth(C, coeff_val)
          st <- sprintf(paste0('A%d -> X%d [label="%.', digits, 'f", penwidth=', penwidth, ']; \n'),
                        j, i, coeff_val)
          scr <- paste0(scr, st)
        }
      }
    }
    # 4.3. Edges from X to Y (Basis contribution)
    max_X <- max(X, na.rm = TRUE)
    for (i in 1:P) {
      for (j in 1:Q) {
        coeff_val <- X[i, j]
        if (coeff_val >= threshold) {
          penwidth <- calculate_penwidth(X, coeff_val)
          st <- sprintf(paste0('X%d -> Y%d [label="%.', digits, 'f", penwidth=', penwidth, ']; \n'),
                        j, i, coeff_val)
          scr <- paste0(scr, st)
        }
      }
    }
  } else if (type == "YA") {
    # 4.1. Define A Nodes (Covariates/Input)
    st <- paste0('subgraph cluster_A{label="', A.title, '" style="rounded"; \n')
    for (j in 1:A_cols_NMF) {
      st <- paste0(st, sprintf('  A%d [label="%s", shape=box]; \n',
                               j,                   # Numerical Node ID
                               A_labels_display[j])) # Display Label
    }
    st <- paste0(st, '}; \n'); scr <- paste0(scr, st)

    # 4.2. Edges from A to Y (via XC)
    max_XC <- max(XC_mat, na.rm = TRUE)
    for (i in 1:P) {
      for (j in 1:A_cols_NMF) {
        coeff_val <- XC_mat[i, j]
        if (coeff_val >= threshold) {
          penwidth <- calculate_penwidth(XC_mat, coeff_val)
          st <- sprintf(paste0('A%d -> Y%d [label="%.', digits, 'f", penwidth=', penwidth, ']; \n'),
                        j, i, coeff_val)
          scr <- paste0(scr, st)
        }
      }
    }
  } else if (type == "YX") {
    # 4.1. Define X Nodes (Latent Variables/Basis)
    st <- paste0('subgraph cluster_X{label="', X.title, '" style="rounded"; \n')
    for (j in 1:Q) {
      st <- paste0(st, sprintf('  X%d [label="%s", shape=ellipse]; \n',
                               j,                   # Numerical Node ID
                               X_labels_display[j])) # Display Label
    }
    st <- paste0(st, '}; \n'); scr <- paste0(scr, st)

    # 4.2. Edges from X to Y (Basis contribution)
    max_X <- max(X, na.rm = TRUE)
    for (i in 1:P) {
      for (j in 1:Q) {
        coeff_val <- X[i, j]
        if (coeff_val >= threshold) {
          penwidth <- calculate_penwidth(X, coeff_val)
          st <- sprintf(paste0('X%d -> Y%d [label="%.', digits, 'f", penwidth=', penwidth, ']; \n'),
                        j, i, coeff_val)
          scr <- paste0(scr, st)
        }
      }
    }
  }
  # --- 5. Graph Finalization ---
  scr <- paste0(scr, "} \n")
  return(scr)
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
  colnames(X) <- NULL
  X
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


#' @title Perform k-fold cross-validation for NMF with kernel covariates
#' @description
#' \code{nmfkc.cv} performs k-fold cross-validation on the model
#' \eqn{Y \approx X C A = X B}, where
#' \itemize{
#' 	\item \eqn{Y(P,N)} is the observation matrix,
#' 	\item \eqn{A(R,N)} is the covariate matrix,
#' 	\item \eqn{X(P,Q)} is the basis matrix with \eqn{Q \le \min(P,N)},
#' 	\item \eqn{C(Q,R)} is the parameter matrix, and
#' 	\item \eqn{B(Q,N)} is the coefficient matrix (\eqn{B = C A}).
#' }
#' Given \eqn{Y} (and optionally \eqn{A}), \eqn{X} and \eqn{C} are estimated,
#' and the predictive performance is assessed by cross-validation.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#' @param Q Rank of the basis matrix \eqn{X}; must satisfy \eqn{Q \le \min(P,N)}.
#' @param ... Additional arguments passed for controlling cross-validation setup and fine-tuning the internal \code{\link{nmfkc}} calls.
#'   These include:
#'   \itemize{
#'     \item \code{Y.weights}: Optional numeric matrix or vector. 0 indicates missing/ignored values.
#'     \item \code{div}: Number of folds (\eqn{k}) for cross-validation (default: 5).
#'     \item \code{seed}: Integer seed for reproducibility of data partitioning (default: 123).
#'     \item \code{shuffle}: Logical. If \code{TRUE} (default), samples are shuffled randomly (Standard CV). If \code{FALSE}, samples are split sequentially (Block CV), recommended for time series.
#'     \item **Arguments passed to \code{\link{nmfkc}}**: \code{gamma}, \code{epsilon}, \code{maxit}, \code{method}, \code{X.restriction}, \code{X.init}, etc.
#'   }
#'
#' @return A list with components:
#' \item{objfunc}{Total objective function value across all folds.}
#' \item{objfunc.block}{Objective function values for each fold.}
#' \item{block}{Vector of block indices (1, …, \code{div}) assigned to each column of \eqn{Y}.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.kernel.beta.cv}}, \code{\link{nmfkc.ar.degree.cv}}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example 1.
#' Y <- matrix(cars$dist,nrow=1)
#' A <- rbind(1,cars$speed)
#' result <- nmfkc.cv(Y,A,Q=1)
#' result$objfunc
#'
#' # Example 2.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' betas <- 25:35/1000
#' objfuncs <- 0*(1:length(betas))
#' for(i in 1:length(betas)){
#'   print(i)
#'   A <- nmfkc.kernel(U,V,beta=betas[i])
#'   result <- nmfkc.cv(Y,A,Q=1,div=10)
#'   objfuncs[i] <- result$objfunc
#' }
#' min(objfuncs)
#' (beta.best <- betas[which.min(objfuncs)])
#' # objective function by beta
#' plot(betas,objfuncs,type="o",log="x")
#' table(result$block) # partition block of cv
#' # fitted curve
#' A <- nmfkc.kernel(U,V,beta=beta.best)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

#' @title Perform k-fold cross-validation for NMF with kernel covariates
#' @description
#' \code{nmfkc.cv} performs k-fold cross-validation on the model
#' \eqn{Y \approx X C A = X B}, where
#' \itemize{
#' 	\item \eqn{Y(P,N)} is the observation matrix,
#' 	\item \eqn{A(R,N)} is the covariate matrix,
#' 	\item \eqn{X(P,Q)} is the basis matrix with \eqn{Q \le \min(P,N)},
#' 	\item \eqn{C(Q,R)} is the parameter matrix, and
#' 	\item \eqn{B(Q,N)} is the coefficient matrix (\eqn{B = C A}).
#' }
#' Given \eqn{Y} (and optionally \eqn{A}), \eqn{X} and \eqn{C} are estimated,
#' and the predictive performance is assessed by cross-validation.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#' @param Q Rank of the basis matrix \eqn{X}; must satisfy \eqn{Q \le \min(P,N)}.
#' @param ... Additional arguments passed for controlling cross-validation setup and fine-tuning the internal \code{\link{nmfkc}} calls.
#'   These include:
#'   \itemize{
#'     \item \code{div}: Number of folds (\eqn{k}) for cross-validation (default: 5).
#'     \item \code{seed}: Integer seed for reproducibility of data partitioning (default: 123).
#'     \item \code{shuffle}: Logical. If \code{TRUE} (default), samples are shuffled randomly (Standard CV). If \code{FALSE}, samples are split sequentially (Block CV), recommended for time series.
#'     \item **Arguments passed to \code{\link{nmfkc}}**: \code{gamma}, \code{epsilon}, \code{maxit}, \code{method}, \code{X.restriction}, \code{X.init}, etc.
#'   }
#'
#' @return A list with components:
#' \item{objfunc}{Total objective function value across all folds.}
#' \item{sigma}{Residual standard error (RMSE). Only available if method="EU". Matches the scale of Y and is consistent with `nmfkc()` output.}
#' \item{objfunc.block}{Objective function values for each fold.}
#' \item{block}{Vector of block indices (1, …, \code{div}) assigned to each column of \eqn{Y}.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.kernel.beta.cv}}, \code{\link{nmfkc.ar.degree.cv}}
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
        epsilon.iter <- abs(newSum-oldSum)/(abs(newSum)+0.1)
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
    A.function <- attr(A, "function")
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
#' \item{objfunc.fold}{List of length equal to Q. Each element contains the MSE values for the k folds.}
#' \item{folds}{A list of length \code{div}, containing the linear indices of held-out elements for each fold (shared across all Q).}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.cv}}
#' @export
nmfkc.ecv <- function(Y, A=NULL, Q=1:3, div=5, seed=123, ...){

  if(!is.matrix(Y)) Y <- as.matrix(Y)
  P <- nrow(Y)
  N <- ncol(Y)

  # 1. Create Folds (Shared across all Q for fair comparison)
  if (!is.null(seed)) set.seed(seed)

  # Indices of observed values (non-NA)
  valid_indices <- which(!is.na(Y))
  n_valid <- length(valid_indices)

  perm_indices <- sample(valid_indices)
  folds <- vector("list", div)

  # Divide indices into chunks
  chunk_size <- n_valid %/% div
  remainder <- n_valid %% div

  start_idx <- 1
  for(k in 1:div){
    current_size <- chunk_size + ifelse(k <= remainder, 1, 0)
    end_idx <- start_idx + current_size - 1
    folds[[k]] <- perm_indices[start_idx:end_idx]
    start_idx <- end_idx + 1
  }

  # Retrieve method to check if RMSE is appropriate
  extra_args <- list(...)
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
    # message(paste0("  Evaluating Q=", q_curr, "..."), appendLF = FALSE)

    objfunc.fold <- numeric(div)

    for(k in 1:div){
      test_idx <- folds[[k]]

      # Create Weights for Training
      weights_train <- matrix(1, nrow=P, ncol=N)
      if(any(is.na(Y))) weights_train[is.na(Y)] <- 0
      weights_train[test_idx] <- 0 # Mask current fold

      # Fit NMF
      fit <- suppressMessages(nmfkc(Y=Y, A=A, Q=q_curr, Y.weights=weights_train, ...))

      # Evaluate on Test Indices ONLY
      pred <- fit$XB
      residuals <- Y[test_idx] - pred[test_idx]

      # MSE Calculation
      objfunc.fold[k] <- mean(residuals^2)
    }

    # Store results for this Q
    result_fold[[i]]  <- objfunc.fold
    result_objfunc[i] <- mean(objfunc.fold)
    result_sigma[i]   <- if(method == "EU") sqrt(result_objfunc[i]) else NA

    # message(" Done.")
  }

  return(list(objfunc = result_objfunc,
              sigma = result_sigma,
              objfunc.fold = result_fold,
              folds = folds))
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

      fit <- suppressMessages(do.call("nmfkc", nmfkc_args))

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

    graphics::lines(results_df$rank, results_df$B.prob.sd.min, col=3, lwd=2)
    legend_txt <- c(legend_txt, "B.prob.sd.min"); legend_col <- c(legend_col, 3); legend_lty <- c(legend_lty, 1)

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

    graphics::legend("bottomright", legend=legend_txt, col=legend_col, lty=legend_lty, lwd=2, bg="white", cex=0.7)
  }

  return(list(rank.best = rank.final, criteria = results_df))
}





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


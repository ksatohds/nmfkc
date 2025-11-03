.onAttach <- function(...) {
  packageStartupMessage("Last update on 03 NOV 2025")
  packageStartupMessage("https://github.com/ksatohds/nmfkc")
}

# internal-utils.R
#' @keywords internal
#' @noRd
.z <- function(x){
  x[is.nan(x)] <- 0
  x[is.infinite(x)] <- 0
  x[x < 0] <- 0
  x
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
  N <- ncol(Y)
  A.columns <- NULL
  for(i in degree:1)A.columns <- rbind(A.columns,i:(i+ncol(Y)-degree-1))
  A <- NULL
  for(i in 1:nrow(A.columns))A <- rbind(A,Y[,A.columns[i,]])
  if(is.null(rownames(Y))) rownames(Y) <- 1:nrow(Y)
  label <- NULL
  for(i in 1:degree)label <- c(label,paste0(rownames(Y),"_",i))
  if(intercept){
    A <- rbind(A,1)
    rownames(A) <- c(label,"(Intercept)")
  }
  Ya <- Y[,A.columns[1,]+1]
  if(!is.matrix(Ya)){
    Ya <- matrix(Ya,nrow=1)
    colnames(Ya) <- colnames(Y)[A.columns[1,]+1]
    rownames(Ya) <- rownames(Y)[1]
  }
  degree.max <- min(ncol(Ya),floor(10*log10(ncol(Ya))))
  list(Y = Ya, A = A, A.columns = A.columns, degree.max = degree.max)
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

nmfkc.ar.DOT <- function(x,degree=1,intercept=FALSE,digits=1,threshold=10^(-digits),rankdir="RL"){
  X <- x$X; C <- round(x$C,digits); D <- min(ncol(C),degree)
  rownames(X) <- gsub(".","",rownames(X),fixed=T)
  colnames(X) <- gsub(".","",colnames(X),fixed=T)
  rownames(C) <- gsub(".","",rownames(C),fixed=T)
  colnames(C) <- gsub(".","",colnames(C),fixed=T)
  Alabels <- unique(gsub("_([0-9]+)","",colnames(C),fixed=F))
  index <-match("(Intercept)", Alabels)
  if(!is.na(index)) Alabels <- Alabels[-index]
  # rankdir=RL # rankdir=TB
  scr <- paste0('digraph XCA {graph [rankdir=',rankdir,' compound=true]; \n')
  # Y
  st <- 'subgraph cluster_Y{label="T" style="rounded"; \n'
  for(j in 1:nrow(X))st <- paste0(st,sprintf('  %s [shape=box]; \n',rownames(X)[j]))
  st <- paste0(st,'}; \n'); scr <- paste0(scr,st)
  # X and element
  st <- 'subgraph cluster_X{label="Latent Variables" style="rounded"; \n'
  for(j in 1:ncol(X))st <- paste0(st,sprintf('  %s [shape=ellipse]; \n',colnames(X)[j]))
  st <- paste0(st,'}; \n'); scr <- paste0(scr,st)
  # edge: X to Y
  for(i in 1:nrow(X))for(j in 1:ncol(X)){
    if(X[i,j]>=threshold){
      st <- sprintf(paste0('%s -> %s [label="%.',digits,'f"]; \n'),
              colnames(X)[j],rownames(X)[i],X[i,j]); scr <- paste0(scr,st)}
  }
  # intercept
  if(intercept==TRUE){
    for(i in 1:ncol(X))
      if(C[i,ncol(C)]>=threshold){
        st <- sprintf(paste0('Const%d [shape=circle label="%.',digits,'f"]; '),i,C[i,ncol(C)])
        st <- paste0(st,sprintf(paste0('Const%d -> %s; \n'),i,colnames(X)[i]))
        scr <- paste0(scr,st)
      }
  }
  # edge: T-k to X
  klist <- NULL; ktoplist <- NULL
  for(k in 1:D){
    Ck <- C[,(k-1)*length(Alabels)+1:length(Alabels)]
    if(is.matrix(Ck)==FALSE){
      Ck <- matrix(Ck,nrow=1)
      colnames(Ck) <- colnames(C)[k]
      rownames(Ck) <- colnames(X)[1]
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
          alabel <- gsub("_([0-9]+)","",colnames(Ck)[j],fixed=F)
          st <- paste0(st,sprintf('  %s [label="%s",shape=box]; \n',
                                  colnames(Ck)[j],alabel))
        }
      }
      st <- paste0(st,"}; \n");scr <- paste0(scr,st)
    }
    for(i in 1:nrow(Ck))for(j in 1:ncol(Ck)){
      if(Ck[i,j]>=threshold){
        st <- sprintf(paste0('%s -> %s [label="%.',digits,'f"]; \n'),
                      colnames(Ck)[j],rownames(Ck)[i],Ck[i,j])
        scr <- paste0(scr,st)
      }
    }
  }
  if(length(klist)>=2){
    for(k in 2:length(klist)){
      st <- sprintf('%s_%d -> %s_%d [ltail=cluster_C%d lhead=cluster_C%d style=invis]; \n',
                    Alabels[ktoplist[k]],klist[k],
                    Alabels[ktoplist[k-1]],klist[k-1],
                          klist[k],klist[k-1])
      scr <- paste0(scr,st)
    }
  }
  scr <- paste0(scr,"} \n")
  return(scr)
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
  objfuncs <- 0*(1:length(degree))
  for(i in 1:length(degree)){
    start.time <- Sys.time()
    message(paste0("degree=",degree[i],"..."),appendLF=FALSE)
    a <- nmfkc.ar(Y=Y,degree=degree[i],intercept=intercept)
    main_args <- list(Y = a$Y,A = a$A,Q = Q)
    all_args <- c(extra_args, main_args)
    result.cv <- do.call("nmfkc.cv", all_args)
    objfuncs[i] <- result.cv$objfunc/ncol(a$Y)
    end.time <- Sys.time()
    diff.time <- difftime(end.time,start.time,units="sec")
    diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                           paste0(round(diff.time/60,1),"min"))
    message(diff.time.st)
  }
  i0 <- which.min(objfuncs)
  best.degree <- degree[i0]
  degree.max <- min(ncol(Y),floor(10*log10(ncol(Y))))
  if(plot){
    plot(degree,objfuncs,type="l",col=2,xlab=paste0("degree (max=",degree.max,")"),ylab="objfunc")
    graphics::points(degree[i0],objfuncs[i0],cex=3,col=2)
    graphics::text(degree,objfuncs,degree)
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
  K
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


#' @title Initialize W (X) matrix using NNDSVD (Vectorized)
#' @description
#'   Internal function to compute the NNDSVD initialization for the basis matrix X.
#'   This version is fully vectorized and removes the for-loop over Q.
#' @param Y Input matrix (P x N)
#' @param Q Rank (number of components)
#' @return X (P x Q) non-negative initial basis matrix
#' @keywords internal
#' @noRd
.nndsvd <- function(Y, Q) {
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
    U_p <- U; U_p[U_p < 0] <- 0
    U_n <- U; U_n[U_n > 0] <- 0; U_n <- -U_n
    V_p <- V; V_p[V_p < 0] <- 0
    V_n <- V; V_n[V_n > 0] <- 0; V_n <- -V_n
    norm_Up <- sqrt(colSums(U_p^2))
    norm_Un <- sqrt(colSums(U_n^2))
    norm_Vp <- sqrt(colSums(V_p^2))
    norm_Vn <- sqrt(colSums(V_n^2))
    Mp <- norm_Up * norm_Vp
    Mn <- norm_Un * norm_Vn
    use_pos <- (Mp > Mn)
    norm_Up_safe <- norm_Up; norm_Up_safe[norm_Up == 0] <- 1e-12
    norm_Un_safe <- norm_Un; norm_Un_safe[norm_Un == 0] <- 1e-12
    W_pos <- sweep(U_p, 2, (D_sqrt * Mp / norm_Up_safe), FUN = "*")
    W_neg <- sweep(U_n, 2, (D_sqrt * Mn / norm_Un_safe), FUN = "*")
    W_combined <- W_pos
    W_combined[, !use_pos] <- W_neg[, !use_pos]
    W[, idx] <- W_combined
  }
  W[is.nan(W)] <- 0
  W[W < 0] <- 0
  return(W)
}


#' @title Optimize NMF with kernel covariates
#' @description
#' \code{nmfkc} fits a nonnegative matrix factorization with kernel covariates
#' under the tri-factorization model \eqn{Y \approx X C A = X B}, where
#' \eqn{Y(P,N)} is the observation matrix, \eqn{A(R,N)} is the covariate matrix,
#' \eqn{X(P,Q)} is the basis matrix (\eqn{Q \le \min(P,N)}), \eqn{C(Q,R)} is the
#' parameter matrix, and \eqn{B(Q,N)=C A} is the coefficient matrix.
#' Given \eqn{Y} and (optionally) \eqn{A}, the algorithm estimates \eqn{X} and \eqn{C}.
#'
#' The estimation is based on minimizing a penalized objective function:
#' \deqn{
#'   J(X,C) =
#'   \begin{cases}
#'     \|Y - XCA\|_F^2 + \gamma\,\mathrm{tr}(C A A^\top C^\top), & \text{for method = "EU"},\\[6pt]
#'     \sum_{p,n}\!\bigl[-Y_{pn}\log(XCA)_{pn} + (XCA)_{pn}\bigr]
#'       + \gamma\,\mathrm{tr}(C A A^\top C^\top), & \text{for method = "KL"}.
#'   \end{cases}
#' }
#' When \code{A = NULL}, the penalty reduces to \eqn{\gamma\,\mathrm{tr}(C C^\top)}.
#' This ridge-type regularization on \eqn{C} (or \eqn{CA}) improves stability
#' and generalization by shrinking coefficient vectors toward zero.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. Default is \code{NULL} (no covariates).
#' @param Q Rank of the basis matrix \eqn{X}; must satisfy \eqn{Q \le \min(P,N)}.
#' @param gamma Nonnegative penalty parameter controlling
#'   the ridge regularization term \eqn{\gamma\,\mathrm{tr}(C A A^\top C^\top)}.
#' @param epsilon Positive convergence tolerance.
#' @param maxit Maximum number of iterations.
#' @param method Objective function: Euclidean distance \code{"EU"} (default) or Kullback–Leibler divergence \code{"KL"}.
#' @param X.restriction Constraint for columns of \eqn{X}:
#'   \code{"colSums"} (default; each column sums to 1),
#'   \code{"colSqSums"} (each column has unit \eqn{\ell_2} norm), or
#'   \code{"totalSum"} (entries sum to 1).
#' @param X.init Method for initializing the basis matrix \eqn{X}.
#'   Default is \code{"kmeans"}.
#'   \code{"nndsvd"} (Nonnegative Double SVD) can also be specified
#'   for a deterministic, often faster, initialization.
#' @param nstart Number of random starts for \code{\link[stats]{kmeans}} when initializing \eqn{X}.
#' @param seed Integer seed passed to \code{\link[base]{set.seed}}.
#' @param prefix Prefix for column names of \eqn{X} and row names of \eqn{B}.
#' @param print.trace Logical. If \code{TRUE}, prints progress every 10 iterations.
#' @param print.dims Logical. If \code{TRUE} (default), prints matrix dimensions and elapsed time.
#' @param save.time Logical. If \code{TRUE} (default), skips some post-computations (e.g., CPCC, silhouette) to save time.
#' @param save.memory Logical. If \code{TRUE}, performs only essential computations (implies \code{save.time = TRUE}) to reduce memory usage.
#' @param fast.calc Logical.
#' If \code{TRUE}, uses an optimized computation mode for large-scale matrices.
#' Specifically, for the Euclidean (\code{"EU"}) method with covariates \code{A},
#' precomputes \eqn{Y A^\top} and \eqn{A A^\top} to eliminate explicit dependence on the sample size \eqn{N}
#' during each iteration.
#' For the Kullback–Leibler (\code{"KL"}) method, automatically applies internal column blocking
#' to avoid constructing large \eqn{P \times N} ratio matrices.
#' If \code{FALSE} (default), runs the reference implementation with the original update rules.
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
#' # Example 1.
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
#' # Example 2.
#' Y <- matrix(cars$dist,nrow=1)
#' A <- rbind(1,cars$speed)
#' result <- nmfkc(Y,A,Q=1)
#' plot(cars$speed,as.vector(Y))
#' lines(cars$speed,as.vector(result$XB),col=2,lwd=2)

nmfkc <- function(Y,A=NULL,Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU",
                  X.restriction="colSums",X.init="kmeans",nstart=1,seed=123,
                  prefix="Basis",print.trace=FALSE,print.dims=TRUE,save.time=TRUE,save.memory=FALSE,fast.calc=FALSE){
  X.restriction <- match.arg(X.restriction, c("colSums", "colSqSums", "totalSum"))
  xnorm <- switch(X.restriction,
                 colSums   = function(X) sweep(X, 2, colSums(X), "/"),
                 colSqSums = function(X) sweep(X, 2, sqrt(colSums(X^2)), "/"),
                 totalSum  = function(X) X / sum(X)
  )
  # simplified silhouette coefficient
  # This internal function computes an approximate version of the silhouette coefficient.
  # Unlike the standard definition (e.g., cluster::silhouette), it does not require
  # pairwise distances among all samples. Instead, it estimates the silhouette by using
  # Euclidean distances between each sample and the mean vectors of clusters.
  # The nearest neighboring cluster is determined from centroid distances rather than
  # all pairwise sample distances.
  # This approximation greatly reduces computational cost (O(N^2) → O(NQ)) and
  # generally preserves the relative tendencies of silhouette values, though
  # the results may differ from the exact definition.
  silhouette.simple <- function(B.prob,B.cluster){
    if(is.matrix(B.prob)){Q <- nrow(B.prob)}else{Q <- 1}
    if(Q==1){
      return(list(cluster=NA,silhouette=NA,
                  silhouette.means=NA,silhouette.mean=NA))
    }else{
      index <-!is.na(B.cluster)
      B.prob <- B.prob[,index]
      B.cluster <- B.cluster[index]
      cluster.means <- matrix(0,nrow=Q,ncol=Q)
      ns <- NULL
      cluster.list <- NULL
      for(q in 1:Q){
        ns <- c(ns,sum(B.cluster==q))
        cluster.list <- c(cluster.list,list(which(B.cluster==q)))
        cluster.means[,q] <- rowMeans(B.prob[,B.cluster==q,drop=F])
      }
      si <- 0*B.cluster
      neighbor.cluster <- 0*B.cluster
      for(q in 1:Q){
        for(i in cluster.list[[q]]){
          di <- colSums((cluster.means-B.prob[,i])^2)
          qn <- ifelse(order(di)[1]==q,order(di)[2],order(di)[1])
          neighbor.cluster[i] <- qn
          if(ns[q]==1){
            si[i] <- 0
          }else{
            ai <- sum(colSums((B.prob[,cluster.list[[q]],drop=F]-B.prob[,i])^2)^0.5)/(ns[q]-1)
            bi <- sum(colSums((B.prob[,cluster.list[[qn]],drop=F]-B.prob[,i])^2)^0.5)/ns[qn]
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
  if(is.null(A)){
    dims <- sprintf("Y(%d,%d)~X(%d,%d)B(%d,%d)",
                    nrow(Y),ncol(Y),nrow(Y),Q,Q,ncol(Y))
  }else{
    dims <- sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)=XB(%d,%d)",
                    nrow(Y),ncol(Y),nrow(Y),Q,Q,nrow(A),nrow(A),ncol(Y),Q,ncol(Y))
  }
  if(print.dims) message(paste0(dims,"..."),appendLF=FALSE)
  start.time <- Sys.time()
  if(is.vector(Y)) Y <- matrix(Y,nrow=1)
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  if(!is.null(A)){
    if(min(A)<0){
      warning("The matrix A should be non-negative.")
      stop()
    }
  }
  if(min(Y)<0){
    warning("The matrix Y should be non-negative.")
    stop()
  }
  is.X.scalar <- FALSE
  if(nrow(Y)>=2){
    if(min(nrow(Y),ncol(Y))>=Q){
      if(ncol(Y)==Q){
        X <- Y
      }else{
        if (!is.null(seed)) {
          set.seed(seed)
        }
        if(X.init=="kmeans"){
          res.kmeans <- stats::kmeans(t(Y),centers=Q,iter.max=maxit,nstart=nstart)
          X <- t(res.kmeans$centers)
        }else{
          X <- .nndsvd(Y,Q)
        }
      }
    }else{
      warning("It does not hold Q<=min(P,N) where dim(Y)=(P,N).")
      stop()
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
    is.X.scalar <- TRUE
  }
  X <- xnorm(X)
  # Initialize C
  if(is.null(A)) C <- matrix(1,nrow=ncol(X),ncol=ncol(Y)) else C <- matrix(1,nrow=ncol(X),ncol=nrow(A))
  # -------- fast.calc precomputations (internal-only) --------
  hasA <- !is.null(A)
  # Auto block size for KL (aim ~32MB per block: ~4e6 doubles). Internal only.
  if (fast.calc && method!="EU") {
    P <- nrow(Y); N <- ncol(Y)
    target_elems <- 4e6
    bsz <- as.integer(max(2000L, min(8000L, floor(target_elems / max(1L,P)))))
    bsz <- min(bsz, N)
    if (print.trace) message(sprintf("KL block size (auto): %d", bsz))
  } # } end if (fast.calc && method!="EU")
  # EU with A: precompute Y%*%t(A) and A%*%t(A) to eliminate N in iterations
  if (hasA) {
    At  <- t(A)                     # (N x R)
    YAt <- Y %*% At                 # (P x R)
    AAt <- A %*% At                 # (R x R)
  }
  if(method=="EU"){
    if (fast.calc) YY_frob2 <- sum(Y*Y)            # constant term for objective
  }else{
    if (hasA) rowSumsA <- rowSums(A)
    rep1ncolY <- rep(1,ncol(Y))
  }
  # ---------------- Main loop ----------------
  epsilon.iter <- Inf  # initialize for post-loop warning safety
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){                                                           # <-- open: main iteration loop

    if (fast.calc) {                                                           # <-- open: FAST path switch

      if (method=="EU") {                                                      # <-- open: FAST EU branch

        if (hasA) {                                                            # <-- open: FAST EU with A
          # Compute XB without forming full B: XB = (X %*% C) %*% A
          XC  <- X %*% C                  # (P x R)
          XB  <- XC %*% A                 # (P x N)

          # Optional progress
          if (print.trace && i %% 10 == 0) message(format(Sys.time(), "%X")," ",i,"...")  # trace

          if(!is.X.scalar){                                                    # <-- open: update X when non-scalar
            # X update: X <- X * ((YAt %*% t(C)) / (X %*% (C AAt C^T)))
            numX <- YAt %*% t(C)                       # (P x Q)
            G    <- C %*% AAt %*% t(C)                 # (Q x Q)
            denX <- X %*% G                            # (P x Q)
            denX <- pmax(denX, 1e-12)                  # protect division
            X    <- X * (numX / denX); X[X<0] <- 0     # multiplicative update + clamp
            X <- xnorm(X)
          } # } end if (!is.X.scalar) (FAST EU with A)

          # C update (multiplicative): C <- C * (numC / denC)
          XtX  <- crossprod(X)                            # (Q x Q)
          numC <- crossprod(X, YAt)                       # (Q x R)
          denC <- (XtX %*% C) %*% AAt                     # (Q x R)
          if (gamma!=0) denC <- denC + gamma*C %*% AAt    # L2 regularization on C
          denC <- pmax(denC, 1e-12)                       # protect division by zero
          C    <- C * (numC / denC); C[C<0] <- 0          # multiplicative update + clamp negatives

          # EU objective without touching large N:
          # ||Y||_F^2 - 2*tr(Y^T X C A) + ||X C A||_F^2 + gamma||C||_F^2
          XC     <- X %*% C
          CtXtXC <- crossprod(C, XtX %*% C)               # (R x R)
          objfunc.iter[i] <- YY_frob2 - 2*sum(YAt * XC) + sum(AAt * CtXtXC) + gamma*sum(diag(C%*% AAt %*% t(C)))

        } else {                                                              # <-- else: FAST EU without A
          B  <- C
          XB <- X %*% B
          if (print.trace && i %% 10 == 0) message(format(Sys.time(), "%X")," ",i,"...")  # trace
          if(!is.X.scalar){                                                    # <-- open: update X when non-scalar
            # Use tcrossprod(B) for efficiency (avoids explicit t(B)%*%B)
            GB   <- tcrossprod(B)                         # (Q x Q)
            numX <- Y %*% t(B)                            # (P x Q)
            denX <- X %*% GB                              # (P x Q)
            denX <- pmax(denX, 1e-12)                     # protect division by zero
            X    <- X * (numX/denX); X[X<0] <- 0          # multiplicative update + clamp negatives
            X <- xnorm(X)
          } # } end if (!is.X.scalar) (FAST EU w/o A)
          # C update
          XtX  <- crossprod(X)
          numC <- crossprod(X, Y)
          denC <- XtX %*% C
          if (gamma!=0) denC <- denC + gamma*C            # L2 regularization on C
          denC <- pmax(denC, 1e-12)                       # protect division by zero
          C    <- C * (numC / denC); C[C<0] <- 0
          objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(C^2)
        } # } end if (hasA) else (FAST EU)
      }else{                                                                  # <-- else: FAST KL branch
        if (hasA) {                                                           # <-- open: FAST KL with A
          # Blocked accumulation to avoid building R=Y/XB as a huge matrix
          XC <- X %*% C                    # (P x R)
          numX <- matrix(0, nrow=nrow(Y), ncol=ncol(X))         # (P x Q) accumulator for X update
          numC_XR <- matrix(0, nrow=ncol(X), ncol=nrow(A))      # (Q x R) accumulator for C update
          for (j in seq(1L, ncol(Y), by=bsz)) {                               # <-- open: KL block loop
            idx   <- j:min(j+bsz-1L, ncol(Y))
            Ablk  <- A[, idx, drop=FALSE]          # (R x b)
            Bblk  <- C %*% Ablk                    # (Q x b)
            XBblk <- XC %*% Ablk                   # (P x b)
            Yblk  <- Y[, idx, drop=FALSE]          # (P x b)
            XBblk <- pmax(XBblk, 1e-12)            # protect division by zero
            Rblk  <- Yblk / XBblk                  # local ratio block
            Rblk[!is.finite(Rblk)] <- 0            # clean up Inf/NaN in block
            # Accumulate numerators for X and C updates
            numX     <- numX + Rblk %*% t(Bblk)              # (P x Q)
            numC_XR  <- numC_XR + crossprod(X, Rblk) %*% t(Ablk)  # (Q x R)
          } # } end for blocks (FAST KL with A)
          # Denominator for X update: rowSums(B) accumulated by blocks
          rsB <- rep(0, ncol(X))    # (Q)
          for (j in seq(1L, ncol(Y), by=bsz)) {                               # <-- open: accumulate rowSums(B)
            idx  <- j:min(j+bsz-1L, ncol(Y))
            Ablk <- A[, idx, drop=FALSE]
            rsB  <- rsB + rowSums(C %*% Ablk)
          } # } end for blocks (rowSums(B))
          # X multiplicative update with blocking results
          X <- sweep(X * pmax(0, numX), 2, pmax(rsB, 1e-12), "/"); X[X<0] <- 0
          X <- xnorm(X)
          # C multiplicative update (blocked numerator, small-matrix denominator)
          cX   <- colSums(X)                           # (Q)
          denC <- outer(cX, rowSums(A))                # (Q x R)
          if (gamma!=0) denC <- denC+2*gamma*(C %*% AAt) # L2 regularization on C in KL (factor 2)
          denC <- pmax(denC, 1e-12)
          C    <- C * (numC_XR / denC); C[C<0] <- 0
          # KL objective (full evaluation; log protected)
          B  <- C %*% A
          XB <- X %*% B
          objfunc.iter[i] <- sum(-Y*.z(log(XB)) + XB) +gamma*sum(diag(C%*% AAt %*% t(C)))
        } else {                                                              # <-- else: FAST KL without A
          B  <- C
          XB <- pmax(X %*% B, 1e-12)

          if (print.trace && i %% 10 == 0) message(format(Sys.time(), "%X")," ",i,"...")  # trace

          if(!is.X.scalar){                                                    # <-- open: update X when non-scalar
            R    <- Y / XB; R[!is.finite(R)] <- 0
            numX <- R %*% t(B)
            rsB  <- rowSums(B)
            X    <- sweep(X * numX, 2, pmax(rsB,1e-12), "/"); X[X<0] <- 0
            X <- xnorm(X)
          } # } end if (!is.X.scalar) (FAST KL w/o A)

          R    <- Y / XB; R[!is.finite(R)] <- 0
          numC <- crossprod(X, R)
          denC <- matrix(colSums(X), nrow=nrow(C), ncol=ncol(C))
          if (gamma!=0) denC <- denC + 2*gamma*C
          denC <- pmax(denC, 1e-12)
          C    <- C * (numC / denC); C[C<0] <- 0
          objfunc.iter[i] <- sum(-Y*.z(log(XB))+XB)+gamma*sum(C^2)
        } # } end if (hasA) else (FAST KL)
      } # } end if (method=="EU") else (FAST KL)

    }else{                                                                       # <-- else: REFERENCE path

      # --- Original (reference) implementation below (kept for compatibility) ---

      if(is.null(A)) B <- C else B <- C %*% A
      XB <- X %*% B
      if(print.trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))

      if(method=="EU"){                                                          # <-- open: REF EU
        if(!is.X.scalar){
          X <- X*.z((Y%*% t(B))/(XB%*%t(B)))
          X <- xnorm(X)
        } # } end if (!is.X.scalar) (REF EU)

        if(is.null(A)) {
          if(gamma!=0){
            C <- C*.z((t(X)%*%Y)/(t(X)%*%XB+gamma*C))
          }else{
            C <- C*.z((t(X)%*%Y)/(t(X)%*%XB))
          }
        } else {
          if(gamma!=0){
            C <- C*.z((t(X)%*% YAt)/(t(X)%*%XB%*%At+gamma*(C %*% AAt)))
          }else{
            C <- C*.z((t(X)%*% YAt)/(t(X)%*%XB%*%At))
          }
        } # } end if (is.null(A)) else (REF EU C-update)

        if(gamma!=0){
          if(is.null(A)) {
            objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(C^2)
          }else{
            objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(diag(C%*% AAt %*% t(C)))
          }
        }else{
          objfunc.iter[i] <- sum((Y-XB)^2)
        }

      }else{                                                                     # <-- else: REF KL
        if(!is.X.scalar){
          X <- X*.z((Y/XB)%*%t(B))/(rep(1,nrow(Y))%o%rowSums(B))
          X <- xnorm(X)
        } # } end if (!is.X.scalar) (REF KL)

        if(is.null(A)) {
          if(gamma!=0){
            C <- C*.z(t(X)%*%.z(Y/XB)/(colSums(X)%o%rep1ncolY+2*gamma*C))
          }else{
            C <- C*.z(t(X)%*%.z(Y/XB)/(colSums(X)%o%rep1ncolY))
          }
        } else {
          if(gamma!=0){
            C <- C*.z(t(X)%*%.z(Y/XB)%*%At/(colSums(X)%o%rowSumsA+2*gamma*(C %*% AAt)))
          }else{
            C <- C*.z(t(X)%*%.z(Y/XB)%*%At/(colSums(X)%o%rowSumsA))
          }
        } # } end if (is.null(A)) else (REF KL C-update)

        if(gamma!=0){
          if(is.null(A)) {
            objfunc.iter[i] <- sum(-Y*.z(log(XB))+XB)+gamma*sum(C^2)
          }else{
            objfunc.iter[i] <- sum(-Y*.z(log(XB))+XB)+gamma*sum(diag(C%*% AAt %*% t(C)))
          }
        }else{
          objfunc.iter[i] <- sum(-Y*.z(log(XB))+XB)
        }
      } # } end if (method=="EU") else (REF)
    } # } end if (fast.calc) else (REFERENCE)

    # Convergence check (shared)
    if(i>=10){
      #epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/max(abs(objfunc.iter[i]),1e-10)
      if(epsilon.iter <= abs(epsilon)){
        objfunc.iter <- objfunc.iter[10:i]
        break
      } # } end if (converged)
    } # } end if (i >= 10)
  } # } end for (main iteration loop)

  if(is.null(A)) B <- C else B <- C %*% A
  XB <- X %*% B
  if(method=="EU"){
    if(is.null(A)) {
      objfunc <- sum((Y-XB)^2)+gamma*sum(C^2)
    }else{
      objfunc <- sum((Y-XB)^2)+gamma*sum(diag(C%*% AAt %*% t(C)))
    }
  }else{
    if(is.null(A)) {
      objfunc <- sum(-Y*.z(log(XB))+XB)+gamma*sum(C^2)
    }else{
      objfunc <- sum(-Y*.z(log(XB))+XB)+gamma*sum(diag(C%*% AAt %*% t(C)))
    }
  }
  if(ncol(X)>1 & sum(rowSums(X)==1)==nrow(X)){
    index <- order(matrix(1:nrow(X)/nrow(X),nrow=1) %*% X)
    X <- X[,index]
    B <- B[index,]
    C <- C[index,]
  }
  rownames(C) <- paste0(prefix,1:nrow(C))
  rownames(X) <- rownames(Y)
  colnames(X) <- paste0(prefix,1:ncol(X))
  rownames(B) <- paste0(prefix,1:nrow(B))
  colnames(B) <- colnames(Y)
  if(X.restriction=="totalSum"){
    if(is.null(A)) nparam <- prod(dim(X))-1+prod(dim(B)) else nparam <- prod(dim(X))-1+prod(dim(A))
  }else{
    if(is.null(A)) nparam <- (nrow(X)-1)*ncol(X)+prod(dim(B)) else nparam <- (nrow(X)-1)*ncol(X)+prod(dim(A))
  }
  if(method=="EU"){
    ICp <- log(objfunc/prod(dim(Y)))+Q*sum(dim(Y))/prod(dim(Y))*log(prod(dim(Y))/sum(dim(Y)))
    sigma2 <- sum((Y-XB)^2)/prod(dim(Y))
    AIC <- prod(dim(Y))*log(sigma2)+2*nparam
    BIC <- prod(dim(Y))*log(sigma2)+nparam*log(ncol(Y))
  }else{
    ICp <- NA
    AIC <- 2*sum(-Y*.z(log(XB))+XB)+2*nparam
    BIC <- 2*sum(-Y*.z(log(XB))+XB)+nparam*log(ncol(Y))
  }
  if(save.memory==FALSE){
    r2 <- stats::cor(as.vector(XB),as.vector(Y))^2
    sigma <- stats::sd(as.vector(Y)-as.vector(XB))
    B.prob <- t(.z(t(B)/colSums(B)))
    B.prob.sd.min <- min(apply(B.prob,1,stats::sd))
    B.cluster <- apply(B.prob,2,which.max)
    B.cluster[colSums(B.prob)==0] <- NA
    rownames(B.prob) <- paste0(prefix,1:nrow(B))
    colnames(B.prob) <- colnames(Y)
    colnames(XB) <- colnames(Y)
    X.prob <- .z(X/rowSums(X))
    X.cluster <- apply(X.prob,1,which.max)
    X.cluster[rowSums(X.prob)==0] <- NA
    if(save.time){
      silhouette <- NA
      CPCC <- NA
      dist.cor <- NA
    }else{
      silhouette <- silhouette.simple(B.prob,B.cluster)
      dist.cor <- stats::cor(as.vector(stats::dist(t(Y))),as.vector(stats::dist(t(B))))
      if(Q>=2){
        M <- t(B.prob) %*% B.prob
        h.dist <- as.matrix(stats::cophenetic(stats::hclust(stats::as.dist(1-M))))
        up <- upper.tri(M)
        CPCC <- stats::cor(h.dist[up],(1-M)[up])
      }else{
        CPCC <- NA
      }
    }
  }else{
    r2 <- NA
    sigma <- NA
    B.prob <- NA
    B.prob.sd.min <- NA
    B.cluster <- NA
    XB <- NA
    X.prob <- NA
    X.cluster <- NA
    silhouette <- NA
    CPCC <- NA
    dist.cor <- NA
  }
  if(epsilon.iter > abs(epsilon)) warning(paste0(
    "maximum iterations (",maxit,
    ") reached and the optimization hasn't converged yet."))
  end.time <- Sys.time()
  diff.time <- difftime(end.time,start.time,units="sec")
  diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                         paste0(round(diff.time/60,1),"min"))
  if(print.dims) message(diff.time.st)
  result <- list(call=match.call(),dims=dims,runtime=diff.time.st,
                 X=X,B=B,XB=XB,C=C,
                 B.prob=B.prob,B.cluster=B.cluster,
                 X.prob=X.prob,X.cluster=X.cluster,
                 objfunc=objfunc,objfunc.iter=objfunc.iter,r.squared=r2,sigma=sigma,
                 criterion=list(B.prob.sd.min=B.prob.sd.min,
                                ICp=ICp,AIC=AIC,BIC=BIC,
                                silhouette=silhouette,CPCC=CPCC,dist.cor=dist.cor))
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
  extra_args <- list(...)
  ans <- list()
  ans$call <- object$call
  ans$dims <- object$dims
  ans$runtime <- object$runtime
  ans$objfunc.iter.length <- length(object$objfunc.iter)
  ans$r.squared <- object$r.squared
  ans$sigma <- object$sigma
  if (!is.null(object$X) && is.matrix(object$X)) {
    ans$X.dist <- summary(as.vector(object$X))
  } else {
    ans$X.dist <- NULL
  }
  if (!is.null(object$B.prob) && is.matrix(object$B.prob)) {
    ans$B.prob.dist <- summary(as.vector(object$B.prob))
  } else {
    ans$B.prob.dist <- NULL
  }
  class(ans) <- "summary.nmfkc"
  return(ans)
}


#' @export
print.summary.nmfkc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("dims: ", x$dims, "\n\n")

  cat("Runtime:", x$runtime, "\n")
  cat("Number of iterations:", x$objfunc.iter.length, "\n")
  cat("Multiple R-squared:", format(x$r.squared, digits = digits), "\n")
  cat("Residual standard error:", format(x$sigma, digits = digits), "\n\n")

  if (!is.null(x$X.dist)) {
    cat("Distribution of X:\n")
    print(x$X.dist, digits = digits)
  }
  if (!is.null(x$B.prob.dist)) {
    cat("Distribution of B.prob:\n")
    print(x$B.prob.dist, digits = digits)
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
      B.prob <- sweep(B, 2, colSums(B), "/")
      B.prob <- .z(B.prob)
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
#'   \item \eqn{Y(P,N)} is the observation matrix,
#'   \item \eqn{A(R,N)} is the covariate matrix,
#'   \item \eqn{X(P,Q)} is the basis matrix with \eqn{Q \le \min(P,N)},
#'   \item \eqn{C(Q,R)} is the parameter matrix, and
#'   \item \eqn{B(Q,N)} is the coefficient matrix (\eqn{B = C A}).
#' }
#' Given \eqn{Y} (and optionally \eqn{A}), \eqn{X} and \eqn{C} are estimated,
#' and the predictive performance is assessed by cross-validation.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#' @param Q Rank of the basis matrix \eqn{X}; must satisfy \eqn{Q \le \min(P,N)}.
#' @param div Number of folds (\eqn{k}) for cross-validation (default: 5).
#' @param seed Integer seed for reproducibility, passed to \code{\link{set.seed}}.
#' @param ... Additional arguments passed to \code{\link{nmfkc}}.
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

nmfkc.cv <- function(Y,A=NULL,Q=2,div=5,seed=123,...){
  extra_args <- list(...)
  gamma <- if (!is.null(extra_args$gamma)) extra_args$gamma else 0
  epsilon <- if (!is.null(extra_args$epsilon)) extra_args$epsilon else 1e-4
  maxit   <- if (!is.null(extra_args$maxit))   extra_args$maxit   else 5000
  method  <- if (!is.null(extra_args$method))  extra_args$method  else "EU"
  is.identity.matrix <- function(A, tol = 1e-12) {
    if (nrow(A) != ncol(A)) return(FALSE)
    isTRUE(all.equal(A, diag(nrow(A)), tolerance = tol))
  }
  optimize.B.from.Y <- function(result,Y,gamma,epsilon,maxit,method){
    X <- result$X
    C <- matrix(1,nrow=ncol(X),ncol=ncol(Y))
    oldSum <- 0
    epsilon.iter <- Inf
    for(l in 1:maxit){
      B <- C
      XB <- X %*% B
      if(method=="EU"){
        C <- C*.z((t(X)%*%Y)/(t(X)%*%XB+gamma*C))
      }else{
        C <- C*(t(X)%*%.z(Y/XB)/(colSums(X)%o%rep(1,ncol(Y))+2*gamma*C))
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
    colnames(B) <- colnames(Y)
    colnames(XB) <- colnames(Y)
    if(epsilon.iter > abs(epsilon)) warning(paste0(
      "maximum iterations (",maxit,
      ") reached and the optimization hasn't converged yet."))
    return(list(B=B,XB=XB))
  }
  if(is.null(A)){
    is_identity <- TRUE
    is_symmetric.matrix <- TRUE
  }else{
    is_identity <- is.identity.matrix(A)
    is_symmetric.matrix <- isSymmetric(A, tol=1e-12)
  }
  n <- ncol(Y)
  remainder <- n %% div
  division <- n %/% div
  block <- 0*(1:n)
  set.seed(seed)
  index <- sample(1:n,n,replace=F)
  for(i in 1:(div-1)){
    plus <- ifelse(i<=remainder,1,0)
    j <- index[1:(division+plus)]
    block[j] <- i
    index <- index[-(1:(division+plus))]
  }
  block[index] <- div
  index <- block
  objfunc.block <- 0*(1:div)
  for(j in 1:div){
    Y_j <- Y[,index!=j]
    Yj <- Y[,index==j]
    if(is_identity){
      A_j <- NULL
    }else{
      if(is_symmetric.matrix){
        A_j <- A[index!=j,index!=j] # kernel matrix or identity matrix
      }else{
        A_j <- A[,index!=j] # ordinary design matrix
      }
    }
    nmfkc_args <- c(
      list(...),
      list(Y = Y_j, A = A_j, Q = Q, seed=NULL,
           print.trace = FALSE,
           print.dims = FALSE,
           save.time = TRUE,
           save.memory = TRUE)
    )
    res_j <- do.call("nmfkc", nmfkc_args)
    if(is_identity){
      resj <- optimize.B.from.Y(res_j,Yj,gamma,epsilon,maxit,method)
      XBj <- resj$XB
    }else{
      if(is_symmetric.matrix){
        Aj <- A[index!=j,index==j]
      }else{
        Aj <- A[,index==j]
      }
      XBj <- res_j$X %*% res_j$C %*% Aj
    }
    if(method=="EU"){
      objfunc.block[j] <- sum((Yj-XBj)^2)
    }else{
      objfunc.block[j] <- sum(-Yj*.z(log(XBj))+XBj)
    }
  }
  objfunc <- sum(objfunc.block)/n
  return(list(objfunc=objfunc,objfunc.block=objfunc.block,block=index))
}


#' @title Rank selection diagnostics with graphical output
#' @description
#' \code{nmfkc.rank} provides diagnostic criteria for selecting the rank (\eqn{Q})
#' in NMF with kernel covariates. Several model selection measures are computed
#' (e.g., R-squared, information criteria, silhouette, CPCC, ARI), and results
#' can be visualized in a plot. This function is still under development.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#' @param Q A vector of candidate ranks to be evaluated.
#' @param plot Logical. If \code{TRUE} (default), draws a plot of the diagnostic criteria.
#' @param ... Additional arguments passed to \code{\link{nmfkc}}.
#'
#' @return A data frame containing rank values (\code{Q}) and the corresponding diagnostic criteria:
#' \itemize{
#'   \item \code{r.squared}: Coefficient of determination (\eqn{R^2}).
#'   \item \code{ICp}, \code{AIC}, \code{BIC}: Information criteria.
#'   \item \code{B.prob.sd.min}: Minimum standard deviation of \code{B.prob}.
#'   \item \code{ARI}: Adjusted Rand Index relative to the previous rank.
#'   \item \code{silhouette}: Mean silhouette score (if computed).
#'   \item \code{CPCC}: Cophenetic correlation coefficient (if computed).
#'   \item \code{dist.cor}: Distance between Y and B (if computed).
#' }
#' @seealso \code{\link{nmfkc}}
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
#' nmfkc.rank(Y,Q=2:4)

nmfkc.rank <- function(Y,A=NULL,Q=1:2,plot=TRUE,...){

  extra_args <- list(...)

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
    Q = Q,
    r.squared = numeric(num_q),
    ICp = numeric(num_q),
    AIC = numeric(num_q),
    BIC = numeric(num_q),
    B.prob.sd.min = numeric(num_q),
    ARI = numeric(num_q),
    silhouette = numeric(num_q),
    CPCC = numeric(num_q),
    dist.cor = numeric(num_q)
  )

  cluster.old <- NULL
  for(q_idx in 1:num_q){
    current_Q <- Q[q_idx]
    save_time_for_this_run <- if (is.null(extra_args$save.time)) TRUE else extra_args$save.time
    extra_args$save.memory <- NULL
    nmfkc_args <- c(
      list(Y = Y,
           A = A,
           Q = current_Q,
           save.memory = FALSE
      ),
      extra_args
    )
    result <- do.call("nmfkc", nmfkc_args)
    results_df$r.squared[q_idx] <- result$r.squared
    results_df$ICp[q_idx] <- result$criterion$ICp
    results_df$AIC[q_idx] <- result$criterion$AIC
    results_df$BIC[q_idx] <- result$criterion$BIC
    results_df$B.prob.sd.min[q_idx] <- result$criterion$B.prob.sd.min

    if(save_time_for_this_run){
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

  if(plot){
    plot(results_df$Q, results_df$r.squared, type="l", col=2, xlab="Rank", ylab="Criterion", ylim=c(0,1), lwd=3)
    graphics::text(results_df$Q, results_df$r.squared, results_df$Q)

    legend_text <- "r.squared"
    legend_cols <- 2

    graphics::lines(results_df$Q, results_df$B.prob.sd.min, col=3, lwd=3)
    graphics::text(results_df$Q, results_df$B.prob.sd.min, results_df$Q)
    legend_text <- c(legend_text, "B.prob.sd.min")
    legend_cols <- c(legend_cols, 3)

    if (any(!is.na(results_df$ARI))) {
      graphics::lines(results_df$Q, results_df$ARI, col=4, lwd=3)
      graphics::text(results_df$Q, results_df$ARI, results_df$Q)
      legend_text <- c(legend_text, "ARI vs Q-1")
      legend_cols <- c(legend_cols, 4)
    }

    if (any(!is.na(results_df$silhouette))) {
      graphics::lines(results_df$Q, results_df$silhouette, col=7, lwd=3)
      graphics::text(results_df$Q, results_df$silhouette, results_df$Q)
      legend_text <- c(legend_text, "silhouette")
      legend_cols <- c(legend_cols, 7)
    }

    if (any(!is.na(results_df$CPCC))) {
      graphics::lines(results_df$Q, results_df$CPCC, col=6, lwd=3)
      graphics::text(results_df$Q, results_df$CPCC, results_df$Q)
      legend_text <- c(legend_text, "CPCC")
      legend_cols <- c(legend_cols, 6)
    }

    if (any(!is.na(results_df$dist.cor))) {
      graphics::lines(results_df$Q, results_df$dist.cor, col=5, lwd=3)
      graphics::text(results_df$Q, results_df$dist.cor, results_df$Q)
      legend_text <- c(legend_text, "dist.cor")
      legend_cols <- c(legend_cols, 5)
    }

    graphics::legend("bottomleft", legend=legend_text, fill=legend_cols, bg="white")
  }

  invisible(results_df)
}

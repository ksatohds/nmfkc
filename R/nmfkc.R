.onAttach <- function(...) {
  packageStartupMessage("Last update on 5 SEP 2025")
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
  list(Y=Ya,A=A,A.columns=A.columns)
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
#' @param div Number of partitions for cross-validation, corresponding to \eqn{k} in k-fold CV.
#' @param seed Integer seed for reproducibility, passed to \code{set.seed}.
#' @param plot Logical. If TRUE (default), a plot of the objective function values is drawn.
#' @param ... Additional arguments passed to \code{nmfkc.cv}.
#'
#' @return A list with components:
#' \item{degree}{The lag order that minimizes the cross-validation objective function.}
#' \item{degree.max}{Maximum recommended lag order, computed as \eqn{10 \log_{10}(N)}
#'   following the \code{ar} function in the \pkg{stats} package.}
#' \item{objfunc}{Objective function values for each candidate lag order.}
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

nmfkc.ar.degree.cv <- function(Y,Q=2,degree=1:2,intercept=T,div=5,seed=123,plot=TRUE,...){
  objfuncs <- 0*(1:length(degree))
  for(i in 1:length(degree)){
    start.time <- Sys.time()
    message(paste0("degree=",degree[i],"..."),appendLF=FALSE)
    a <- nmfkc.ar(Y=Y,degree=degree[i],intercept=intercept)
    result.cv <- nmfkc.cv(Y=a$Y,A=a$A,Q=Q,div=div,seed=seed,...)
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
#' @param beta Kernel parameter. Default is 0.5. For the \code{"Periodic"} kernel, specify as \code{c(beta1, beta2)}.
#' @param kernel Kernel function to use. Default is \code{"Gaussian"}. Options are \code{"Gaussian"}, \code{"Exponential"}, \code{"Periodic"}, \code{"Linear"}, \code{"NormalizedLinear"}, and \code{"Polynomial"}.
#' @param degree Degree parameter for the \code{"Polynomial"} kernel. Default is 2.
#'
#' @return Kernel matrix \eqn{A(N,M)}.
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

nmfkc.kernel <- function(U, V = NULL, beta=0.5,
                         kernel = c("Gaussian","Exponential","Periodic",
                                    "Linear","NormalizedLinear","Polynomial"),
                         degree=2){
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
              Gaussian = exp(-beta * D2),
              Exponential = {
                d <- sqrt(D2)
                exp(-beta * d)
              },
              Periodic = {
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
              Polynomial = (G + beta)^degree
  )
  if (min(K) < 0) {
    warning("The constructed matrix is not non-negative.")
  }
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
#'   varying the scale parameter \eqn{\sigma=d_{med}} by factors
#'   \eqn{2^{-2},\,2^{-1},\,1,\,2^{1},\,2^{2}} and converting each to
#'   \eqn{\beta=1/(2\sigma^2)}. The values are ordered from small to large.}
#' \item{dist_median}{the median nearest-neighbor distance}
#' \item{block_size_used}{actual block size used in computation}
#' @details
#' The function computes all pairwise squared distances between columns of
#' \eqn{U}, excludes self-distances, and takes the median of the nearest-neighbor
#' distances (after square root). This median is then used to set \eqn{\beta}.
#' @seealso \code{\link{nmfkc.kernel}} for creating kernel matrices from covariates.
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
  beta_candidates <- 1 / (2 *(d_med*c(4,2,0,2,4))^2)
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
#' @param kernel Kernel function to use. Default is \code{"Gaussian"}. See \code{\link{nmfkc.kernel}} for available options.
#' @param degree Degree parameter for the \code{"Polynomial"} kernel. Default is 2.
#' @param div Number of partitions for cross-validation, corresponding to \eqn{k} in k-fold CV.
#' @param seed Integer seed for reproducibility, passed to \code{set.seed}.
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

nmfkc.kernel.beta.cv <- function(Y,Q=2,U,V=NULL,beta=NULL,
                                 kernel="Gaussian",degree=2,div=5,seed=123,plot=TRUE,...){
  if(is.null(beta)){
    if(is.null(V)) V <- U
    result.beta <- nmfkc.kernel.beta.nearest.med(V)
    beta.med <- result.beta$beta
    beta <- beta.med*10^(-3:1)
  }
  objfuncs <- 0*(1:length(beta))
  for(i in 1:length(beta)){
    start.time <- Sys.time()
    message(paste0("beta=",beta[i],"..."),appendLF=FALSE)
    A <- nmfkc.kernel(U=U,V=V,beta=beta[i],kernel=kernel,degree=degree)
    result <- nmfkc.cv(Y=Y,A=A,Q=Q,div=div,seed=seed,...)
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


#' @title Optimize NMF with kernel covariates
#' @description
#' \code{nmfkc} fits a nonnegative matrix factorization with kernel covariates
#' under the tri-factorization model \eqn{Y \approx X C A = X B}, where
#' \eqn{Y(P,N)} is the observation matrix, \eqn{A(R,N)} is the covariate matrix,
#' \eqn{X(P,Q)} is the basis matrix (\eqn{Q \le \min(P,N)}), \eqn{C(Q,R)} is the
#' parameter matrix, and \eqn{B(Q,N)=C A} is the coefficient matrix.
#' Given \eqn{Y} and (optionally) \eqn{A}, the algorithm estimates \eqn{X} and \eqn{C}.
#'
#' @param Y Observation matrix.
#' @param A Covariate matrix. Default is \code{NULL} (no covariates).
#' @param Q Rank of the basis matrix \eqn{X}; must satisfy \eqn{Q \le \min(P,N)}.
#' @param gamma Nonnegative penalty parameter for the parameter matrix \eqn{C}.
#' @param epsilon Positive convergence tolerance.
#' @param maxit Maximum number of iterations.
#' @param method Objective function: Euclidean distance \code{"EU"} (default) or Kullback–Leibler divergence \code{"KL"}.
#' @param X.restriction Constraint for columns of \eqn{X}:
#'   \code{"colSums"} (default; each column sums to 1),
#'   \code{"colSqSums"} (each column has unit \eqn{\ell_2} norm), or
#'   \code{"totalSum"} (entries sum to 1).
#' @param nstart Number of random starts for \code{\link[stats]{kmeans}} when initializing \eqn{X}.
#' @param seed Integer seed passed to \code{\link[base]{set.seed}}.
#' @param prefix Prefix for column names of \eqn{X} and row names of \eqn{B}.
#' @param print.trace Logical. If \code{TRUE}, prints progress every 10 iterations.
#' @param print.dims Logical. If \code{TRUE} (default), prints matrix dimensions and elapsed time.
#' @param save.time Logical. If \code{TRUE} (default), skips some post-computations (e.g., CPCC, silhouette) to save time.
#' @param save.memory Logical. If \code{TRUE}, performs only essential computations (implies \code{save.time = TRUE}) to reduce memory usage.
#'
#' @return A list with components:
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
#' \item{criterion}{A list of selection criteria, including \code{ICp}, \code{CPCC}, \code{silhouette}, \code{AIC}, and \code{BIC}.}
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
                  X.restriction="colSums",nstart=1,seed=123,
                  prefix="Basis",print.trace=FALSE,print.dims=TRUE,save.time=TRUE,save.memory=FALSE){
  z <- function(x){
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
    return(x)
  }
  mysilhouette <- function(B.prob,B.cluster){
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
  if(print.dims)if(is.null(A)){
    message(
      sprintf("Y(%d,%d)~X(%d,%d)B(%d,%d)...",
              nrow(Y),ncol(Y),nrow(Y),Q,Q,ncol(Y)),appendLF=FALSE)
  }else{
    message(
      sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)=XB(%d,%d)...",
              nrow(Y),ncol(Y),nrow(Y),Q,Q,nrow(A),nrow(A),ncol(Y),Q,ncol(Y)),appendLF=FALSE)
  }
  start.time <- Sys.time()
  set.seed(seed)
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
        res.kmeans <- stats::kmeans(t(Y),centers=Q,iter.max=maxit,nstart=nstart)
        X <- t(res.kmeans$centers)
      }
    }else{
      warning("It does not hold Q<=min(P,N) where dim(Y)=(P,N).")
      stop()
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
    is.X.scalar <- TRUE
  }
  if(!is.X.scalar){
    if(X.restriction=="colSums"){
      X <- t(t(X)/colSums(X))
    }else if(X.restriction=="colSqSums"){
      X <- t(t(X)/colSums(X^2)^0.5)
    }else if(X.restriction=="totalSum") X <- X/sum(X)
  }
  if(is.null(A)) C <- matrix(1,nrow=ncol(X),ncol=ncol(Y)) else C <- matrix(1,nrow=ncol(X),ncol=nrow(A))
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){
    if(is.null(A)) B <- C else B <- C %*% A
    XB <- X %*% B
    if(print.trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))
    if(method=="EU"){
      if(!is.X.scalar){
        X <- X*z((Y%*% t(B))/(XB%*%t(B)))
        if(X.restriction=="colSums"){
          X <- t(t(X)/colSums(X))
        }else if(X.restriction=="colSqSums"){
          X <- t(t(X)/colSums(X^2)^0.5)
        }else if(X.restriction=="totalSum") X <- X/sum(X)
      }
      if(is.null(A)) C <- C*z((t(X)%*%Y)/(t(X)%*%XB+gamma*C)) else
        C <- C*z((t(X)%*%Y%*%t(A))/(t(X)%*%XB%*%t(A)+gamma*C))
      objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(C^2)
    }else{
      if(!is.X.scalar){
        X <- X*z((Y/XB)%*%t(B))/(rep(1,nrow(Y))%o%rowSums(B))
        if(X.restriction=="colSums"){
          X <- t(t(X)/colSums(X))
        }else if(X.restriction=="colSqSums"){
          X <- t(t(X)/colSums(X^2)^0.5)
        }else if(X.restriction=="totalSum") X <- X/sum(X)
      }
      if(is.null(A)) C <- C*z(t(X)%*%z(Y/XB)/(colSums(X)%o%rep(1,ncol(Y))+2*gamma*C)) else
        C <- C*z(t(X)%*%z(Y/XB)%*%t(A)/(colSums(X)%o%rowSums(A)+2*gamma*C))
      objfunc.iter[i] <- sum(-Y*z(log(XB))+XB)+gamma*sum(C^2)
    }
    if(i>=10){
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
      if(epsilon.iter <= abs(epsilon)){
        objfunc.iter <- objfunc.iter[10:i]
        break
      }
    }
  }
  if(is.null(A)) B <- C else B <- C %*% A
  XB <- X %*% B
  if(method=="EU"){
    objfunc <- sum((Y-XB)^2)+gamma*sum(C^2)
  }else{
    objfunc <- sum(-Y*z(log(XB))+XB)+gamma*sum(C^2)
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
  r2 <- stats::cor(as.vector(XB),as.vector(Y))^2
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
    AIC <- 2*sum(-Y*z(log(XB))+XB)+2*nparam
    BIC <- 2*sum(-Y*z(log(XB))+XB)+nparam*log(ncol(Y))
  }
  if(save.memory==FALSE){
    B.prob <- t(z(t(B)/colSums(B)))
    B.prob.sd.min <- min(apply(B.prob,1,stats::sd))
    B.cluster <- apply(B.prob,2,which.max)
    B.cluster[colSums(B.prob)==0] <- NA
    rownames(B.prob) <- paste0(prefix,1:nrow(B))
    colnames(B.prob) <- colnames(Y)
    colnames(XB) <- colnames(Y)
    X.prob <- z(X/rowSums(X))
    X.cluster <- apply(X.prob,1,which.max)
    X.cluster[rowSums(X.prob)==0] <- NA
    if(save.time){
      silhouette <- NA
      CPCC <- NA
    }else{
      silhouette <- mysilhouette(B.prob,B.cluster)
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
    B.prob <- NA
    B.prob.sd.min <- NA
    B.cluster <- NA
    XB <- NA
    X.prob <- NA
    X.cluster <- NA
    silhouette <- NA
    CPCC <- NA
  }
  if(epsilon.iter > abs(epsilon)) warning(paste0(
    "maximum iterations (",maxit,
    ") reached and the optimization hasn't converged yet."))
  end.time <- Sys.time()
  diff.time <- difftime(end.time,start.time,units="sec")
  diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                         paste0(round(diff.time/60,1),"min"))
  if(print.dims) message(diff.time.st)
  result <- list(X=X,B=B,XB=XB,C=C,
                 B.prob=B.prob,B.cluster=B.cluster,
                 X.prob=X.prob,X.cluster=X.cluster,
                 objfunc=objfunc,objfunc.iter=objfunc.iter,r.squared=r2,
                 criterion=list(B.prob.sd.min=B.prob.sd.min,ICp=ICp,AIC=AIC,BIC=BIC,silhouette=silhouette,CPCC=CPCC))
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
  plot(x$objfunc.iter,xlab="iter",ylab="objfunc",main=paste0("r.squared=",round(x$r.squared,3)),...)
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
#' @export
predict.nmfkc <- function(x,newA=NULL,type="response"){
  z <- function(x){
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
    return(x)
  }
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
      B.prob <- t(z(t(B)/colSums(B)))
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
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- nmfkc.class(iris$Species)
#' Y[,1:6]
nmfkc.class <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  unix <- levels(x)
  X <- matrix(0,nrow=length(unix),ncol=length(x))
  rownames(X) <- unix
  for(j in 1:length(unix)) X[j,] <- ifelse(x==unix[j],1,0)
  result <- X
  return(result)}


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
  y <- 0*x
  for(j in 1:ncol(x))y[,j] <- (x[,j]-r[1,j])/(r[2,j]-r[1,j])
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
  y <- 0 * x
  for (j in 1:ncol(x)) {
    y[, j] <- x[, j] * (r[2, j] - r[1, j]) + r[1, j]
  }
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
  arglist=list(...)
  gamma <- ifelse("gamma" %in% names(arglist),arglist$gamma,0)
  epsilon <- ifelse("epsilon" %in% names(arglist),arglist$epsilon,1e-4)
  maxit <- ifelse("maxit" %in% names(arglist),arglist$maxit,5000)
  method <- ifelse("method" %in% names(arglist),arglist$method,"EU")
  X.restriction <- ifelse("X.restriction" %in% names(arglist),arglist$X.restriction,"colSums")
  nstart <- ifelse("nstart" %in% names(arglist),arglist$nstart,1)
  prefix <- ifelse("prefix" %in% names(arglist),arglist$prefix,"Basis")
  print.trace <- ifelse("print.trace" %in% names(arglist),arglist$print.trace,FALSE)
  print.dims <- ifelse("print.dims" %in% names(arglist),arglist$print.dims,FALSE)
  save.time <- TRUE
  save.memory <- TRUE
  if(is.null(A)) A <- diag(ncol(Y))
  is.symmetric.matrix <- function(A){
    result <- FALSE
    if(nrow(A)==ncol(A)){
      result <- sum(abs(t(A)-A))==0
    }
    return(result)
  }
  is.identity.matrix <- function(A){
    result <- FALSE
    if(nrow(A)==ncol(A)&min(A)==0&max(A)==1){
      if(prod(diag(A))==1&sum(A-diag(nrow(A)))==0) result <- TRUE
    }
    return(result)
  }
  z <- function(x){
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
    return(x)
  }
  optimize.B.from.Y <- function(result,Y,gamma,epsilon,maxit,method){
    X <- result$X
    C <- matrix(1,nrow=ncol(X),ncol=ncol(Y))
    oldSum <- 0
    for(l in 1:maxit){
      B <- C
      XB <- X %*% B
      if(method=="EU"){
        C <- C*z((t(X)%*%Y)/(t(X)%*%XB+gamma*C))
      }else{
        C <- C*(t(X)%*%z(Y/XB)/(colSums(X)%o%rep(1,ncol(Y))+2*gamma*C))
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
  is.identity <- is.identity.matrix(A)
  is.symmetric.matrix <- is.symmetric.matrix(A)
  n <- ncol(Y)
  (remainder <- n %% div)
  (division <- n %/% div)
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
    if(is.symmetric.matrix){
      A_j <- A[index!=j,index!=j] # kernel matrix or identity matrix
    }else{
      A_j <- A[,index!=j] # ordinary design matrix
    }
    res_j <- nmfkc(Y_j,A_j,Q,gamma,epsilon,maxit,method,X.restriction,nstart,seed,prefix,print.trace,print.dims,save.time,save.memory)
    if(is.identity){
      resj <- optimize.B.from.Y(res_j,Yj,gamma,epsilon,maxit,method)
      XBj <- resj$XB
    }else{
      if(is.symmetric.matrix){
        Aj <- A[index!=j,index==j]
      }else{
        Aj <- A[,index==j]
      }
      XBj <- res_j$X %*% res_j$C %*% Aj
    }
    if(method=="EU"){
      objfunc.block[j] <- sum((Yj-XBj)^2)+gamma*sum(res_j$C^2)
    }else{
      objfunc.block[j] <- sum(-Yj*z(log(XBj))+XBj)+gamma*sum(res_j$C^2)
    }
  }
  objfunc <- sum(objfunc.block)
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
#' }
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

nmfkc.rank <- function(Y,A=NULL,Q=2:min(5,ncol(Y),nrow(Y)),plot=TRUE,...){
  arglist=list(...)
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
  gamma <- ifelse("gamma" %in% names(arglist),arglist$gamma,0)
  epsilon <- ifelse("epsilon" %in% names(arglist),arglist$epsilon,1e-4)
  maxit <- ifelse("maxit" %in% names(arglist),arglist$maxit,5000)
  method <- ifelse("method" %in% names(arglist),arglist$method,"EU")
  X.restriction <- ifelse("X.restriction" %in% names(arglist),arglist$X.restriction,"colSums")
  seed <- ifelse("seed" %in% names(arglist),arglist$nstart,123)
  nstart <- ifelse("nstart" %in% names(arglist),arglist$nstart,1)
  prefix <- ifelse("prefix" %in% names(arglist),arglist$prefix,"Basis")
  print.trace <- ifelse("print.trace" %in% names(arglist),arglist$print.trace,FALSE)
  print.dims <- ifelse("print.dims" %in% names(arglist),arglist$print.dims,TRUE)
  save.time <- ifelse("save.time" %in% names(arglist),arglist$save.time,TRUE)
  save.memory <- FALSE
  r.squared <- 0*Q; names(r.squared) <- Q
  ICp <- 0*Q; names(ICp) <- Q
  AIC <- 0*Q; names(AIC) <- Q
  BIC <- 0*Q; names(BIC) <- Q
  silhouette <- 0*Q; names(silhouette) <- Q
  CPCC <- 0*Q; names(CPCC) <- Q
  B.prob.sd.min <- 0*Q; names(B.prob.sd.min) <- Q
  ARI <- 0*Q; names(ARI) <- Q
  for(q in 1:length(Q)){
    if(save.time){
      result <- nmfkc(Y,A,Q=Q[q],gamma,epsilon,maxit,method,X.restriction,nstart,seed,prefix,print.trace,print.dims,save.time=T,save.memory)
      CPCC[q] <- NA
      silhouette[q] <- NA
    }else{
      result <- nmfkc(Y,A,Q=Q[q],gamma,epsilon,maxit,method,X.restriction,nstart,seed,prefix,print.trace,print.dims,save.time=F,save.memory)
      CPCC[q] <- result$criterion$CPCC
      silhouette[q] <- result$criterion$silhouette$silhouette.mean
    }
    r.squared[q] <- result$r.squared
    ICp[q] <- result$criterion$ICp
    AIC[q] <- result$criterion$AIC
    BIC[q] <- result$criterion$BIC
    if(q==1){
      ARI[q] <- NA
      cluster.old <- result$B.cluster
    }else{
      df <- data.frame(old=cluster.old,new=result$B.cluster)
      df <- df[stats::complete.cases(df),]
      f <- table(df$old,df$new)
      ARI[q] <- AdjustedRandIndex(f)$ARI
      cluster.old <- result$B.cluster
    }
    B.prob.sd.min[q] <- result$criterion$B.prob.sd.min
  }
  if(plot){
    # Criterion
    plot(Q,r.squared,type="l",col=2,xlab="Rank",ylab="Criterion",ylim=c(0,1),lwd=3)
    graphics::text(Q,r.squared,Q)
    legend <- "r.squared"
    fill <- 2
    graphics::lines(Q,B.prob.sd.min,col=3,lwd=3)
    graphics::text(Q,B.prob.sd.min,Q)
    legend <- c(legend,"B.prob.sd.min")
    fill <- c(fill,3)
    graphics::lines(Q[-1],ARI[-1],col=4,lwd=3)
    graphics::text(Q[-1],ARI[-1],Q[-1])
    legend <- c(legend,"ARI for Q-1")
    fill <- c(fill,4)
    if(!save.time){
      graphics::lines(Q,silhouette,col=7,lwd=3)
      graphics::text(Q,silhouette,Q)
      legend <- c(legend,"silhouette")
      fill <- c(fill,7)
      graphics::lines(Q,CPCC,col=6,lwd=3)
      graphics::text(Q,CPCC,Q)
      legend <- c(legend,"CPCC")
      fill <- c(fill,6)
    }
    graphics::legend("bottomleft",legend=legend,fill=fill,bg=NULL)
  }
  invisible(data.frame(Q,r.squared,ICp,AIC,BIC,B.prob.sd.min,ARI,silhouette,CPCC))
}


.onAttach <- function(...) {
  packageStartupMessage("Last update on 6th Feb 2024")
  packageStartupMessage("https://github.com/ksatohds/nmfkc")
}

#' @title Creating kernel matrix from covariates
#' @description \code{nmfkc.kernel} create kernel matrix A(R,N) from covariate matrix U(R,N)
#' @param U covariate matrix U(K,N)=(u_1,...,u_N) each row should be normalized in advance
#' @param V covariate matrix V(K,M)=(v_1,...,v_N) usually used for prediction, and the default value is U.
#' @param beta parameter of kernel matrix A of which element is defined by exp(-beta*|u_n-v_m|^2)
#' @return kernel matrix A(N,M)
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' U <- matrix(1:3,nrow=1,ncol=3)
#' print(U)
#' A <- nmfkc.kernel(U,beta=1)
#' print(A)
#' print(log(A))

nmfkc.kernel <- function(U,V=U,beta){
  kernel <- function(m){
    vm <- t(rep(1,ncol(U)) %o% V[,m])
    exp(-beta*colSums((U-vm)^2))}
  A <- NULL
  for(m in 1:ncol(V)) A <- cbind(A,kernel(m))
  return(A)
}

#' @title Optimizing NMF (Non-negative Matrix Factorization) with kernel covariates
#' @description \code{nmkcfreg} The goal of the package is to perform NMF (Non-negative Matrix Factorization) with kernel covariates described by Y(P,N)~X(P,Q)C(Q,R)A(R,N)
#'  where observation matrix Y(P,N),
#'  covariate matrix A(R,N),
#'  basis matrix X(P,Q) and Q<=min(P,N)
#'  and coefficient matrix C(Q,R).
#'  Note that Y(N,P) and A(R,N) are known, and X(P,Q) and C(Q,R) are unknown.
#' @param Y observation matrix
#' @param A covariate matrix. Without covariate, identity matrix is used.
#' Or matrix A(R,N) having N columns can be accepted.
#' kernel matrix A(N,N) can be created by nmfkc.kernel function.
#' @param Q rank of basis matrix and Q<=min(P,N)
#' @param gamma penalty parameter for C(Q,R) in objective function
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method default objective function is Euclid distance "EU", otherwise Kullback–Leibler divergence "KL"
#' @param X.column The default is X.column="sum" and the column sum of basis matrix is 1, and it is interpreted as probability. The column of basis matrix is unit vector when X.column="squared".
#' @param print.trace display current iteration every 10 times if print.trace=TRUE
#' @param print.dims display dimensions of matrix sizes if  print.dim=TRUE. The default is set by  print.dim=FALSE.
#' @return X(P,Q): basis matrix. The column sum depends on X.column.
#' @return B(Q,N): B(Q,N)=C(Q,R)A(R,N) regression coefficient matrix
#' @return B.prob(Q,N): probability matrix whose column sum is 1
#' for soft clustering based on regression coefficient matrix B(Q,N).
#' @return B.cluster: the number of the basis that takes the maximum value of each column of B.prob(Q,N)
#' for hard clustering
#' @return XB(P,N): XB(P,N)=X(P,Q)B(Q,N) prediction matrix or fitted values for observation matrix Y
#' @return C(Q,R): parameter matrix
#' @return objfunc: last objective function
#' @return objfunc.iter: objective function at each iteration
#' @return r.squared: coefficient of determination R^2, squared correlation between Y(P,N) and XB(P,N)
#' @export
#' @examples
#' library(nmfkc)
#' Y <- t(iris[,-5])
#' result <- nmfkc(Y,Q=2)
#' # visualization of some results
#' plot(result$objfunc.iter) # convergence
#'
#' # goodness of fit
#' plot(as.vector(result$XB),as.vector(Y),
#' main=paste0("r.squared=",round(result$r.squared,3)))
#' abline(a=0,b=1,col=2)
#'
#' # dimension reduction based on regression coefficient B
#' plot(t(result$B))

nmfkc <- function(Y,A=diag(ncol(Y)),Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU",X.column="sum",print.trace=FALSE,print.dims=TRUE){
  is.identity.matrix <- function(A){
    result <- FALSE
    if(nrow(A)==ncol(A)&min(A)==0&max(A)==1){
      if(prod(diag(A))==1&sum(A-diag(nrow(A)))==0) result <- TRUE
    }
    return(result)
  }
  set.seed(123)
  if(is.vector(Y)) Y <- t(as.matrix(Y))
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  if(min(A)<0) warning("Minimum value of A is negative. It should be a non-negative matrix.")
  if(min(Y)<0) warning("Minimum value of Y is negative. It should be a non-negative matrix.")
  if(nrow(Y)>=2 & sum(colSums(Y)==0)>0) warning("There is a column of Y of which elements are all zero. Some elements should be positive.")
  if(nrow(Y)>=2 & sum(rowSums(Y)==0)>0) warning("There is a row of Y of which elements are all zero. Some elements should be positive.")
  if(nrow(Y)>=2){
    if(min(nrow(Y),ncol(Y))>=Q){
      if(ncol(Y)==Q){
        X <- Y
      }else{
        res.kmeans <- stats::kmeans(t(Y),centers=Q,iter.max=maxit)
        X <- t(res.kmeans$centers)
      }
    }else{
      warning("It does not hold Q<=min(P,N) where dim(Y)=(P,N).")
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
  }
  if(X.column=="sum") X <- t(t(X)/colSums(X)) else X <- t(t(X)/colSums(X^2)^0.5)
  C <- matrix(1,nrow=ncol(X),ncol=nrow(A))
  B <- C %*% A
  XB <- X %*% B
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){
    if(print.trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))
    if(method=="EU"){
      X <- X * ((Y %*% t(B)) / (XB %*% t(B)))
      if(X.column=="sum") X <- t(t(X)/colSums(X)) else X <- t(t(X)/colSums(X^2)^0.5)
      C <- C*((t(X)%*%Y%*%t(A))/(t(X)%*%XB%*%t(A)+gamma*C))
      B <- C %*% A
      XB <- X %*% B
      objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(C^2)
    }else{
      X <- t(t(X*(Y/XB)%*%t(B))/rowSums(B))
      if(X.column=="sum") X <- t(t(X)/colSums(X)) else X <- t(t(X)/colSums(X^2)^0.5)
      C0 <- t(X)%*%(Y/XB)%*%t(A)
      C <- C*(C0/(colSums(X)%o%rowSums(A)+2*gamma*C))
      B <- C %*% A
      XB <- X %*% B
      objfunc.iter[i] <- sum(-Y*log(XB)+XB)+gamma*sum(C^2)
    }
    if(i>=10){
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
      if(epsilon.iter <= epsilon){
        objfunc.iter <- objfunc.iter[1:i]
        break
      }
    }
  }
  if(method=="EU"){
    objfunc <- sum((Y-XB)^2)+gamma*sum(C^2)
    #aic <- ncol(Y)*log(2*pi*sum((Y-XB)^2)/ncol(Y))+ncol(Y)+2*(ncol(X)*(nrow(X)-1)+ncol(C)*nrow(C))
  }else{
    objfunc <- sum(-Y*log(XB)+XB)+gamma*sum(C^2)
  }
  r2 <- stats::cor(as.vector(XB),as.vector(Y))^2
  colnames(B) <- colnames(Y)
  colnames(XB) <- colnames(Y)
  B.prob <- t(t(B)/colSums(B))
  B.cluster <- apply(B.prob,2,which.max)
  colnames(B.prob) <- colnames(Y)
  if(epsilon.iter > epsilon) warning(paste0(
    "maximum iterations (",maxit,
    ") reached and the optimization hasn't converged yet."))
  if(print.dims)if(is.identity.matrix(A)){
    packageStartupMessage(
      sprintf("Y(%d,%d)~X(%d,%d)B(%d,%d)",
              nrow(Y),ncol(Y),nrow(Y),Q,Q,ncol(Y)))
  }else{
    packageStartupMessage(
      sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)=XB(%d,%d)",
              nrow(Y),ncol(Y),nrow(Y),Q,Q,nrow(A),nrow(A),ncol(Y),Q,ncol(Y)))
  }
  return(list(X=X,B=B,B.prob=B.prob,B.cluster=B.cluster,XB=XB,C=C,
              objfunc=objfunc,objfunc.iter=objfunc.iter,r.squared=r2))
}

#' @title Performing k-fold cross validation on NMF (Non-negative Matrix Factorization) with kernel covariates
#' @description \code{nmfkc.cv} apply cross validation method for k-partitioned columns of Y(P,N)
#' @param Y observation matrix
#' @param A covariate matrix. Without covariate, identity matrix is used.
#' Or matrix A(R,N) having N columns can be accepted.
#' kernel matrix A(N,N) can be created by nmfkc.kernel function.
#' @param Q rank of basis matrix and Q<=min(P,N)
#' @param gamma penalty parameter for C(Q,R) in objective function
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method default objective function is Euclid distance "EU", otherwise Kullback–Leibler divergence "KL"
#' @param X.column The default is X.column="sum" and the column sum of basis matrix is 1, and it is interpreted as probability. The column of basis matrix is unit vector when X.column="squared".
#' @param div number of partition usually described as "k" of k-fold
#' @param seed integer used as argument in set.seed function
#'  which controls the reproducibility of the partition.
#' @return objfunc: last objective function
#' @return objfunc.block: objective function at each block
#' @return block: partition block index (1,...,div) assigned to each column of Y
#' @return r.squared: coefficient of determination R^2, squared correlation between Y(P,N) and XB(P,N)
#' @export
#' @examples
#' library(nmfkc)
#' Y <- t(iris[,-5])
#' result <- nmfkc.cv(Y,Q=2)
#' table(result$block)
#' result$objfunc

nmfkc.cv <- function(Y,A=diag(ncol(Y)),Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU",X.column="sum",div=5,seed=123){
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
  optimize.B.from.Y <- function(result,Y,gamma,epsilon,maxit,method){
    X <- result$X
    C <- matrix(1,nrow=ncol(X),ncol=ncol(Y))
    A <- diag(ncol(Y))
    oldSum <- 0
    for(l in 1:maxit){
      XB <- X %*% C
      if(method=="EU"){
        C <- C*((t(X)%*%Y)/(t(X)%*%XB+gamma*C))
      }else{
        C0 <- t(X)%*%(Y/XB)%*%t(A)
        C <- C*(C0/(colSums(X)%o%rowSums(A)+2*gamma*C))
      }
      newSum <- sum(C)
      if(l>=10){
        epsilon.iter <- abs(newSum-oldSum)/(abs(newSum)+0.1)
        if(epsilon.iter <= epsilon) break
      }
      oldSum <- sum(C)
    }
    B <- C
    XB <- X %*% B
    colnames(B) <- colnames(Y)
    colnames(XB) <- colnames(Y)
    colnames(C) <- rownames(A)
    rownames(C) <- colnames(X)
    if(epsilon.iter > epsilon) warning(paste0(
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
  Yvec <- NULL
  XBvec <- NULL
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
    res_j <- nmfkc(Y_j,A_j,Q,gamma,epsilon,maxit,method,X.column,print.dims=FALSE)
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
      objfunc.block[j] <- sum(-Yj*log(XBj)+XBj)+gamma*sum(res_j$C^2)
    }
    Yvec <- c(Yvec,as.vector(Yj))
    XBvec <- c(XBvec,as.vector(XBj))
  }
  objfunc <- sum(objfunc.block)
  r2 <- stats::cor(XBvec,Yvec)^2
  return(list(objfunc=objfunc,objfunc.block=objfunc.block,block=index,r.squared=r2))
}
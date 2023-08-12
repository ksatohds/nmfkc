#' @title Creating kernel matrix A from covariate matrix U
#' @description \code{create.kernel} create kernel matrix A from covariate matrix U
#' @param U covariate matrix U: KxN, each row should be normalized
#' @param beta parameter of kernel matrix A of which element is defined by exp(-beta*|u_i-u_j|^2)
#' @return kernel matrix A: NxN
#' @export
#' @examples
#' # A <- create.kernel(U,beta=5)

create.kernel <- function(U,beta){
  kernel <- function(n){
    un <- t(rep(1,ncol(U)) %o% U[,n])
    exp(-beta*colSums((U-un)^2))}
  A <- NULL
  for(n in 1:ncol(U)) A <- rbind(A,kernel(n))
  return(A)
}

#' @title Optimizing X and C on NMF:Y~XCA where Y and A are given
#' @description \code{nmfreg} optimize basis matrix X whose column sum is 1 and coefficient matrix C
#' @param Y observation matrix Y:PxN
#' @param A kernel matrix A:NxN
#' @param Q rank of basis matrix X:PxQ
#' @param gamma penalty parameter for C:QxN where objective function:tr(Y-YHAT)'(Y-YHAT)+gamma*trC'C for method="EU" and sum(-Y*log(YHAT)+YHAT)+gamma*sum(C^2) for method="KL"
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method default objective function uses Euclid distance "EU" and Kullback–Leibler divergence is used by method="KL".
#' @return X:NxQ whose column sum is 1
#' @return C:QxN parameter matrix
#' @return B=XC: QxN coefficient matrix
#' @return P:PxN, calculated from B as P <- t(t(B)/colSums(B)) and probability matrix whose column sum is 1
#' @return objfunc: last objective function
#' @return objfunc.iter: objective function at each iteration
#' @return r.squared: squared correlation, i.e., coefficient of determination R^2 between Y and YHAT
#' @export
#' @examples
#' # library(fda)
#' # data(CanadianWeather)
#' # d <- CanadianWeather$dailyAv[,,1]
#' # Y <- d-min(d)
#' # A <- diag(ncol(Y))
#' # library(nmfreg)
#' # result <- nmfreg(Y,A,Q=2)
#' # result$r.squared

nmfreg <- function(Y,A,Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU"){
  set.seed(123)
  res <- stats::kmeans(t(Y),centers=Q)
  X <- t(res$centers)
  X <- t(t(X)/colSums(X))
  C <- matrix(stats::rnorm(ncol(X)*nrow(A),mean=2,sd=0.3),
              nrow=ncol(X),ncol=nrow(A))
  B <- C %*% A
  YHAT <- X %*% B
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){
    if(method=="EU"){
      X <- X * ((Y %*% t(B)) / (YHAT %*% t(B)))
      X <- t(t(X)/colSums(X))
      C <- C*((t(X)%*%Y%*%t(A))/(t(X)%*%YHAT%*%t(A)+gamma*C))
      B <- C %*% A
      YHAT <- X %*% B
      objfunc.iter[i] <- sum((Y-YHAT)^2)+gamma*sum(C^2)
    }else{
      X <- t(t(X*(Y/YHAT)%*%t(B))/rowSums(B))
      X <- t(t(X)/colSums(X))
      C0 <- t(X)%*%(Y/YHAT)%*%t(A)
      C <- C*(C0/(colSums(X)%o%rowSums(A)+2*gamma*C))
      B <- C %*% A
      YHAT <- X %*% B
      objfunc.iter[i] <- sum(-Y*log(YHAT)+YHAT)+gamma*sum(C^2)
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
    objfunc <- sum((Y-YHAT)^2)+gamma*sum(C^2)
  }else{
    objfunc <- sum(-Y*log(YHAT)+YHAT)+gamma*sum(C^2)
  }
  r2 <- stats::cor(as.vector(YHAT),as.vector(Y))^2
  P <- t(t(B)/colSums(B))
  if(epsilon.iter > epsilon) print(paste0(
    "maximum iterations (",maxit,
    ") reached and the optimization hasn't converged yet."))
  return(list(X=X,C=C,B=B,YHAT=YHAT,P=P,
              objfunc=objfunc,objfunc.iter=objfunc.iter,r.squared=r2))
}

#' @title Performing k-fold cross validation
#' @description \code{nmfreg.cv} apply cv method for k-pertitioned columns of Y:PxN
#' @param Y observation matrix Y:PxN
#' @param A kernel matrix A:NxN
#' @param Q rank of basis matrix X:P*Q
#' @param gamma penalty parameter for C:QxN, where objective function:tr(Y-YHAT)'(Y-YHAT)+gamma*trC'C
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param div number of partition, usually referred to as "k"
#' @param seed random seed which decides the partition at random
#' @param method default objective function uses Euclid distance "EU" and Kullback–Leibler divergence is used by method="KL".
#' @return objfunc: last objective function
#' @return objfunc.block: objective function at each block
#' @return block: partition block index {1,...,div} assigned to each column of Y
#' @export
#' @examples
#' # library(fda)
#' # data(CanadianWeather)
#' # d <- CanadianWeather$dailyAv[,,1]
#' # Y <- d-min(d)
#' # A <- diag(ncol(Y))
#' # library(nmfreg)
#' # result.cv <- nmfreg.cv(Y,A,Q=2)

nmfreg.cv <- function(Y,A,Q,gamma=0,epsilon=1e-4,maxit=5000,div=5,seed=123,method="EU"){
  is.identity.matrix <- function(A){
    result <- FALSE
    if(nrow(A)==ncol(A)&min(A)==0&max(A)==1){
      if(prod(diag(A))==1&sum(A-diag(nrow(A)))==0) result <- TRUE
    }
    return(result)
  }
  is.identity <- is.identity.matrix(A)
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
    A_j <- A[index!=j,index!=j]
    res <- nmfreg(Y_j,A_j,Q,gamma,epsilon,maxit)
    X_j <- res$X
    C_j <- res$C
    if(is.identity){
      C_j <- matrix(stats::rnorm(ncol(X_j)*ncol(Yj),mean=2,sd=0.3),
                    nrow=ncol(X_j),ncol=ncol(Yj))
      oldSum <- 0
      for(l in 1:maxit){
        YHATj <- X_j %*% C_j
        C_j <- C_j*((t(X_j)%*%Yj)/(t(X_j)%*%YHATj+gamma*C_j))
        newSum <- sum(C_j)
        if(l>=10){
          epsilon.iter <- abs(newSum-oldSum)/(abs(newSum)+0.1)
          if(epsilon.iter <= epsilon) break
        }
        oldSum <- sum(C_j)
      }
      YHATj <- X_j %*% C_j
    }else{
      YHATj <- X_j %*% C_j %*% A[index!=j,index==j] # Aのサイズに注意！
    }
    objfunc.block[j] <- sum((Yj-YHATj)^2)+gamma*sum(C_j^2)
  }
  objfunc <- sum(objfunc.block)
  return(list(objfunc=objfunc,objfunc.block=objfunc.block,block=index))
}

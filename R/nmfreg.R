#' @title Creating kernel matrix A from covariate matrix U
#' @description \code{create.kernel} create kernel matrix A from covariate matrix U
#' @param U covariate matrix U: KxN, each row should be normalized
#' @param beta parameter of kernel matrix A of which element is defined by exp(-beta*|u_i-u_j|^2)
#' @return kernel matrix A: NxN
#' @export
#' @examples
#' # A <- create.kernel(U,beta=10)

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
#' @param gamma penalty parameter for C:QxN, where objective function:tr(Y-YHAT)'(Y-YHAT)+gamma*trC'C
#' @param iter number of iterations
#' @return X:NxQ whose column sum is 1
#' @return C:QxN parameter matrix
#' @return B=XC: QxN coefficient matrix
#' @return P:PxN, calculated from B as P <- t(t(B)/colSums(B)) and probability matrix whose column sum is 1
#' @return objective.function: last objective function
#' @return objective.function.iter: objective function at each iteration
#' @return r.squared: squared correlation, i.e., coefficient of determination R^2 between Y and YHAT
#' @export
#' @examples
#' # library(fda)
#' # data(CanadianWeather)
#' # d <- CanadianWeather$dailyAv[,,1]
#' # Y <- d-min(d)
#' # A <- diag(ncol(Y))
#' #
#' # library(nmfreg)
#' # result <- nmfreg(Y,A,Q=2)
#' # result$r.squared

nmfreg <- function(Y,A,Q=2,gamma=0,iter=500){
  set.seed(123)
  res <- stats::kmeans(t(Y),centers=Q)
  X <- t(res$centers)
  X <- t(t(X)/colSums(X))
  C <- matrix(stats::rnorm(ncol(X)*nrow(A),mean=2,sd=0.3),
               nrow=ncol(X),ncol=nrow(A))
  B <- C %*% A
  YHAT <- X %*% B
  objective.function.iter <- 0*(1:iter)
  for(i in 1:iter){
    X <- X * ((Y %*% t(B)) / (YHAT %*% t(B)))
    X <- t(t(X)/colSums(X))
    C <- C*((t(X)%*%Y%*%t(A))/(t(X)%*%YHAT%*%t(A)+gamma*C))
    B <- C %*% A
    YHAT <- X %*% B
    objective.function.iter[i] <- sum((Y-YHAT)^2)+gamma*sum(C^2)
  }
  P <- t(t(B)/colSums(B))
  objective.function <- sum((Y-YHAT)^2)+gamma*sum(C^2)
  r2 <- stats::cor(as.vector(YHAT),as.vector(Y))^2
  return(list(X=X,C=C,B=B,YHAT=YHAT,P=P,
              objective.function=objective.function,
              objective.function.iter=objective.function.iter,
              r.squared=r2))
}

#' @title Checking if kernel matrix is identity matrix
#' @description \code{is.identity.matrix} check if kernel matrix A is identity matrix
#' @param A kernel matrix A:NxN
#' @return logical value
#' @export
#' @examples
#' # if(is.identity.matrix(A)) print("A is identity matrix") else print("A is not identity matrix")

is.identity.matrix <- function(A){
  result <- FALSE
  if(min(A)==0&max(A)==1){
    result <- TRUE
    for(i in 1:nrow(A))
      for(j in 1:ncol(A))
        if(i==j){
          if(A[i,j]!=1) result <- FALSE
        }else{
          if(A[i,j]!=0) result <- FALSE
        }
  }
  return(result)
}

#' @title Finding best seed for nmfreg.cv
#' @description \code{nmfreg.cv.seed} find best seed for k-fold cross validation
#' @param n number of columns of Y: PxN
#' @param div number of partition, usually referred to as "k"
#' @param iter number of iterations
#' @return seed: best seed which minimizes sd between group index
#' @return group: partition index between 1 and div for each column of Y
#' @export
#' @examples
#' # best.seed  <- nmfreg.cv.seed(n=ncol(Y),div=10)
#' # result.cv <- nmfreg.cv(Y,A,seed=best.seed)

nmfreg.cv.seed <- function(n,div=5,iter=500){
  sds <- 0*(1:iter)
  for(i in 1:iter){
    set.seed(i)
    index <- sample(1:div,n,replace=T)
    f <- table(index)
    sds[i] <- stats::sd(f)
  }
  return(which.min(sds))
}

#' @title Performing k-fold cross validation
#' @description \code{nmfreg.cv} apply cv method for k-pertitioned columns of Y:PxN
#' @param Y observation matrix Y:PxN
#' @param A kernel matrix A:NxN
#' @param Q rank of basis matrix X:P*Q
#' @param gamma penalty parameter for C:QxN, where objective function:tr(Y-YHAT)'(Y-YHAT)+gamma*trC'C
#' @param iter number of iterations
#' @param div number of partition, usually referred to as "k"
#' @param seed random seed which decides the partition at random
#' @return X:NxQ, C:QxN, B=XC, YHAT=XB
#' @return objective.function: last objective function
#' @return objective.function.iter: objective function at each iteration
#' @return group: partition index between 1 and div for each column of Y
#' @export
#' @examples
#' # best.seed  <- nmfreg.cv.seed(n=ncol(Y),div=10)
#' # result.cv <- nmfreg.cv(Y,A,seed=best.seed)

nmfreg.cv <- function(Y,A,Q,gamma=0,iter=500,div=5,seed=123){
  set.seed(seed)
  index <- sample(1:div,ncol(Y),replace=T)
  objective.function.partition <- 0*(1:div)
  for(j in 1:div){
    Y_j <- Y[,index!=j]
    Yj <- Y[,index==j]
    A_j <- A[index!=j,index!=j]
    res <- nmfreg(Y_j,A_j,Q,gamma,iter)
    X_j <- res$X
    C_j <- res$C
    if(is.identity.matrix(A)){
      C_j <- matrix(stats::rnorm(ncol(X_j)*ncol(Yj),mean=2,sd=0.3),
                    nrow=ncol(X_j),ncol=ncol(Yj))
      for(l in 1:iter){
        YHATj <- X_j %*% C_j
        C_j <- C_j*((t(X_j)%*%Yj)/(t(X_j)%*%YHATj+gamma*C_j))
      }
      YHATj <- X_j %*% C_j
    }else{
      YHATj <- X_j %*% C_j %*% A[index!=j,index==j] # Aのサイズに注意！
    }
    objective.function.partition[j] <- sum((Yj-YHATj)^2)+gamma*sum(C_j^2)
  }
  objective.function <- sum(objective.function.partition)
  return(list(objective.function=objective.function,
              objective.function.partition=objective.function.partition,
              group=index))
}

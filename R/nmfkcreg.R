.onAttach <- function(...) {
  packageStartupMessage("This package is a development version.")
  packageStartupMessage("https://github.com/ksatohds/nmfkcreg")
}

#' @title Creating kernel matrix from covariates
#' @description \code{create.kernel} create kernel matrix A(N,N) from covariate matrix U(R,N)
#' @param U covariate matrix U(K,N)=(u_1,...,u_N) each row should be normalized in advance
#' @param beta parameter of kernel matrix A of which element is defined by exp(-beta*|u_i-u_j|^2)
#' @return kernel matrix A(N,N)
#' @export
#' @examples
#' # install.packages("devtools")
#' # devtools::install_github("ksatohds/nmfkcreg")
#' U <- matrix(1:3,nrow=1,ncol=3)
#' print(U)
#' A <- create.kernel(U,beta=1)
#' print(A)
#' print(log(A))

create.kernel <- function(U,beta){
  kernel <- function(n){
    un <- t(rep(1,ncol(U)) %o% U[,n])
    exp(-beta*colSums((U-un)^2))}
  A <- NULL
  for(n in 1:ncol(U)) A <- cbind(A,kernel(n))
  return(A)
}

#' @title Optimizing NMF (Non-negative Matrix Factorization) with kernel covariates regression
#' @description \code{nmkcfreg} The goal of the package is to perform NMF (Non-negative Matrix Factorization) with kernel covariates regression described by Y(P,N)~X(P,Q)C(Q,N)A(N,N)
#'  where observation matrix Y(P,N),
#'  kernel matrix A(N,N) with parameter beta,
#'  basis matrix X(P,Q) whose column sum is 1 and Q<=min{P,N}
#'  and coefficient matrix C(Q,N).
#'  Note that Y(N,P) and A(N,N) are known, and X(P,Q) and C(Q,N) are unknown.
#' @param Y observation matrix
#' @param A kernel matrix. Without covariate, identity matrix is used. Or matrix A(R,N) having N columns can be accepted.
#' @param Q rank of basis matrix and Q<=min{P,N}
#' @param gamma penalty parameter for C(Q,N) in objective function
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method default objective function is Euclid distance "EU", otherwise Kullback–Leibler divergence "KL"
#' @param trace display current iteration every 10 times if trace=TRUE
#' @return X(P,Q): basis matrix whose column sum is 1
#' @return C(Q,N): parameter matrix
#' @return B(Q,N): B(Q,N)=C(Q,N)A(N,N) regression coefficient matrix
#' @return XB(P,N): XB(P,N)=X(P,Q)B(Q,N) prediction matrix or fitted values for observation matrix Y
#' @return P(Q,N): probability matrix whose column sum is 1
#' for soft clustering based on regression coefficient matrix B(Q,N).
#' @return cluster: the number of the basis that takes the maximum value of each column of P(Q,N)
#' for hard clustering
#' @return objfunc: last objective function
#' @return objfunc.iter: objective function at each iteration
#' @return r.squared: coefficient of determination R^2, squared correlation between Y(P,N) and XB(P,N)
#' @export
#' @examples
#' library(nmfkcreg)
#' Y <- t(iris[,-5])
#' result <- nmfkcreg(Y,Q=2)
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

nmfkcreg <- function(Y,A=diag(ncol(Y)),Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU",trace=FALSE){
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
        res.kmeans <- stats::kmeans(t(Y),centers=Q)
        X <- t(res.kmeans$centers)
      }
    }else{
      warning("It does not hold Q<=min{P,N} where dim(Y)=(P,N).")
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
  }
  X <- t(t(X)/colSums(X))
  C <- matrix(1,nrow=ncol(X),ncol=nrow(A))
  B <- C %*% A
  XB <- X %*% B
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){
    if(trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))
    if(method=="EU"){
      X <- X * ((Y %*% t(B)) / (XB %*% t(B)))
      X <- t(t(X)/colSums(X))
      C <- C*((t(X)%*%Y%*%t(A))/(t(X)%*%XB%*%t(A)+gamma*C))
      B <- C %*% A
      XB <- X %*% B
      objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(C^2)
    }else{
      X <- t(t(X*(Y/XB)%*%t(B))/rowSums(B))
      X <- t(t(X)/colSums(X))
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
  }else{
    objfunc <- sum(-Y*log(XB)+XB)+gamma*sum(C^2)
  }
  r2 <- stats::cor(as.vector(XB),as.vector(Y))^2
  colnames(B) <- colnames(Y)
  colnames(XB) <- colnames(Y)
  P <- t(t(B)/colSums(B))
  cluster <- apply(P,2,which.max)
  colnames(P) <- colnames(Y)
  if(epsilon.iter > epsilon) warning(paste0(
    "maximum iterations (",maxit,
    ") reached and the optimization hasn't converged yet."))
  return(list(X=X,C=C,B=B,XB=XB,P=P,cluster=cluster,
              objfunc=objfunc,objfunc.iter=objfunc.iter,r.squared=r2))
}

#' @title Performing k-fold cross validation on NMF (Non-negative Matrix Factorization) with kernel covariates regression
#' @description \code{nmfkcreg.cv} apply cross validation method for k-partitioned columns of Y(P,N)
#' @param Y observation matrix
#' @param A kernel matrix. Without covariates, identity matrix is used.  Or matrix A(R,N) having N columns can be accepted.
#' @param Q rank of basis matrix and Q<=min{P,N}
#' @param gamma penalty parameter for C(Q,N) in objective function
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method default objective function is Euclid distance "EU", otherwise Kullback–Leibler divergence "KL"
#' @param div number of partition usually described as "k" of k-fold
#' @param seed integer used as argument in set.seed function
#'  which controls the reproducibility of the partition.
#' @return objfunc: last objective function
#' @return objfunc.block: objective function at each block
#' @return block: partition block index {1,...,div} assigned to each column of Y
#' @export
#' @examples
#' library(nmfkcreg)
#' Y <- t(iris[,-5])
#' result <- nmfkcreg.cv(Y,Q=2)
#' table(result$block)
#' result$objfunc

nmfkcreg.cv <- function(Y,A=diag(ncol(Y)),Q=2,gamma=0,epsilon=1e-4,maxit=5000,div=5,seed=123,method="EU"){
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
    res <- nmfkcreg(Y_j,A_j,Q,gamma,epsilon,maxit,method)
    X_j <- res$X
    C_j <- res$C
    if(is.identity){
      A_j <- diag(ncol(Yj))
      C_j <- matrix(1,nrow=ncol(X_j),ncol=ncol(Yj))
      oldSum <- 0
      for(l in 1:maxit){
        XBj <- X_j %*% C_j
        if(method=="EU"){
          C_j <- C_j*((t(X_j)%*%Yj)/(t(X_j)%*%XBj+gamma*C_j))
        }else{
          C0_j <- t(X_j)%*%(Yj/XBj)%*%t(A_j)
          C_j <- C_j*(C0_j/(colSums(X_j)%o%rowSums(A_j)+2*gamma*C_j))
        }
        newSum <- sum(C_j)
        if(l>=10){
          epsilon.iter <- abs(newSum-oldSum)/(abs(newSum)+0.1)
          if(epsilon.iter <= epsilon) break
        }
        oldSum <- sum(C_j)
      }
      XBj <- X_j %*% C_j
    }else{
      XBj <- X_j %*% C_j %*% A[index!=j,index==j] # Aのサイズに注意！
    }
    if(method=="EU"){
      objfunc.block[j] <- sum((Yj-XBj)^2)+gamma*sum(C_j^2)
    }else{
      objfunc.block[j] <- sum(-Yj*log(XBj)+XBj)+gamma*sum(C_j^2)
    }
  }
  objfunc <- sum(objfunc.block)
  return(list(objfunc=objfunc,objfunc.block=objfunc.block,block=index))
}

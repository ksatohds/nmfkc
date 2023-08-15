#' @title Creating kernel matrix A from covariate matrix U
#' @description \code{create.kernel} create kernel matrix A from covariate matrix U
#' @param U covariate matrix U(K,N)=(u_1,...,u_N) each row should be normalized in advance
#' @param beta parameter of kernel matrix A of which element is defined by exp(-beta*|u_i-u_j|^2)
#' @return kernel matrix A(N,N)
#' @export
#' @examples
#' #library(nmfkcreg)
#' #library(fda)
#' #data(CanadianWeather)
#' #d <- CanadianWeather$dailyAv[,,1]
#' #Y <- d-min(d)
#' #u0 <- CanadianWeather$coordinates[,2:1]
#' #u0[,1] = -u0[,1]
#' #u = t(u0)
#' #umin <- apply(u,1,min)
#' #umax <- apply(u,1,max)
#' # install.packages("devtools")
#' # devtools::install_github("ksatohds/nmfkcreg")
#' # A <- create.kernel(U,beta=5)

create.kernel <- function(U,beta){
  kernel <- function(n){
    un <- t(rep(1,ncol(U)) %o% U[,n])
    exp(-beta*colSums((U-un)^2))}
  A <- NULL
  for(n in 1:ncol(U)) A <- rbind(A,kernel(n))
  return(A)
}

#' @title Optimizing X and C on NMF (Non-negative Matrix Factorization) kernel covariate regression model Y~XCA where Y and A are given
#' @description \code{nmkcfreg} The goal of nmfkcreg is to perform NMF (Non-negative Matrix Factorization) regression model described by Y~XCA
#'  where observation matrix Y(P,N),
#'  kernel matrix A(N,N) with parameter beta,
#'  basis matrix X(P,Q) whose column sum is 1 and Q<=min{P,N}
#'  and coefficient matrix C(Q,N).
#'  Note that Y and A are known, and X and C are unknown.
#' @param Y observation matrix
#' @param A kernel matrix. Without covariate, use identity matrix A=diag(ncol(Y)). Or matrix A(R,N) having N columns can be accepted.
#' @param Q rank of basis matrix and Q<=P
#' @param gamma penalty parameter for C:QxN where
#' objective function:tr(Y-YHAT)'(Y-YHAT)+gamma*trC'C for method="EU"
#' and sum(-Y*log(YHAT)+YHAT)+gamma*sum(C^2) for method="KL"
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method default objective function is Euclid distance "EU", otherwise Kullback–Leibler divergence "KL"
#' @param trace display current iteration every 10 times if trace=TRUE
#' @return X(N,Q) whose column sum is 1
#' @return C(Q,N) parameter matrix
#' @return B, B=XC regression coefficient matrix
#' @return YHAT(P,N), YHAT=XB=XCA prediction matrix or fitted values for observation matrix Y
#' @return P(P,N) probability matrix
#' for soft clustering based on regression coeffient matrix B.
#' It is given by P <- t(t(B)/colSums(B)) whose column sum is 1.
#' @return objfunc: last objective function
#' @return objfunc.iter: objective function at each iteration
#' @return r.squared: coefficient of determination R^2, squared correlation between Y and YHAT
#' @export
#' @examples
#' # library(fda)
#' # data(CanadianWeather)
#' # d <- CanadianWeather$dailyAv[,,1]
#' # Y <- d-min(d)
#' # A <- diag(ncol(Y))
#' # install.packages("devtools")
#' # devtools::install_github("ksatohds/nmfkcreg")
#' # library(nmfkcreg)
#' # result <- nmfkcreg(Y,A,Q=2)
#' # result$r.squared

nmfkcreg <- function(Y,A,Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU",trace=FALSE){
  set.seed(123)
  if(nrow(Y)>=2){
    if(nrow(Y)>=Q){
      if(nrow(Y)==Q){
        print("Q==P")
        X <- Y
      }else{
        print("Q<P")
        res.kmeans <- stats::kmeans(t(Y),centers=Q)
        X <- t(res.kmeans$centers)
      }
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
  }
  X <- t(t(X)/colSums(X))
  C <- matrix(stats::rnorm(ncol(X)*nrow(A),mean=2,sd=0.3),
              nrow=ncol(X),ncol=nrow(A))
  B <- C %*% A
  YHAT <- X %*% B
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){
    if(trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))
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

#' @title Performing k-fold cross validation on NMF (Non-negative Matrix Factorization) kernel covariate regression model
#' @description \code{nmfkcreg.cv} apply cv method for k-partitioned columns of Y on NMF (Non-negative Matrix Factorization) regression model
#' @param Y observation matrix
#' @param A kernel matrix. Without covariate, use identity matrix A=diag(ncol(Y)).  Or matrix A(R,N) having N columns can be accepted.
#' @param Q rank of basis matrix and Q<=min{P,N}
#' @param gamma penalty parameter for C:QxN where
#' objective function:tr(Y-YHAT)'(Y-YHAT)+gamma*trC'C for method="EU"
#' and sum(-Y*log(YHAT)+YHAT)+gamma*sum(C^2) for method="KL"
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
#' #library(nmfkcreg)
#' #library(fda)
#' #data(CanadianWeather)
#' #d <- CanadianWeather$dailyAv[,,1]
#' #Y <- d-min(d)
#' #u0 <- CanadianWeather$coordinates[,2:1]
#' #u0[,1] = -u0[,1]
#' #u = t(u0)
#' #umin <- apply(u,1,min)
#' #umax <- apply(u,1,max)
#' # install.packages("devtools")
#' # devtools::install_github("ksatohds/nmfkcreg")
#' # A <- create.kernel(U,beta=5)
#' # result.cv <- nmfkcreg.cv(Y,A,Q=2)

nmfkcreg.cv <- function(Y,A,Q,gamma=0,epsilon=1e-4,maxit=5000,div=5,seed=123,method="EU"){
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
      C_j <- matrix(stats::rnorm(ncol(X_j)*ncol(Yj),mean=2,sd=0.3),
                    nrow=ncol(X_j),ncol=ncol(Yj))
      oldSum <- 0
      for(l in 1:maxit){
        YHATj <- X_j %*% C_j
        if(method=="EU"){
          C_j <- C_j*((t(X_j)%*%Yj)/(t(X_j)%*%YHATj+gamma*C_j))
        }else{
          C0_j <- t(X_j)%*%(Yj/YHATj)%*%t(A_j)
          C_j <- C_j*(C0_j/(colSums(X_j)%o%rowSums(A_j)+2*gamma*C_j))
        }
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
    if(method=="EU"){
      objfunc.block[j] <- sum((Yj-YHATj)^2)+gamma*sum(C_j^2)
    }else{
      objfunc.block[j] <- sum(-Yj*log(YHATj)+YHATj)+gamma*sum(C_j^2)
    }
  }
  objfunc <- sum(objfunc.block)
  return(list(objfunc=objfunc,objfunc.block=objfunc.block,block=index))
}

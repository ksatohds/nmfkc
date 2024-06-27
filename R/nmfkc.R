.onAttach <- function(...) {
  packageStartupMessage("Last update on 25th Jun 2024")
  packageStartupMessage("https://github.com/ksatohds/nmfkc")
}


#' @title Creating kernel matrix from covariates
#' @description \code{nmfkc.kernel} create kernel matrix from covariate matrix
#' @param U covariate matrix U(K,N)=(u_1,...,u_N) each row might be normalized in advance
#' @param V covariate matrix V(K,M)=(v_1,...,v_M) usually used for prediction, and the default value is U.
#' @param method The default kernel function is Gaussian kernel. For other functions, check by typing "nmfkc.kernel".
#' @param beta The default parameter of kernel function is 0.5.
#' @return kernel matrix A(N,M)
#' @export
#' @source Satoh, K. (2024) Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv preprint arXiv:2403.05359. \url{https://arxiv.org/abs/2403.05359}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example 1.
#' U <- matrix(1:3,nrow=1,ncol=3)
#' print(U)
#' A <- nmfkc.kernel(U,beta=1)
#' print(A)
#' print(log(A))
#'
#' # Example 2.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' A <- nmfkc.kernel(U,V,beta=0.031) # see, nmfkc.cv
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.kernel <- function(U,V=U,method="Gaussian",beta=0.5){
  kernel <- function(m){
    vm <- t(rep(1,ncol(U)) %o% V[,m])
    d <- colSums((U-vm)^2)^0.5
    k <- 0
    if(method=="Gaussian") k <- exp(-beta*d^2) # Gaussian
    if(method=="Exponential") k <- exp(-beta*d)
    if(method=="Periodic") k <- exp(-beta[1]*sin(beta[2]*d)^2)
    if(method=="Linear") k <- t(U) %*% V[,m]
    if(method=="NormalizedLinear") k <- diag(1/colSums(U^2)^0.5) %*% t(U) %*% V[,m]/sum(V[,m]^2)^0.5
    return(k)}
  A <- NULL; for(m in 1:ncol(V)) A <- cbind(A,kernel(m))
  if(min(A)<0){
    warning("The constructed matrix is not non-negative.")
  }
  return(A)
}


#' @title Optimizing NMF (Non-negative Matrix Factorization) with Kernel Covariate
#' @description \code{nmfkc} The goal of the package is to perform NMF (Non-negative Matrix Factorization) with Kernel Covariate described by Y~XCA=XB
#'  where observation matrix Y(P,N),
#'  covariate matrix A(R,N),
#'  basis matrix X(P,Q) and Q<=min(P,N),
#'  parameter matrix C(Q,R)
#'  and coefficient matrix B(Q,N).
#'  Note that Y(N,P) and A(R,N) are given, and X(P,Q) and C(Q,R) are unknown.
#' @param Y observation matrix
#' @param A covariate matrix. Without covariate, identity matrix is used.
#' @param Q rank of basis matrix and Q<=min(P,N)
#' @param gamma penalty parameter for parameter matrix C
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method The default objective function is Euclid distance "EU", otherwise Kullback–Leibler divergence "KL"
#' @param X.column The default is X.column="sum" and the column sum of basis matrix is 1, and it is interpreted as probability. The column of basis matrix is unit vector when X.column="squared".
#' @param nstart The default is one. It is the "nstart" option of "kmeans" function used for the initial values of basis matrix.
#' @param print.trace display current iteration every 10 times if print.trace=TRUE
#' @param print.dims display dimensions of matrix sizes if print.dim=TRUE. The default is set by  print.dim=FALSE.
#' @return X: basis matrix. The column sum depends on X.column.
#' @return B: coefficient matrix, B=CA
#' @return B.prob: probability matrix for soft clustering based on coefficient matrix B. Those column sum is 1.
#' @return B.cluster: the number of the basis that takes the maximum value of each column of B.prob for hard clustering
#' @return XB: fitted values for observation matrix Y
#' @return C: parameter matrix
#' @return objfunc: last objective function
#' @return objfunc.iter: objective function at each iteration
#' @return r.squared: coefficient of determination R^2, squared correlation between Y and XB
#' @export
#' @source Satoh, K. (2024) Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv preprint arXiv:2403.05359. \url{https://arxiv.org/abs/2403.05359}
#' @references Ding, C., Li, T., Peng, W. and Park, H. (2006) Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining, 126-135. \url{https://doi.org/10.1145/1150402.1150420}
#' @references Potthoff, R.F. and Roy, S.N. (1964). A generalized multivariate analysis of variance model useful especially for growth curve problems. Biometrika, 51, 313-326. \url{https://doi.org/10.2307/2334137}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example 1.
#' library(nmfkc)
#' Y <- t(iris[,-5])
#' Q <- 2
#' result <- nmfkc(Y,Q=Q)
#' # visualization of some results
#' plot(result$objfunc.iter) # convergence
#' # goodness of fit
#' plot(as.vector(result$XB),as.vector(Y),
#' main=paste0("r.squared=",round(result$r.squared,3)))
#' abline(a=0,b=1,col=2)
#' # dimension reduction based on coefficient matrix B
#' plot(t(result$B))
#' # soft clustering based on coefficient matrix B
#' plot(t(result$B),col=as.numeric(iris$Species))
#' legend("topright",legend=1:Q,fill=1:Q+1)
#' stars(t(result$B.prob),locations=t(result$B),scale=FALSE,
#'   draw.segments=TRUE,col.segments=1:Q+1,len=0.2,add=TRUE)
#'
#' # Example 2.
#' Y <- matrix(cars$dist,nrow=1)
#' A <- rbind(1,cars$speed)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(A[2,]),as.vector(Y))
#' lines(as.vector(A[2,]),as.vector(result$XB),col=2,lwd=2)

nmfkc <- function(Y,A=diag(ncol(Y)),Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU",X.column="sum",nstart=1,print.trace=FALSE,print.dims=TRUE){
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
  set.seed(123)
  if(is.vector(Y)) Y <- t(as.matrix(Y))
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  if(min(A)<0){
    warning("The matrix A should be non-negative.")
    stop()
  }
  if(min(Y)<0){
    warning("The matrix Y should be non-negative.")
    stop()
  }
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
  }
  if(X.column=="sum") X <- t(t(X)/colSums(X)) else X <- t(t(X)/colSums(X^2)^0.5)
  C <- matrix(1,nrow=ncol(X),ncol=nrow(A))
  B <- C %*% A
  XB <- X %*% B
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){
    if(print.trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))
    if(method=="EU"){
      X <- X*z((Y%*%t(B))/(XB%*%t(B)))
      if(X.column=="sum") X <- t(t(X)/colSums(X)) else X <- t(t(X)/colSums(X^2)^0.5)
      C <- C*z((t(X)%*%Y%*%t(A))/(t(X)%*%XB%*%t(A)+gamma*C))
      B <- C %*% A
      XB <- X %*% B
      objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(C^2)
    }else{
      X <- t(t(X*z(Y/XB)%*%t(B))/rowSums(B))
      if(X.column=="sum") X <- t(t(X)/colSums(X)) else X <- t(t(X)/colSums(X^2)^0.5)
      C0 <- t(X)%*%z(Y/XB)%*%t(A)
      C <- C*(C0/(colSums(X)%o%rowSums(A)+2*gamma*C))
      B <- C %*% A
      XB <- X %*% B
      objfunc.iter[i] <- sum(-Y*z(log(XB))+XB)+gamma*sum(C^2)
    }
    if(i>=10){
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
      if(epsilon.iter <= abs(epsilon)){
        objfunc.iter <- objfunc.iter[1:i]
        break
      }
    }
  }
  if(method=="EU"){
    objfunc <- sum((Y-XB)^2)+gamma*sum(C^2)
  }else{
    objfunc <- sum(-Y*z(log(XB))+XB)+gamma*sum(C^2)
  }
  r2 <- stats::cor(as.vector(XB),as.vector(Y))^2
  colnames(B) <- colnames(Y)
  colnames(XB) <- colnames(Y)
  B.prob <- t(z(t(B)/colSums(B)))
  B.cluster <- apply(B.prob,2,which.max)
  B.cluster[colSums(B.prob)==0] <- NA
  colnames(B.prob) <- colnames(Y)
  if(epsilon.iter > abs(epsilon)) warning(paste0(
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


#' @title Performing k-fold cross validation on NMF (Non-negative Matrix Factorization) with Kernel Covariate
#' @description \code{nmfkc.cv} apply cross validation method for k-partitioned columns of Y~XCA=XB
#'  where observation matrix Y(P,N),
#'  covariate matrix A(R,N),
#'  basis matrix X(P,Q) and Q<=min(P,N),
#'  parameter matrix C(Q,R)
#'  and coefficient matrix B(Q,N).
#'  Note that Y(N,P) and A(R,N) are given, and X(P,Q) and C(Q,R) are unknown.
#' @param Y observation matrix
#' @param A covariate matrix. Without covariate, identity matrix is used.
#' @param Q rank of basis matrix and Q<=min(P,N) where Y(P,N)
#' @param div number of partition usually described as "k" of k-fold
#' @param seed integer used as argument in set.seed function
#' @param ... arguments to be passed to nmfkc function.
#'  which controls the reproducibility of the partition.
#' @return objfunc: last objective function
#' @return objfunc.block: objective function at each block
#' @return block: partition block index (1,...,div) assigned to each column of Y
#' @return r.squared: coefficient of determination R^2, squared correlation between Y and XB
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
#' betas <- 20:40/1000
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

nmfkc.cv <- function(Y,A=diag(ncol(Y)),Q=2,div=5,seed=123,...){
  arglist=list(...)
  gamma <- ifelse("gamma" %in% names(arglist),arglist$gamma,0)
  epsilon <- ifelse("epsilon" %in% names(arglist),arglist$epsilon,1e-4)
  maxit <- ifelse("maxit" %in% names(arglist),arglist$maxit,5000)
  method <- ifelse("method" %in% names(arglist),arglist$method,"EU")
  X.column <- ifelse("X.column" %in% names(arglist),arglist$X.column,"sum")
  nstart <- ifelse("nstart" %in% names(arglist),arglist$nstart,1)
  print.trace <- ifelse("print.trace" %in% names(arglist),arglist$print.trace,FALSE)
  print.dims <- ifelse("print.dims" %in% names(arglist),arglist$print.dims,FALSE)
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
    A <- diag(ncol(Y))
    oldSum <- 0
    for(l in 1:maxit){
      XB <- X %*% C
      if(method=="EU"){
        C <- C*z((t(X)%*%Y)/(t(X)%*%XB+gamma*C))
      }else{
        C0 <- t(X)%*%z(Y/XB)%*%t(A)
        C <- C*z(C0/(colSums(X)%o%rowSums(A)+2*gamma*C))
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
    colnames(C) <- rownames(A)
    rownames(C) <- colnames(X)
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
    res_j <- nmfkc(Y_j,A_j,Q,gamma,epsilon,maxit,method,X.column,nstart,print.trace,print.dims)
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
    Yvec <- c(Yvec,as.vector(Yj))
    XBvec <- c(XBvec,as.vector(XBj))
  }
  objfunc <- sum(objfunc.block)
  r.squared <- stats::cor(XBvec,Yvec)^2
  return(list(objfunc=objfunc,objfunc.block=objfunc.block,block=index,r.squared=r.squared))
}


#' @title Rank selection diagnostics using figure
#' @description \code{nmfkc.rank} provides rank selection diagnostics. The method is still under development.
#' @param Y observation matrix
#' @param A covariate matrix. Without covariate, identity matrix is used.
#' @param Q vector of ranks to be diagnosed.
#' @param draw.figure draw a diagram for diagnosis
#' @param hclust.method option of hclust for calculating Cophenetic distances
#' @param ... arguments to be passed to nmfkc function.
#' @return r.squared
#' @return CPCC: Cophenetic correlation coefficient based on B.prob
#' @export
#' @references Brunet, J.P., Tamayo, P., Golub, T.R., Mesirov, J.P. (2004) Metagenes and molecular pattern discovery using matrix factorization. Proc. Natl. Acad. Sci. USA 2004, 101, 4164–4169. \url{https://doi.org/10.1073/pnas.0308531101}
#' @references Punera, K. and Ghosh, J. (2008). CONSENSUS-BASED ENSEMBLES OF SOFT CLUSTERINGS. Applied Artificial Intelligence, 22(7–8), 780–810. \url{https://doi.org/10.1080/08839510802170546}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' library(nmfkc)
#' Y <- t(iris[,-5])
#' nmfkc.rank(Y,Q=2:4)

nmfkc.rank <- function(Y,A=diag(ncol(Y)),Q=2:min(5,ncol(Y),nrow(Y)),draw.figure=TRUE,hclust.method="average",...){
  arglist=list(...)
  gamma <- ifelse("gamma" %in% names(arglist),arglist$gamma,0)
  epsilon <- ifelse("epsilon" %in% names(arglist),arglist$epsilon,1e-4)
  maxit <- ifelse("maxit" %in% names(arglist),arglist$maxit,5000)
  method <- ifelse("method" %in% names(arglist),arglist$method,"EU")
  X.column <- ifelse("X.column" %in% names(arglist),arglist$X.column,"sum")
  nstart <- ifelse("nstart" %in% names(arglist),arglist$nstart,1)
  print.trace <- ifelse("print.trace" %in% names(arglist),arglist$print.trace,FALSE)
  print.dims <- ifelse("print.dims" %in% names(arglist),arglist$print.dims,TRUE)
  myrmultinom <- function(prob.mat) stats::rmultinom(n=1,size=1,prob=prob.mat)
  r.squared <- 0*Q
  names(r.squared) <- Q
  correlation <- 0*Q
  names(correlation) <- Q
  for(q in 1:length(Q)){
    result <- nmfkc(Y,A,Q=Q[q],gamma,epsilon,maxit,method,X.column,nstart,print.trace,print.dims)
    r.squared[q] <- result$r.squared
    M <- t(result$B.prob) %*% result$B.prob
    up <- upper.tri(M)
    h.dist <- as.matrix(stats::cophenetic(stats::hclust(stats::as.dist(1-M),method=hclust.method)))
    correlation[q] <- stats::cor(h.dist[up],(1-M)[up])
  }
  if(draw.figure){
    graphics::par(mar=c(5,4,4,4)+0.1)
    plot(Q,correlation,type="l",col=2,axes=F,xlab="Rank",ylab="")
    graphics::text(Q,correlation,Q)
    graphics::axis(side=1,at=Q)
    graphics::axis(side=2,col=2,col.axis=2)
    graphics::mtext(side=2,text="Cophenetic correlation coefficient (CPCC)",line=2.5,col=2)
    graphics::par(new=T)
    plot(Q,r.squared,type="l",col=4,axes=F,xlab="",ylab="")
    graphics::text(Q,r.squared,Q)
    graphics::axis(side=4,col=4,col.axis=4)
    graphics::mtext(side=4,text="R.squared",line=2.5,col=4)
    graphics::box()
  }
  invisible(list(Q=Q,r.squared=r.squared,CPCC=correlation))
}
